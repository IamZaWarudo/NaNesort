"""Bateman half-life fitter - local Streamlit interface.

Run it with:

    streamlit run app.py

It opens in your browser at http://localhost:8501 and talks to nothing but your
own machine. The fit itself lives in bateman_fit.py, which is the same code the
web tool runs - verified to produce bit-identical results.
"""

import io
import json
import math
import os
import shutil
import tempfile
import zipfile

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

import grid_scan
from bateman_fit import run_fit

st.set_page_config(page_title="Bateman Half-Life Fitter", page_icon="⏱", layout="wide")

MODE_HELP = ("free = fitted with no constraint;  fixed = held at the value you type;  "
             "constrained = fitted but pulled towards your value by a Gaussian penalty "
             "of width sigma, so its uncertainty propagates into the parent's syst error.")

DEFAULT_CHAIN = pd.DataFrame([
    {"label": "44S",  "t12": 118.2, "t12Sigma": 1.6, "fixed": "constrained",
     "eff": 1.0, "effFixed": "fixed", "parent": 0, "br": 100.0, "brSigma": 0.0},
    {"label": "44Cl", "t12": 700.0, "t12Sigma": 70.0, "fixed": "fixed",
     "eff": 1.0, "effFixed": "free", "parent": 1, "br": 80.0, "brSigma": 0.0},
    {"label": "43Cl", "t12": 3130.0, "t12Sigma": 90.0, "fixed": "fixed",
     "eff": 1.0, "effFixed": "fixed", "parent": 1, "br": 20.0, "brSigma": 5.0},
])

# the engine wants "prior"; "constrained" is the word the interface uses
MODE_OUT = {"free": "free", "fixed": "fixed", "constrained": "prior"}
MODE_IN = {v: k for k, v in MODE_OUT.items()}


# --------------------------------------------------------------------- helpers
def fmt(value, err=None, stat=None, syst=None, unit=""):
    if value is None or not np.isfinite(value):
        return "—"
    s = f"{value:.6g}"
    if stat is not None and syst is not None and np.isfinite(stat) and np.isfinite(syst):
        s += f" ± {stat:.3g} stat ± {syst:.3g} syst"
    elif err is not None and np.isfinite(err):
        s += f" ± {err:.3g}"
    return s + (f" {unit}" if unit else "")


def load_config(cfg):
    """Push a saved settings file into the controls, so the sidebar and the chain
    table come up showing exactly what the file says. Keys the file omits keep
    whatever they are now."""
    def put(widget_key, cfg_key, cast=None, default=None):
        if cfg_key not in cfg:
            return
        v = cfg[cfg_key]
        if v is None and default is not None:
            v = default
        st.session_state[widget_key] = cast(v) if (cast and v is not None) else v

    put("w_unit", "timeUnit")
    put("w_edge", "edgeRef")
    put("w_rebin", "rebin", int)
    put("w_crop", "cropToWindow", bool)
    put("w_tmin", "tMin", float)
    put("w_tmax", "tMax", float)
    put("w_nbins", "nBins", int)
    put("w_binmode", "binMode")
    put("w_fitbkg", "fitBackground", bool)
    put("w_bkgmode", "bkgMode")
    put("w_autostart", "autoStart", bool)
    put("w_bkgrate", "bkgRate", float)
    put("w_bkgt12", "bkgT12", float, default=0.0)
    put("w_bkgt12fixed", "bkgT12Fixed")
    put("w_bkgt12sigma", "bkgT12Sigma", float, default=0.0)
    put("w_tailfit", "bkgTailFit", bool)
    put("w_taillo", "bkgTailMin", float)
    put("w_tailhi", "bkgTailMax", float)
    put("w_tie", "tieSiblings", bool)
    put("w_capeff", "capEff", bool)
    put("w_tag", "tagOutputs", bool)
    put("w_scan", "scanParent", bool)
    put("w_scanlo", "scanMin", float)
    put("w_scanhi", "scanMax", float)
    put("w_scann", "scanPoints", int)

    if cfg.get("useDemoData"):
        st.session_state["w_source"] = "Demo data"

    if cfg.get("chain"):
        rows = []
        for m in cfg["chain"]:
            rows.append({
                "label": str(m.get("label", "")),
                "t12": float(m.get("t12", 0.0)),
                "t12Sigma": float(m.get("t12Sigma") or 0.0),
                # the file says "prior"; the table says "constrained"
                "fixed": MODE_IN.get(str(m.get("fixed", "free")), "free"),
                "eff": float(m.get("eff", 1.0)),
                "effFixed": str(m.get("effFixed", "fixed")),
                "parent": int(m.get("parent") or 0),
                "br": float(m.get("br", 100.0)),
                "brSigma": float(m.get("brSigma") or 0.0),
            })
        st.session_state.chain = pd.DataFrame(rows)
        # drop the editor's pending edits, or they reapply on top of the new rows
        st.session_state.pop("chain_editor", None)


def levels(text, cast=float):
    """'2, 5, 10' or '2 5 10' -> [2.0, 5.0, 10.0]. Blank means 'derive them'."""
    out = []
    for tok in str(text or "").replace(",", " ").split():
        try:
            out.append(cast(float(tok)))
        except (TypeError, ValueError):
            pass
    return out


def bkg_levels(text):
    """Like levels(), but 'flat' (or 'none') survives as the no-slope level.

    '100000, 400000, flat' -> [100000.0, 400000.0, 'flat']. Blank means 'derive them'.
    """
    out = []
    for tok in str(text or "").replace(",", " ").split():
        if tok.strip().lower() in ("flat", "none", "inf"):
            out.append("flat")
            continue
        try:
            out.append(float(tok))
        except (TypeError, ValueError):
            pass
    return out


def walk_outputs(outdir):
    """Every file under outdir, as (path-relative-to-outdir, full path), sorted.

    Recursive on purpose: the grid scan writes a systematics/ subfolder, and a
    flat os.listdir would hand that directory to open() and raise IsADirectoryError.
    """
    out = []
    for root, _dirs, files in os.walk(outdir):
        for fn in files:
            full = os.path.join(root, fn)
            out.append((os.path.relpath(full, outdir), full))
    return sorted(out)


def zip_outputs(outdir):
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for rel, full in walk_outputs(outdir):
            z.write(full, rel)
    return buf.getvalue()


# -------------------------------------------------------------------- controls
st.title("Bateman Half-Life Fitter")
st.caption("Local build · the fit engine is identical to the web tool · nothing leaves this machine")

with st.sidebar:
    st.header("Settings file")
    cfg_up = st.file_uploader("Load a settings JSON, e.g. config_44S.json",
                              type=["json"], key="w_cfg",
                              help="The same file run_fit.py takes. Loading it fills in "
                                   "the chain table and every control below.")
    if cfg_up is not None and st.button("Load these settings", use_container_width=True):
        try:
            load_config(json.loads(cfg_up.getvalue().decode("utf-8")))
            st.success(f"Loaded {cfg_up.name}")
            st.rerun()
        except Exception as exc:
            st.error(f"Could not read that file: {exc}")

    st.header("Data")
    source = st.radio("Spectrum", ["Upload a file", "Path on disk", "Demo data"],
                      index=1, label_visibility="collapsed", key="w_source")
    upload = path = None
    if source == "Upload a file":
        upload = st.file_uploader("CSV: dt,counts  ·  or t_low,t_high,counts  ·  or an event list",
                                  type=["csv", "txt", "dat", "tsv"])
    elif source == "Path on disk":
        path = st.text_input("File path", "test_44S.csv", key="w_path")

    unit = st.selectbox("Time unit of everything you type", ["ms", "s", "us", "ns"],
                        index=0, key="w_unit")

    st.header("Reading the file")
    edge_ref = st.selectbox(
        "Two-column times are", ["center", "left"], index=0, key="w_edge",
        format_func=lambda v: "bin centres" if v == "center" else "left edges",
        help="How to turn a single time column into bin edges.")
    rebin = st.number_input("Rebin by", 1, 1000, 1, step=1, key="w_rebin",
                            help="Sum groups of this many channels before fitting.")
    crop = st.checkbox("Crop to the fit window below", False, key="w_crop")
    c1, c2 = st.columns(2)
    t_min = c1.number_input(f"Fit from ({unit})", value=2.0, format="%g", key="w_tmin")
    t_max = c2.number_input(f"Fit to ({unit})", value=5000.0, format="%g", key="w_tmax")
    st.caption("Window and binning below apply only to event lists; a pre-binned "
               "file keeps its own channels unless you tick crop.")
    c1, c2 = st.columns(2)
    n_bins = c1.number_input("Bins", 4, 4096, 64, step=1, key="w_nbins")
    bin_mode = c2.selectbox("Bin spacing", ["linear", "log"], index=0, key="w_binmode")

    st.header("Background")
    fit_bkg = st.checkbox("Fit the background level", True, key="w_fitbkg")
    bkg_mode = st.selectbox("Background shape", ["flat", "decaying"], index=1, key="w_bkgmode",
                            help="A decaying pedestal is one more chain component with "
                                 "its own half-life.")
    auto_start = st.checkbox("Estimate start values from the data", True, key="w_autostart")
    bkg_rate = st.number_input("Background start (counts/s)", value=8.0, format="%g",
                               disabled=auto_start, key="w_bkgrate")
    bkg_t12 = st.number_input(f"Background T½ start ({unit}), 0 = very long",
                              value=0.0, format="%g", disabled=bkg_mode == "flat",
                              key="w_bkgt12")
    bkg_t12_fixed = st.selectbox("Background T½", ["free", "fixed", "constrained"], index=0,
                                 disabled=bkg_mode == "flat", key="w_bkgt12fixed",
                                 help="constrained floats the slope inside a Gaussian prior of the "
                                      "width below - for a weak gate whose own tail cannot measure "
                                      "the pedestal slope, carrying a value and uncertainty "
                                      "measured on a high-statistics spectrum.")
    bkg_t12_sigma = st.number_input(f"Background T½ uncertainty ({unit}), 1σ",
                                    value=0.0, format="%g",
                                    disabled=(bkg_mode == "flat"
                                              or bkg_t12_fixed != "constrained"),
                                    key="w_bkgt12sigma")
    tail_fit = st.checkbox("Take the background slope from the tail", True,
                           disabled=bkg_mode == "flat", key="w_tailfit",
                           help="Fit b·exp(−λt) to the tail bins alone, then hold that "
                                "half-life, so the pedestal's slope cannot trade against "
                                "the parent's.")
    c1, c2 = st.columns(2)
    tail_lo = c1.number_input(f"Tail from ({unit})", value=3500.0, format="%g",
                              disabled=not tail_fit, key="w_taillo")
    tail_hi = c2.number_input(f"Tail to ({unit})", value=5000.0, format="%g",
                              disabled=not tail_fit, key="w_tailhi")

    st.header("Options")
    tie_sib = st.checkbox("Tie sibling branches to 100 %", True, key="w_tie",
                          help="The dominant branch out of each parent is set to "
                               "100 % minus the others, so only the minor branch is fitted.")
    cap_eff = st.checkbox("Cap efficiencies at 1", False, key="w_capeff")
    tag_out = st.checkbox("Name output files after the parent", True, key="w_tag")

    st.header("Parent half-life scan")
    scan = st.checkbox("Scan the parent's half-life", False, key="w_scan",
                       help="Hold T½(parent) across a grid, re-minimise everything "
                            "else, and report the cost in Cash.")
    c1, c2, c3 = st.columns(3)
    scan_lo = c1.number_input(f"from ({unit})", value=0.0, format="%g",
                              disabled=not scan, key="w_scanlo")
    scan_hi = c2.number_input(f"to ({unit})", value=0.0, format="%g",
                              disabled=not scan, key="w_scanhi")
    scan_n = c3.number_input("points", 5, 41, 15, step=1, disabled=not scan, key="w_scann")

    st.header("Systematic grid scan")
    st.caption(
        "Re-runs this same fit once per analysis choice you could defensibly have made "
        "differently, throws out the runs that no longer describe the data, and reduces "
        "what survives to one ΔT½ per axis and a quadrature total. **Run grid scan** is "
        "next to **Run fit** below.")
    g_mode = st.selectbox(
        "Topology", ["oat", "factorial"], index=0, key="w_gmode",
        format_func=lambda v: ("one axis at a time" if v == "oat" else "full factorial"),
        help="A tornado chart and a quadrature budget are by construction one-axis-at-a-"
             "time about a reference. Factorial looks for interactions between choices, "
             "but its per-axis numbers are main effects, not clean attributions.")

    with st.expander("Axes", expanded=False):
        st.caption("Leave a level list blank to derive it from the reference window.")
        g_ws = st.checkbox("Fit window start", True, key="w_gws")
        g_ws_v = st.text_input(f"tMin levels ({unit})", "", key="w_gwsv",
                               placeholder=", ".join(f"{v:g}" for v in
                                                     [t_min, t_min * 2.5, t_min * 5, t_min * 10]))
        g_we = st.checkbox("Fit window end", True, key="w_gwe")
        g_we_v = st.text_input(f"tMax levels ({unit})", "", key="w_gwev",
                               placeholder=", ".join(f"{v:g}" for v in
                                                     [t_max * 0.2, t_max * 0.4,
                                                      t_max * 0.7, t_max]))
        st.caption("Both window axes force **crop to window** on. Without it tMin and tMax "
                   "do nothing to an already-binned file and the axis reports a spurious zero.")
        g_rb = st.checkbox("Rebin factor", True, key="w_grb")
        g_rb_v = st.text_input("Rebin levels", "1, 2, 4, 10", key="w_grbv")
        g_bm = st.checkbox("Bin spacing (linear / log)", True, key="w_gbm",
                           help="Legitimately reads zero on a file that arrives already binned.")
        g_bg = st.checkbox("Background treatment", True, key="w_gbg",
                           help="flat floated · decaying with the tail-fitted slope · "
                                "decaying with a free slope")
        g_dt = st.checkbox("Daughter half-lives", True, key="w_gdt",
                           help="literature prior · fixed · free. The prior level is skipped, "
                                "and reported as skipped, on any row with no σ on its value.")
        g_bt = st.checkbox("Background half-life", False, key="w_gbt",
                           help="Explicit pedestal half-lives, each held fixed while only its "
                                "amplitude floats, plus a flat level. Use it instead of "
                                "Background treatment once a tail fit has told you which "
                                "slopes are allowed: the three flat/decaying labels charge "
                                "you for slopes the tail has already excluded.")
        g_bt_v = st.text_input(f"T½(bkg) levels ({unit}, or 'flat')", "", key="w_gbtv",
                               placeholder=", ".join(
                                   [f"{t_max * k:g}" for k in (20, 40, 80, 200)] + ["flat"]))
        g_ce = st.checkbox("Efficiency cap", True, key="w_gce")

    with st.expander("Admission rule", expanded=False):
        g_chi = st.number_input("Reject above χ²/ndf", 1.0, 100.0, 3.0, step=0.5,
                                format="%g", key="w_gchi")
        g_p = st.number_input("Reject below p", 0.0, 0.5, 0.001, step=0.001,
                              format="%g", key="w_gp")
        g_strict = st.checkbox("Also require the minimiser's convergence flag", False,
                               key="w_gstrict",
                               help="Off by default: on this likelihood the flag fires on the "
                                    "near-degenerate ε/N₀ direction even at a good minimum, and "
                                    "it threw out fits with χ²/ndf = 1.03. It is recorded on "
                                    "every row either way.")
        g_max = st.number_input("Max runs (factorial cap)", 1, 5000, 200, step=10,
                                key="w_gmax", disabled=g_mode != "factorial")
        st.caption("Rejected runs stay in the CSV with their reason — nothing disappears.")

    # filled in below, once params and grid_spec exist, with the exact run list
    grid_plan_slot = st.empty()

    st.header("Output")
    outdir = st.text_input("Output folder", "out", key="w_outdir")

# ------------------------------------------------------------------ chain table
st.subheader("Decay chain")
st.caption(
    "One row per member, the parent first. **Fed by** is the row number that feeds this "
    "one (0 = the row above), so a β⁻n daughter and a β⁻ daughter both point at the same "
    "parent and are parallel branches, not successive generations. **Branch σ** turns that "
    "branching into a fitted parameter with a Gaussian prior — which is how a poorly known "
    "Pn ends up in the parent's systematic error."
)

if "chain" not in st.session_state:
    st.session_state.chain = DEFAULT_CHAIN.copy()

chain_df = st.data_editor(
    st.session_state.chain,
    num_rows="dynamic",
    use_container_width=True,
    key="chain_editor",
    column_config={
        "label": st.column_config.TextColumn("Nuclide", width="small"),
        "t12": st.column_config.NumberColumn(f"T½ ({unit})", format="%g"),
        "t12Sigma": st.column_config.NumberColumn(f"T½ σ ({unit})", format="%g"),
        "fixed": st.column_config.SelectboxColumn(
            "T½ mode", options=["free", "fixed", "constrained"], help=MODE_HELP),
        "eff": st.column_config.NumberColumn("Efficiency", format="%g"),
        "effFixed": st.column_config.SelectboxColumn("ε mode", options=["free", "fixed"]),
        "parent": st.column_config.NumberColumn("Fed by #", format="%d", min_value=0),
        "br": st.column_config.NumberColumn("Branch %", format="%g"),
        "brSigma": st.column_config.NumberColumn("Branch σ %", format="%g"),
    },
)

chain = []
for _, r in chain_df.iterrows():
    if not str(r.get("label", "")).strip():
        continue
    chain.append({
        "label": str(r["label"]).strip(),
        "t12": float(r["t12"]),
        "t12Sigma": float(r.get("t12Sigma") or 0.0),
        "fixed": MODE_OUT.get(str(r.get("fixed", "free")), "free"),
        "eff": float(r.get("eff") if r.get("eff") is not None else 1.0),
        "effFixed": str(r.get("effFixed", "fixed")),
        "parent": int(r.get("parent") or 0),
        "br": float(r.get("br") if r.get("br") is not None else 100.0),
        "brSigma": float(r.get("brSigma") or 0.0),
    })

params = {
    "chain": chain, "timeUnit": unit,
    "useDemoData": source == "Demo data",
    "tMin": t_min, "tMax": t_max, "cropToWindow": crop,
    "nBins": int(n_bins), "binMode": bin_mode, "edgeRef": edge_ref, "rebin": int(rebin),
    "fitBackground": fit_bkg, "bkgRate": bkg_rate, "autoStart": auto_start,
    "bkgMode": bkg_mode, "bkgT12": (bkg_t12 or None), "bkgT12Fixed": bkg_t12_fixed,
    "bkgT12Sigma": (bkg_t12_sigma or None),
    "bkgTailFit": tail_fit and bkg_mode == "decaying",
    "bkgTailMin": tail_lo, "bkgTailMax": tail_hi,
    "capEff": cap_eff, "tieSiblings": tie_sib, "tagOutputs": tag_out,
    "scanParent": scan, "scanMin": scan_lo, "scanMax": scan_hi,
    "scanPoints": int(scan_n),
}

grid_spec = {
    "mode": g_mode,
    "maxRedChi2": float(g_chi),
    "minPValue": float(g_p),
    "maxRuns": int(g_max),
    "strictConvergence": bool(g_strict),
    "axes": {
        "window_start": {"enabled": g_ws, "values": levels(g_ws_v) or "auto"},
        "window_end": {"enabled": g_we, "values": levels(g_we_v) or "auto"},
        "rebin": {"enabled": g_rb, "values": levels(g_rb_v, int) or [1, 2, 4, 10]},
        "bin_mode": {"enabled": g_bm, "values": ["linear", "log"]},
        "background": {"enabled": g_bg,
                       "values": ["flat_floated", "decay_tail_slope", "decay_free_slope"]},
        "daughter_t12": {"enabled": g_dt, "values": ["prior", "fixed", "free"]},
        "bkg_t12": {"enabled": g_bt, "values": bkg_levels(g_bt_v) or "auto"},
        "cap_eff": {"enabled": g_ce, "values": [False, True]},
    },
}

def grid_plan():
    """What Run grid scan would do. (plan, None) or (None, message to show instead)."""
    try:
        import grid_scan
    except ImportError as exc:
        return None, (f"Copy **grid_scan.py** into this folder to see the run list "
                      f"and to use the scan ({exc}).")
    if not hasattr(grid_scan, "plan_preview"):
        return None, ("The grid_scan.py beside app.py predates the run-list preview. "
                      "Copy the current one over it.")
    try:
        return grid_scan.plan_preview(params, grid_spec), None
    except Exception as exc:                 # never leave the slot silently blank
        return None, f"Could not work out the run list: {type(exc).__name__}: {exc}"


plan, plan_msg = grid_plan()
with grid_plan_slot.container():
    if plan is None:
        st.warning(plan_msg)
    else:
        st.markdown(f"**This scan: {plan['nRuns']} fits** "
                    f"({'one axis at a time' if plan['mode'] == 'oat' else 'full factorial'})")
        if plan["nCapped"]:
            st.warning(f"{plan['nCapped']} further combinations are above the cap and "
                       f"will not be run.")
        if not plan["axes"]:
            st.info("Every axis is switched off — the scan would only fit the reference.")
        for ax in plan["axes"]:
            # " · ", not ", " — several level names contain a comma of their own
            st.markdown(f"- **{ax['title']}** — {' · '.join(ax['levels']) or '—'}")
            for s in ax["skipped"]:
                st.caption(f"    skipped: {s}")
        st.caption("The reference fit is the first of those runs, and it is the "
                   "configuration set above — so its result also fills the other tabs.")

b1, b2, b3 = st.columns([2, 1, 1])
run = b1.button("Run fit", type="primary", use_container_width=True)
run_grid = b2.button("Run grid scan", use_container_width=True,
                     help="One fit per analysis choice, then the tornado and the budget. "
                          "Results appear in the Systematics tab.")
b3.download_button(
    "Save these settings", json.dumps(params, indent=2),
    file_name=f"config_{chain[0]['label'] if chain else 'chain'}.json",
    mime="application/json", use_container_width=True,
    help="Writes exactly what is on screen now. Reload it here, or hand it to run_fit.py.")

# ------------------------------------------------------------------------- run
def resolve_spectrum():
    """The spectrum argument, or an error string to show instead."""
    if source == "Upload a file":
        if upload is None:
            return None, "Choose a file, or switch the picker to demo data."
        return upload, None
    if source == "Path on disk":
        if not os.path.exists(str(path)):
            return None, (f"No file at {path!r}. Give a path relative to where you launched "
                          f"Streamlit, or an absolute one.")
        return path, None
    return None, None


if run:
    spec, err = resolve_spectrum()
    if err:
        st.error(err)
        st.stop()

    st.session_state.chain = chain_df.copy()
    with st.spinner("Minimising…"):
        st.session_state.res = run_fit(params, outdir=outdir, spectrum=spec)
        st.session_state.outdir = outdir
        st.session_state.params = params

if run_grid:
    spec, err = resolve_spectrum()
    if err:
        st.error(err)
        st.stop()

    try:
        import grid_scan
    except ImportError as exc:
        st.error(f"The grid scan needs grid_scan.py beside app.py ({exc}).")
        st.stop()

    # Every grid point re-reads the spectrum. A path can be read any number of times;
    # an uploaded file handle cannot - it is exhausted after the first fit, and every
    # later point would silently fall back to demo data. So snapshot it to bytes.
    blob = None
    if spec is not None and not isinstance(spec, (str, bytes, os.PathLike)):
        try:
            spec.seek(0)
        except Exception:
            pass
        blob = spec.read()
        if isinstance(blob, str):
            blob = blob.encode("utf-8")

    st.session_state.chain = chain_df.copy()
    gdir = os.path.join(outdir, "systematics")
    scratch = tempfile.mkdtemp(prefix="bateman_grid_")
    captured = {}
    status = st.status("Scanning…", expanded=True)
    line = status.empty()
    n = [0]

    def fitfn(p, run_id):
        n[0] += 1
        # the reference IS this configuration, so let it write the normal output
        # folder and keep its result - that saves refitting it for the other tabs
        ref = run_id == "ref"
        d = outdir if ref else os.path.join(scratch, str(run_id).replace(":", "_"))
        try:
            r = run_fit(p, outdir=d,
                        spectrum=(io.BytesIO(blob) if blob is not None else spec))
            if ref:
                captured["ref"] = r
            return r
        except Exception as exc:                 # a dead point must not kill the scan
            return {"error": f"{type(exc).__name__}: {exc}"}
        finally:
            if not ref:
                shutil.rmtree(d, ignore_errors=True)

    try:
        grid = grid_scan.run_grid(
            params, grid_spec, fitfn,
            progress=lambda msg: line.write(f"fit {n[0] + 1} · {msg}"))
    finally:
        shutil.rmtree(scratch, ignore_errors=True)

    if not grid.get("error"):
        grid["written"] = grid_scan.write_outputs(grid, gdir, figure=True)
        grid["summaryTable"] = grid_scan.summary_table(grid)
    status.update(label=(grid.get("error") or
                         f"{grid['nAdmitted']} of {grid['nRuns']} runs admitted"),
                  state="error" if grid.get("error") else "complete", expanded=False)

    st.session_state.grid = grid
    st.session_state.griddir = gdir
    st.session_state.params = params
    if captured.get("ref") is not None and "error" not in captured["ref"]:
        st.session_state.res = captured["ref"]
        st.session_state.outdir = outdir

# --------------------------------------------------------------------- results
res = st.session_state.get("res")
if res is None:
    st.info("Set up the chain above and press **Run fit**. "
            "Pick *Demo data* in the sidebar to try it with a synthetic spectrum first.")
    st.stop()

if "error" in res:
    st.error(res["error"])
    with st.expander("Traceback"):
        st.code(res.get("traceback", ""))
    st.stop()

UF, unit_name = res["unitFactor"], res["timeUnit"]
sp = res["spectrum"]
x = np.array(sp["center"]) / UF
counts = np.array(sp["counts"])
mu = np.array(sp["mu"])
bkg = np.array(sp["bkg"])
st_ = res["stats"]

k1, k2, k3, k4 = st.columns(4)
p = res["fit"][0]
k1.metric(f"T½ ({res['parentLabel']})",
          f"{p['t12']/UF:.5g} {unit_name}",
          f"± {p['sigma']/UF:.3g}" if p.get("sigma") else None)
k2.metric("N₀ correlated implants", f"{res['N0']:.5g}",
          f"± {res['N0Err']:.3g}" if res.get("N0Err") else None)
k3.metric("A₀ at t = 0", f"{res['A0']:.5g} /s",
          f"± {res['A0Err']:.3g}" if res.get("A0Err") else None)
k4.metric("χ²/ndf", f"{st_['redChi2']:.4f}", f"p = {st_['pValue']:.3g}")

tabs = st.tabs(["Fit", "Residuals", "Parameters", "Profiles", "Parent scan",
                "Correlation", "Notes", "Systematics", "Files"])

# ------------------------------------------------------------------------ fit
with tabs[0]:
    c1, c2, c3 = st.columns([1, 1, 2])
    subtract = c1.checkbox("Subtract background", False)
    logy = c2.checkbox("Log y", True)

    fig = go.Figure()
    y_d = counts - bkg if subtract else counts
    fig.add_trace(go.Scatter(
        x=x, y=y_d, mode="markers", name="data",
        marker=dict(size=3, color="#555"),
        error_y=dict(type="data", array=np.sqrt(np.maximum(counts, 0)),
                     visible=len(x) < 300, width=0, thickness=0.6)))
    fig.add_trace(go.Scatter(
        x=x, y=(mu - bkg if subtract else mu), mode="lines", name="total model",
        line=dict(color="#c0392b", width=2.5)))
    if not subtract:
        fig.add_trace(go.Scatter(x=x, y=bkg, mode="lines", name="background",
                                 line=dict(color="#7f8c8d", width=1.5, dash="dash")))
    for r, comp in zip(res["fit"], sp["comps"]):
        fig.add_trace(go.Scatter(x=x, y=np.array(comp), mode="lines", name=r["label"],
                                 line=dict(width=1.4)))
    fig.update_layout(
        height=520, hovermode="x unified", margin=dict(l=10, r=10, t=30, b=10),
        xaxis_title=f"implant–decay time difference ({unit_name})",
        yaxis_title="counts / channel",
        legend=dict(orientation="h", y=1.08))
    if logy:
        fig.update_yaxes(type="log")
    if res["binMode"] == "log":
        fig.update_xaxes(type="log")
    st.plotly_chart(fig, use_container_width=True)

    if res.get("flatMembers"):
        st.warning("Flat over this window, so degenerate with the pedestal: "
                   + ", ".join(res["flatMembers"]))

# ------------------------------------------------------------------ residuals
with tabs[1]:
    kind = st.radio("Residual", ["pearson", "deviance"], horizontal=True,
                    format_func=lambda v: "Pearson (data−fit)/σ" if v == "pearson"
                    else "Signed deviance")
    r = np.array(res["residuals"][kind])
    fig = go.Figure()
    fig.add_trace(go.Bar(x=x, y=r, marker_color=np.where(r >= 0, "#3884ff", "#e65046"),
                         name=kind))
    for lev, col in ((0, "#888"), (1, "#ccc"), (-1, "#ccc"), (3, "#eee"), (-3, "#eee")):
        fig.add_hline(y=lev, line_color=col, line_width=1)
    fig.update_layout(height=380, margin=dict(l=10, r=10, t=30, b=10),
                      xaxis_title=f"time difference ({unit_name})", yaxis_title="σ",
                      showlegend=False)
    if res["binMode"] == "log":
        fig.update_xaxes(type="log")
    st.plotly_chart(fig, use_container_width=True)
    st.caption(f"Cash = {st_['cash']:.1f} · χ² = {st_['chi2']:.1f} / {st_['ndf']} "
               f"· {st_['nBins']} bins, {st_['nParams']} free parameters · "
               f"{st_['totalCounts']:.6g} counts in the window")

# ----------------------------------------------------------------- parameters
with tabs[2]:
    def unc(r):
        """± stat ± syst when the split exists, otherwise the single total."""
        if not r.get("sigma"):
            return "—"
        stat, syst = r.get("sigmaStat"), r.get("sigmaSyst")
        if stat and syst:
            return f"± {stat/UF:.3g} stat ± {syst/UF:.3g} syst"
        return f"± {r['sigma']/UF:.3g}"

    rows = []
    for r in res["fit"]:
        rows.append({
            "nuclide": r["label"],
            f"T½ ({unit_name})": r["t12"] / UF,
            "T½ uncertainty": unc(r),
            "T½ status": "fixed" if r["fixed"] else ("constrained" if r["constrained"] else "free"),
            "profile ±1σ": (f"{r['loBound']/UF:.5g} – {r['hiBound']/UF:.5g}"
                            if r.get("loBound") and r.get("hiBound") else "—"),
            "fed by": r.get("fedBy") or "—",
            "branch %": (f"{r['branchPercent']:.4g}" if r.get("branchPercent") is not None else "—"),
            "branch ±": (f"{r['branchErr']:.3g}" if r.get("branchErr") else "—"),
            "branch status": ("free (Pn prior)" if r["branchFree"]
                              else "tied (100 % − siblings)" if r.get("branchTied")
                              else "fixed"),
            "efficiency": f"{r['eff']:.4g}" + (f" ± {r['effErr']:.2g}" if r.get("effErr") else ""),
            "ε status": "free" if r["effFree"] else "fixed",
            "fitted counts": f"{r['counts']:.6g}",
        })
    st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

    st.markdown("**Derived quantities**")
    st.dataframe(pd.DataFrame([
        {"quantity": f"N₀ (correlated {res['parentLabel']} implants)",
         "value": fmt(res["N0"], res.get("N0Err"), res.get("N0ErrStat"), res.get("N0ErrSyst"))},
        {"quantity": "A₀ (initial activity)",
         "value": fmt(res["A0"], res.get("A0Err"), res.get("A0ErrStat"), res.get("A0ErrSyst"),
                      "decays/s")},
        {"quantity": f"A₀ per {res['binWidth']/UF:.4g} {unit_name} channel",
         "value": fmt(res.get("A0PerBin"), res.get("A0PerBinErr"), unit="counts")},
        {"quantity": "background amplitude at t = 0",
         "value": fmt(res.get("bkg"), res.get("bkgErr"), unit="counts/s")},
        {"quantity": "background T½",
         "value": (fmt(res["bkgT12"] / UF,
                       (res["bkgT12Err"] / UF) if res.get("bkgT12Err") else None, unit=unit_name)
                   + ("  (lower limit — flat over this window)" if res.get("bkgT12Plateau") else "")
                   + (f"  (constrained to {res['bkgT12Prior']/UF:.6g} ± "
                      f"{res['bkgT12PriorSigma']/UF:.6g} {unit_name})"
                      if res.get("bkgT12Constrained") else ""))
                  if res.get("bkgT12") else "—"},
        {"quantity": "ρ(ln T½, ln N₀)",
         "value": f"{res['A0Rho']:+.4f}" if res.get("A0Rho") is not None else "—"},
    ]), use_container_width=True, hide_index=True)

# -------------------------------------------------------------------- profiles
with tabs[3]:
    if not res.get("profiles"):
        st.info("No profile scans — every half-life is fixed.")
    else:
        fig = go.Figure()
        for pr in res["profiles"]:
            fig.add_trace(go.Scatter(x=np.array(pr["t12Grid"]) / UF, y=pr["deltaC"],
                                     mode="lines+markers", name=pr["label"]))
        fig.add_hline(y=1, line_dash="dash", line_color="#888",
                      annotation_text="ΔC = 1  (1σ)")
        fig.add_hline(y=9, line_dash="dot", line_color="#bbb",
                      annotation_text="ΔC = 9  (3σ)")
        fig.update_layout(height=440, margin=dict(l=10, r=10, t=30, b=10),
                          xaxis_title=f"T½ ({unit_name})", yaxis_title="ΔCash",
                          yaxis_range=[0, 12])
        st.plotly_chart(fig, use_container_width=True)
        st.caption("Where the curve crosses ΔC = 1 is the asymmetric 1σ interval. "
                   "With millions of counts per channel ΔC is dominated by model "
                   "mismatch, so divide an apparent significance by √(χ²/ndf) "
                   f"= {math.sqrt(st_['redChi2']):.3f} before believing it.")

# ------------------------------------------------------------------ parent scan
with tabs[4]:
    sc = res.get("parentScan")
    if not sc:
        st.info("Tick **Scan the parent's half-life** in the sidebar and refit.")
    else:
        pts = sc["points"]
        t = np.array([q["t12"] for q in pts]) / UF
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=t, y=[q["deltaC"] for q in pts],
                                 mode="lines+markers", name="ΔCash",
                                 line=dict(color="#c0392b")))
        fig.add_hline(y=1, line_dash="dash", line_color="#888", annotation_text="1σ")
        fig.add_hline(y=9, line_dash="dot", line_color="#bbb", annotation_text="3σ")
        fig.add_vline(x=sc["best"] / UF, line_color="#3884ff",
                      annotation_text=f"best {sc['best']/UF:.4g}")
        fig.update_layout(height=400, margin=dict(l=10, r=10, t=30, b=10),
                          xaxis_title=f"T½ ({sc['label']}) held at ({unit_name})",
                          yaxis_title="ΔCash")
        st.plotly_chart(fig, use_container_width=True)

        st.dataframe(pd.DataFrame([{
            f"T½ ({unit_name})": q["t12"] / UF, "ΔCash": q["deltaC"],
            "N₀": q["N0"], "A₀ (/s)": q["A0"], "background (/s)": q["bkg"],
            **{f"branch {k} %": v for k, v in q["branch"].items()},
            **{f"ε {k}": v for k, v in q["eff"].items()},
        } for q in pts]), use_container_width=True, hide_index=True)

# ------------------------------------------------------------------ correlation
with tabs[5]:
    cv = res.get("covariance")
    if not cv or len(cv.get("names", [])) < 2:
        st.info("Fewer than two free parameters, or the Hessian was singular.")
    else:
        short = [nm.replace("T12_", "T½ ").replace("N0_", "N₀ ")
                   .replace("eff_", "ε ").replace("branch_", "br ")
                   .replace("_", " ") for nm in cv["names"]]
        M = np.array([[np.nan if v is None else v for v in row] for row in cv["corr"]])
        np.fill_diagonal(M, 1.0)
        fig = go.Figure(go.Heatmap(
            z=M, x=short, y=short, zmin=-1, zmax=1, colorscale="RdBu", reversescale=True,
            text=[[f"{v:+.3f}" if np.isfinite(v) else "—" for v in row] for row in M],
            texttemplate="%{text}", textfont=dict(size=11),
            colorbar=dict(title="ρ", thickness=14)))
        fig.update_layout(height=110 + 46 * len(short),
                          margin=dict(l=10, r=10, t=30, b=10),
                          yaxis=dict(autorange="reversed"))
        st.plotly_chart(fig, use_container_width=True)
        st.markdown(
            "Blue is positive, red negative, intensity is |ρ|. Every error quoted "
            "elsewhere is a diagonal element of this matrix; these are the off-diagonals. "
            "A derived quantity has variance **σ²(f) = wᵀ·cov·w**, with **w** the partial "
            "derivatives — adding errors in quadrature instead assumes every off-diagonal "
            "cell is zero. The full matrix, absolute and log-space, is in `covariance.csv`.")

# ---------------------------------------------------------------------- notes
with tabs[6]:
    for note in res["notes"]:
        st.markdown(f"- {note}")

# ----------------------------------------------------------------- systematics
with tabs[7]:
    grid = st.session_state.get("grid")
    if grid is None:
        st.info("Press **Run grid scan** above. It re-fits this configuration once per "
                "analysis choice — window edges, binning, background treatment, daughter "
                "half-lives, the efficiency cap — and turns the spread into a budget.")
    elif grid.get("error"):
        st.error(grid["error"])
    else:
        gunit = grid["timeUnit"]
        ref = grid["reference"]
        gd = st.session_state.get("griddir", "out/systematics")

        # The same axes do not cost the three quantities the same, so each gets its
        # own budget; pick which one is on display.
        _labels = {"t12": "T½", "N0": "N₀", "A0": "A₀"}
        _avail = [r for r in grid_scan.RESPONSES
                  if (ref.get("t12") if r == "t12" else ref.get(r)) is not None]
        resp = st.radio("Budget for", _avail, horizontal=True,
                        format_func=lambda r: _labels[r], key="grid_resp")
        rlab = _labels[resp]
        runit = grid_scan.response_unit(resp, gunit)
        rref = ref["t12"] if resp == "t12" else ref[resp]
        rstat = (ref.get("sigma") if resp == "t12"
                 else (ref.get(f"{resp}_stat") or ref.get(f"{resp}_sigma")))
        rsyst = grid.get(grid_scan.SIGMA_KEY[resp]) or 0.0
        dkey = f"delta_{resp}"

        m1, m2, m3 = st.columns(3)
        m1.metric(f"Reference {rlab} ({ref.get('parent') or ''})",
                  f"{rref:.5g} {runit}",
                  f"± {rstat:.3g} stat" if rstat else None)
        m2.metric("Quadrature systematic", f"± {rsyst:.3g} {runit}",
                  f"{100.0 * rsyst / abs(rref):.2f} %" if rref else None)
        m3.metric("Runs admitted", f"{grid['nAdmitted']} / {grid['nRuns']}",
                  f"mode: {grid['mode']}")
        st.markdown(f"### {rlab} = {rref:.6g} ± {rstat or 0:.3g} (stat) "
                    f"± {rsyst:.3g} (syst) {runit}")

        bars = []
        for b in grid["budget"]:
            if b.get(dkey) is None:
                continue
            vals = [r[resp] for r in grid["rows"]
                    if r.get("axis") == b["axis"] and r.get("admitted")
                    and r.get(resp) is not None]
            lo = (min(vals) - rref) if vals else -b[dkey]
            hi = (max(vals) - rref) if vals else b[dkey]
            bars.append((b["axis_label"], min(lo, 0.0), max(hi, 0.0), b[dkey],
                         b.get(f"worst_level_{resp}") or b.get("worst_level", "")))
        bars.sort(key=lambda z: z[3])

        if bars:
            fig = go.Figure()
            fig.add_trace(go.Bar(
                y=[z[0] for z in bars], x=[z[2] - z[1] for z in bars],
                base=[z[1] for z in bars], orientation="h", name="admitted envelope",
                marker=dict(color="#8fb4dd", line=dict(color="#2d5f92", width=1)),
                customdata=[[z[3], z[4]] for z in bars],
                hovertemplate="%{y}<br>Δ" + rlab + " = %{customdata[0]:.4g} " + runit +
                              "<br>worst: %{customdata[1]}<extra></extra>"))
            if rstat:
                fig.add_vrect(x0=-rstat, x1=rstat, fillcolor="#d9534f",
                              opacity=0.10, line_width=0,
                              annotation_text=f"±{rstat:.3g} stat",
                              annotation_position="top left")
            if rsyst:
                for s in (-rsyst, rsyst):
                    fig.add_vline(x=s, line_dash="dash", line_color="#2d5f92", line_width=1.2)
            fig.add_vline(x=0, line_color="#111", line_width=1.2)
            fig.update_layout(
                height=140 + 40 * len(bars), showlegend=False,
                margin=dict(l=10, r=10, t=40, b=10),
                xaxis_title=f"deviation of {rlab} from the reference fit ({runit})",
                title=dict(text=f"{ref.get('parent') or 'parent'}: "
                                f"{rlab} = {rref:.4g} {runit} (reference)", font=dict(size=13)))
            st.plotly_chart(fig, use_container_width=True)

        st.markdown("**Budget**")
        # all three deviations side by side, so an axis that is cheap for the yield
        # and ruinous for the half-life shows that in one glance
        def _pct(d, v):
            return (100.0 * d / abs(v)) if (d is not None and v) else None
        st.dataframe(pd.DataFrame([{
            "source": b["axis_label"], f"ΔT½ ({gunit})": b.get("delta_t12"),
            "ΔN₀ (%)": _pct(b.get("delta_N0"), ref.get("N0")),
            "ΔA₀ (%)": _pct(b.get("delta_A0"), ref.get("A0")),
            "runs admitted": b.get("n_admitted"),
            "worst level": b.get(f"worst_level_{resp}") or b.get("worst_level", ""),
            f"envelope ({runit})": (f"{b['lo_' + resp]:.6g} – {b['hi_' + resp]:.6g}"
                                    if b.get(f"lo_{resp}") is not None else "—"),
        } for b in grid["budget"]] + [{
            "source": "quadrature total", f"ΔT½ ({gunit})": grid.get("sigmaSyst"),
            "ΔN₀ (%)": _pct(grid.get("sigmaSystN0"), ref.get("N0")),
            "ΔA₀ (%)": _pct(grid.get("sigmaSystA0"), ref.get("A0")),
            "runs admitted": grid["nAdmitted"], "worst level": "",
            f"envelope ({runit})": "",
        }]), use_container_width=True, hide_index=True)

        with st.expander(f"Every grid point ({grid['nRuns']} runs)"):
            st.dataframe(pd.DataFrame(grid["rows"]), use_container_width=True,
                         hide_index=True)

        st.markdown("**Notes**")
        for note in grid.get("notes", []):
            st.markdown(f"- {note}")

        with st.expander("Summary table (paste into the paper)"):
            st.code(grid.get("summaryTable", ""))

        _stem = "systematic_tornado" if resp == "t12" else f"systematic_tornado_{resp}"
        png = os.path.join(gd, _stem + ".png")
        if os.path.exists(png):
            st.markdown("**Publication figure**")
            st.image(png)
        st.caption(f"Written to `{os.path.abspath(gd)}`")
        _files = ["systematic_grid_results.csv", "systematic_budget.csv"]
        for r in grid_scan.RESPONSES:
            s = "systematic_tornado" if r == "t12" else f"systematic_tornado_{r}"
            _files += [s + ".pdf", s + ".png"]
        _files = [f for f in _files if os.path.exists(os.path.join(gd, f))]
        cols = st.columns(4)
        for i, fn in enumerate(_files):
            full = os.path.join(gd, fn)
            with open(full, "rb") as fh:
                cols[i % 4].download_button(fn, fh.read(), file_name=fn,
                                            key=f"gdl_{fn}", use_container_width=True)

# ---------------------------------------------------------------------- files
with tabs[8]:
    od = st.session_state.get("outdir", "out")
    if not os.path.isdir(od):
        st.warning(f"No output folder at {od!r} yet.")
        st.stop()
    st.caption(f"Written to `{os.path.abspath(od)}`")
    for rel, full in walk_outputs(od):
        c1, c2 = st.columns([4, 1])
        c1.write(f"**{rel}** · {os.path.getsize(full):,} bytes")
        with open(full, "rb") as fh:
            c2.download_button("Download", fh.read(),
                               file_name=os.path.basename(rel), key=f"dl_{rel}")
    st.download_button("Download everything as a zip", zip_outputs(od),
                       file_name=f"{res['parentLabel'].strip() or 'fit'}_outputs.zip",
                       type="primary")
    with st.expander("Settings used (paste into a config file for run_fit.py)"):
        st.code(json.dumps(st.session_state.get("params", {}), indent=2), language="json")
