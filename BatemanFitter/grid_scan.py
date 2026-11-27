#!/usr/bin/env python3
"""Automated systematic parameter grid scans for the Bateman fitter.

A half-life's systematic uncertainty is not a property of one fit - it is the
spread of the answer over the analysis choices you could defensibly have made
instead.  This module makes that spread reproducible: it re-runs the same Cash
likelihood fit over a declared grid of analysis choices, throws out the runs
that did not converge or that no longer describe the data, and reduces what is
left to one number per axis and a quadrature total.

It does this for three quantities at once - the parent half-life T1/2, the
correlated-implant yield N0, and the initial activity A0 - because the same
axes do not cost them the same.  Over a window many half-lives long N0 is close
to a counted quantity, the excess over the pedestal, while T1/2 is the *shape*
of that excess; a background-slope assumption can therefore move the half-life
by 10 % and leave the yield at the per-mille level.  Each gets its own tornado
and its own quadrature total, and the one that belongs to a number is the one
to quote beside it.

The engine is untouched - every grid point goes through the same
``bateman_fit.run_fit`` the single fit uses, so a grid entry and a hand-typed
configuration give bit-identical numbers.

Two modes:

``oat``  (default, and what a tornado plot means)
    One axis is varied at a time from the reference configuration.  N runs for
    N declared levels, and the deviation of each run is attributable to exactly
    one choice.  This is the mode whose quadrature total is a systematic budget.

``factorial``
    The full Cartesian product of the enabled axes.  Use it to look for
    interactions between choices - it costs the product of the level counts and
    its per-axis numbers are main effects, not clean attributions, so the
    quadrature total from a factorial scan is a rough indicator rather than a
    budget.  Capped by ``maxRuns``; anything dropped is reported, never silent.

Admission rule
    A run is admitted only if the minimiser converged, the fitted parent
    half-life is finite and positive, chi2/ndf <= ``maxRedChi2`` and the Pearson
    p-value >= ``minPValue``.  A rejected run still appears in the CSV with the
    reason recorded - a configuration that will not fit is itself a result.

The window axes force ``cropToWindow`` on.  On an already-binned file ``tMin``
and ``tMax`` are otherwise ignored, so a window scan without the crop silently
returns the same fit at every level and reports a zero systematic.
"""

import itertools
import json
import math
import os

# ─────────────────────────────────────────────────────────────── grid spec ──

DEFAULT_SPEC = {
    "mode": "oat",                     # oat | factorial
    "maxRedChi2": 3.0,
    "minPValue": 0.001,
    "maxRuns": 200,
    "strictConvergence": False,        # see admit() - the flag is unreliable here
    "scanParentInGrid": False,         # the parent scan is slow and per-point useless
    "axes": {
        # "auto" derives the levels from the reference window, so the numbers are
        # always in the file's own time unit. A 2 - 5000 ms reference gives
        # 2/5/10/20 ms and 1000/2000/3500/5000 ms; a 2 - 30 s one scales with it.
        "window_start": {"enabled": True, "values": "auto"},
        "window_end":   {"enabled": True, "values": "auto"},
        "rebin":        {"enabled": True, "values": [1, 2, 4, 10]},
        "bin_mode":     {"enabled": True, "values": ["linear", "log"]},
        "background":   {"enabled": True,
                         "values": ["flat_floated", "decay_tail_slope", "decay_free_slope"]},
        "daughter_t12": {"enabled": True, "values": ["prior", "fixed", "free"]},
        "cap_eff":      {"enabled": True, "values": [False, True]},
        # Explicit pedestal half-lives, in the file's own unit, each held fixed while
        # only the amplitude floats - so the tornado reads DeltaT12 as a function of the
        # assumed slope rather than of the three flat/decaying labels above. Off by
        # default: switch it on instead of "background" once a tail fit has told you
        # which slopes are actually allowed. "auto" = 20/40/80/200 x tMax, and flat.
        "bkg_t12":      {"enabled": False, "values": "auto"},
    },
}

AXIS_TITLE = {
    "window_start": "Fit window start",
    "window_end":   "Fit window end",
    "rebin":        "Rebin factor",
    "bin_mode":     "Bin spacing",
    "background":   "Background treatment",
    "bkg_t12":      "Background half-life",
    "daughter_t12": "Daughter half-lives",
    "cap_eff":      "Efficiency cap",
}

# The quantities the budget is built for.  T1/2 is the measurement, but the same
# axes move the yield by a very different amount: on a window many half-lives
# long N0 is close to a counted quantity - the excess over the pedestal - while
# T1/2 is the *shape* of that excess, so a slope assumption that ruins one can
# leave the other alone.  Budgeting all three is the only way to see that.
RESPONSES = ("t12", "N0", "A0")
RESPONSE_TITLE = {"t12": "T1/2", "N0": "N0", "A0": "A0"}
# key each response's quadrature total lands under in the returned dict; "t12"
# keeps the historical name so existing readers of ``sigmaSyst`` still work
SIGMA_KEY = {"t12": "sigmaSyst", "N0": "sigmaSystN0", "A0": "sigmaSystA0"}


def response_unit(resp, unit):
    """The unit string for a response, given the file's own time unit."""
    if resp == "t12":
        return unit
    if resp == "N0":
        return "counts"
    return "counts/s"          # A0 is a rate per second whatever the file's unit is


LEVEL_TITLE = {
    "flat_floated":     "flat, floated",
    "decay_tail_slope": "decaying, slope from tail",
    "decay_free_slope": "decaying, free slope",
    "prior":            "literature prior",
    "fixed":            "fixed",
    "free":             "free",
    True:               "capped at 1",
    False:              "uncapped",
}


def merge_spec(user_spec):
    """Overlay a user grid spec on the defaults, axis by axis."""
    spec = json.loads(json.dumps(DEFAULT_SPEC))
    if not user_spec:
        return spec
    for k, v in user_spec.items():
        if k == "axes" and isinstance(v, dict):
            for ax, cfg in v.items():
                base = spec["axes"].setdefault(ax, {"enabled": True, "values": []})
                if isinstance(cfg, dict):
                    base.update(cfg)
                elif isinstance(cfg, list):
                    base["values"] = cfg
                elif isinstance(cfg, bool):
                    base["enabled"] = cfg
        else:
            spec[k] = v
    return spec


# ───────────────────────────────────────────────────── applying one variant ──

def auto_values(axis, params):
    """Window levels derived from the reference window, in the file's own unit."""
    t0 = float(params.get("tMin") or 0)
    t1 = float(params.get("tMax") or 0)
    if axis == "window_start":
        if t0 <= 0:
            return []
        out = [t0 * k for k in (1, 2.5, 5, 10)]
        return [round(x, 9) for x in out if not t1 or x < 0.5 * t1]
    if axis == "window_end":
        if t1 <= 0:
            return []
        out = [t1 * k for k in (0.2, 0.4, 0.7, 1.0)]
        return [round(x, 9) for x in out if x > t0]
    if axis == "bkg_t12":
        # a pedestal only looks like a pedestal if it decays far more slowly than the
        # window is long; below ~20x tMax it is a component, not a background
        if t1 <= 0:
            return ["flat"]
        return [round(t1 * k, 9) for k in (20, 40, 80, 200)] + ["flat"]
    return []


def resolve_values(axis, cfg, params):
    v = cfg.get("values")
    if v == "auto" or not v:
        return auto_values(axis, params)
    return list(v)


def _daughter_rows(params):
    return [i for i in range(1, len(params.get("chain") or []))]


def apply_level(params, axis, value):
    """Return (new_params, note_or_None, skip_reason_or_None).

    Never mutates the configuration passed in.
    """
    p = json.loads(json.dumps(params))
    note = None

    if axis == "window_start":
        p["tMin"] = float(value)
        p["cropToWindow"] = True
        note = "cropToWindow forced on - tMin is ignored on a pre-binned file otherwise"
        if float(value) >= float(p.get("tMax", 0) or 0):
            return None, None, "window start is at or past the window end"

    elif axis == "window_end":
        p["tMax"] = float(value)
        p["cropToWindow"] = True
        note = "cropToWindow forced on - tMax is ignored on a pre-binned file otherwise"
        if float(value) <= float(p.get("tMin", 0) or 0):
            return None, None, "window end is at or before the window start"

    elif axis == "rebin":
        p["rebin"] = max(1, int(value))

    elif axis == "bin_mode":
        p["binMode"] = str(value)

    elif axis == "background":
        p["fitBackground"] = True
        if value == "flat_floated":
            p["bkgMode"], p["bkgTailFit"] = "flat", False
        elif value == "decay_tail_slope":
            p["bkgMode"], p["bkgTailFit"] = "decaying", True
        elif value == "decay_free_slope":
            p["bkgMode"], p["bkgTailFit"] = "decaying", False
            p["bkgT12Fixed"] = "free"
        else:
            return None, None, f"unknown background treatment {value!r}"

    elif axis == "bkg_t12":
        p["fitBackground"] = True
        if isinstance(value, str) and value.strip().lower() in ("flat", "none", "inf"):
            p["bkgMode"], p["bkgTailFit"] = "flat", False
            p["bkgT12"] = None
        else:
            try:
                t12b = float(value)
            except (TypeError, ValueError):
                return None, None, f"{value!r} is not a half-life or 'flat'"
            if not t12b > 0:
                return None, None, "a background half-life has to be positive"
            # tailFit off and Fixed on: the slope is the level being scanned, so nothing
            # may re-estimate it from the tail or float it away from the value asked for
            p["bkgMode"], p["bkgTailFit"] = "decaying", False
            p["bkgT12"], p["bkgT12Fixed"] = t12b, "fixed"
            note = ("background half-life held at this level - only the pedestal "
                    "amplitude floats")

    elif axis == "daughter_t12":
        rows = _daughter_rows(p)
        if not rows:
            return None, None, "no daughters in the chain"
        if value == "prior":
            bad = [i for i in rows if not float(p["chain"][i].get("t12Sigma") or 0) > 0]
            if bad:
                return None, None, ("no literature uncertainty typed for row(s) "
                                    + ", ".join(str(i + 1) for i in bad))
        for i in rows:
            p["chain"][i]["fixed"] = str(value)

    elif axis == "cap_eff":
        p["capEff"] = bool(value)

    else:
        return None, None, f"unknown axis {axis!r}"

    return p, note, None


def level_label(axis, value):
    if axis in ("window_start", "window_end"):
        return f"{value:g}"
    if axis == "bkg_t12":
        if isinstance(value, str):
            return value.strip().lower()
        try:
            return f"{float(value):g}"
        except (TypeError, ValueError):
            return str(value)
    if axis == "rebin":
        return f"x{int(value)}"
    return LEVEL_TITLE.get(value, str(value))


def plan_preview(params, spec):
    """What a scan with this spec would run, worked out without running it.

    Same resolution path as run_grid - auto levels derived from the reference
    window, levels that cannot be applied (a prior with no sigma, say) marked
    as skipped - so the preview and the scan cannot disagree.
    """
    spec = merge_spec(spec)
    base = json.loads(json.dumps(params or {}))
    mode = str(spec.get("mode", "oat")).lower()

    axes = []
    for a, c in spec["axes"].items():
        if not c.get("enabled", True):
            continue
        vals = resolve_values(a, c, base)
        if not vals:
            continue
        levels, skipped = [], []
        for v in vals:
            p, _note, why = apply_level(base, a, v)
            if p is None:
                skipped.append(f"{level_label(a, v)} — {why}")
            else:
                levels.append(level_label(a, v))
        if levels or skipped:
            axes.append({"axis": a, "title": AXIS_TITLE.get(a, a),
                         "levels": levels, "skipped": skipped, "n": len(levels)})

    if mode == "factorial":
        n = 1
        for ax in axes:
            n *= max(1, ax["n"])
        cap = int(spec.get("maxRuns", 200))
        runs, capped = min(n, cap), max(0, n - cap)
    else:
        runs, capped = sum(ax["n"] for ax in axes), 0
    return {"mode": mode, "axes": axes, "nRuns": runs + 1, "nCapped": capped}


# ──────────────────────────────────────────────────── admission + recording ──

def admit(res, spec):
    """(admitted, reason). A run that is not admitted still goes in the CSV.

    The test is the physics one - does this fit still describe the data - not the
    minimiser's self-report. ``converged`` here is L-BFGS-B's ``success`` flag,
    and on this likelihood it reports failure routinely while sitting exactly at
    the minimum: the efficiency/N0 direction is near-degenerate, so the line
    search terminates abnormally with a converged answer. Admitting on that flag
    threw out runs with chi2/ndf = 1.03 and a stable half-life. It is recorded on
    every row as ``converged`` so you can see which fits it flagged, and
    ``strictConvergence`` restores it as a hard cut if you want it.
    """
    if not res or res.get("error"):
        return False, f"fit failed: {(res or {}).get('error', 'no result')}"
    if spec.get("strictConvergence", False) and not res.get("converged", False):
        return False, "minimiser reported non-convergence (strictConvergence is on)"
    fit = res.get("fit") or []
    if not fit:
        return False, "no fitted members"
    t12 = fit[0].get("t12")
    if t12 is None or not math.isfinite(float(t12)) or float(t12) <= 0:
        return False, "parent half-life is not finite and positive"
    st = res.get("stats") or {}
    rc, pv = st.get("redChi2"), st.get("pValue")
    if rc is None or not math.isfinite(float(rc)):
        return False, "chi2/ndf unavailable"
    if float(rc) > float(spec["maxRedChi2"]):
        return False, f"chi2/ndf {float(rc):.3f} > {float(spec['maxRedChi2']):g}"
    if pv is not None and math.isfinite(float(pv)) and float(pv) < float(spec["minPValue"]):
        return False, f"p = {float(pv):.3g} < {float(spec['minPValue']):g}"
    return True, ""


def record(res, run_id, axis, value, note):
    """One flat row per grid point, in the file's own time unit."""
    row = {"run": run_id, "axis": axis, "axis_label": AXIS_TITLE.get(axis, axis),
           "level": level_label(axis, value) if axis != "reference" else "reference",
           "note": note or ""}
    if not res or res.get("error"):
        row["error"] = (res or {}).get("error", "no result")
        return row

    UF = res.get("unitFactor") or 1.0
    unit = res.get("timeUnit") or "s"
    fit = res["fit"]
    p0 = fit[0]
    st = res.get("stats") or {}

    row.update({
        "time_unit": unit,
        "parent": p0.get("label"),
        "t12": (p0["t12"] / UF) if p0.get("t12") is not None else None,
        "t12_sigma": (p0["sigma"] / UF) if p0.get("sigma") is not None else None,
        "t12_stat": (p0["sigmaStat"] / UF) if p0.get("sigmaStat") is not None else None,
        "t12_syst": (p0["sigmaSyst"] / UF) if p0.get("sigmaSyst") is not None else None,
        "t12_profile_lo": (p0["loBound"] / UF) if p0.get("loBound") is not None else None,
        "t12_profile_hi": (p0["hiBound"] / UF) if p0.get("hiBound") is not None else None,
        "chi2_per_ndf": st.get("redChi2"),
        "p_value": st.get("pValue"),
        "cash": st.get("cash"),
        "ndf": st.get("ndf"),
        "n_bins": st.get("nBins"),
        "n_params": st.get("nParams"),
        "N0": res.get("N0"),
        "N0_sigma": res.get("N0Err"),
        "N0_stat": res.get("N0ErrStat"),
        "A0": res.get("A0"),
        "A0_sigma": res.get("A0Err"),
        "A0_stat": res.get("A0ErrStat"),
        "bkg_rate": res.get("bkg"),
        "bkg_rate_sigma": res.get("bkgErr"),
        "bkg_t12": (res["bkgT12"] / UF) if res.get("bkgT12") is not None else None,
        "bkg_t12_plateau": res.get("bkgT12Plateau"),
        "converged": res.get("converged"),
    })
    # the profile interval is asymmetric and that asymmetry is real
    if row["t12"] is not None:
        if row["t12_profile_lo"] is not None:
            row["t12_minus"] = row["t12"] - row["t12_profile_lo"]
        if row["t12_profile_hi"] is not None:
            row["t12_plus"] = row["t12_profile_hi"] - row["t12"]
    for i, f in enumerate(fit):
        if i:
            row[f"eff_{f.get('label')}"] = f.get("eff")
            if f.get("branchPercent") is not None:
                row[f"branch_pct_{f.get('label')}"] = f.get("branchPercent")
    return row


# ───────────────────────────────────────────────────────────── the scan ──

def _key(params, axis_levels):
    return json.dumps([axis_levels], sort_keys=True)


def run_grid(params, spec, fitfn, progress=None):
    """Execute the grid. ``fitfn(params, run_id) -> res`` does one fit.

    Returns the full result dict; writing files is the caller's job.
    """
    spec = merge_spec(spec)
    mode = str(spec.get("mode", "oat")).lower()
    base = json.loads(json.dumps(params))
    if not spec.get("scanParentInGrid", False):
        base["scanParent"] = False

    axes = []
    for a, c in spec["axes"].items():
        if not c.get("enabled", True):
            continue
        vals = resolve_values(a, c, base)
        if vals:
            axes.append((a, dict(c, values=vals)))
    spec["resolvedAxes"] = {a: c["values"] for a, c in axes}

    notes, rows, skipped = [], [], []

    def _say(msg):
        if progress:
            progress(msg)

    # ── reference ────────────────────────────────────────────────────────
    _say("reference configuration")
    ref_res = fitfn(base, "ref")
    ref_row = record(ref_res, "ref", "reference", None, "the configuration as you set it")
    ref_ok, ref_reason = admit(ref_res, spec)
    ref_row["admitted"], ref_row["admit_reason"] = ref_ok, ref_reason
    for resp in RESPONSES:
        ref_row[f"delta_{resp}"] = 0.0
    rows.append(ref_row)
    if not ref_ok:
        return {"error": f"The reference configuration itself is not admissible ({ref_reason}). "
                         "Fix the single fit before scanning systematics around it.",
                "rows": rows, "spec": spec}
    ref_t12 = ref_row["t12"]
    ref_val = {resp: ref_row.get(resp) for resp in RESPONSES}
    unit = ref_row.get("time_unit", "s")

    # ── build the run list ───────────────────────────────────────────────
    plan = []                                   # (run_id, {axis: value}, params, note)
    if mode == "factorial":
        levels = [[(a, v) for v in c["values"]] for a, c in axes]
        combos = list(itertools.product(*levels)) if levels else []
        cap = int(spec.get("maxRuns", 200))
        dropped = max(0, len(combos) - cap)
        if dropped:
            notes.append(f"Factorial grid asked for {len(combos)} combinations; "
                         f"capped at maxRuns = {cap}, so {dropped} were not run. "
                         "The budget below therefore covers only the first "
                         f"{cap} combinations in axis order.")
            combos = combos[:cap]
        for k, combo in enumerate(combos):
            p, ns, skip = base, [], None
            for a, v in combo:
                p2, note, why = apply_level(p, a, v)
                if p2 is None:
                    skip = f"{AXIS_TITLE.get(a, a)} = {level_label(a, v)}: {why}"
                    break
                p, _ = p2, ns.append(note) if note else None
            if skip:
                skipped.append(skip)
                continue
            plan.append((f"f{k:03d}", dict(combo), p,
                         "; ".join(sorted(set(x for x in ns if x)))))
    else:
        for a, c in axes:
            for j, v in enumerate(c["values"]):
                p, note, why = apply_level(base, a, v)
                if p is None:
                    skipped.append(f"{AXIS_TITLE.get(a, a)} = {level_label(a, v)}: {why}")
                    continue
                plan.append((f"{a}:{j}", {a: v}, p, note))

    # ── execute ──────────────────────────────────────────────────────────
    for run_id, combo, p, note in plan:
        # a variant that reproduces the reference exactly is not worth a refit
        if json.dumps(p, sort_keys=True) == json.dumps(base, sort_keys=True):
            r = dict(ref_row)
            r.update({"run": run_id, "note": (note or "") + " (identical to the reference)"})
            if mode != "factorial":
                a, v = next(iter(combo.items()))
                r["axis"], r["axis_label"] = a, AXIS_TITLE.get(a, a)
                r["level"] = level_label(a, v)
            for resp in RESPONSES:
                r[f"delta_{resp}"] = 0.0
            rows.append(r)
            continue
        _say(run_id)
        res = fitfn(p, run_id)
        if mode == "factorial":
            axis_key = "factorial"
            lvl = ", ".join(f"{AXIS_TITLE.get(a, a)}={level_label(a, v)}"
                            for a, v in combo.items())
            row = record(res, run_id, axis_key, None, note)
            row["level"] = lvl
            for a, v in combo.items():
                row[f"axis__{a}"] = level_label(a, v)
        else:
            a, v = next(iter(combo.items()))
            row = record(res, run_id, a, v, note)
        ok, why = admit(res, spec)
        row["admitted"], row["admit_reason"] = ok, why
        for resp in RESPONSES:
            got, base_v = row.get(resp), ref_val.get(resp)
            row[f"delta_{resp}"] = ((got - base_v)
                                    if (ok and got is not None and base_v is not None)
                                    else None)
        rows.append(row)

    # ── synthesise the budget ────────────────────────────────────────────
    budget, total_sq = [], {resp: 0.0 for resp in RESPONSES}
    if mode == "factorial":
        for a, _c in axes:
            vals = [r["t12"] for r in rows
                    if r.get("admitted") and r.get(f"axis__{a}") is not None
                    and r.get("t12") is not None]
            if not vals:
                continue
            entry = {"axis": a, "axis_label": AXIS_TITLE.get(a, a),
                     "n_admitted": len(vals), "worst_level": "(main effect across levels)"}
            for resp in RESPONSES:
                # main effect: the spread of the answer across this axis's levels
                per = {}
                for r in rows:
                    if (r.get("admitted") and r.get(f"axis__{a}") is not None
                            and r.get(resp) is not None):
                        per.setdefault(r[f"axis__{a}"], []).append(r[resp])
                if not per:
                    entry[f"delta_{resp}"] = None
                    continue
                meds = [sorted(v)[len(v) // 2] for v in per.values()]
                d = (max(meds) - min(meds)) if len(meds) > 1 else 0.0
                entry[f"delta_{resp}"] = d
                entry[f"worst_level_{resp}"] = "(main effect across levels)"
                total_sq[resp] += d * d
            budget.append(entry)
        adm = [r["t12"] for r in rows if r.get("admitted") and r.get("t12") is not None]
        notes.append("Factorial mode: per-axis numbers are main effects (the spread of the "
                     "per-level medians), not clean one-at-a-time attributions, so treat the "
                     "quadrature total as an indicator. The full envelope over all admitted "
                     f"runs is {min(adm):.4g} - {max(adm):.4g} {unit}." if adm else "")
    else:
        for a, _c in axes:
            mine = [r for r in rows if r.get("axis") == a and r.get("admitted")
                    and r.get("delta_t12") is not None]
            if not mine:
                rej = [r for r in rows if r.get("axis") == a]
                if rej:
                    notes.append(f"{AXIS_TITLE.get(a, a)}: every level was rejected, so this "
                                 "axis contributes nothing to the budget. That is a gap in the "
                                 "budget, not a zero systematic.")
                continue
            entry = {"axis": a, "axis_label": AXIS_TITLE.get(a, a), "n_admitted": len(mine)}
            for resp in RESPONSES:
                got = [r for r in rows if r.get("axis") == a and r.get("admitted")
                       and r.get(f"delta_{resp}") is not None]
                if not got:
                    entry[f"delta_{resp}"] = None
                    continue
                worst = max(got, key=lambda r: abs(r[f"delta_{resp}"]))
                d = abs(worst[f"delta_{resp}"])
                entry[f"delta_{resp}"] = d
                entry[f"worst_level_{resp}"] = worst["level"]
                entry[f"lo_{resp}"] = min(r[resp] for r in got)
                entry[f"hi_{resp}"] = max(r[resp] for r in got)
                total_sq[resp] += d * d
            # the half-life keeps the unsuffixed names: it is the measurement, and
            # every reader written before N0 and A0 joined the budget looks for these
            entry["worst_level"] = entry.get("worst_level_t12", "")
            entry["lo"], entry["hi"] = entry.get("lo_t12"), entry.get("hi_t12")
            budget.append(entry)

    budget.sort(key=lambda b: -(b.get("delta_t12") or 0.0))
    sigma = {resp: math.sqrt(total_sq[resp]) for resp in RESPONSES}
    sigma_syst = sigma["t12"]

    n_flag = sum(1 for r in rows if r.get("admitted") and r.get("converged") is False)
    if n_flag and not spec.get("strictConvergence", False):
        notes.append(f"{n_flag} admitted runs carry the minimiser's non-convergence flag. On "
                     "this likelihood that flag fires on the near-degenerate efficiency/N0 "
                     "direction even at a good minimum, so it is recorded rather than used as "
                     "a cut; check their chi2/ndf in the table before trusting them, or set "
                     "strictConvergence to exclude them.")

    n_rej = sum(1 for r in rows if r.get("admitted") is False)
    if n_rej:
        notes.append(f"{n_rej} of {len(rows)} runs were not admitted and are excluded from the "
                     "budget; each carries its reason in the results table.")
    if skipped:
        notes.append("Levels that could not be built: " + "; ".join(skipped))
    notes.append("The window axes force 'crop to window' on. Without it tMin and tMax do "
                 "nothing to an already-binned file and the axis would report a spurious zero.")

    rel = [f"{RESPONSE_TITLE[resp]} {100.0 * sigma[resp] / abs(ref_val[resp]):.2f} %"
           for resp in RESPONSES if ref_val.get(resp)]
    if rel:
        notes.append("Relative systematic by quantity: " + ", ".join(rel) + ". The same axes "
                     "do not cost the three the same. Over a window many half-lives long N0 "
                     "is close to the counted excess over the pedestal, while T1/2 is the "
                     "shape of that excess, so an assumption that ruins the half-life can "
                     "leave the yield almost untouched - and A0 = N0 ln2 / T1/2 then carries "
                     "whatever the half-life does. Quote each from the budget that belongs "
                     "to it rather than assuming one number covers all three.")

    return {
        "mode": mode, "spec": spec, "rows": rows, "budget": budget,
        "reference": {"t12": ref_t12,
                      "sigma": ref_row.get("t12_sigma"),
                      "stat": ref_row.get("t12_stat"),
                      "syst": ref_row.get("t12_syst"),
                      "chi2_per_ndf": ref_row.get("chi2_per_ndf"),
                      "p_value": ref_row.get("p_value"),
                      "parent": ref_row.get("parent"),
                      "N0": ref_row.get("N0"),
                      "N0_sigma": ref_row.get("N0_sigma"),
                      "N0_stat": ref_row.get("N0_stat"),
                      "A0": ref_row.get("A0"),
                      "A0_sigma": ref_row.get("A0_sigma"),
                      "A0_stat": ref_row.get("A0_stat")},
        "sigmaSyst": sigma_syst,
        "sigmaSystN0": sigma["N0"], "sigmaSystA0": sigma["A0"],
        "timeUnit": unit,
        "nRuns": len(rows), "nAdmitted": sum(1 for r in rows if r.get("admitted")),
        "notes": [n for n in notes if n],
    }


# ─────────────────────────────────────────────────────────────── outputs ──

CSV_LEAD = ["run", "axis", "axis_label", "level", "admitted", "admit_reason",
            "t12", "t12_sigma", "t12_stat", "t12_syst",
            "t12_profile_lo", "t12_profile_hi", "t12_minus", "t12_plus", "delta_t12",
            "chi2_per_ndf", "p_value", "cash", "ndf", "n_bins", "n_params",
            "N0", "N0_sigma", "N0_stat", "delta_N0",
            "A0", "A0_sigma", "A0_stat", "delta_A0",
            "bkg_rate", "bkg_rate_sigma", "bkg_t12", "bkg_t12_plateau",
            "converged", "time_unit", "parent", "note", "error"]


def write_outputs(grid, outdir, figure=True):
    """Write systematic_grid_results.csv, systematic_budget.csv and the tornado."""
    import pandas as pd
    os.makedirs(outdir, exist_ok=True)
    written = []

    rows = grid.get("rows") or []
    df = pd.DataFrame(rows)
    cols = [c for c in CSV_LEAD if c in df.columns] + \
           [c for c in df.columns if c not in CSV_LEAD]
    path = os.path.join(outdir, "systematic_grid_results.csv")
    df[cols].to_csv(path, index=False)
    written.append(path)

    unit = grid.get("timeUnit", "s")
    bud = list(grid.get("budget") or [])
    bdf = pd.DataFrame(bud)
    if len(bdf):
        bdf = bdf.rename(columns={"delta_t12": f"delta_T12_{unit}",
                                  "delta_N0": "delta_N0_counts",
                                  "delta_A0": "delta_A0_per_s"})
    total = {"axis": "TOTAL", "axis_label": "quadrature total",
             f"delta_T12_{unit}": grid.get("sigmaSyst"),
             "delta_N0_counts": grid.get("sigmaSystN0"),
             "delta_A0_per_s": grid.get("sigmaSystA0"),
             "n_admitted": grid.get("nAdmitted"), "worst_level": ""}
    bdf = pd.concat([bdf, pd.DataFrame([total])], ignore_index=True)
    bpath = os.path.join(outdir, "systematic_budget.csv")
    bdf.to_csv(bpath, index=False)
    written.append(bpath)

    with open(os.path.join(outdir, "grid_results.json"), "w") as f:
        json.dump(_clean(grid), f, allow_nan=False)
    written.append(os.path.join(outdir, "grid_results.json"))

    if figure:
        for resp in RESPONSES:
            stem = "systematic_tornado" if resp == "t12" else f"systematic_tornado_{resp}"
            written += tornado_figure(grid, outdir, stem=stem, response=resp)
    return written


def _clean(o):
    if isinstance(o, dict):
        return {k: _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if isinstance(o, bool):
        return o
    if isinstance(o, float):
        return o if math.isfinite(o) else None
    try:
        import numpy as np
        if isinstance(o, np.floating):
            return float(o) if math.isfinite(float(o)) else None
        if isinstance(o, np.integer):
            return int(o)
        if isinstance(o, np.bool_):
            return bool(o)
    except Exception:
        pass
    return o


def tornado_figure(grid, outdir, stem="systematic_tornado", response="t12"):
    """Publication-grade horizontal tornado, written as both PDF and PNG.

    ``response`` picks the quantity: "t12", "N0" or "A0". One figure each, because
    the three have different units and, more to the point, different shapes.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    tunit = grid.get("timeUnit", "s")
    unit = response_unit(response, tunit)
    refs = grid.get("reference") or {}
    ref = refs.get("t12") if response == "t12" else refs.get(response)
    dkey = f"delta_{response}"
    bud = [b for b in (grid.get("budget") or []) if b.get(dkey) is not None]
    if ref is None or not bud:
        return []
    rows = grid.get("rows") or []

    # for each axis: the admitted envelope around the reference
    bars = []
    for b in bud:
        a = b["axis"]
        if a == "factorial":
            continue
        vals = [r[response] for r in rows
                if r.get("axis") == a and r.get("admitted") and r.get(response) is not None]
        if not vals:
            vals = [ref - b[dkey], ref + b[dkey]]
        bars.append((b["axis_label"], min(vals) - ref, max(vals) - ref, b[dkey]))
    if not bars:
        return []
    bars.sort(key=lambda x: x[3])

    h = max(2.2, 0.46 * len(bars) + 1.5)
    fig, ax = plt.subplots(figsize=(7.0, h))
    y = range(len(bars))
    for i, (_lab, lo, hi, _d) in enumerate(bars):
        lo, hi = min(lo, 0.0), max(hi, 0.0)
        ax.barh(i, hi - lo, left=lo, height=0.6,
                color="#8fb4dd", edgecolor="#2d5f92", linewidth=0.8, zorder=3)
    ax.axvline(0.0, color="#111111", linewidth=1.1, zorder=4)

    stat = refs.get("sigma") if response == "t12" else refs.get(f"{response}_sigma")
    if stat:
        ax.axvspan(-stat, stat, color="#d9534f", alpha=0.10, zorder=1)
        ax.axvline(-stat, color="#d9534f", lw=0.8, ls=":", zorder=2)
        ax.axvline(stat, color="#d9534f", lw=0.8, ls=":", zorder=2,
                   label=f"statistical  ±{stat:.3g} {unit}")

    syst = grid.get(SIGMA_KEY[response])
    if syst:
        ax.axvline(-syst, color="#2d5f92", lw=1.0, ls="--", zorder=2)
        ax.axvline(syst, color="#2d5f92", lw=1.0, ls="--", zorder=2,
                   label=f"quadrature systematic  ±{syst:.3g} {unit}")

    ax.set_yticks(list(y))
    ax.set_yticklabels([b[0] for b in bars], fontsize=9)
    rlab = {"t12": "T½", "N0": "N₀", "A0": "A₀"}[response]
    ax.set_xlabel(f"Deviation of {rlab} from the reference fit  ({unit})", fontsize=10)
    parent = (grid.get("reference") or {}).get("parent") or "parent"
    ax.set_title(f"{parent}:  {rlab} = {ref:.4g} {unit}   (reference configuration)",
                 fontsize=10.5, loc="left")
    ax.grid(axis="x", color="#e4e4e4", lw=0.7, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    if stat or syst:
        ax.legend(fontsize=8, frameon=False, loc="lower right")
    fig.tight_layout()

    out = []
    for ext in ("pdf", "png"):
        p = os.path.join(outdir, f"{stem}.{ext}")
        fig.savefig(p, dpi=200 if ext == "png" else None, bbox_inches="tight")
        out.append(p)
    plt.close(fig)
    return out


def summary_table(grid):
    """The budget as the plain-text table a paper's systematics section wants."""
    unit = grid.get("timeUnit", "s")
    ref = grid.get("reference") or {}
    lines = []
    t12 = ref.get("t12")
    lines.append(f"  reference  T1/2 = {t12:.6g} {unit}"
                 + (f"  +/-{ref['sigma']:.4g} (stat)" if ref.get("sigma") else "")
                 + (f"   chi2/ndf = {ref['chi2_per_ndf']:.4f}" if ref.get("chi2_per_ndf") else ""))
    n0, a0 = ref.get("N0"), ref.get("A0")
    if n0:
        lines.append(f"             N0   = {n0:.6g} counts"
                     + (f"  +/-{ref['N0_sigma']:.4g} (fit)" if ref.get("N0_sigma") else ""))
    if a0:
        lines.append(f"             A0   = {a0:.6g} counts/s"
                     + (f"  +/-{ref['A0_sigma']:.4g} (fit)" if ref.get("A0_sigma") else ""))
    lines.append("")
    # T1/2 in its own unit; the yield and rate as percentages, which is how they read
    lines.append(f"  {'source':<26} {'delta T1/2 ('+unit+')':>18} {'dN0 (%)':>9} "
                 f"{'dA0 (%)':>9}  {'runs':>5}  worst level (T1/2)")
    lines.append("  " + "-" * 96)

    def _pct(d, v):
        return f"{100.0 * d / abs(v):>9.3g}" if (d is not None and v) else f"{'-':>9}"

    for b in grid.get("budget") or []:
        d12 = b.get("delta_t12")
        lines.append(f"  {b['axis_label']:<26} "
                     + (f"{d12:>18.4g}" if d12 is not None else f"{'-':>18}")
                     + f" {_pct(b.get('delta_N0'), n0)} {_pct(b.get('delta_A0'), a0)}  "
                     + f"{b.get('n_admitted', 0):>5}  {b.get('worst_level', '')}")
    lines.append("  " + "-" * 96)
    lines.append(f"  {'quadrature total':<26} {grid.get('sigmaSyst', 0.0):>18.4g}"
                 f" {_pct(grid.get('sigmaSystN0'), n0)} {_pct(grid.get('sigmaSystA0'), a0)}")
    lines.append("")
    if t12 is not None:
        lines.append(f"  T1/2 = {t12:.6g} +/-{ref.get('sigma') or 0:.4g} (stat) "
                     f"+/-{grid.get('sigmaSyst', 0.0):.4g} (syst) {unit}")
    if n0 and grid.get("sigmaSystN0") is not None:
        lines.append(f"  N0   = {n0:.6g} +/-{ref.get('N0_stat') or ref.get('N0_sigma') or 0:.4g}"
                     f" (stat) +/-{grid['sigmaSystN0']:.4g} (syst) counts")
    if a0 and grid.get("sigmaSystA0") is not None:
        lines.append(f"  A0   = {a0:.6g} +/-{ref.get('A0_stat') or ref.get('A0_sigma') or 0:.4g}"
                     f" (stat) +/-{grid['sigmaSystA0']:.4g} (syst) counts/s")
    return "\n".join(lines)
