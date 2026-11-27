#!/usr/bin/env python3
"""Command-line Bateman fitter - no GUI, no Streamlit, nothing but numpy/scipy/pandas.

This is the version to use on a locked-down analysis node.

    python run_fit.py config_44S.json test_44S.csv -o out_44S

It writes the same files the interface writes:

    <parent>_fit_results.csv        half-lives, errors, branchings, efficiencies
    <parent>_derived_quantities.csv N0, A0, background, rho(lnT12, lnN0)
    <parent>_covariance.csv         the four stacked matrices
    <parent>_spectrum_fit.csv       every curve, for drawing the paper figure
    <parent>_parent_scan.csv        only if scanParent is true
    results.json                    everything above, machine-readable

Add --plot to also write a PNG of the fit and residuals (needs matplotlib).

With --grid it instead runs a systematic parameter scan around that configuration:

    python run_fit.py config_44S.json test_44S.csv -o out_44S_grid --grid

Every analysis choice you could defensibly have made differently is re-fitted, the
runs that no longer describe the data are thrown out, and what survives is reduced
to one deviation per axis and a quadrature total - for the half-life, for the
correlated-implant yield N0 and for the initial activity A0, since the same axes
do not cost the three the same.  Written as

    systematic_grid_results.csv     one row per grid point, every recorded quantity
    systematic_budget.csv           delta T1/2, delta N0 and delta A0 per axis,
                                    plus the quadrature total for each
    systematic_tornado.pdf / .png   the publication figure (T1/2)
    systematic_tornado_N0.*         the same figure for the yield
    systematic_tornado_A0.*         the same figure for the initial activity
    grid_results.json               everything above, machine-readable

--grid-spec points at a JSON file overriding the axes; --grid-mode factorial runs
the full Cartesian product instead of one-axis-at-a-time.
"""

import argparse
import json
import math
import os
import sys

from bateman_fit import run_fit


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("config", help="JSON settings file, e.g. config_44S.json")
    ap.add_argument("spectrum", nargs="?", default=None,
                    help="data file (dt,counts). Omit to use the built-in demo spectrum.")
    ap.add_argument("-o", "--outdir", default="out", help="output directory (default: out)")
    ap.add_argument("--plot", action="store_true", help="also write a PNG of the fit")
    ap.add_argument("--subtract", action="store_true",
                    help="with --plot: show the background-subtracted spectrum, which is "
                         "the only readable view when the pedestal dominates")
    ap.add_argument("--set", action="append", default=[], metavar="KEY=VALUE",
                    help="override one config key, e.g. --set tieSiblings=false. Repeatable.")
    ap.add_argument("--grid", action="store_true",
                    help="run a systematic parameter grid scan around this configuration "
                         "and write the tornado figure and budget instead of a single fit")
    ap.add_argument("--grid-spec", default=None, metavar="FILE",
                    help="with --grid: JSON file overriding the scan axes, admission "
                         "thresholds or mode. Omitted keys keep their defaults.")
    ap.add_argument("--grid-mode", default=None, choices=["oat", "factorial"],
                    help="with --grid: one-axis-at-a-time (default) or the full product")
    ap.add_argument("--grid-max-runs", type=int, default=None, metavar="N",
                    help="with --grid --grid-mode factorial: hard cap on the run count")
    args = ap.parse_args()

    with open(args.config) as f:
        params = json.load(f)

    for kv in args.set:
        if "=" not in kv:
            sys.exit(f"--set needs KEY=VALUE, got {kv!r}")
        k, v = kv.split("=", 1)
        try:
            params[k] = json.loads(v)          # numbers, true/false, null
        except json.JSONDecodeError:
            params[k] = v                      # plain string
        print(f"override: {k} = {params[k]!r}")

    if args.grid:
        return run_grid_cli(params, args)

    res = run_fit(params, outdir=args.outdir, spectrum=args.spectrum)

    if "error" in res:
        print("\nFIT FAILED:", res["error"], file=sys.stderr)
        print(res.get("traceback", ""), file=sys.stderr)
        return 1

    unit = res["timeUnit"]
    UF = res["unitFactor"]

    print()
    print(f"  {'nuclide':<16} {'T1/2':>12}  {'error':<28} status")
    print("  " + "-" * 74)
    for r in res["fit"]:
        val = r["t12"] / UF
        if r["sigma"] is None:
            err = "(held)"
        elif r["sigmaStat"] is not None and r["sigmaSyst"] is not None:
            err = (f"+/-{r['sigmaStat']/UF:.4g} stat "
                   f"+/-{r['sigmaSyst']/UF:.4g} syst")
        else:
            err = f"+/-{r['sigma']/UF:.4g}"
        status = "fixed" if r["fixed"] else ("constrained" if r["constrained"] else "free")
        if r["loBound"] and r["hiBound"]:
            status += f"  [profile {r['loBound']/UF:.4g}, {r['hiBound']/UF:.4g}]"
        print(f"  {r['label']:<16} {val:>12.6g}  {err:<28} {status}")

    st = res["stats"]
    print()
    print(f"  N0 = {res['N0']:.6g} correlated {res['parentLabel']} implants"
          + (f"  +/-{res['N0Err']:.4g}" if res.get("N0Err") else ""))
    print(f"  A0 = {res['A0']:.6g} decays/s at t = 0"
          + (f"  +/-{res['A0Err']:.4g}" if res.get("A0Err") else ""))
    print(f"  background {res['bkg']:.6g} counts/s"
          + (f", T1/2(bkg) = {res['bkgT12']/UF:.6g} {unit}" if res.get("bkgT12") else ""))
    print(f"  Cash = {st['cash']:.1f}   chi2/ndf = {st['redChi2']:.4f} "
          f"({st['chi2']:.1f}/{st['ndf']})   p = {st['pValue']:.4g}")
    print(f"  converged: {res['converged']}")

    print("\n  Notes")
    for note in res["notes"]:
        print("   -", note)

    if args.plot:
        make_plot(res, args.outdir, subtract=args.subtract)

    print(f"\n  Files written to {os.path.abspath(args.outdir)}:")
    for fn in sorted(os.listdir(args.outdir)):
        print("   -", fn)

    return 0


def run_grid_cli(params, args):
    """--grid: scan the analysis choices, synthesise the budget, draw the tornado."""
    import shutil
    import tempfile
    import time

    import grid_scan

    spec = {}
    if args.grid_spec:
        with open(args.grid_spec) as f:
            spec = json.load(f)
    if args.grid_mode:
        spec["mode"] = args.grid_mode
    if args.grid_max_runs:
        spec["maxRuns"] = args.grid_max_runs

    os.makedirs(args.outdir, exist_ok=True)
    scratch = tempfile.mkdtemp(prefix="bateman_grid_")
    t0 = time.time()
    n = [0]

    def fitfn(p, run_id):
        n[0] += 1
        d = os.path.join(scratch, str(run_id).replace(":", "_"))
        try:
            return run_fit(p, outdir=d, spectrum=args.spectrum)
        except Exception as e:                       # a dead point must not kill the scan
            return {"error": f"{type(e).__name__}: {e}"}
        finally:
            shutil.rmtree(d, ignore_errors=True)

    def progress(msg):
        print(f"  [{n[0]:3d}] {msg}", flush=True)

    try:
        grid = grid_scan.run_grid(params, spec, fitfn, progress=progress)
    finally:
        shutil.rmtree(scratch, ignore_errors=True)

    if grid.get("error"):
        print("\nGRID FAILED:", grid["error"], file=sys.stderr)
        return 1

    written = grid_scan.write_outputs(grid, args.outdir, figure=True)

    print()
    print(f"  {grid['nAdmitted']} of {grid['nRuns']} runs admitted "
          f"in {time.time() - t0:.1f} s  (mode: {grid['mode']})")
    print()
    print(grid_scan.summary_table(grid))
    if grid.get("notes"):
        print("  Notes")
        for note in grid["notes"]:
            print("   -", note)
    print(f"\n  Files written to {os.path.abspath(args.outdir)}:")
    for p in written:
        print("   -", os.path.basename(p))
    return 0


def make_plot(res, outdir, subtract=False):
    """Fit and residuals, the same two panels the interface shows."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("\n  --plot needs matplotlib: pip install matplotlib", file=sys.stderr)
        return

    sp, UF, unit = res["spectrum"], res["unitFactor"], res["timeUnit"]
    import numpy as np
    x = np.array(sp["center"]) / UF
    counts = np.array(sp["counts"])
    mu = np.array(sp["mu"])
    bkg = np.array(sp["bkg"])

    fig, (ax, axr) = plt.subplots(
        2, 1, figsize=(8, 7), sharex=True,
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.06})

    y_d = counts - bkg if subtract else counts
    ax.errorbar(x, y_d, yerr=np.sqrt(np.maximum(counts, 0)), fmt=".", ms=2,
                color="0.35", elinewidth=0.5, label="data", zorder=1)
    ax.plot(x, mu - bkg if subtract else mu, "-", lw=1.8, color="#c0392b",
            label="chain model" if subtract else "total model", zorder=5)
    if not subtract:
        ax.plot(x, bkg, "--", lw=1.2, color="#7f8c8d", label="background", zorder=3)
    for lab, comp in zip([r["label"] for r in res["fit"]], sp["comps"]):
        ax.plot(x, comp, "-", lw=1.0, alpha=0.9, label=lab, zorder=4)
    ax.set_yscale("log")
    if subtract:
        pos = y_d[y_d > 0]
        if pos.size:
            ax.set_ylim(max(pos.min() * 0.5, 1e-3), y_d.max() * 3)
    ax.set_ylabel("counts / channel" + (" (background subtracted)" if subtract else ""))
    ax.legend(fontsize=8, ncol=2)
    ax.set_title(f"{res['parentLabel']}  T1/2 = {res['fit'][0]['t12']/UF:.4g} {unit}"
                 f"   chi2/ndf = {res['stats']['redChi2']:.3f}")

    pear = np.array(res["residuals"]["pearson"])
    axr.axhline(0, color="0.6", lw=0.8)
    axr.plot(x, pear, ".", ms=2, color="#2c3e50")
    axr.set_ylabel("(data-fit)/sigma")
    axr.set_xlabel(f"implant-decay time difference ({unit})")

    path = os.path.join(outdir, "fit.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(main())
