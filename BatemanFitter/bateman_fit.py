"""Bateman half-life fitter - the fit engine, lifted verbatim from the
Academia Co-scientist mini-app so it runs on your own machine.

Nothing in this file talks to a network. It takes a parameter dictionary and a
spectrum file, and writes the same CSVs the web tool writes.

    from bateman_fit import run_fit
    res = run_fit(params, outdir="out", spectrum="test_44S.csv")

`params` keys are documented in README.md and exemplified in config_44S.json.
"""

import os, re, json, math, traceback

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import chi2 as chi2dist

LN2 = math.log(2.0)


def run_fit(params, outdir, spectrum=None, verbose=False):
    """Fit a decay chain to an implant-decay time-difference spectrum.

    params    dict of settings (see README / config_44S.json)
    outdir    directory the CSVs and results.json are written to
    spectrum  path to the data file, or a file-like object. None (or
              params["useDemoData"] = True) runs on synthetic data instead.
    returns   the full results dictionary, identical to results.json
    """
    params = dict(params)
    if spectrum is not None and "useDemoData" not in params:
        params["useDemoData"] = False   # a file was handed in, so use it
    os.makedirs(outdir, exist_ok=True)
    res = {}

    try:
        chain = params.get("chain") or []
        if not chain:
            raise ValueError("Decay chain is empty - add at least the parent nuclide.")

        n = len(chain)
        labels = [str(c.get("label") or f"member {i+1}") for i, c in enumerate(chain)]
        t12_guess = np.array([max(float(c.get("t12") or 1.0), 1e-9) for c in chain], float)
        eff = np.array([float(c.get("eff") if c.get("eff") is not None else 1.0) for c in chain], float)
        # three modes per half-life: free, fixed, or "prior" - free but pulled towards the
        # literature value by a Gaussian penalty of the width the user quotes
        modes = [str(c.get("fixed", "free")).lower() for c in chain]
        t12_sig_raw = np.array([float(c.get("t12Sigma") or 0.0) for c in chain], float)
        free = [i for i, m in enumerate(modes) if m != "fixed"]
        eff_free = [i for i, c in enumerate(chain)
                    if str(c.get("effFixed", "fixed")).lower() == "free"]

        # every time the user types is in this unit; internally everything is seconds
        UF = {"s": 1.0, "ms": 1e-3, "us": 1e-6, "ns": 1e-9}.get(str(params.get("timeUnit", "s")), 1.0)
        unit_name = str(params.get("timeUnit", "s"))
        t12_guess = t12_guess * UF
        t12_sig = t12_sig_raw * UF

        use_demo = bool(params.get("useDemoData", True))
        t_min = float(params.get("tMin", 0.002)) * UF
        t_max = float(params.get("tMax", 30.0)) * UF
        n_bins = int(params.get("nBins", 64))
        bin_mode = str(params.get("binMode", "log"))
        edge_ref = str(params.get("edgeRef", "center"))
        rebin = max(int(params.get("rebin", 1)), 1)
        fit_bkg = bool(params.get("fitBackground", True))
        bkg0 = max(float(params.get("bkgRate", 8.0)), 1e-8)
        auto_start = bool(params.get("autoStart", True))
        bkg_mode = str(params.get("bkgMode", "flat"))        # flat | decaying
        # the background is treated as one more chain component: its own half-life,
        # with its own start value and free/fixed switch
        _bt = params.get("bkgT12", None)
        bkg_t12_user = (float(_bt) * UF) if (_bt is not None and float(_bt) > 0) else None
        _bkg_t12_mode = str(params.get("bkgT12Fixed", "free")).lower()
        bkg_t12_free = _bkg_t12_mode != "fixed"
        # "constrained" is the third option, and the one to use on a weak gate whose own
        # tail cannot measure the pedestal slope: the slope still floats, but a Gaussian
        # prior holds it to a value measured elsewhere - typically the tail fit of a
        # high-statistics spectrum - so the prior's width propagates into every error
        # instead of being either ignored (fixed) or unconstrained (free). Same machinery
        # as a constrained daughter half-life; it needs both a value and a width, and
        # falls back to plain "free" if either is missing.
        _bt_sig = params.get("bkgT12Sigma", None)
        bkg_t12_sigma = (float(_bt_sig) * UF) if (_bt_sig is not None and float(_bt_sig) > 0) else None
        bkg_t12_constrained = (_bkg_t12_mode == "constrained"
                               and bkg_t12_user is not None and bkg_t12_sigma is not None)

        notes = []

        # every output file carries the parent's name, so a folder holding a 44S run next to
        # a 43S run stays legible and nothing overwrites anything. results.json is the one
        # exception - the interface reads it by a fixed name.
        _tag = ""
        if bool(params.get("tagOutputs", True)):
            _tag = re.sub(r"[^0-9A-Za-z._-]+", "_", labels[0]).strip("_")

        def OUT(name):
            return os.path.join(outdir, f"{_tag}_{name}" if _tag else name)

        # the chain need not be a single line. A beta-delayed neutron branch feeds a
        # *different* nuclide directly from the same parent, so it is a parallel branch,
        # not the next generation. Each member therefore names which earlier member feeds
        # it and with what branching fraction; the default is the row above at 100%.
        par, brf, br_sig = [], [], []
        for i, c in enumerate(chain):
            try:
                p = int(c.get("parent") or 0)
            except Exception:
                p = 0
            p = (p - 1) if p > 0 else (i - 1)          # typed value is a 1-based row number
            if i == 0:
                p = -1
            elif not (0 <= p < i):
                raise ValueError(
                    f"{labels[i]} (row {i+1}) is set to be fed by row {p+1}. A member can only be fed "
                    f"by an earlier row - feeding runs one way down the chain.")
            par.append(p)
            try:
                f_ = float(c.get("br", 100.0))
            except Exception:
                f_ = 100.0
            brf.append(1.0 if i == 0 else min(max(f_ / 100.0, 0.0), 1.0))
            # Pn values are usually the worst-known number in the chain (25-50% relative is
            # typical), so a branching may carry its own uncertainty, in percentage points
            try:
                s_ = float(c.get("brSigma") or 0.0)
            except Exception:
                s_ = 0.0
            br_sig.append(0.0 if i == 0 else max(s_ / 100.0, 0.0))

        # each member's own path back to the parent
        paths = []
        for i in range(n):
            p_, q = [i], i
            while par[q] >= 0:
                q = par[q]
                p_.append(q)
            p_.reverse()
            paths.append(p_)

        def path_weights(bf):
            """Product of the branching fractions along each member's path back to row 1.

            Recomputed on every likelihood evaluation, because a branching with a quoted
            uncertainty is a fitted parameter - and one uncertain Pn propagates into every
            descendant that is fed through it.
            """
            w = []
            for i in range(n):
                m_ = 1.0
                for q2 in paths[i][1:]:
                    m_ *= bf[q2]
                w.append(m_)
            return w

        # Branches out of one nuclide are not independent numbers: a beta-delayed neutron
        # emitter has one measured quantity, Pn, and the beta- branch is whatever is left.
        # Fitted separately they become free scale factors that absorb whatever the shape
        # cannot explain, and their sum wanders away from 100%. So one child of each
        # multi-branch parent - the dominant one - is made the complement of the others.
        tie_sib = bool(params.get("tieSiblings", True))
        sib_of = {}
        for i in range(1, n):
            sib_of.setdefault(par[i], []).append(i)
        comp = {}                      # complement child -> its siblings
        if tie_sib:
            for p_, kids in sib_of.items():
                if len(kids) < 2:
                    continue
                c_ = max(kids, key=lambda k: brf[k])
                comp[c_] = [k for k in kids if k != c_]
                new = min(max(1.0 - sum(brf[k] for k in comp[c_]), 1e-9), 1.0)
                if abs(new - brf[c_]) > 1e-9:
                    notes.append(
                        f"{labels[c_]} is the dominant branch out of {labels[p_]}, so it is tied to "
                        f"100% minus its siblings: {100.0*brf[c_]:.4g}% -> {100.0*new:.4g}%. Only the "
                        f"minor branch(es) are fitted; the sum is 100% by construction.")
                else:
                    notes.append(
                        f"{labels[c_]} is tied to 100% minus the other branch(es) out of {labels[p_]}, "
                        f"so the branchings out of {labels[p_]} always sum to 100% and only the minor "
                        f"branch is a fitted parameter.")
                brf[c_] = new

        brmul = path_weights(brf)

        branched = any(par[i] != i - 1 for i in range(1, n)) or \
                   any(abs(brf[i] - 1.0) > 1e-9 for i in range(1, n))
        if branched:
            notes.append("Chain topology as fitted: " + "; ".join(
                f"{labels[i]} fed by {labels[par[i]]} at {100.0*brf[i]:.4g}%" for i in range(1, n))
                + ". Two members fed by the same ancestor are parallel branches, each computed from "
                  "its own path back to the parent - not successive generations.")
            for i in range(1, n):
                sib = [k for k in range(1, n) if par[k] == par[i]]
                if len(sib) > 1 and i == sib[0]:
                    tot = sum(brf[k] for k in sib)
                    if abs(tot - 1.0) > 1e-6:
                        notes.append(
                            f"The branches out of {labels[par[i]]} sum to {100.0*tot:.4g}%, not 100%. "
                            f"That is legitimate if you are only following some of them, but the "
                            f"missing fraction is then absorbed by N0 and the efficiencies.")

        if 0 in eff_free:
            # eff[0] multiplies every term exactly as N0 does - freeing both is singular
            eff_free = [i for i in eff_free if i != 0]
            notes.append("The parent's efficiency is exactly degenerate with the number of correlated "
                         "implants N0, so it is held at the value you typed and N0 carries the overall scale.")
        if not free and not eff_free:
            raise ValueError("Every half-life and efficiency is fixed - free at least one parameter.")

        # a branching with a quoted uncertainty becomes a fitted parameter carrying a Gaussian
        # prior, exactly like a constrained half-life. Nothing in the data pins it down on its
        # own - brf[i] multiplies member i's whole term just as eff[i] does - so the prior is
        # what makes it identifiable, and its width is what propagates into the parent.
        br_free = [i for i in range(1, n)
                   if br_sig[i] > 0 and brf[i] > 0 and i not in comp]
        for i in comp:
            if br_sig[i] > 0:
                notes.append(f"The uncertainty you typed on {labels[i]}'s branching is not used: it is "
                             f"tied to 100% minus its siblings, so its error follows theirs.")
        _clash = [i for i in br_free if i in eff_free]
        if _clash:
            eff_free = [i for i in eff_free if i not in _clash]
            notes.append(
                "Branching and efficiency multiply the same term, so freeing both for "
                + ", ".join(labels[i] for i in _clash) +
                " is exactly degenerate. The efficiency is held at the value you typed and the "
                "branching carries that scale, within its quoted error.")
        for i in range(1, n):
            if br_sig[i] > 0 and brf[i] <= 0:
                notes.append(f"{labels[i]} has a branching uncertainty but a branching of 0%, "
                             f"so there is nothing to vary - it stays at zero.")

        # echo back exactly how the typed start values were interpreted, so a unit
        # mix-up is visible rather than buried in the fit
        # members with a quoted uncertainty are constrained, not hard-fixed: free
        # parameters carrying a Gaussian penalty, so their uncertainty propagates into
        # every other parameter instead of being silently ignored
        prior = {i: float(t12_sig[i]) for i, m in enumerate(modes)
                 if m == "prior" and t12_sig[i] > 0}
        for i, m in enumerate(modes):
            if m == "prior" and i not in prior:
                notes.append(f"{labels[i]} is set to constrained but its uncertainty is zero, "
                             f"so it is fitted freely instead.")

        def _mode_word(i):
            if i not in free:
                return "held"
            if i in prior:
                return f"constrained to +/-{t12_sig[i]/UF:.4g} {unit_name}"
            return "free"

        notes.append("Start half-lives as read: " + "; ".join(
            f"{labels[i]} {t12_guess[i]/UF:.6g} {unit_name} = {t12_guess[i]:.6g} s "
            f"({_mode_word(i)})" for i in range(n)) + ".")

        # ------------------------------------------------------------------ data
        upload = spectrum
        lo = hi = counts = None
        input_kind = "events"

        if use_demo or not upload:
            if not use_demo:
                raise ValueError("No spectrum file selected. Upload one or switch the picker to demo data.")
            rng = np.random.default_rng(20240917)
            true_t12 = np.array([0.043, 0.226, 0.945])
            true_eff = np.array([1.00, 0.85, 0.70])
            n_par, b_true, win = 1500, 8.0, 30.0
            tau = true_t12 / LN2
            t_cum = np.zeros(n_par)
            parts = []
            for i in range(3):
                t_cum = t_cum + rng.exponential(tau[i], n_par)
                parts.append(t_cum[rng.random(n_par) < true_eff[i]])
            parts.append(rng.uniform(0.0, win, rng.poisson(b_true * win)))
            dt = np.concatenate(parts)
            dt = dt[dt > 0]
            res["dataSource"] = "demo"
            notes.append("Synthetic Monte Carlo: 1500 implanted parents, T-half = 43 / 226 / 945 ms, "
                         "detection efficiency 1.00 / 0.85 / 0.70, flat random-correlation background 8 s^-1.")
        else:
            def _read(**kw):
                """Read the file, trying the sniffing parser first.

                A path can be opened again; an uploaded file handle cannot - it is
                left at EOF, and every read after the first comes back empty. Since
                a headerless file is read twice (below), rewind before each attempt.
                """
                for attempt in ({"sep": None, "engine": "python"}, {}):
                    if hasattr(upload, "seek"):
                        try:
                            upload.seek(0)
                        except Exception:
                            pass
                    try:
                        return pd.read_csv(upload, **dict(attempt, **kw))
                    except Exception:
                        if attempt:
                            continue
                        raise

            df = _read()
            df.columns = [str(c).strip().lower() for c in df.columns]

            def _isnum(s):
                try:
                    float(str(s).strip()); return True
                except Exception:
                    return False

            if len(df.columns) and all(_isnum(c) for c in df.columns):
                # no header row - pandas ate the first data line as column names
                df = _read(header=None)
                nc = df.shape[1]
                if nc == 1:
                    df.columns = ["dt"]
                elif nc == 2:
                    df.columns = ["dt", "counts"]
                elif nc == 3:
                    df.columns = ["t_low", "t_high", "counts"]
                else:
                    df.columns = [f"c{i}" for i in range(nc)]
                notes.append(f"No header row in the file - read {nc} unnamed column(s) as "
                             f"{', '.join(map(str, df.columns))} ({len(df)} rows).")

            res["dataSource"] = "upload"
            cols = list(df.columns)
            time_cand = [c for c in ("dt_s", "dt_ms", "dt_us", "dt_ns", "dt", "delta_t", "deltat",
                                     "time", "t", "x", "channel", "bin") if c in cols]
            cnt_cand = [c for c in ("counts", "count", "n", "yield", "freq", "frequency",
                                    "y", "entries", "events") if c in cols]
            binned_cols = {"t_low", "t_high", "counts"}

            if binned_cols.issubset(set(cols)):
                input_kind = "binned"
                lo = df["t_low"].to_numpy(float) * UF
                hi = df["t_high"].to_numpy(float) * UF
                counts = df["counts"].to_numpy(float)
                keep = (hi > lo) & np.isfinite(counts)
                lo, hi, counts = lo[keep], hi[keep], counts[keep]
                notes.append(f"Pre-binned histogram read from file: {len(counts)} bins "
                             f"({lo.min()/UF:.4g}-{hi.max()/UF:.4g} {unit_name}).")

            elif cnt_cand and (time_cand or len(cols) >= 2):
                # two-column tick histogram: one time column + one counts column
                input_kind = "binned"
                ccol = cnt_cand[0]
                tcol = time_cand[0] if time_cand else [c for c in cols if c != ccol][0]
                tc = pd.to_numeric(df[tcol], errors="coerce").to_numpy(float)
                cc = pd.to_numeric(df[ccol], errors="coerce").to_numpy(float)
                keep = np.isfinite(tc) & np.isfinite(cc)
                tc, cc = tc[keep] * UF, cc[keep]
                if tc.size < 4:
                    raise ValueError(f"Only {tc.size} usable rows in columns '{tcol}'/'{ccol}'.")
                order = np.argsort(tc)
                tc, cc = tc[order], cc[order]
                if edge_ref == "left":
                    step = np.diff(tc)
                    step = np.append(step, step[-1])
                    lo, hi = tc, tc + step
                else:                                   # times are bin centres (default)
                    mid = 0.5 * (tc[1:] + tc[:-1])
                    lo = np.concatenate(([tc[0] - (mid[0] - tc[0])], mid))
                    hi = np.concatenate((mid, [tc[-1] + (tc[-1] - mid[-1])]))
                lo = np.maximum(lo, 0.0)
                counts = cc
                tick = float(np.median(hi - lo))
                notes.append(f"Two-column histogram read: time '{tcol}', contents '{ccol}' - "
                             f"{counts.size} channels of ~{tick/UF:.4g} {unit_name} "
                             f"spanning {lo[0]/UF:.4g}-{hi[-1]/UF:.4g} {unit_name}.")
            else:
                col = time_cand[0] if time_cand else cols[0]
                dt = pd.to_numeric(df[col], errors="coerce").to_numpy(float) * UF
                dt = dt[np.isfinite(dt) & (dt > 0)]
                if dt.size == 0:
                    raise ValueError(f"Column '{col}' contained no positive time differences.")
                notes.append(f"Event list read from column '{col}': {dt.size} implant-decay time differences.")

        if counts is not None and rebin > 1:
            m = (len(counts) // rebin) * rebin
            if m >= rebin * 2:
                lo2, hi2 = lo[:m].reshape(-1, rebin), hi[:m].reshape(-1, rebin)
                lo, hi = lo2[:, 0], hi2[:, -1]
                counts = counts[:m].reshape(-1, rebin).sum(axis=1)
                notes.append(f"Rebinned by a factor of {rebin} -> {len(counts)} bins.")

        if counts is not None and bool(params.get("cropToWindow", False)):
            sel = (hi > t_min) & (lo < t_max)
            if sel.sum() >= 4:
                if int((~sel).sum()):
                    notes.append(f"{int((~sel).sum())} bins fall outside the "
                                 f"{t_min/UF:.4g}-{t_max/UF:.4g} {unit_name} fit window and are excluded.")
                lo, hi, counts = lo[sel], hi[sel], counts[sel]
            else:
                notes.append(f"Fit window {t_min/UF:.4g}-{t_max/UF:.4g} {unit_name} keeps fewer than 4 bins - "
                             f"using the file's full range instead.")

        if counts is None:
            if bin_mode == "log":
                if t_min <= 0:
                    t_min = max(float(np.min(dt)) * 0.5, 1e-6)
                edges = np.logspace(math.log10(t_min), math.log10(t_max), n_bins + 1)
            else:
                edges = np.linspace(max(t_min, 0.0), t_max, n_bins + 1)
            counts, _ = np.histogram(dt, bins=edges)
            counts = counts.astype(float)
            lo, hi = edges[:-1], edges[1:]
            n_out = int(np.sum((dt < lo[0]) | (dt > hi[-1])))
            if n_out:
                notes.append(f"{n_out} events fall outside the {lo[0]:.4g}-{hi[-1]:.4g} s fit window and are excluded.")

        width = hi - lo
        n_b = len(counts)
        if n_b < 4:
            raise ValueError("Fewer than 4 bins in the fit window - widen the range or reduce the bin count.")

        # ------------------------------------------ starting amplitudes from data
        # a typed background rate is meaningless when a channel holds millions of
        # counts, so by default read the flat level off the tail of the spectrum
        span = float(hi[-1] - lo[0])
        rate = counts / np.maximum(width, 1e-300)
        n_tail = max(n_b // 4, 3)
        b_data = float(np.median(rate[-n_tail:]))
        excess = float(np.sum(np.maximum(counts - b_data * width, 0.0)))
        N0_data = max(excess / max(float(eff[0]), 1e-3), 1.0)
        if auto_start:
            bkg_start = max(b_data, 1e-12)
            N0_start = N0_data
            notes.append(f"Starting values estimated from the data: flat level {b_data:.4g} counts s^-1 "
                         f"(median of the last {n_tail} bins), decay amplitude {N0_data:.4g} correlated parents.")
        else:
            bkg_start = bkg0
            N0_start = max(float(counts.sum()), 1.0)
        if not fit_bkg:
            bkg_start = bkg0
        # background half-life: use the value the user typed, else "very long" (20x the
        # spectrum span), which is numerically the same as flat until the fit pulls it in
        t12b_start = bkg_t12_user if bkg_t12_user is not None else max(20.0 * span, 1e-6)

        # -------------------------------------------------------- Bateman model
        def components(t12, bw=None):
            """Detected counts per bin for each member, per implanted parent.

            The general Bateman solution applies to a *linear* sequence, so each member is
            evaluated along its own path back to the parent, with the branching fractions on
            that path multiplied in. For an unbranched chain this is identical to the usual
            formula; with branches, two members fed by the same ancestor are parallel, and
            neither appears in the other's path.
            """
            l = LN2 / np.asarray(t12, float)
            for a in range(n):                      # lift exact degeneracies
                for b_ in range(a + 1, n):
                    if abs(l[a] - l[b_]) < 1e-9 * max(abs(l[a]), 1e-12):
                        l[b_] = l[b_] * (1.0 + 1e-6 * (b_ + 1))
            E = (np.exp(-np.outer(lo, l)) - np.exp(-np.outer(hi, l))) / l   # [bin, member]
            out = np.zeros((len(lo), n))
            for i in range(n):
                p = paths[i]
                m_ = len(p)
                lp = l[p]
                pl = float(np.prod(lp[:-1])) if m_ > 1 else 1.0   # decay constants of the ancestors
                acc = np.zeros(len(lo))
                for k in range(m_):
                    den = 1.0
                    for j in range(m_):
                        if j != k:
                            den *= (lp[j] - lp[k])
                    acc += (pl / den) * E[:, p[k]]
                w_ = brmul[i] if bw is None else bw[i]
                out[:, i] = w_ * l[i] * acc         # integral of lambda_i N_i over each bin
            return out

        def _exp(v):                                # overflow-safe exp in log-space
            return math.exp(min(max(float(v), -700.0), 700.0))

        # background: flat, or a slowly decaying level (long-lived activity / beam-off
        # decay). Both are integrated exactly over each bin.
        # shape: is the background allowed a slope at all?  decay_bkg: is that slope
        # a free parameter, or held at the half-life the user typed?
        shape_decay = bkg_mode == "decaying"
        decay_bkg = shape_decay and bkg_t12_free
        t12b_fixed = t12b_start if (shape_decay and not bkg_t12_free) else None
        if shape_decay and not bkg_t12_free:
            notes.append(f"Background half-life held at {t12b_start/UF:.6g} {unit_name}.")
        # "constrained" asked for but not usable: say so rather than silently floating
        bkg_t12_constrained = bkg_t12_constrained and shape_decay and decay_bkg
        if _bkg_t12_mode == "constrained" and not bkg_t12_constrained:
            notes.append(
                "A constrained background half-life was asked for but not applied: it needs a "
                "decaying background shape, a background half-life and a positive uncertainty "
                "on it. The background half-life was left free.")

        def bkg_bins(b, t12b):
            if t12b is None:
                return b * width
            lb = LN2 / max(float(t12b), 1e-300)
            if lb * span < 1e-8:                    # indistinguishable from flat
                return b * width
            return b * np.exp(-lb * lo) * (-np.expm1(-lb * width)) / lb

        # the tail of the spectrum, where the chain has died away, is almost pure background:
        # fit b*exp(-lambda_b t) to those bins alone and then hold that half-life in the full
        # fit, so the pedestal's slope cannot trade against the parent's half-life
        tail_fit = None
        if shape_decay and bool(params.get("bkgTailFit", False)):
            tl = float(params.get("bkgTailMin", 0) or 0.0) * UF
            th = float(params.get("bkgTailMax", 0) or 0.0) * UF
            if not (th > tl):
                tl, th = lo[0] + 0.75 * span, hi[-1]
            m = (lo >= tl - 1e-12) & (hi <= th + 1e-12)
            if int(m.sum()) < 3:
                notes.append(f"Background tail window {tl/UF:.6g}-{th/UF:.6g} {unit_name} holds fewer than "
                             f"3 bins, so the slope was not estimated from it.")
            else:
                cl, wl, ll = counts[m], width[m], lo[m]

                def tail_cash(z):
                    bb, tb = _exp(z[0]), _exp(z[1])
                    lb = LN2 / max(tb, 1e-300)
                    mm = (bb * np.exp(-lb * ll) * (-np.expm1(-lb * wl)) / lb
                          if lb * span > 1e-8 else bb * wl)
                    mm = np.maximum(mm, 1e-12)
                    v = 2.0 * float(np.sum(mm - cl + np.where(
                        cl > 0, cl * np.log(np.maximum(cl, 1e-12) / mm), 0.0)))
                    return v if np.isfinite(v) else 1e12

                z0 = np.array([math.log(max(float(np.median(cl / wl)), 1e-12)),
                               math.log(max(t12b_start, span))])
                zb = minimize(tail_cash, z0, method="Nelder-Mead",
                              options={"maxiter": 3000, "xatol": 1e-7, "fatol": 1e-7}).x
                b_tail, t12b_tail = _exp(zb[0]), _exp(zb[1])
                if bkg_t12_constrained:
                    notes.append("The tail fit takes precedence over the constrained background "
                                 "half-life: the slope came from this spectrum's own tail, so the "
                                 "external prior was not applied.")
                    bkg_t12_constrained = False
                decay_bkg = False                      # slope now comes from the tail, not the fit
                t12b_fixed = t12b_tail
                t12b_start = t12b_tail
                bkg_start = b_tail
                tail_fit = {"tLo": float(tl), "tHi": float(th), "nBins": int(m.sum()),
                            "t12": float(t12b_tail), "b": float(b_tail)}
                _chain_tail = float(np.sum(np.maximum(cl - b_tail * wl, 0.0)))
                notes.append(
                    f"Background slope taken from the {tl/UF:.6g}-{th/UF:.6g} {unit_name} tail "
                    f"({int(m.sum())} bins): T12(bkg) = {t12b_tail/UF:.6g} {unit_name}, amplitude "
                    f"{b_tail:.6g} counts s^-1. That half-life is then held in the full fit; only the "
                    f"amplitude still floats. Counts in the tail above that pedestal: {_chain_tail:.4g} "
                    f"({100.0 * _chain_tail / max(float(cl.sum()), 1.0):.3f}% of the window) - if that is "
                    f"not small, the chain has not died away yet and the window is too early.")

        # parameter layout:
        # [ln T12 (free members)] [ln N0] [ln b] [ln T12_bkg] [ln eff (free)] [ln branching (free)]
        jN = len(free)
        j_b = jN + 1 if fit_bkg else None
        j_tb = (j_b + 1) if decay_bkg else None
        _base = jN + 1 + (1 if fit_bkg else 0) + (1 if decay_bkg else 0)
        j_eff = {i: _base + m for m, i in enumerate(eff_free)}
        _base2 = _base + len(eff_free)
        j_br = {i: _base2 + m for m, i in enumerate(br_free)}

        def unpack(x):
            t12 = t12_guess.copy()
            for j, i in enumerate(free):
                t12[i] = max(_exp(x[j]), 1e-300)
            N0 = _exp(x[jN])
            b = _exp(x[j_b]) if fit_bkg else bkg0
            t12b = _exp(x[j_tb]) if decay_bkg else t12b_fixed
            ef = eff.copy()
            for i, j in j_eff.items():
                ef[i] = _exp(x[j])
            bf = list(brf)
            for i, j in j_br.items():
                bf[i] = min(max(_exp(x[j]), 1e-9), 1.0)
            for c_, others in comp.items():          # dominant branch = 100% - the rest
                bf[c_] = min(max(1.0 - sum(bf[k] for k in others), 1e-9), 1.0)
            return t12, N0, b, t12b, ef, bf

        def mu_of(x):
            t12, N0, b, t12b, ef, bf = unpack(x)
            per = N0 * ef[None, :] * components(t12, path_weights(bf))
            bk = bkg_bins(b, t12b)
            return np.maximum(per.sum(axis=1) + bk, 1e-12), per, b, N0, t12, t12b, bk, ef, bf

        # positions in the parameter vector that carry a Gaussian prior: constrained
        # half-lives, and branchings with a quoted Pn uncertainty. Both enter the Cash
        # statistic the same way and both are struck out of the unpenalised Hessian, so
        # their contribution lands in the syst half of every error below.
        pen_j = {j: (t12_guess[i], prior[i]) for j, i in enumerate(free) if i in prior}
        for i, j in j_br.items():
            pen_j[j] = (brf[i], br_sig[i])
        # a constrained pedestal slope joins them on exactly the same footing: penalised in
        # the Cash statistic, struck out of the unpenalised Hessian, so what its prior width
        # costs the half-life and N0 shows up in their syst column rather than vanishing
        if bkg_t12_constrained and j_tb is not None:
            pen_j[j_tb] = (bkg_t12_user, bkg_t12_sigma)
            notes.append(
                f"Background half-life constrained to {bkg_t12_user/UF:.6g} +/-{bkg_t12_sigma/UF:.6g} "
                f"{unit_name} by a Gaussian prior - it still floats, but within that width, and the "
                f"width is propagated into the syst half of every error below. Use this on a spectrum "
                f"whose own tail is too thin to measure the slope, with a value and uncertainty taken "
                f"from one that can.")

        def cash_stat(x):
            mu = mu_of(x)[0]
            term = np.where(counts > 0, counts * np.log(np.maximum(counts, 1e-12) / mu), 0.0)
            v = 2.0 * float(np.sum(mu - counts + term))
            return v if np.isfinite(v) else 1e12

        def penalty(x):
            # (T - T_lit)^2 / sigma^2 added to C: the standard way to let a measured
            # daughter half-life float within its own error bar
            p = 0.0
            for j, (mu_t, s_t) in pen_j.items():
                p += ((_exp(x[j]) - mu_t) / s_t) ** 2
            return p

        def cash(x):
            return cash_stat(x) + penalty(x) if pen_j else cash_stat(x)

        x0 = [math.log(t12_guess[i]) for i in free] + [math.log(N0_start)]
        bnds = [(math.log(1e-9), math.log(1e9))] * len(free) + [(math.log(1e-6), math.log(1e18))]
        if fit_bkg:
            x0.append(math.log(bkg_start))
            bnds.append((math.log(1e-12), math.log(1e18)))
        if decay_bkg:
            x0.append(math.log(t12b_start))
            # wide enough to reach a background that decays within the window as well as
            # one that is effectively flat, and always wide enough to hold the start value
            bnds.append((math.log(min(span * 0.01, t12b_start * 0.01, 1e-6)),
                         math.log(max(span * 1e9, t12b_start * 100.0))))
        # an efficiency is the fraction of that member's decays that land in this spectrum, so
        # physically it cannot exceed 1. Left unbounded it is the parameter that quietly absorbs
        # anything the model is missing, which is informative - hence the switch rather than a
        # hard-wired cap.
        cap_eff = bool(params.get("capEff", False))
        eff_hi_bound = 1.0 if cap_eff else 10.0
        for i in eff_free:
            x0.append(math.log(min(max(float(eff[i]), 1e-4), eff_hi_bound)))
            bnds.append((math.log(1e-6), math.log(eff_hi_bound)))
        # branchings start at the quoted Pn and are bounded by physics: a fraction between
        # 0 and 1. The Gaussian prior, not the bound, is what holds them near the literature value.
        for i in br_free:
            x0.append(math.log(max(brf[i], 1e-6)))
            bnds.append((math.log(1e-6), math.log(1.0)))
        x0 = np.array(x0, float)

        def _solve(start):
            y = minimize(cash, start, method="Nelder-Mead",
                         options={"maxiter": 4000, "xatol": 1e-6, "fatol": 1e-6}).x
            r = minimize(cash, y, method="L-BFGS-B", bounds=bnds)
            return r, float(r.fun)

        # the background half-life start value decides which basin the minimiser lands in:
        # start it too short and the chain swallows the pedestal (background amplitude runs
        # to zero, N0 explodes, the parent runs off to seconds). Try several and keep the
        # best Cash rather than trusting one start.
        starts = [x0]
        if decay_bkg:
            b_lo, b_hi = bnds[j_tb]
            for mult in (1.0, 20.0, 200.0, 5000.0):
                cand = x0.copy()
                cand[j_tb] = min(max(math.log(mult * span), b_lo), b_hi)
                starts.append(cand)
            if bkg_t12_user is not None and bkg_t12_user < span:
                notes.append(f"The background half-life you typed ({bkg_t12_user/UF:.6g} {unit_name}) is "
                             f"shorter than the {span/UF:.6g} {unit_name} spectrum, so it decays away "
                             f"inside the window and behaves like another chain member - a classic route "
                             f"into a degenerate solution.")

        best = None
        for st in starts:
            try:
                r, c = _solve(st)
            except Exception:
                continue
            if best is None or c < best[1]:
                best = (r, c)
        if best is None:
            raise RuntimeError("The minimiser failed from every starting point.")
        fit, Cmin = best[0], best[1]
        xh = fit.x
        if len(starts) > 1:
            notes.append(f"Fit repeated from {len(starts)} background half-life starting points; "
                         f"the one with the lowest Cash ({Cmin:.1f}) is reported.")

        # -------------------------------------------- covariance from Hessian(C)
        def hessian(f, x, h=1e-4):
            m = len(x); H = np.zeros((m, m)); f0 = f(x)
            for i in range(m):
                for j in range(i, m):
                    if i == j:
                        xp = x.copy(); xp[i] += h
                        xm = x.copy(); xm[i] -= h
                        H[i, i] = (f(xp) - 2 * f0 + f(xm)) / h ** 2
                    else:
                        a = x.copy(); a[i] += h; a[j] += h
                        b_ = x.copy(); b_[i] += h; b_[j] -= h
                        c_ = x.copy(); c_[i] -= h; c_[j] += h
                        d_ = x.copy(); d_[i] -= h; d_[j] -= h
                        H[i, j] = H[j, i] = (f(a) - f(b_) - f(c_) + f(d_)) / (4 * h ** 2)
            return H

        try:
            cov = 2.0 * np.linalg.inv(hessian(cash, xh))   # C = 2*NLL  ->  cov = 2 H_C^-1
            sig = np.sqrt(np.clip(np.diag(cov), 0.0, None))
        except Exception:
            cov = None
            sig = np.full(len(xh), np.nan)
            notes.append("Hessian was singular - symmetric errors unavailable; profile-likelihood bounds still reported.")

        # how much of each error is the daughter's uncertainty?  invert the unpenalised
        # Hessian with the constrained parameters struck out: that is the error the fit
        # would have if the literature half-lives were exact. The rest, in quadrature,
        # is what their quoted uncertainty contributes.
        sig_stat = np.full(len(xh), np.nan)
        cov_stat = np.full((len(xh), len(xh)), np.nan)
        if pen_j:
            try:
                Hn = hessian(cash_stat, xh)
                keep_k = [k for k in range(len(xh)) if k not in pen_j]
                cs = 2.0 * np.linalg.inv(Hn[np.ix_(keep_k, keep_k)])
                d = np.sqrt(np.clip(np.diag(cs), 0.0, None))
                for a, k in enumerate(keep_k):
                    sig_stat[k] = d[a]
                cov_stat[np.ix_(keep_k, keep_k)] = cs
            except Exception:
                pass

        def split_err(j, value):
            """(total, statistical, systematic) absolute errors for parameter slot j."""
            tot = float(value * sig[j]) if np.isfinite(sig[j]) else None
            st = float(value * sig_stat[j]) if np.isfinite(sig_stat[j]) else None
            sy = None
            if tot is not None and st is not None:
                sy = float(math.sqrt(max(tot * tot - st * st, 0.0)))
            return tot, st, sy

        def comb_err(weights, value):
            """Error on a product of fitted parameters, e.g. A0 = ln2/T12 * N0.

            Every parameter is fitted in log space, so a product/ratio is a linear
            combination there: w = +1 for a factor, -1 for a divisor.  The variance is
            w^T cov w, which keeps the T12-N0 correlation - ignoring it would badly
            over-state the error on A0, because T12 and N0 move together (a longer
            half-life is paid for by more implants) and the correlation cancels in the ratio.
            """
            def one(C):
                if C is None:
                    return None
                w = np.zeros(len(xh))
                for j, wt in weights.items():
                    w[j] = wt
                M = np.asarray(C, float).copy()
                M[~np.isfinite(M)] = 0.0      # struck-out (constrained) slots carry NaN; w is 0 there
                v = float(w @ M @ w)
                return float(value * math.sqrt(v)) if np.isfinite(v) and v >= 0 else None
            tot = one(cov)
            st = one(cov_stat) if pen_j else None
            if st is not None and not np.isfinite(st):
                st = None
            sy = None
            if tot is not None and st is not None:
                sy = float(math.sqrt(max(tot * tot - st * st, 0.0)))
            return tot, st, sy

        # ------------------------------------------- profile likelihood (dC = 1)
        profiles = []
        prof_bounds = {}
        n_grid, n_it = (21, 800) if n_b <= 256 else (13, 250)   # keep fine-binned files responsive
        for j, i in enumerate(free[:4]):
            s = sig[j] if np.isfinite(sig[j]) and 0 < sig[j] < 5.0 else 0.3
            grid = np.clip(np.linspace(xh[j] - 3.0 * s, xh[j] + 3.0 * s, n_grid), -700.0, 700.0)
            others = [k for k in range(len(xh)) if k != j]
            dC = []
            for g in grid:
                def f2(y, g=g, j=j, others=others):
                    xx = xh.copy(); xx[j] = g
                    for a, k in enumerate(others):
                        xx[k] = y[a]
                    return cash(xx)
                r2 = minimize(f2, xh[others], method="Nelder-Mead",
                              options={"maxiter": n_it, "xatol": 1e-5, "fatol": 1e-5}) if others else None
                dC.append((float(r2.fun) if r2 is not None else cash(np.array([g]))) - Cmin)
            dC = np.array(dC)

            def cross(side):
                idx = range(len(grid) - 1) if side > 0 else range(len(grid) - 1, 0, -1)
                for k in idx:
                    k2 = k + 1 if side > 0 else k - 1
                    if (grid[k] - xh[j]) * side >= 0 and (dC[k] - 1.0) * (dC[k2] - 1.0) < 0:
                        f = (1.0 - dC[k]) / (dC[k2] - dC[k])
                        return math.exp(grid[k] + f * (grid[k2] - grid[k]))
                return None

            prof_bounds[i] = (cross(-1), cross(+1))
            profiles.append({"label": labels[i],
                             "t12Grid": np.exp(grid).tolist(),
                             "deltaC": dC.tolist()})

        # ------------------------------------------------- parent half-life scan
        # A gamma-gated fit measures one nuclide in isolation; this ungated spectrum measures
        # it through everything sitting on top of it. When the two disagree, the question is
        # not which minimiser is better but how much the likelihood actually cares - and what
        # the rest of the model does when the parent is forced to the gated value. So: hold
        # T12(parent) at each point of a grid, re-minimise everything else, and report the
        # cost in Cash together with the Pn, N0 and pedestal the fit needs at that point.
        scan = None
        if bool(params.get("scanParent", False)):
            if 0 not in free:
                notes.append("Parent half-life scan skipped: the parent's half-life is fixed, so "
                             "there is nothing to scan. Set it to free or constrained.")
            else:
                j0 = free.index(0)
                t12_hat = unpack(xh)[0]                 # the output section's t12_fit, computed early
                s_lo = float(params.get("scanMin", 0) or 0.0) * UF
                s_hi = float(params.get("scanMax", 0) or 0.0) * UF
                if not (s_hi > s_lo > 0):
                    _s = sig[j0] if np.isfinite(sig[j0]) and 0 < sig[j0] < 2.0 else 0.15
                    s_lo = float(t12_hat[0]) * math.exp(-4.0 * _s)
                    s_hi = float(t12_hat[0]) * math.exp(+4.0 * _s)
                npts = int(min(max(int(params.get("scanPoints", 15) or 15), 5), 41))
                g0 = np.linspace(math.log(s_lo), math.log(s_hi), npts)
                oth = [k for k in range(len(xh)) if k != j0]
                pts = []
                for g in g0:
                    xx = xh.copy()
                    xx[j0] = g
                    if oth:
                        def f3(y, g=g, oth=oth):
                            z = xh.copy(); z[j0] = g
                            for a, k in enumerate(oth):
                                z[k] = y[a]
                            return cash(z)
                        r3 = minimize(f3, xh[oth], method="Nelder-Mead",
                                      options={"maxiter": 2000, "xatol": 1e-6, "fatol": 1e-6})
                        for a, k in enumerate(oth):
                            xx[k] = r3.x[a]
                    _t12, _N0, _b, _t12b, _ef, _bf = unpack(xx)
                    _lam = LN2 / max(float(_t12[0]), 1e-300)
                    pts.append({
                        "t12": float(_t12[0]),
                        "deltaC": float(cash(xx) - Cmin),
                        "cash": float(cash_stat(xx)),
                        "N0": float(_N0),
                        "A0": float(_lam * _N0),
                        "bkg": float(_b),
                        "branch": {labels[i]: 100.0 * float(_bf[i]) for i in br_free},
                        "eff": {labels[i]: float(_ef[i]) for i in eff_free},
                    })
                scan = {"label": labels[0], "points": pts,
                        "best": float(t12_hat[0]), "cmin": float(Cmin)}

                # how far can the parent be pushed before the ungated data objects?  dC = 1 is
                # 1 sigma, dC = 9 is 3 sigma, on the same profile-likelihood scale as the errors
                def _cross_at(level):
                    out = []
                    for k in range(len(pts) - 1):
                        a_, b2 = pts[k]["deltaC"] - level, pts[k + 1]["deltaC"] - level
                        if a_ * b2 < 0:
                            f_ = a_ / (a_ - b2)
                            out.append(pts[k]["t12"] + f_ * (pts[k + 1]["t12"] - pts[k]["t12"]))
                    return out
                _c1, _c9 = _cross_at(1.0), _cross_at(9.0)
                scan["sigma1"] = _c1
                scan["sigma3"] = _c9
                _pn_swing = ""
                if br_free:
                    _lo_pt, _hi_pt = pts[0], pts[-1]
                    _pn_swing = " Across the scan the fitted branching moves " + "; ".join(
                        f"{labels[i]}: {_lo_pt['branch'][labels[i]]:.3g}% -> "
                        f"{_hi_pt['branch'][labels[i]]:.3g}%" for i in br_free) + "."
                notes.append(
                    f"Parent scan: T12({labels[0]}) held across {s_lo/UF:.4g}-{s_hi/UF:.4g} {unit_name} "
                    f"in {npts} steps, everything else re-minimised at each point. "
                    + (f"dC = 1 (1 sigma) at " + ", ".join(f"{v/UF:.4g}" for v in _c1) + f" {unit_name}; "
                       if _c1 else "")
                    + (f"dC = 9 (3 sigma) at " + ", ".join(f"{v/UF:.4g}" for v in _c9) + f" {unit_name}. "
                       if _c9 else "")
                    + "Read dC as the price the ungated spectrum charges for a half-life: below about 1 "
                      "the two values are equally good fits to these data, and the difference is not a "
                      "disagreement."
                    + _pn_swing)
                pd.DataFrame([{
                    f"t12_{unit_name}": p["t12"] / UF, "t12_s": p["t12"],
                    "delta_cash": p["deltaC"], "cash": p["cash"],
                    "N0": p["N0"], "A0_per_s": p["A0"], "background_per_s": p["bkg"],
                    **{f"branch_{k}_percent": v for k, v in p["branch"].items()},
                    **{f"eff_{k}": v for k, v in p["eff"].items()},
                } for p in pts]).to_csv(OUT("parent_scan.csv"), index=False)

        # ---------------------------------------------------------------- output
        mu, per, b_fit, N0, t12_fit, t12b_fit, bk_fit, eff_fit, brf_fit = mu_of(xh)
        brmul_fit = path_weights(brf_fit)
        k_par = (len(free) + 1 + (1 if fit_bkg else 0) + (1 if decay_bkg else 0)
                 + len(eff_free) + len(br_free))
        pearson = (counts - mu) / np.sqrt(mu)
        dev_term = 2.0 * (mu - counts + np.where(counts > 0, counts * np.log(np.maximum(counts, 1e-12) / mu), 0.0))
        deviance = np.sign(counts - mu) * np.sqrt(np.clip(dev_term, 0.0, None))
        chi2v = float(np.sum(pearson ** 2))
        ndf = max(n_b - k_par, 1)

        # ------------------------------------------------------ full covariance matrix
        # Every error quoted above is a diagonal element of this. Propagating a derived
        # quantity - a branching that divides counts by N0, an integrated activity - with
        # errors added in quadrature assumes the off-diagonals are zero, and here they are
        # emphatically not: rho(lnT12, lnN0) alone runs above +0.9. The correct variance is
        # w^T cov w with w the vector of partial derivatives, which needs this matrix.
        pnames, pvals = [], []
        for i in free:
            pnames.append(f"T12_{labels[i]}"); pvals.append(float(t12_fit[i]))
        pnames.append(f"N0_{labels[0]}"); pvals.append(float(N0))
        if fit_bkg:
            pnames.append("background_amplitude"); pvals.append(float(b_fit))
        if decay_bkg:
            pnames.append("T12_background"); pvals.append(float(t12b_fit))
        for i in eff_free:
            pnames.append(f"eff_{labels[i]}"); pvals.append(float(eff_fit[i]))
        for i in br_free:
            pnames.append(f"branch_{labels[i]}"); pvals.append(float(brf_fit[i]))
        pv = np.asarray(pvals, float)

        cov_block = None
        if cov is not None and len(pnames) == cov.shape[0]:
            # the fit runs in log space, so cov is the covariance of ln(parameter). The
            # relative form is what most propagation formulae want; the absolute form is
            # cov_ln scaled by the values themselves, and is written out too.
            outer = np.outer(pv, pv)
            d = np.sqrt(np.clip(np.diag(cov), 0.0, None))
            with np.errstate(divide="ignore", invalid="ignore"):
                corr = cov / np.outer(d, d)
            corr = np.where(np.isfinite(corr), corr, np.nan)
            mats = [("covariance_absolute", cov * outer),
                    ("covariance_log", cov),
                    ("correlation", corr)]
            if pen_j and np.isfinite(cov_stat).any():
                mats.insert(1, ("covariance_absolute_stat_only", cov_stat * outer))
            blocks = []
            for kind, M in mats:
                for a, nm in enumerate(pnames):
                    row = {"matrix": kind, "parameter": nm, "value": pvals[a],
                           "unit": ("s" if nm.startswith("T12") else
                                    "parents" if nm.startswith("N0") else
                                    "counts/s" if nm == "background_amplitude" else
                                    "fraction" if nm.startswith("branch") else
                                    "ratio" if nm.startswith("eff") else "")}
                    row.update({nm2: float(M[a, b_]) for b_, nm2 in enumerate(pnames)})
                    blocks.append(row)
            pd.DataFrame(blocks).to_csv(OUT("covariance.csv"), index=False)
            cov_block = {"names": pnames, "values": pvals,
                         "corr": [[None if not np.isfinite(corr[a, b_]) else float(corr[a, b_])
                                   for b_ in range(len(pnames))] for a in range(len(pnames))],
                         "cov": [[float(cov[a, b_] * pv[a] * pv[b_]) for b_ in range(len(pnames))]
                                 for a in range(len(pnames))]}
            _worst = ""
            _off = [(abs(corr[a, b_]), pnames[a], pnames[b_], corr[a, b_])
                    for a in range(len(pnames)) for b_ in range(a + 1, len(pnames))
                    if np.isfinite(corr[a, b_])]
            if _off:
                _off.sort(reverse=True)
                _worst = (" The strongest off-diagonal term is "
                          f"rho({_off[0][1]}, {_off[0][2]}) = {_off[0][3]:+.3f}"
                          + (f", then rho({_off[1][1]}, {_off[1][2]}) = {_off[1][3]:+.3f}"
                             if len(_off) > 1 else "") + ".")
            notes.append(
                f"Covariance written to {os.path.basename(OUT('covariance.csv'))} for all "
                f"{len(pnames)} free parameters, as the absolute matrix (units: parameter x "
                f"parameter, half-lives in seconds), the log-space matrix the fit actually uses, "
                f"and the correlation matrix"
                + (", plus the stat-only absolute matrix with the constrained parameters removed"
                   if pen_j else "")
                + f". Use sigma^2(f) = w^T cov w with w the partial derivatives of your derived "
                  f"quantity; adding errors in quadrature instead assumes every off-diagonal is "
                  f"zero." + _worst)

        rows = []
        for i in range(n):
            j = free.index(i) if i in free else None
            s_t = s_stat = s_syst = None
            if j is not None:
                s_t, s_stat, s_syst = split_err(j, t12_fit[i])
            plo, phi = prof_bounds.get(i, (None, None))
            # a branching that floats within its Pn error gets the same stat/syst split as a
            # constrained half-life, reported in percentage points
            _bre = (None, None, None)
            if i in j_br:
                _bre = tuple((None if v is None else 100.0 * v)
                             for v in split_err(j_br[i], brf_fit[i]))
            elif i in comp:
                # the tied branch is 1 - sum(siblings), so its error is the siblings' errors
                # added in quadrature - the tie makes it move exactly opposite to them
                def _q(sel):
                    vs = [split_err(j_br[k], brf_fit[k])[sel] for k in comp[i] if k in j_br]
                    vs = [v for v in vs if v is not None and np.isfinite(v)]
                    return (100.0 * math.sqrt(sum(v * v for v in vs))) if vs else None
                _bre = (_q(0), _q(1), _q(2))
            rows.append({
                "label": labels[i], "t12": float(t12_fit[i]), "sigma": s_t,
                "sigmaStat": s_stat, "sigmaSyst": s_syst,
                "constrained": i in prior,
                "priorSigma": float(prior[i]) if i in prior else None,
                "loBound": plo, "hiBound": phi, "fixed": i not in free,
                "fedBy": (labels[par[i]] if par[i] >= 0 else None),
                "fedByRow": (par[i] + 1) if par[i] >= 0 else None,
                "branchPercent": (100.0 * brf_fit[i]) if i > 0 else None,
                "branchStart": (100.0 * brf[i]) if i > 0 else None,
                "branchFree": i in br_free,
                "branchTied": i in comp,
                "branchErr": _bre[0], "branchErrStat": _bre[1], "branchErrSyst": _bre[2],
                "pathWeight": float(brmul_fit[i]),
                "eff": float(eff_fit[i]), "effFree": i in eff_free,
                "effErr": (float(eff_fit[i] * sig[j_eff[i]])
                           if (i in eff_free and np.isfinite(sig[j_eff[i]])) else None),
                "counts": float(per[:, i].sum()),
            })

        # a background half-life much longer than the window is indistinguishable from a flat
        # pedestal: the likelihood is flat along that direction, so the fitted value is only a
        # lower limit even though the half-lives and N0 sitting on top of it are perfectly stable
        bkg_plateau = False
        if decay_bkg and t12b_fit is not None and not bkg_t12_constrained:
            # with a prior on it the slope is no longer a lower limit, so the plateau
            # warning below would be actively misleading - skip it in that case
            rel = float(sig[j_tb]) if np.isfinite(sig[j_tb]) else np.inf
            bkg_plateau = (t12b_fit > 5.0 * span) or (rel > 0.5)
            if bkg_plateau:
                notes.append(
                    f"The fitted background half-life ({t12b_fit/UF:.4g} {unit_name}) is long compared "
                    f"with the {span/UF:.4g} {unit_name} fit window, so the background is effectively "
                    "flat here and this number is a lower limit, not a measurement - the likelihood "
                    "barely changes as it grows. That is the expected outcome when both background "
                    "parameters float, and it does not harm the chain: the half-lives and N0 are the "
                    "ones the data constrains.")

        _n0e = split_err(jN, N0)

        # what a decay curve actually measures is the initial activity A0 = lambda1 * N0:
        # the height of the parent's exponential at t = 0.  N0 is that height divided by the
        # fitted decay constant, so it inherits the half-life's error - and because T12 and N0
        # are strongly correlated (a longer half-life is paid for by more implants, leaving the
        # t = 0 height unchanged), A0 is far better determined than either factor on its own.
        lam1 = LN2 / float(t12_fit[0])
        A0 = lam1 * float(N0)
        _w = {jN: 1.0}
        if 0 in free:
            _w[free.index(0)] = -1.0                       # A0 ~ N0 / T12
        _a0e = comb_err(_w, A0)
        rho_tn = None
        if cov is not None and 0 in free:
            _jt = free.index(0)
            _d = math.sqrt(max(cov[_jt, _jt] * cov[jN, jN], 0.0))
            if _d > 0:
                rho_tn = float(cov[_jt, jN] / _d)
        notes.append(
            f"Initial activity A0 = ln2/T12 x N0 = {A0:.6g} decays s^-1 at t = 0 "
            f"({A0 * UF:.6g} per {unit_name}), for the parent alone." +
            (f" T12 and N0 are correlated with rho = {rho_tn:+.3f}, which is why A0 carries a "
             f"smaller relative error than N0 does." if rho_tn is not None else "") +
            f" The rate you would actually count is A0 x eff({labels[0]}) = "
            f"{A0 * float(eff_fit[0]):.6g} s^-1, with the parent efficiency held at "
            f"{float(eff_fit[0]):.4g}.")

        # an efficiency is a pure scale factor on that member's Bateman shape, so it only means
        # "fraction detected" as long as the shape is distinguishable from everything else in the
        # fit. A member whose half-life is much longer than the window contributes an essentially
        # flat curve, which is degenerate with the background pedestal - its efficiency then stops
        # measuring detection and starts absorbing background, and can run past 1.
        med_w = float(np.median(width))
        flat_members = [i for i in range(n) if float(t12_fit[i]) > 2.0 * span]
        hot_eff = [i for i in eff_free if float(eff_fit[i]) > 1.0]
        for i in flat_members:
            frac = float(per[:, i].sum()) / max(float(counts.sum()), 1.0)
            notes.append(
                f"{labels[i]} has T12 = {t12_fit[i]/UF:.6g} {unit_name}, longer than the "
                f"{span/UF:.6g} {unit_name} window, so its decay curve is nearly flat across the "
                f"spectrum and is degenerate with the background pedestal. It carries "
                f"{per[:, i].sum():.4g} counts ({100.0*frac:.4f}% of the window)" +
                (f", and its fitted efficiency {eff_fit[i]:.3g} is absorbing background rather than "
                 f"measuring detection - do not read it as a detection fraction."
                 if i in hot_eff or float(eff_fit[i]) > 1.0 else
                 ", so its efficiency is barely constrained by these data.") +
                " Holding that efficiency fixed, or shortening the chain, removes the ambiguity.")
        for i in hot_eff:
            if i in flat_members:
                continue
            notes.append(
                f"Fitted efficiency for {labels[i]} is {eff_fit[i]:.3g}, above 1. Efficiencies here are "
                f"relative to the parent's, which is held at {float(eff_fit[0]):.4g}: a value above 1 means "
                f"this member deposits more counts per parent decay than the parent does, which happens "
                f"when its branch feeds more than one detected transition, when the parent's efficiency "
                f"is set too low, or when the component is soaking up mismodelled background. "
                f"Tick 'cap efficiencies at 1' to forbid it and see how much the fit degrades.")
        if cap_eff and any(float(eff_fit[i]) > 0.999 for i in eff_free):
            notes.append("One or more efficiencies are sitting on the upper bound of 1, so their errors "
                         "are not meaningful - the fit wants to go higher than physics allows.")

        if pen_j:
            notes.append(
                "Constrained half-lives and branchings are fitted, not frozen: each is pulled towards "
                "the value you typed by its own uncertainty. Every error below is therefore quoted as "
                "+/-stat +/-syst, where syst is the part contributed by those literature uncertainties "
                "(total^2 = stat^2 + syst^2)." +
                (f" {len(j_br)} branching(s) float under a Pn prior, so the parent's syst error now "
                 f"includes what a poorly known Pn costs it." if j_br else ""))

        res.update({
            "fit": rows,
            "N0": float(N0),
            "N0Err": _n0e[0], "N0ErrStat": _n0e[1], "N0ErrSyst": _n0e[2],
            "A0": float(A0),
            "A0Err": _a0e[0], "A0ErrStat": _a0e[1], "A0ErrSyst": _a0e[2],
            "A0Detected": float(A0 * float(eff_fit[0])),
            "A0Rho": rho_tn,
            # A0 is a rate per second; the plot's y-axis is counts per channel. Export both so
            # the number in the table can be checked directly against the height of the curve.
            "binWidth": med_w,
            "A0PerBin": float(A0 * med_w),
            "A0PerBinErr": (float(_a0e[0] * med_w) if _a0e[0] is not None else None),
            "tailFit": tail_fit,
            "capEff": bool(cap_eff),
            "flatMembers": [labels[i] for i in flat_members],
            "hotEff": [labels[i] for i in hot_eff],
            "parentLabel": labels[0],
            "parentEff": float(eff_fit[0]),
            "branched": bool(branched),
            "topology": [{"row": i + 1, "label": labels[i],
                          "fedByRow": (par[i] + 1) if par[i] >= 0 else None,
                          "fedBy": labels[par[i]] if par[i] >= 0 else None,
                          "branchPercent": (100.0 * brf_fit[i]) if i > 0 else None,
                          "branchStart": (100.0 * brf[i]) if i > 0 else None,
                          "branchFree": i in br_free,
                          "branchTied": i in comp,
                          "pathWeight": float(brmul_fit[i]),
                          "path": [labels[q] for q in paths[i]]} for i in range(n)],
            "bkg": float(b_fit),
            "bkgErr": float(b_fit * sig[j_b]) if (fit_bkg and np.isfinite(sig[j_b])) else None,
            "bkgT12": float(t12b_fit) if t12b_fit is not None else None,
            "bkgT12Err": (float(t12b_fit * sig[j_tb])
                          if (decay_bkg and np.isfinite(sig[j_tb])) else None),
            "bkgT12Free": bool(decay_bkg),
            "bkgT12Constrained": bool(bkg_t12_constrained),
            "bkgT12Prior": (float(bkg_t12_user) if bkg_t12_constrained else None),
            "bkgT12PriorSigma": (float(bkg_t12_sigma) if bkg_t12_constrained else None),
            "bkgT12Plateau": bool(bkg_plateau),
            "bkgCounts": float(np.asarray(bk_fit, float).sum()),
            "fitBackground": fit_bkg,
            "converged": bool(fit.success),
            "stats": {"cash": cash_stat(xh), "cashPenalised": Cmin, "chi2": chi2v, "ndf": ndf, "redChi2": chi2v / ndf,
                      "pValue": float(chi2dist.sf(chi2v, ndf)), "nBins": n_b, "nParams": k_par,
                      "totalCounts": float(counts.sum()), "modelCounts": float(mu.sum())},
            "spectrum": {"lo": lo.tolist(), "hi": hi.tolist(),
                         "center": (np.where(lo > 0, np.sqrt(np.maximum(lo, 1e-300) * hi), 0.5 * (lo + hi))
                                    if bin_mode == "log" else 0.5 * (lo + hi)).tolist(),
                         "width": width.tolist(), "counts": counts.tolist(), "mu": mu.tolist(),
                         "comps": [per[:, i].tolist() for i in range(n)],
                         "bkg": np.asarray(bk_fit, float).tolist()},
            "residuals": {"pearson": np.where(np.isfinite(pearson), pearson, 0.0).tolist(),
                          "deviance": np.where(np.isfinite(deviance), deviance, 0.0).tolist()},
            "profiles": profiles,
            "parentScan": scan,
            "covariance": cov_block,
            "inputKind": input_kind,
            "binMode": bin_mode,
            "timeUnit": unit_name,
            "unitFactor": UF,
            "bkgMode": "decaying" if shape_decay else "flat",
            "notes": notes,
        })

        pd.DataFrame([{"nuclide": r["label"], "half_life_s": r["t12"],
                       "sigma_s": r["sigma"], "sigma_stat_s": r["sigmaStat"], "sigma_syst_s": r["sigmaSyst"],
                       "profile_lo_s": r["loBound"], "profile_hi_s": r["hiBound"],
                       "status": "fixed" if r["fixed"] else ("constrained" if r["constrained"] else "free"),
                       "prior_sigma_s": r["priorSigma"],
                       "fed_by": r["fedBy"], "branch_percent": r["branchPercent"],
                       "branch_percent_sigma": r["branchErr"],
                       "branch_percent_sigma_stat": r["branchErrStat"],
                       "branch_percent_sigma_syst": r["branchErrSyst"],
                       "branch_start_percent": r["branchStart"],
                       "branch_status": ("free (Pn prior)" if r["branchFree"]
                                         else "tied (100% - siblings)" if r.get("branchTied")
                                         else "fixed"),
                       "path_weight": r["pathWeight"],
                       "efficiency": r["eff"], "efficiency_sigma": r["effErr"],
                       "efficiency_status": "free" if r["effFree"] else "fixed",
                       "fitted_counts": r["counts"]} for r in rows]
                     ).to_csv(OUT("fit_results.csv"), index=False)
        pd.DataFrame([
            {"quantity": f"N0 (correlated {labels[0]} implants)", "value": float(N0), "unit": "parents",
             "sigma": _n0e[0], "sigma_stat": _n0e[1], "sigma_syst": _n0e[2]},
            {"quantity": f"A0 (initial {labels[0]} activity)", "value": float(A0), "unit": "decays/s",
             "sigma": _a0e[0], "sigma_stat": _a0e[1], "sigma_syst": _a0e[2]},
            {"quantity": f"A0 per {med_w/UF:.6g} {unit_name} channel (height of the curve at dt=0)",
             "value": float(A0 * med_w), "unit": "counts/channel",
             "sigma": (_a0e[0] * med_w) if _a0e[0] is not None else None,
             "sigma_stat": (_a0e[1] * med_w) if _a0e[1] is not None else None,
             "sigma_syst": (_a0e[2] * med_w) if _a0e[2] is not None else None},
            {"quantity": f"A0 x efficiency (detected initial rate)",
             "value": float(A0 * float(eff_fit[0])), "unit": "counts/s",
             "sigma": (_a0e[0] * float(eff_fit[0])) if _a0e[0] is not None else None,
             "sigma_stat": (_a0e[1] * float(eff_fit[0])) if _a0e[1] is not None else None,
             "sigma_syst": (_a0e[2] * float(eff_fit[0])) if _a0e[2] is not None else None},
            {"quantity": "background amplitude at t=0", "value": float(b_fit), "unit": "counts/s",
             "sigma": (float(b_fit * sig[j_b]) if (fit_bkg and np.isfinite(sig[j_b])) else None),
             "sigma_stat": None, "sigma_syst": None},
            {"quantity": "correlation rho(lnT12, lnN0)", "value": rho_tn, "unit": "",
             "sigma": None, "sigma_stat": None, "sigma_syst": None},
        ]).to_csv(OUT("derived_quantities.csv"), index=False)

        # this file is what a paper figure gets drawn from, so it carries every curve the
        # interface plots, in the unit the spectrum was written in as well as in seconds,
        # plus the Poisson error bar and the background-subtracted view.
        _bkf = np.asarray(bk_fit, float)
        # geometric bin centre on a log axis, arithmetic otherwise - the same x the plot uses
        ctr = (np.where(lo > 0, np.sqrt(np.maximum(lo, 1e-300) * hi), 0.5 * (lo + hi))
               if bin_mode == "log" else 0.5 * (lo + hi))
        _sp = {}
        if unit_name != "s":
            _sp.update({f"t_low_{unit_name}": lo / UF, f"t_high_{unit_name}": hi / UF,
                        f"t_center_{unit_name}": ctr / UF, f"bin_width_{unit_name}": width / UF})
        _sp.update({"t_low_s": lo, "t_high_s": hi, "t_center_s": ctr, "bin_width_s": width,
               "counts": counts, "counts_err": np.sqrt(np.maximum(counts, 0.0)),
               "model_total": mu, "model_background": _bkf})
        for i in range(n):
            _sp[f"model_{labels[i]}"] = per[:, i]
        _sp["model_chain_total"] = mu - _bkf
        _sp["counts_bkg_subtracted"] = counts - _bkf
        _sp["residual_pearson"] = pearson
        _sp["residual_deviance"] = deviance
        pd.DataFrame(_sp).to_csv(OUT("spectrum_fit.csv"), index=False)

    except Exception as e:
        res["error"] = f"{type(e).__name__}: {e}"
        res["traceback"] = traceback.format_exc()

    def _clean(o):
        if isinstance(o, dict):
            return {k: _clean(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return [_clean(v) for v in o]
        if isinstance(o, (np.floating, float)):
            return float(o) if math.isfinite(float(o)) else None
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.bool_,)):
            return bool(o)
        return o

    with open(os.path.join(outdir, "results.json"), "w") as f:
        json.dump(_clean(res), f, allow_nan=False)
    if verbose: print("done:", res.get("error") or f"{len(res.get('fit', []))} chain members fitted")

    return res
