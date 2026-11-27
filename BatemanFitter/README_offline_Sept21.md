# Bateman Half-Life Fitter — offline version

The fit engine from the web tool, packaged to run on your own machine. No network
access at any point: the data file never leaves the computer you run this on.

This is what you use for ⁴¹P / ⁴²P / ⁴³P / ⁴⁴P. The ⁴⁴S spectrum is already
uploaded to the web tool, so there is no reason not to keep refining the
procedure there — and anything you settle on there transfers here unchanged,
because it is the same code.

---

## 1. Install

Python 3.10 or newer.

```bash
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

That pulls numpy, scipy, pandas (the engine), streamlit and plotly (the
interface), and matplotlib (the publication tornado figure, and `--plot`).

## 2. Run the graphical version

```bash
streamlit run app.py
```

It opens a browser tab on `http://localhost:8501`. The browser is only the
display surface — the process is on your machine and serves only to it.

Everything the web tool has is here: the chain table, the binning and window
controls, background handling, the systematic grid scan, the nine result tabs
(Fit, Residuals, Parameters, Profiles, Parent scan, Correlation, Notes,
Systematics, Files), and download buttons for every CSV.

Upload your spectrum with the file box in the sidebar, or type a path if the
file is large.

**Settings files.** The top of the sidebar has a *Settings file* box: drop
`config_44S.json` in it, press **Load these settings**, and the chain table and
every control fill in from the file. **Save these settings**, beside the Run
button, writes the current screen back out as a JSON named after the parent —
so once ⁴¹P is set up the way you want, you have a file you can reload next
week or hand to `run_fit.py`. It is the same format either way.

You do not have to load anything to start: the app opens with the ⁴⁴S
configuration already in place.

## 3. Run without a GUI

For a locked-down analysis node, or for scripting a batch of isotopes:

```bash
python run_fit.py config_44S.json test_44S.csv -o out_44S
```

| flag | meaning |
|---|---|
| `-o DIR` | output directory (default `out`) |
| `--plot` | also write `fit.png` (needs matplotlib) |
| `--subtract` | with `--plot`: subtract the background, the only readable view when the pedestal dominates — as it does in the ⁴⁴S data |
| `--set KEY=VALUE` | override one config key, repeatable, e.g. `--set tieSiblings=false --set nBins=128` |
| `--grid` | run the systematic grid scan instead of a single fit (section 4) |
| `--grid-spec FILE` | JSON overriding the scan axes or thresholds; omitted keys keep their defaults |
| `--grid-mode oat\|factorial` | one axis at a time (default) or the full product |
| `--grid-max-runs N` | hard cap in factorial mode |

It prints the half-life table, N₀, A₀, background, Cash/χ², and all the notes,
then lists the files it wrote.

## 4. The systematic grid scan

A half-life's systematic uncertainty is not a property of one fit — it is the
spread of the answer over the analysis choices you could defensibly have made
instead. `grid_scan.py` makes that spread reproducible, and both the GUI and the
CLI drive the same module, so the button and the command line give the same
budget.

```bash
python run_fit.py config_44S.json test_44S.csv -o out_44S --grid
```

In the GUI: set the fit up the way you want it, open **Systematic grid scan** in
the sidebar, press **Run grid scan** beside **Run fit**, and read the
**Systematics** tab. Your current configuration is the *reference* — everything
is measured as a deviation from it — and the reference fit also populates the
other tabs, so you do not need to run the single fit first.

The sidebar prints the run list before you press anything — total number of fits,
every axis, and the levels it resolved for each, including the ones it will skip
and why. With the default settings and a 2 – 5000 ms window that reads **23 fits**.
Switch an axis off, or type your own levels, and the list updates immediately.

**The axes**, each varied one at a time about the reference:

| axis | levels | note |
|---|---|---|
| Fit window start | tMin values; blank derives `t₀, 2.5t₀, 5t₀, 10t₀` | forces `cropToWindow` on |
| Fit window end | tMax values; blank derives `0.2, 0.4, 0.7, 1.0 × t₁` | forces `cropToWindow` on |
| Rebin factor | 1, 2, 4, 10 | coarser bins raise per-bin counts, so mild mismodelling shows up as a worse χ²/ndf |
| Bin spacing | linear, log | legitimately reads 0 on an already-binned file |
| Background | flat floated · decaying, slope from tail · decaying, free slope | |
| Daughter T½ | literature prior · fixed · free | the prior level is skipped, and reported as skipped, on a row with no σ |
| Background T½ | explicit pedestal half-lives, plus `flat`; blank derives `20, 40, 80, 200 × t₁` | **off by default** — see below |
| Efficiency cap | off, on | |

**Background treatment vs background half-life.** The two background axes ask
different questions and you normally want one of them, not both. *Background
treatment* varies the three labels — flat, slope fitted from the tail, slope
free — which is the right first pass when you have no idea what the pedestal is
doing. But once a tail fit has told you, those labels overcharge you: they still
include slopes the tail has already excluded, and on a spectrum where the
pedestal outweighs the signal that single axis can dominate the whole budget.

*Background half-life* replaces them with the numbers you choose. Each level is
a pedestal half-life in the file's own unit, held **fixed** while only its
amplitude floats (`bkgTailFit` is forced off, so nothing re-estimates the slope
behind your back), and the word `flat` is accepted as a level of its own. The
tornado then reads ΔT½ *as a function of the assumed slope*, so you can take the
systematic over the range the tail fit actually allows — its value ± 1σ — rather
than over every slope in principle. Blank levels derive `20, 40, 80, 200 × tMax`;
on a 2 – 5000 ms window that is 100, 200, 400, 1000 s, plus flat. A level at or
below zero is reported as skipped, never silently dropped.

**Changing what is scanned.** In the GUI, each axis has a checkbox, and the four
numeric axes take your own levels — type `2, 5, 10, 20` (commas or spaces) into
**tMin levels**, **tMax levels**, **Rebin levels** or **T½(bkg) levels**
(where `flat` is also a valid entry), or leave them blank to derive them. The
four categorical axes are on/off only.

For finer control use a spec file, which overlays the defaults axis by axis and
is the only way to take a *subset* of a categorical axis:

```json
{
  "mode": "oat",
  "maxRedChi2": 2.0,
  "minPValue": 0.001,
  "strictConvergence": false,
  "maxRuns": 200,
  "axes": {
    "window_start": {"enabled": true, "values": [2, 4, 8]},
    "window_end":   {"enabled": true, "values": [2000, 3500, 5000]},
    "rebin":        {"enabled": true, "values": [1, 2, 5]},
    "bin_mode":     {"enabled": false},
    "background":   {"values": ["flat_floated", "decay_tail_slope"]},
    "daughter_t12": {"values": ["prior", "free"]},
    "bkg_t12":      {"enabled": true, "values": [350000, 407000, 470000, "flat"]},
    "cap_eff":      {"enabled": false}
  }
}
```

```bash
python run_fit.py config_44S.json test_44S.csv -o out_44S --grid --grid-spec myspec.json
```

The eight axis *names* are fixed: each is a branch in `apply_level()` that knows
how to rewrite the parameters, and an unrecognised name is reported as skipped
rather than ignored. Adding a ninth axis — a free vs constrained parent, the
tail-fit range, a branching ratio within its σ — means adding a branch there and
a label in `AXIS_TITLE`.

**Why the window axes force the crop.** On a pre-binned file `tMin`/`tMax` do
nothing unless `cropToWindow` is on. Scan the window without it and every level
returns the identical fit and the axis reports a systematic of exactly zero —
which is how the first ⁴¹P window-edge number came out 0.000 ms. The grid sets
the flag itself so that cannot happen.

**Admission.** A run enters the budget only if χ²/ndf ≤ 3 and p ≥ 0.001 (both
adjustable). Rejected runs stay in the CSV with their reason — a configuration
that will not fit is itself a result. The minimiser's own convergence flag is
recorded on every row but is **not** a cut by default: on this likelihood it
fires on the near-degenerate ε/N₀ direction even at a good minimum, and it threw
out fits with χ²/ndf = 1.03 and a stable half-life. `strictConvergence` restores
it if you want it.

**Synthesis.** Per axis, ΔT½ = the largest |T½ − T½_ref| over its admitted runs;
the total is σ_syst = √(Σ ΔT½ᵢ²). An axis whose every level was rejected is
reported as a gap in the budget, not as a zero.

**`--grid-mode factorial`** runs the full Cartesian product instead. Use it to
look for interactions between choices: its per-axis numbers are main effects
(the spread of the per-level medians), not clean attributions, so its quadrature
total is an indicator and not a budget. It is capped by `maxRuns`, and whatever
the cap drops is reported.

Written to `<outdir>/systematics` from the GUI, and to `-o DIR` from the CLI:

- `systematic_grid_results.csv` — one row per grid point: T½, asymmetric profile interval, stat/syst split, χ²/ndf, Pearson p, Cash, ndf, N₀, A₀, background terms, per-daughter ε and branching, the convergence flag, admitted and why
- `systematic_budget.csv` — ΔT½ per axis plus the quadrature `TOTAL`
- `systematic_tornado.pdf` / `.png` — the publication figure: admitted envelope per axis, shaded ±stat band, dashed ±σ_syst
- `grid_results.json` — all of it, machine-readable

On ⁴⁴S with `config_44S.json`, the full scan admits 21 of 23 runs and gives
daughters 3.17, window end 2.41, background 1.92, window start 0.91, capEff
0.16, rebin 0.02, bin spacing 0 ms → **σ_syst = 4.52 ms**, i.e.
T½ = 109.394 ± 0.466 (stat) ± 4.52 (syst) ms. That is the grid's own budget
around the reference in `config_44S.json`; it is not a reproduction of the
robustness numbers in `44S_half_life_summary.md` §3, which were measured around
a different reference configuration.

## 5. The data file

Two columns, `dt,counts`, with or without a header. `dt` is the implant–decay
time difference in whatever `timeUnit` says; `counts` is the channel content.
Your ⁴⁴S file is 1666 channels of ~3 ms spanning 2–5000 ms, headerless — that
form is read correctly with no configuration.

A list of individual decay times (one column) also works: it is histogrammed
into `nBins`.

## 6. The chain table

One row per chain member. Row 0 is the parent — the nuclide whose half-life you
are measuring.

| column | meaning |
|---|---|
| **Nuclide** | label, used in output filenames and plots |
| **T½** | half-life in `timeUnit` |
| **T½ σ** | its literature uncertainty, used only when the mode is `constrained` |
| **T½ mode** | `free` / `fixed` / `constrained` — see below |
| **Efficiency** | counts detected per decay, *relative to the parent's* |
| **ε mode** | `free` or `fixed` |
| **Fed by #** | 1-based row number of the parent that feeds this one; `0` for the root |
| **Branch %** | percentage of that parent's decays that produce this member |
| **Branch σ %** | its uncertainty; non-zero makes it a constrained (fitted, pulled) parameter |

**T½ modes.** `free` floats with no constraint. `fixed` holds the typed value
exactly. `constrained` fits it but pulls it towards the typed value with a
Gaussian penalty of width `T½ σ` — which is what makes the stat/syst split
meaningful: `total² = stat² + syst²`, where syst is the part the literature
uncertainties contribute.

**Branching topology.** "Fed by #" makes a tree, not a ladder. Two members that
name the same parent are *parallel branches* — each one's population is the
classic chain solution along its own path back to the root, multiplied by the
branching fractions on that path. ⁴⁴Cl and ⁴³Cl both fed by ⁴⁴S is a β⁻ / β⁻n
split, not two generations.

**Tie siblings** (on by default) makes the dominant child of each multi-branch
parent the complement `1 − Σ(siblings)`, so branchings out of one parent always
sum to 100 % and only the minor branch is a fitted parameter. Its error is the
siblings' errors in quadrature.

## 7. Full config key reference

`config_44S.json` reproduces your last web run and is the template to copy.

| key | meaning |
|---|---|
| `timeUnit` | `"s"`, `"ms"`, `"us"`, `"ns"` — the unit of every time you type |
| `chain` | list of member dicts: `label`, `t12`, `t12Sigma`, `fixed` (`"free"`/`"fixed"`/`"prior"`), `eff`, `effFixed`, `parent`, `br`, `brSigma`. Note the JSON word for constrained is **`"prior"`**; the interface shows it as "constrained". |
| `useDemoData` | `true` runs on a synthetic spectrum — useful to check the install |
| `tMin`, `tMax` | fit window |
| `cropToWindow` | drop channels outside the window instead of just excluding them |
| `nBins`, `binMode` | binning for unbinned input; `"linear"` or `"log"` |
| `edgeRef` | `"center"` or `"left"` — what the `dt` column means for a pre-binned file |
| `rebin` | combine N adjacent channels |
| `fitBackground` | float a background component |
| `bkgRate` | starting level, counts/s |
| `bkgMode` | `"flat"` or `"decaying"` |
| `bkgT12`, `bkgT12Fixed` | background half-life and whether it floats |
| `bkgTailFit`, `bkgTailMin`, `bkgTailMax` | take the background slope from a late tail, then hold it in the full fit and float only the amplitude |
| `autoStart` | estimate starting values from the data |
| `capEff` | forbid efficiencies above 1 |
| `tieSiblings` | see above |
| `tagOutputs` | prefix output filenames with the parent label |
| `scanParent`, `scanMin`, `scanMax`, `scanPoints` | Cash-vs-T½ scan of the parent |

## 8. Output files

Written to the output directory, prefixed with the parent label when
`tagOutputs` is on:

- `<parent>_fit_results.csv` — half-lives, ±stat ±syst, profile bounds, branchings, efficiencies
- `<parent>_derived_quantities.csv` — N₀, A₀, background, ρ(lnT½, lnN₀)
- `<parent>_covariance.csv` — four stacked matrices: absolute, log-space (what the fit uses), correlation, and stat-only
- `<parent>_spectrum_fit.csv` — every curve, bin by bin, for drawing the paper figure yourself
- `<parent>_parent_scan.csv` — only when `scanParent` is true
- `results.json` — everything above, machine-readable
- `fit.png` — only with `--plot`
- the grid scan's own files are listed in section 4

**Reading the covariance.** The fit works in log space, so the log-space block
is the covariance of `ln(parameter)`; the absolute block is that multiplied by
`outer(p, p)`. For a derived quantity use `σ²(f) = wᵀ·cov·w` with `w` the
partial derivatives. Adding errors in quadrature instead assumes every
off-diagonal is zero — and for this chain ρ(ε_44Cl, br_43Cl) = +0.98, so it
is not.

## 9. Is it the same code?

Yes, and this was checked rather than assumed. `bateman_fit.py` is the web
tool's kernel cell lifted verbatim — the only changes are mechanical: the
parameters arrive as a dict instead of a JSON string, the spectrum arrives as a
path instead of a global, and the function returns its result.

The check: the original notebook cell was executed directly, `run_fit()` was
run on the same config and the same ⁴⁴S file, and the two `results.json` files
were compared field by field to 1 part in 10⁹. They are identical.

For the record, that run gives T½(⁴⁴S) = 109.394 ms ± 0.466 (constrained,
profile 108.9–109.9), with ⁴⁴Cl at 700 ms and ⁴³Cl at 3130 ms held,
N₀ = 1.9256 × 10⁷, χ²/ndf = 1.0303, p = 0.19. If you run the CLI on that file
with `config_44S.json` and get those numbers, your install is correct.

## 10. Two things worth remembering

**The ΔC scale.** With ~4.5 M counts per channel, ΔC is dominated by model
mismatch rather than by statistics. Divide any apparent significance by
√(χ²/ndf) before quoting it.

**Efficiency above 1.** The fitted ε(⁴⁴Cl) ≈ 1.39 means that member deposits
more counts per parent decay than the parent does. That is physical if its
branch feeds more than one detected transition, but it is also what a component
does when it soaks up mismodelled background. `capEff` forbids it, and how much
the fit degrades tells you which case you are in.
