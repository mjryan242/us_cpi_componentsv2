# DECISIONS.md

One entry per methodological choice: what was chosen, what alternatives were
rejected and why, and a confidence level. Add an entry whenever a design
choice is made or reconsidered.

---

## Connectedness estimator: pseudo-quantile $R^2$ (Genizi/nearPD), via `use_R2Q.R`

**Chosen:** `R2ConnectednessQ2()` — quantile-correlation input to the Genizi
(1993) order-invariant $R^2$ decomposition, with Schäfer–Strimmer shrinkage
(PSD-targeted override) and `Matrix::nearPD` repair of the quantile-
correlation matrix when it isn't PSD (mainly at $\tau=0.9$).

**Rejected:**
- **Ando, Greenwood-Nimmo & Shin (2022)-style FEVD-based quantile
  connectedness.** Propagates non-zero-mean quantile-regression residuals
  through a generalised FEVD as if they were zero-mean shocks; inflates
  measured connectedness in both tails (documented in Shahzad et al.).
- **Eigenvalue-clipping estimator** (the project's earlier, pre-`use_R2Q.R`
  approach, used in the now-deprecated `us_cpi_bootstrap.R`). Repairs
  non-PSD matrices per-block by clipping negative eigenvalues rather than a
  whole-matrix `nearPD` repair. Agrees with the current estimator at
  $\tau\le0.7$; diverges only in the extreme tail. Superseded, not
  currently used anywhere in the live pipeline.

**Confidence:** High. This is the estimator the whole current paper is built
on; changing it would require re-running every driver script.

---

## Headline system: six BLS CPI subcomponents, no oil

**Chosen:** Food, Gasoline, Household energy, Core goods, Shelter, Core
services ex-shelter ($w=0.60$ shelter weight). `us_cpi_6cpi_driver.R` and
siblings.

**Rejected:** An earlier 7-series system (the same six groups + WTI crude
oil). Scripts for this variant (`*_6group_oil_*`) still exist in the wider
project this repo was extracted from but are superseded — dropping oil
simplified the story and removed a series whose variation is dominated by a
market outside the CPI system entirely (see the appendix note on gasoline's
weak within-system connectedness).

**Confidence:** High — this is the paper's current stated headline spec.

---

## H5/H6 policy regime classifier: Markov-switching (2-state, AR(1)), not a fixed threshold

**Chosen:** Fit `MSwM::msmFit()` to the connectedness series itself
(AR(1), 2 states, all parameters switching); classify each period into the
higher-mean-connectedness state via smoothed probabilities; interact that
state indicator with the inflation gap in an inertial Taylor rule.

**Rejected:** A Hansen-style fixed threshold on the connectedness level
(grid search over candidate $\gamma^*$, minimise split SSR). Both were built
and estimated; the Markov-switching version is preferred because the regime
is treated as latent and probabilistic rather than a single hard cut, and
because the fixed-threshold version's regimes are extremely persistent
(rolling-window connectedness is itself near-integrated), which risked
conflating "regime" with a long historical block.

**Confidence:** Medium. Neither classifier produces state-dependent
responsiveness that is robust across window length (120 vs 60 month) or
across the matched/mismatched inflation-measure choice — see the paper's own
"suggestive, not robust" framing in the H5/H6 results section.

---

## Expectations-connectedness regime: total (bidirectional) index, not directional

**Chosen:** The regime classifier for H6 uses the *total* connectedness
index of the {expectations, inflation} system (contemporaneous + lagged, both
directions summed), matching the H5 (six-CPI) construction.

**Rejected:** A directional variant using only the lagged
inflation→expectations spillover (the "adaptive channel" element of the
connectedness matrix, verified against a synthetic DGP to confirm the
row/column convention). Built and run (`ms_dir_exp12_*`,`ms_dir_expmm_*`,
never wired into the paper) — gives a *noisier*, less stable regime
interaction than the total index, so it does not strengthen H6 and was left
as an unused robustness check rather than promoted into the paper.

**Confidence:** Medium — a total, bidirectional index is the more natural
match to "spiral risk" (expectations and inflation feeding each other), but
this was a judgment call, not a decisive empirical result either way.

---

## Inflation measure for the expectations system: matched (12-month), not mismatched (m/m)

**Chosen:** Pair the Michigan 1-year-ahead expectation with 12-month
(year-over-year) CPI inflation ("matched" — same horizon as the survey
question).

**Rejected:** Pairing with month-on-month CPI inflation ("mismatched" — a
horizon mismatch between a 12-month-ahead expectation and a 1-month
realised rate). Both versions are kept in the paper (Panel A / Panel B of
the expectations table) precisely because the result is sensitive to this
choice: the one significant H6 result (matched, $\tau=0.7$, 120-month) does
not survive under the mismatched measure.

**Confidence:** High on the *preference* for the matched measure
(conceptually correct horizon alignment); Medium on what that implies for
H6, since the result doesn't survive other robustness cuts either
(60-month window, directional-only spillover, extended 1960- sample).

---

## Two-month lag limit → quarterly-frequency robustness

**Chosen:** Main spec uses `nlag=2` on **monthly** data (spillovers
propagate over 2 months); as a robustness check, connectedness is
re-estimated on **quarterly**-aggregated data with `nlag=2` (spillovers now
propagate over ~6 months), classifying whole quarters into regimes and
broadcasting the regime label to each quarter's three months for the
(still-monthly) Taylor-rule estimation.

**Rejected:** Simply raising `nlag` on monthly data. Additional lags rapidly
exhaust the degrees of freedom in the short rolling windows (worst in the
six-variable system and the 60-month window) and degrade the estimates —
this is why `nlag=2` is the ceiling at monthly frequency, and why the
quarterly-aggregation route was used to extend the effective horizon instead
of just increasing `nlag` directly.

**Confidence:** High on the mechanism (why nlag can't just increase);
Medium on the substantive finding it produces — quarterly connectedness
*strengthens* H5 (six of eight $\tau\in\{0.5,0.7\}$ interaction estimates
turn significant, vs. two of twelve at monthly frequency), which is a
genuine, only partially understood sensitivity to the propagation horizon,
not yet fully reconciled with the monthly null.

---

## Bootstrap design: fixed-design quantile-wild, not recursive

**Chosen:** Fixed-design wild bootstrap (Caporin–Bonaccolto–Shahzad
residual construction): bootstrap values are `fitted + wild residual`
using the *actual* lags, no feedback loop.

**Rejected:** Recursive QVAR simulation (feed bootstrap draws back in as
lags for the next period). Explodes at $\tau=0.9$ — the tail quantile-VAR is
non-stationary, so a recursive path diverges (~$10^{24}$) and every draw
collapses to a constant.

**Confidence:** High — this is a numerical-stability issue, not a
close judgment call.

---

## `us_cpi_bootstrap_r2q.R`: adapted from 7-series+oil to six-CPI, `B` left at 10

**Chosen:** The inherited bootstrap script targeted the old 7-series+oil
system; it was edited to use the same six no-oil columns as
`us_cpi_6cpi_driver.R` (oil merge removed), and its hardcoded `setwd()` to
the old project path was removed. Verified: a `B=10` smoke-test run
reproduces the *exact* base TCI values in `results_6cpi/TCI_by_quantile.csv`
before any bootstrap resampling is applied, confirming the series-list edit
is correct.

**Not yet done:** `B` is still `10` (a ~7-second smoke test). A real run
needs `B<-1000` (or `BOOT_B=1000`), and the resulting CIs are not yet wired
into `make_tables_hyp.R`/`paper/tables/h1.tex` — `h1.tex`'s own notes still
say "bootstrap confidence intervals for the six-CPI system are to be added."

**Confidence:** High that the adaptation is correct (verified against the
cached base numbers); the actual CIs are simply not computed yet at
publication-quality `B`.

---

## Paper's `\input`/`\includegraphics`/`\addbibresource` paths: root-relative, not paper/-relative

**Chosen:** All relative paths in `main (1).tex` are now written relative to
the **repo root** (e.g. `\input{paper/tables/h1}`,
`\includegraphics{paper/TCI_decomp.pdf}`,
`\addbibresource{paper/references.bib}`). Compile from the repo root with
`pdflatex -output-directory=paper "paper/main (1).tex"` and
`biber --output-directory=paper "main (1)"`.

**Rejected:** Paths relative to `paper/` (the original form, e.g.
`\input{tables/h1}`), compiled by `cd paper && pdflatex "main (1).tex"`.
This worked locally but breaks under Overleaf's GitHub-sync integration:
Overleaf always resolves relative paths against the **project root** (the
top of the synced repo), regardless of which subfolder the main `.tex` file
lives in — it has no equivalent of "cd into paper/ first." Verified the fix
by recompiling from the repo root locally (not `cd paper`) and confirming an
identical PDF (242,816 vs 242,810 bytes; same content, only the resolution
mode differs) — this is the same resolution behaviour Overleaf will use.

**Confidence:** High — this is a documented Overleaf behaviour (project
root is always the sync root), not a judgment call, and the fix was
verified by reproducing Overleaf's exact resolution mode locally.
