# LOG.md

A plain-language, chronological narrative of the project's progress — for a
coauthor catching up, not a commit log. Newest entries at the bottom.

---

### 2026-07-11 — Repo extracted, made reproducible, moved out of OneDrive

The working paper (`main (1).tex`) had accumulated a lot of scaffolding
across a much larger, messier project directory — dozens of superseded
driver-script variants (with-oil vs no-oil, 4/5/6/7-series, LASSO vs Genizi,
monthly vs quarterly), most no longer relevant to the current spec. The goal
this round was to pull out *just* what's needed to reproduce every table and
figure the paper currently uses, verify it actually works, and set the
project up somewhere safe to put under version control.

**What we found tracing the dependencies:** the H1–H3 (six-CPI broadening,
episode mechanism), H2 (wage-price), and H6/H6net (expectations) tables all
funnel through `paper/make_tables_hyp.R`, which reads cached CSVs from
`results_*/` folders written by individual driver scripts. The appendix
connectedness matrices and the H5 quarterly sequencing table are self-
contained — their driver scripts write LaTeX fragments directly. The Markov-
switching regime tables (H5/H6 policy section) are newer work and don't
follow that pattern at all — all three regime scripts dump everything into
one flat `regime_tex/` folder, and nothing was actually copying that output
into the subfolders (`paper/tables/regime_tex/`, `.../regime_texQtr/`) the
paper's `\input` paths expect. In fact, checking the timestamps, it looks
like **the paper had never successfully compiled** with the regime section
in place before this — that folder split had just never been done.

Along the way we caught two real bugs, not just missing plumbing:
1. The deployed `msq_exp12_120.tex` (quarterly expectations regime table)
   was the wrong version — a 3-column short-sample table, when the paper's
   prose and footnote describe a 4-column table with an "Extended (1969-)"
   column from the longer quarterly Michigan survey. The correct version
   comes from a different, later script (`us_cpi_msq_exp12_addcol.R`) that
   has to run *after* the main quarterly regime script to overwrite it.
2. `references.bib` was missing 7 citation keys the paper actually cites
   (Powell speeches, Hamilton 1989, Wu-Xia, Clarida-Gali-Gertler, Goodfriend-
   King, Bryan-Cecchetti) — the paper would have compiled with `[?]` markers
   or biber errors without these.

We ran the whole pipeline end-to-end (not just assembled files and hoped):
six-CPI H1 numbers reproduced byte-identical to the existing cached results;
all ~84 Markov-switching fits across the regime scripts converged; the paper
compiled to a real 22-page PDF. Needed one small LaTeX preamble fix
(`natbib=true` on the `biblatex` package) — the paper mixes `\citet`/`\citep`
commands with `style=apa`, which doesn't define `\citet` on its own. Left
three pre-existing `\label`/`\ref` mismatches in the paper's own prose
unfixed and flagged instead, since those are content decisions, not build
plumbing.

Separately, adapted the connectedness bootstrap script
(`us_cpi_bootstrap_r2q.R`), which had been written for the old 7-series
system, to the current six-CPI system, and fixed a hardcoded `setwd()` that
pointed at the old project path (would have broken outright here). Smoke-
tested at `B=10` — the base TCI values before resampling matched the cached
H1 numbers exactly, confirming the series-list swap was done correctly.
`B` is still set to 10 (seconds to run); a real run needs `B=1000` and isn't
wired into the paper's tables yet.

Finally, moved the whole thing from
`OneDrive - The University of Waikato\Desktop\connectedness\...` (where
continuous OneDrive sync of a `.git/` folder risks corruption) to a plain
local path, `C:\Users\mryan\us_cpi_componentsv2`, ahead of setting it up
under git and pushing to GitHub.

**Open questions for the next session:**
- Wire the bootstrap CIs into `h1.tex` once `B` is bumped to 1000.
- Decide whether to fix the `\label`/`\ref` mismatches in `main (1).tex`.
- The quarterly-connectedness robustness (H5) shows a genuinely different
  result from the monthly spec (see `DECISIONS.md`) — worth reconciling in
  the paper's discussion rather than leaving as an unexplained sensitivity.

### 2026-07-11 (later) — Repo pushed to GitHub; fixed paths for Overleaf sync

Local-only git steps run first (`init`, `add -A`, `commit`, `branch -M
main`), then the user created the GitHub repo and pushed. That closed out
the last open item from the earlier entry.

Next: the user plans to link this GitHub repo to Overleaf (Claude Code →
GitHub → Overleaf workflow) and asked whether `main (1).tex`'s `\input`
paths would need editing first. They were right to ask — Overleaf resolves
relative paths against the *project root* it syncs, not against the folder
containing the main `.tex` file, so the existing paths (written assuming
`cd paper && pdflatex ...` locally) would have broken as soon as Overleaf
tried to compile. Rewrote every `\input`, `\includegraphics`, and
`\addbibresource` path in `main (1).tex` to be root-relative (prefixed with
`paper/`), and updated the local compile recipe in `README.md`/`CLAUDE.md`
to match (`pdflatex -output-directory=paper "paper/main (1).tex"` run from
the repo root, not `cd paper` first). Verified by actually recompiling from
the repo root rather than assuming the fix was right — biber needed the
same treatment (`--output-directory=paper`, run from root) since it hit the
identical resolution issue on the first attempt (`Cannot find
'paper/references.bib'` when run from inside `paper/`). Final PDF is
byte-count-identical in substance to the pre-fix version (242,816 vs
242,810 bytes), confirming only the path-resolution mode changed.

**Open questions for the next session:**
- Wire the bootstrap CIs into `h1.tex` once `B` is bumped to 1000.
- Decide whether to fix the `\label`/`\ref` mismatches in `main (1).tex`.
- The quarterly-connectedness robustness (H5) shows a genuinely different
  result from the monthly spec (see `DECISIONS.md`) — worth reconciling in
  the paper's discussion rather than leaving as an unexplained sensitivity.
- Link the repo to Overleaf via GitHub integration and do a live test
  compile there (the local proxy is verified, but Overleaf itself hasn't
  been tried yet).
- Commit and push this round of changes (paths fix + updated
  README/CLAUDE.md/DECISIONS.md/LOG.md).

**Also caught while committing:** `paper/references.bib` had been reduced
from 20 entries to 5 sometime between the initial commit and now (not by
me — I only touched `main (1).tex` this round; git history confirms the
initial commit had the full file). Flagged it before pushing rather than
assuming; the user confirmed it's intentional (they're manually reconciling
their own bibliography) and asked for a reference sheet instead of a
restore. Wrote `paper/MISSING_REFERENCES.md` listing the 11 cited-but-
missing keys with ready-to-paste BibTeX (pulled from the initial commit,
not reconstructed from memory), plus a flag on an existing typo'd key
(`ClaridaGaliGertler200`, missing the trailing `0`, so it won't match the
paper's `\citep{ClaridaGaliGertler2000}` even though a near-identical entry
exists). Did not touch `references.bib` or recompile — left both for the
user to finish. Until the bib is reconciled, the committed
`paper/main (1).pdf` reflects a compile from *before* the trim (all
citations resolved) and will go stale relative to source if recompiled now.

### 2026-07-11 (later still) — Introduction strengthened for JMCB

Pulled the user's `revised introduction` commit (`3e79d64`) from GitHub —
they'd been editing the intro in Overleaf, which fixed an earlier duplicated-
H1/H2 bug I'd flagged. Reviewed the intro against JMCB expectations and made
a round of edits: (i) reordered so the closest-prior-work + three
contributions paragraph sits *before* the six hypotheses rather than after,
merging it with the "what we do" paragraph to remove a duplicated method
explanation; (ii) wrote the previously-missing results preview (two
paragraphs, real numbers: 11.8->35.2 broadening, lagged 23.3 vs contemp 11.9,
wage-price 25.2 vs 3.6, expectations NET +29.1 -> -20.2/-6.7, the hedged
H5/H6 line, and the causal caveat); (iii) fixed the H6 run-on and four
typos; (iv) inserted a `\todo{[lit review here]}` placeholder in the
related-literature slot with a commented list of suggested cites by theme
(breadth/core, 2021-22 surge, Great Inflation, expectations anchoring,
wage-price, disaggregated inflation/transmission, connectedness methods).
The working title also changed (in the user's Overleaf) to "When prices move
together: connectedness among CPI components, wages, and expectations across
the inflation distribution and its consequences" — flagged that the "and its
consequences" tail leans on the weak H5/H6 result.

Committed the source edits and pushed. Did **not** commit the recompiled PDF:
with `references.bib` still trimmed to 5 entries (the user's intentional,
in-progress bibliography reconciliation), a fresh compile shows `[?]` for the
11 missing keys, so pushing that PDF would be misleading. The repo keeps the
last good PDF; Overleaf regenerates on its own. Open: bib reconciliation
(see `paper/MISSING_REFERENCES.md`), write the lit review, and the abstract
is still empty.

### 2026-07-13 — Methodology section reviewed vs the companion methods paper

Read the private companion methods paper (Shahzad-Ryan-Gabauer, in
`additional_materials/` — gitignored, never pushed) and reviewed the CPI
paper's methodology against it. The core (Genizi decomposition, quantile-
correlation substitution, stacking, TCI = TCI^C + TCI^L) matches the companion
paper exactly. Implemented the gaps: (i) added the Choi-Shin quantile-
correlation formula; (ii) added a paragraph describing the PSD regularisation
actually used (Schafer-Strimmer PSD-targeted shrinkage + Higham nearPD repair),
with schafer2005 and higham2002 added to references.bib; (iii) fixed the TCI^L
formula to sum over j != k (own-lags excluded), matching both the "off-diagonal"
prose and the code (`od()` zeroes the diagonal); (iv) unified FROM/TO/NET to a
single index k (row-k sum vs column-k sum), which also resolves the draft's two
bracketed queries -- both answered "correct as written": TCI = avg FROM = avg TO
(row-total = column-total / K, holds despite asymmetry), and NET_k = TO_k - FROM_k
is properly indexed on the same series k.

**Did NOT add the requested "own-lags dropped in the rolling estimates"
footnote.** Verified in `use_R2Q.R` (lines 63-66) that the headline estimator
`R2ConnectednessQ2` explicitly *refuses* to drop own-lags (`stop("drop_own_lags
is not supported...")`) and instead secures PSD via shrinkage + nearPD. Dropping
own-lags is the *companion paper's* rolling-window approach, not this paper's, so
the footnote would misdescribe the method. Flagged for the user to decide:
keep the accurate shrinkage/nearPD description (done), or switch the code to the
own-lag-dropping (clipping) estimator to match the companion paper (a code
change, not a footnote).

Also added `additional_materials/` to `.gitignore`.

### 2026-07-13 (later) — Matched-horizon robustness for H4/h6; decided on mismatched throughout

Ran a new standalone check (`us_cpi_h6_matched_robustness.R`) re-estimating table h6's
{MICH, CPI} directional decomposition with YoY (4-quarter) CPI inflation — matched to
the Michigan survey's 1-year horizon — alongside the QoQ version the paper currently
uses ("mismatched"). The QoQ panel exactly reproduced the published h6 numbers first,
as a check that the new script's data construction was right, before trusting the YoY
comparison. Result: the net lead/follow sign pattern (expectations lead in Core,
follow in Great Inflation/COVID) survives the horizon change, but the headline GI
$\tau=0.9$ lagged-dominance claim (inflation→expectations $50.7$ vs reverse $27.5$)
does not — under YoY the two are roughly equal ($37.1$ vs $39.8$). The COVID
$\tau=0.9$ matched cell is also degenerate (shrinkage saturates with only ~25
quarters of post-2020 data).

Based on this, the user decided to adopt the **mismatched (short-horizon) measure as
the paper's headline throughout** — wherever realised inflation is paired with
expectations, use m/m (monthly) or QoQ (quarterly), not YoY/12-month — reversing an
earlier documented preference for the matched measure. This is mainly an internal-
consistency move (every other system in the paper — six CPI subcomponents, wages —
already uses short-horizon changes) plus a mechanical argument (overlapping YoY
changes share most of their underlying months with neighbouring observations, which
inflates lagged connectedness and blurs the contemporaneous/lagged and directional
attributions H3 and H4 rely on). Logged in `DECISIONS.md`, with the earlier "matched"
entry marked superseded rather than deleted.

**Consequences still to implement in the paper (not yet done):** the H5/H6 policy
tables' headline expectations panel needs to switch from the matched variant
(`ms_expM_*`, `msq_exp12_*`) to the mismatched one (`ms_expX_*`, `msq_expmm_*`); the
"we prefer the annual measure" policy-section prose needs rewriting since it now says
the opposite of the working choice; and under the mismatched measure the one
significant H6-policy result (matched, $\tau=0.7$, 120-month) is no longer headline —
H6-policy reads as a clean null throughout, which should be reflected in the
hypothesis framing. The h6/h6net (H4) prose should also gain a footnote noting the
horizon mismatch, citing the new matched-YoY results as robustness rather than as an
equally-weighted alternative panel.

**Not done in this pass:** touching `main (1).tex` itself. The user pushed a large
"post first readthrough" revision to the paper directly from Overleaf in parallel
(commit `b095edc`, pulled in) that already rewrote the methodology section
substantially — including their own version of the shrinkage/PSD paragraph, so an
earlier locally-drafted "smaller"→"larger" correction to that paragraph was discarded
as superseded rather than reapplied. The mismatched-measure consequences above are
documented but still need to be implemented in the paper text and table selection.

### 2026-07-17 — GitHub issues #1/#2/#4/#5/#6 + all bracketed draft notes resolved

One coordinated pass over the open GitHub issues (except #3, bootstrapping,
deferred by decision) and every bracketed `[...]` note the user left in the
draft. New computations first, then a single paper edit, then a full compile.

**New analysis.** (i) The episode-mechanism (era) scripts now run at tau=0.7
as well as 0.9, monthly and quarterly, and also save the calm-core base
levels (issue #6). The tau=0.7 results strengthen H3: quarterly, the Great
Inflation adds 25.6 points of connectedness (73% lagged) while COVID adds
essentially nothing. (ii) A monthly analogue of the expectations table
(issue #1): monthly MICH starts Jan 1978 (verified — the draft's "[check?]"),
so the Great Inflation appears only as a caveated 1978-82 sample; the
monthly pattern matches the quarterly one (expectations lead in the calm
core, NET +24.1 at tau=0.9; follow during COVID). (iii) The "taus and
levels" robustness (issue #4), both designs: matching episodes at the same
inflation level rather than the same tau. Key number: at the common 8.7%
annualised level, the GI marginal is +12.2 vs COVID's +1.1 — the episode
difference is not a level artefact, which supports the anchoring reading.
(iv) A one-in/one-out marginal construction for the expectations system,
matching h3's design (the draft asked whether this was worth adding).

**Table changes.** h3 now has tau sub-panels; h6 gained a NET column (the
separate h6net table is retired, folded in) and its samples reordered to
Full, Core, Great Inflation, COVID; the replication guard passed (existing
h6 cells reproduce exactly). Three new appendix tables: h6_monthly,
h6_marginal, tau_levels.

**Paper edits.** endfloat added (issue #5) — all floats now collect at the
back with "[Table X about here.]" markers (19 in the compiled PDF). The
shelter-separation sentence expanded with verified references (issue #2):
Genesove 2003 REStat; Gallin-Verbrugge 2019 JET; Adams-Loewenstein-Montag-
Verbrugge 2024 AER:Insights; Bolhuis-Cramer-Summers 2022 Rev. Finance. The
H2 monthly/quarterly direction flip is rationalised as a lag-window effect
(wage-setting adjusts at multi-month horizons; verified in the CSVs: monthly
CPI->Wages 16.2 vs 11.5, quarterly flips to 23.8 vs 35.1). The H3 quarterly
caveat now cites the actual base levels (76.1 overall / 51.4 lagged at
tau=0.9 — the "[is this true??]" was true) AND the 26-quarter COVID window.
The H4 footnote on annual-vs-quarterly inflation now explains the overlap-
induced moving-average component properly. The policy section's Panel B
(msq_expmm_120, mismatched quarterly) is wired in with corrected panel
headings and notes, and the superseded matched-measure commented blocks were
deleted.

**A correction to our own records:** DECISIONS.md previously claimed the
mismatched-measure H6 policy interaction was "insignificant throughout".
Reading the actual fragments, that is wrong: the monthly system is null, but
the quarterly mismatched system has delta significant at tau=0.5 (1993***,
2003***) and tau=0.7 (2003**). The paper now frames H6 as "suggestive at
best", parallel to H5, and DECISIONS.md carries the correction.

Compile: clean (3 passes + biber), no undefined references, 36 pages with
the float pages at the back. The only bracketed note remaining sits inside
an already-commented-out draft line.
