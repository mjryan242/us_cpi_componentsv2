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
- `git init` this new location and push to GitHub (not yet done).
