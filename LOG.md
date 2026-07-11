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
