# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A research project reproducing and extending the paper *"Comparing two
periods of high inflation: the Great Inflation and the 2020s inflation"*
(`paper/main (1).tex`). It studies whether the dependence structure among US
CPI subcomponents, wages, and inflation expectations changes when inflation
is high, using a pseudo-quantile $R^2$ connectedness estimator, and tests
whether that connectedness feeds into the Fed's Taylor-rule reaction function
(H5/H6, via a Markov-switching regime classifier). See `README.md` for the
full run order — this file covers things the README doesn't.

## Commands

All R scripts use **bare relative paths** (`source("use_R2Q.R")`,
`read.csv("MICH.csv")`, `dir.create("results_...")`) and must be run **from
this repo root**, not from `paper/` or any subfolder — except
`paper/make_tables_hyp.R`, which is itself invoked from the repo root but
writes into `paper/tables/`.

```sh
Rscript us_cpi_6cpi_driver.R          # rebuild any single driver
Rscript paper/make_tables_hyp.R       # rebuild h1/h2/h3/h5/h6/h6net from cached results_*/ CSVs
```

There is no test suite, linter, or build system — this is a research
pipeline, not a package. "Correctness" is checked by comparing a script's
output CSV/tex against the numbers already cited in `paper/main (1).tex`'s
prose (e.g. `results_6cpi/TCI_by_quantile.csv` should reproduce `11.808`,
`35.216` etc. — see `README.md`'s H1 section). When re-running a script,
diff its output against the existing file in `results_*/` or `paper/tables/`
before assuming a change in the numbers reflects your edit rather than a
FRED vintage update or a bug.

Compile the paper **from the repo root**, not from `paper/` — every path in
`main (1).tex` (`\input`, `\includegraphics`, `\addbibresource`) is written
relative to the repo root (e.g. `\input{paper/tables/h1}`), because that is
how Overleaf resolves relative paths when this repo is synced via its
GitHub integration (project root = repo root, regardless of which subfolder
the main file lives in). See `README.md`'s compile section for the exact
`-output-directory=paper` invocation. Note the literal space and `(1)` in
the filename — always quote it.

## Architecture

**Estimator layer** (never edit without understanding the whole project
depends on it): `use_R2Q.R` sources `ConnectednessTools/R/*.R` (a vendored,
trimmed copy of S.J.H. Shahzad's package — just the 4 files actually used,
not a full clone) and exposes `R2ConnectednessQ2()`, which wraps
`R2QConnectedness` with a PSD-targeted Schäfer–Strimmer shrinkage override.
Every driver script calls this same function, so the whole project shares
one estimator. `R2Q_lasso_dir.R` is a separate, LASSO-based pseudo-$R^2$
metric used only for the H5 goods→services sequencing robustness (not on the
same scale as the Genizi/`R2ConnectednessQ2` metric — don't compare them
directly). `quarterly_utils.R` provides the monthly→quarterly aggregation
convention used everywhere quarterly results appear: **aggregate the price
index to a quarterly period-average, then take the % change** (never
aggregate monthly % changes directly); survey rates (MICH) take the
quarterly mean of the level, no differencing.

**Driver scripts → cached CSVs → table assembler.** Each `us_cpi_*.R` driver
computes one hypothesis's connectedness numbers and writes a CSV to a
`results_*/` folder (e.g. `us_cpi_6cpi_driver.R` → `results_6cpi/`). Three
scripts (`us_cpi_6cpi_matrices.R`/`_q.R`, `us_cpi_seq_q.R`) instead write
`.tex` fragments directly into `paper/tables/`. Everything else funnels
through `paper/make_tables_hyp.R`, which reads the `results_*/` CSVs and
formats `paper/tables/{h1,h2,h3,h5,h6,h6net}.tex` (H5 and H6 are partly
recomputed live inside this script rather than cached — see its source).
This two-layer structure means: to change a number in the paper, find the
driver that owns that `results_*/` folder, not `make_tables_hyp.R` itself.

**The "no oil" pivot.** The paper's headline system is six BLS CPI
subcomponents with core-services measured *ex-shelter* (`W <- 0.60` shelter
weight) — no crude oil. Scripts named `*_6cpi_*` are the current, correct
variant. If you ever see or are asked to reuse a `*_6group_oil_*` or
7-series script elsewhere in the wider project this repo was extracted from,
that is the **superseded, pre-pivot** spec — do not mix its outputs with the
current tables without translating the series list first (see how
`us_cpi_bootstrap_r2q.R` was adapted, in git history / `DECISIONS.md`).

**Markov-switching regime scripts** (`us_cpi_regime_markov_switch_MR.R`,
`us_cpi_regime_ms_quarterly.R`, `us_cpi_msq_exp12_addcol.R`) all write into a
single flat `regime_tex/` folder at the repo root, using filename prefixes
(`ms_*` = monthly connectedness regime, `msq_*` = quarterly connectedness
regime) to disambiguate — they do **not** write directly to
`paper/tables/`. The paper expects them split into two subfolders,
`paper/tables/regime_tex/` (ms_*) and `paper/tables/regime_texQtr/` (msq_*);
copying them there is a manual "glue" step (see `README.md` step 8), not
something the scripts do themselves.

**Known gotcha:** `us_cpi_regime_ms_quarterly.R` writes a 3-column
`msq_exp12_120.tex` (short Michigan-survey sample, monthly MICH from 1978
aggregated to quarters). The paper's actual table has 4 columns, including
an "Extended (1969-)" column using the native quarterly Michigan survey
(`MICH_QTR.csv`, from 1960). `us_cpi_msq_exp12_addcol.R` must run **after**
`us_cpi_regime_ms_quarterly.R` to overwrite the short version with the
correct one — order matters, see `README.md` step 7.

**Bootstrap.** `us_cpi_bootstrap_r2q.R` computes fixed-design quantile-wild
bootstrap CIs for the H1 headline. `B` is currently `10` (smoke-test only,
seconds to run) — bump to `1000` (or set `BOOT_B` env var) before citing any
CI in the paper. It is not yet wired into `make_tables_hyp.R`/`h1.tex`.

## Project practice: DECISIONS.md and LOG.md

This project keeps two running project-management files at the repo root,
**in addition to** normal code comments and commit messages. Keep both up to
date as you work — do not treat them as optional or backfill them only at
the end of a session.

- **`DECISIONS.md`** — one entry per methodological choice (estimator
  variant, lag length, quantile set, sample split, bootstrap design, series
  inclusion/exclusion, etc.). For each: what was chosen, what alternative(s)
  were rejected and why, and a confidence level (e.g. High/Medium/Low, or a
  short caveat) reflecting how settled the choice is. Add an entry whenever
  a design choice is made or revisited — not just for new work, but when an
  existing choice is reconsidered or reaffirmed.
- **`LOG.md`** — a plain-language, chronological narrative a coauthor could
  read to follow the project's progress without reading code or diffs: what
  was tried, what was found, what changed and why, dead ends, and open
  questions. Write it for a collaborator seeing the update for the first
  time, not as a terse commit log.

If either file doesn't exist yet, create it the first time you make a
decision or complete a unit of work worth logging.
