# Replication package — *Comparing two periods of high inflation*

Code and local data inputs needed to reproduce every table and figure in
`paper/main (1).tex` (six BLS CPI subcomponents, no oil; pseudo-quantile $R^2$
connectedness; H1–H6 including the Markov-switching connectedness-regime
Taylor rule). Nothing here was rewritten — every script is copied unmodified
from the working project; this README only documents run order and two small
"glue" steps (copying script output into the folders the paper's `\input`
paths expect).

## Requirements
- **R** (developed on 4.5.2) with packages: `quantmod, zoo, xts, quantreg,
  MASS, corpcor, Matrix, ConnectednessApproach, sandwich, lmtest, MSwM, dplyr,
  lubridate, tidyr` (scripts auto-install anything missing via
  `requireNamespace`/`install.packages`).
- **Internet** — most series (six CPI subcomponents, headline CPI, FEDFUNDS,
  UNRATE, NROU) are downloaded live from FRED via `quantmod::getSymbols`.
  Local inputs shipped here: `MICH.csv` (Michigan 1-yr expectations, monthly),
  `MICH_QTR.csv` (same, quarterly mean, from 1960 Q1), `AHETPI.csv` (average
  hourly earnings), `WTISPLC.csv` (WTI crude spot), `WuXiaShadowRate.csv`
  (Wu–Xia shadow policy rate).
- **pdfLaTeX + biber** to compile the paper.

## Directory layout
Run everything **from this folder** (`us_cpi_componentsv2/`) — every script
uses bare relative paths (`source("use_R2Q.R")`, `read.csv("MICH.csv")`,
`dir.create("results_...")`), exactly as in the original project, so the
working directory must be this one.

```
us_cpi_componentsv2/
├── use_R2Q.R, quarterly_utils.R, R2Q_lasso_dir.R   # shared estimator/helpers
├── ConnectednessTools/R/*.R                        # Shahzad et al. estimator (4 files)
├── us_cpi_6cpi_driver.R / _q.R                     # H1: six-CPI TCI by quantile (monthly/quarterly)
├── us_cpi_6cpi_era.R / _q.R                        # H3: episode marginal decomposition
├── us_cpi_wage_spiral.R / _q.R                     # H2: wage-CPI connectedness
├── us_cpi_2var_exp_q.R                             # H6net: expectations NET directional
├── us_cpi_6cpi_matrices.R / _q.R                   # Appendix: full pairwise matrices
├── us_cpi_seq_q.R                                  # H5 (quarterly): goods->services sequencing
├── us_cpi_regime_markov_switch_MR.R                # H5/H6 MS Taylor rule (monthly, 120/60mo)
├── us_cpi_regime_ms_quarterly.R                    # H5/H6 MS Taylor rule (quarterly connectedness)
├── us_cpi_msq_exp12_addcol.R                       # rebuilds msq_exp12_120 w/ the "Extended" (1960-) column
├── MICH.csv, MICH_QTR.csv, AHETPI.csv, WTISPLC.csv, WuXiaShadowRate.csv
└── paper/
    ├── main (1).tex, references.bib
    ├── make_tables_hyp.R                           # builds h1,h2,h3,h5,h6,h6net from the CSVs above
    └── tables/                                     # <- all .tex fragments land here
```

## How to run

```sh
R="/path/to/Rscript"

# 1. Six-CPI headline (H1) + episode mechanism (H3), monthly & quarterly
$R us_cpi_6cpi_driver.R        # -> results_6cpi/, plots_6cpi/TCI_decomp.pdf
$R us_cpi_6cpi_driver_q.R      # -> results_6cpi_q/
$R us_cpi_6cpi_era.R           # -> results_6cpi_era/
$R us_cpi_6cpi_era_q.R         # -> results_6cpi_era_q/

# 2. Wage-CPI (H2)
$R us_cpi_wage_spiral.R        # -> results_wage_spiral/
$R us_cpi_wage_spiral_q.R      # -> results_wage_spiral_q/

# 3. Expectations NET directional (H6net)
$R us_cpi_2var_exp_q.R         # -> results_2var_exp_q/

# 4. Appendix connectedness matrices (writes paper/tables/conn_*.tex directly)
$R us_cpi_6cpi_matrices.R
$R us_cpi_6cpi_matrices_q.R    # (companion; not currently \input by main.tex)

# 5. Goods->services sequencing, quarterly (H5_q; writes paper/tables/h5_q.tex directly)
$R us_cpi_seq_q.R

# 6. Assemble h1, h2, h3, h5 (monthly), h6, h6net from the CSVs written in 1-3
$R paper/make_tables_hyp.R     # -> paper/tables/{h1,h2,h3,h5,h6,h6net}.tex

# 6b. (optional) Bootstrap CIs for the H1 headline TCI-by-quantile result.
#     B=10 by default (smoke test) -- bump to B<-1000 for a real run; see
#     the "Bootstrap" section below. Not currently \input by the paper.
$R us_cpi_bootstrap_r2q.R      # -> results_boot_r2q/, plots_boot_r2q/

# 7. Markov-switching connectedness-regime Taylor rule (H5/H6 policy section)
#    Run ms_quarterly BEFORE the addcol fix, since addcol overwrites its
#    msq_exp12_120.tex with the correct 4-column "Extended" (1960-) version.
$R us_cpi_regime_markov_switch_MR.R   # -> regime_tex/ms_{six,expM,expX}_{120,60}.tex
$R us_cpi_regime_ms_quarterly.R       # -> regime_tex/msq_{six,exp12,expmm}_{120,60}.tex
$R us_cpi_msq_exp12_addcol.R          # -> regime_tex/msq_exp12_120.tex (OVERWRITES with 4-col version)

# 8. Glue: main.tex expects the MS-regime fragments split into two
#    subfolders under paper/tables/ (monthly "ms_*" vs quarterly "msq_*"),
#    but the scripts above write everything flat into regime_tex/.
mkdir -p paper/tables/regime_tex paper/tables/regime_texQtr
cp regime_tex/ms_*.tex  paper/tables/regime_tex/
cp regime_tex/msq_*.tex paper/tables/regime_texQtr/

# 9. Glue: the one figure (\includegraphics{TCI_decomp.pdf}, no path prefix)
#    must sit next to main.tex.
cp plots_6cpi/TCI_decomp.pdf paper/

# 10. Compile -- FROM THE REPO ROOT, not from paper/. main.tex's \input,
#     \includegraphics and \addbibresource paths are all written relative to
#     the repo root (e.g. \input{paper/tables/h1}), because that is how
#     Overleaf resolves relative paths when this repo is synced via GitHub
#     integration (project root = repo root, regardless of which subfolder
#     the main .tex file lives in). Compiling from inside paper/ (the old
#     "cd paper && pdflatex ..." recipe) will NOT find these files.
PDFLATEX="/path/to/pdflatex"; BIBER="/path/to/biber"
"$PDFLATEX" -interaction=nonstopmode -output-directory=paper "paper/main (1).tex"
"$BIBER" --output-directory=paper "main (1)"
"$PDFLATEX" -interaction=nonstopmode -output-directory=paper "paper/main (1).tex"
"$PDFLATEX" -interaction=nonstopmode -output-directory=paper "paper/main (1).tex"
```

**Overleaf (GitHub sync):** once this repo is linked as an Overleaf project
via GitHub integration, set `paper/main (1).tex` as the project's "Main
document" — no other configuration is needed; Overleaf's compiler already
resolves relative paths from the project root, matching the recipe above.

## Bootstrap (adapted to the six-CPI system; smoke-tested at B=10)
`us_cpi_bootstrap_r2q.R` computes fixed-design quantile-wild bootstrap CIs
(Caporin–Bonaccolto–Shahzad family) for the headline TCI-by-quantile result,
using the same estimator as everything else (`use_R2Q.R`). It originally
targeted the old 7-series (six BLS groups + WTI oil) headline; it has been
adapted here to the **current six-CPI no-oil system** (same six columns as
`us_cpi_6cpi_driver.R`: Food, Gasoline, HHEnergy, CoreGoods, Shelter,
CoreServ\_xS, oil merge removed), and the hardcoded `setwd()` to the old
project path was dropped (run it from this folder like everything else).

`h1.tex`'s notes still say "bootstrap confidence intervals for the six-CPI
system are to be added" — this script is that piece; it just hasn't been
wired into `make_tables_hyp.R`/`h1.tex` yet.

**`B` is currently set to 10** (a ~7-second smoke test, not a real CI) — the
comment above the config block says exactly where to change it back. A
B=10 run reproduces the exact base TCI values in `results_6cpi/TCI_by_quantile.csv`
(5.55 / 7.48 / 11.81 / 19.70 / 35.22 across $\tau=0.1,\dots,0.9$), confirming the
adaptation targets the right system; set `B <- 1000` (or `BOOT_B=1000`) for
publication-quality CIs (the header estimates several minutes at that B with
the default `NWORK = detectCores()-2` parallel workers).

(A separate, fully deprecated script, `us_cpi_bootstrap.R`, uses a different
eigenvalue-clipping estimator entirely and was not copied in.)

## Notes / known issues in the source scripts (left as-is, not "fixed")
- `paper/make_tables_hyp.R` requires the package `ConnectednessApproach` even
  though it is not actually called by this script's logic (a harmless leftover
  from an earlier estimator; `install.packages("ConnectednessApproach")` will
  satisfy it).
- **`msq_exp12_120.tex` has two versions.** `us_cpi_regime_ms_quarterly.R`
  writes a 3-column version (short Michigan sample, monthly survey from 1978
  aggregated to quarters). The paper's actual table (4 columns: `Extended
  (1969) | Earliest (1987) | 1993 | 2003`, with the `$\dagger$` footnote
  referencing the native quarterly Michigan survey back to 1960) comes from
  `us_cpi_msq_exp12_addcol.R`, which must run *after* step 7 above to
  overwrite the short version.
- `paper/references.bib` was missing 7 citation keys used in `main (1).tex`
  (`Powell2021`, `Powell2022JH`, `hamilton1989new`, `wu2016measuring`,
  `ClaridaGaliGertler2000`, `goodfriend1997new`, `bryan1994measuring`); they
  have been appended (see the bottom of the file, clearly marked) rather than
  regenerated, since these are motivating/methodological citations, not model
  output.
- `tables/regime_tex/{ms,msq}_{six,expM,expX,exp12,expmm}_60.tex` and
  `us_cpi_6cpi_matrices_q.R`'s `conn_*_q.tex` outputs are produced by these
  same scripts but are **not** currently `\input` by the paper (they support
  the 60-month and quarterly-matrix robustness narrative in prose only,
  without an inline table). No action needed; harmless if generated.
