# JM_sam

This repository contains State-Space Assessment Model (SAM) analyses for the
SPRFMO Chilean jack mackerel (*Trachurus murphyi*) assessment work.

The SAM model is used here as an independent diagnostic tool for the H1 jack
mackerel assessment configuration. It is not intended to replace the Joint Jack
Mackerel Model (JJM). Instead, the analyses support review of operating-model
conditioning choices for the jack mackerel MSE process, including:

- candidate selectivity block structures;
- productivity and reference point diagnostics;
- index uncertainty, catchability, and fleet-conflict diagnostics;
- retrospective and leave-one-out behavior.

The repository was split from the `sam/` directory of `JM_SCW_prep` and is now
the canonical home for the SAM workflow, rendered diagnostics, and supporting
data products. The former `JM_SCW_prep/sam` folder should be treated as retired
after migration.

## Contents

- `01_run_assessment.R` runs the H1 SAM assessment workflow.
- `01_stock.R` prepares stock and index inputs.
- `01b_retro_loo_parallel.R` runs retrospective and leave-one-out analyses.
- `02_standard_diagnostics.R` through `06_selectivity_blocks.R` generate the
  main diagnostic summaries.
- `sam_diagnostics.qmd` is the Quarto diagnostics report.
- `diagnostics/` contains generated diagnostic figures and tables.

## Reproducible setup

This repository does not currently ship an `renv.lock`, so install the runtime
dependencies explicitly before running the workflow. The analysis uses R,
Quarto, and the FLR/SAM package stack:

```r
install.packages(c(
  "ggplot2",
  "patchwork",
  "RColorBrewer",
  "quarto",
  "R.utils"
))

install.packages("remotes")
remotes::install_github(c(
  "flr/FLCore",
  "flr/FLSAM"
))
```

Quarto must also be available on `PATH` for command-line rendering. The active
render path is:

```sh
Rscript render_from_output.R
```

To refit the assessment before rendering, run:

```sh
Rscript 01_run_assessment.R
```

## Report

The rendered diagnostics report is available through GitHub Pages:

<https://sprfmo.github.io/JM_sam/>

The source report file is `sam_diagnostics.qmd`, and the rendered HTML is
available as `sam_diagnostics.html`. The site landing page is generated from
`index.qmd` and links to the active diagnostics report plus a frozen
`JM_SCW_prep/sam` render from commit `f5a21f0` for migration provenance.
