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

The repository was split from the `sam/` directory of `JM_SCW_prep` so that the
SAM workflow, rendered diagnostics, and supporting data products can be developed
and versioned independently.

## Contents

- `01_run_assessment.R` runs the H1 SAM assessment workflow.
- `01_stock.R` prepares stock and index inputs.
- `01b_retro_loo_parallel.R` runs retrospective and leave-one-out analyses.
- `02_standard_diagnostics.R` through `06_selectivity_blocks.R` generate the
  main diagnostic summaries.
- `sam_diagnostics.qmd` is the Quarto diagnostics report.
- `diagnostics/` contains generated diagnostic figures and tables.

## Report

The rendered diagnostics report is available through GitHub Pages:

<https://sprfmo.github.io/JM_sam/>

The source report file is `sam_diagnostics.qmd`, and the rendered HTML is
available in both `sam_diagnostics.html` and `index.html`. To refresh the site,
run the assessment and diagnostics scripts, render `sam_diagnostics.qmd` with
Quarto, and copy the rendered report to `index.html`.
