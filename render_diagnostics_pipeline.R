library(FLCore)
library(FLSAM)
library(quarto)
source("generate_diagnostics_assets.R")

quarto::quarto_render("sam_diagnostics.qmd")
