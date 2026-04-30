library(FLCore)
library(FLSAM)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
source("R/sam_helpers.R")

load("output_h1.RData")

if (!all(c("stk_h1", "idx_h1", "ctrl_h1", "sam_h1", "retro_h1", "loo_fits") %in% ls())) {
  stop("output_h1.RData is missing one or more required objects.")
}

dir.create("diagnostics", showWarnings = FALSE)

inject_context(prepare_diagnostics_context(stk_h1, idx_h1, sam_h1))

source("02_standard_diagnostics.R")
source("03_retro_diagnostics.R")
source("04_catchability_diagnostics.R")
source("05_productivity.R")
source("06_selectivity_blocks.R")

make_placeholder_diag <- function(path, label) {
  ext <- tolower(tools::file_ext(path))
  if (ext != "png") {
    return(invisible(NULL))
  }

  png(path, width = 1400, height = 900, res = 150)
  plot.new()
  text(0.5, 0.56, "Diagnostic figure unavailable", cex = 1.6, font = 2)
  text(0.5, 0.46, basename(path), cex = 1.1)
  text(0.5, 0.38, label, cex = 0.9)
  dev.off()
}

qmd_lines <- readLines("sam_diagnostics.qmd", warn = FALSE)
diag_refs <- regmatches(
  qmd_lines,
  gregexpr('show_diag\\("([^"]+)"', qmd_lines, perl = TRUE)
)
diag_refs <- unique(unlist(lapply(diag_refs, function(x) {
  if (length(x) == 0 || identical(x, character(0))) {
    return(character(0))
  }
  sub('show_diag\\("([^"]+)"', "\\1", x, perl = TRUE)
})))

for (ref in diag_refs) {
  path <- file.path("diagnostics", ref)
  if (!file.exists(path)) {
    make_placeholder_diag(path, "Created by generate_diagnostics_assets.R fallback")
  }
}
