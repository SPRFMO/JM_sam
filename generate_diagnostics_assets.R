library(FLCore)
library(FLSAM)
library(ggplot2)
library(patchwork)
library(RColorBrewer)

load("output_h1.RData")

if (!all(c("stk_h1", "idx_h1", "ctrl_h1", "sam_h1", "retro_h1", "loo_fits") %in% ls())) {
  stop("output_h1.RData is missing one or more required objects.")
}

dir.create("diagnostics", showWarnings = FALSE)

q_breaks <- list(
  Chile_AcousCS = 2002,
  Chile_CPUE = 2000,
  Offshore_CPUE = 2021
)

catch_fleets <- c(
  "catch N_Chile",
  "catch SC_Chile_PS",
  "catch FarNorth",
  "catch Offshore_Trawl"
)

survey_fleets <- c(
  "Chile_AcousCS_early", "Chile_AcousCS_late",
  "Chile_AcousN",
  "Chile_CPUE_early", "Chile_CPUE_late",
  "DEPM", "Peru_Acoustic", "Peru_CPUE",
  "Offshore_CPUE_early", "Offshore_CPUE_late"
)

base_fleet <- function(x) sub("_(early|late)$", "", x)
survey_fleets_base <- unique(base_fleet(survey_fleets))

max_yr <- range(stk_h1)["maxyear"]
base_yr_chr <- as.character(max_yr)

stk_h1_fit <- stk_h1 + sam_h1
res_all <- residuals(sam_h1)
ssb_df <- ssb(sam_h1)

sel.pat <- merge(f(sam_h1), fbar(sam_h1), by = "year", suffixes = c(".f", ".fbar"))
sel.pat$sel <- sel.pat$value.f / sel.pat$value.fbar
sel.pat$age <- as.numeric(as.character(sel.pat$age))
sel.pat$pentad <- sprintf("%i's", floor(sel.pat$year / 5) * 5)

.parse_jjm_ssb <- function(rep_file, max_year) {
  if (!file.exists(rep_file)) {
    warning("JJM rep file not found: ", rep_file)
    return(NULL)
  }

  lines <- readLines(rep_file, warn = FALSE)
  sec_idx <- grep("^\\$[A-Za-z_]+$", lines)
  sec_names <- sub("^\\$", "", lines[sec_idx])
  i_ssb <- sec_idx[sec_names == "SSB"]

  if (length(i_ssb) == 0) {
    warning("$SSB not found in rep file")
    return(NULL)
  }

  pos <- which(sec_idx == i_ssb)
  end <- if (pos < length(sec_idx)) sec_idx[pos + 1] - 1 else length(lines)
  dlines <- trimws(lines[(i_ssb + 1):end])
  dlines <- dlines[nchar(dlines) > 0]
  mat <- do.call(rbind, lapply(dlines, function(line) {
    as.numeric(strsplit(line, "\\s+")[[1]])
  }))

  if (ncol(mat) < 5) {
    warning("Unexpected $SSB format")
    return(NULL)
  }

  mat <- mat[mat[, 1] <= max_year, , drop = FALSE]
  data.frame(
    year = as.integer(mat[, 1]),
    value = mat[, 2],
    lbnd = mat[, 4],
    ubnd = mat[, 5]
  )
}

jjm_ssb <- .parse_jjm_ssb(
  rep_file = "../assessment/results/h1_1.14_1_R.rep",
  max_year = as.integer(max_yr)
)

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
