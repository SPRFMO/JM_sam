library(FLCore)
library(FLSAM)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
source("R/sam_helpers.R")

skip_flsam_osa_residuals()

# ============================================================
# Load FLStock and FLIndices objects from RDA files
# ============================================================

# H1 - Single stock hypothesis
load("h1_1.14_2025.rda")   # loads: stk (FLStock), idx (FLIndices)
stk_h1 <- stk
idx_h1 <- idx
rm(stk, idx)

idx_h1 <- FLIndices(lapply(idx_h1, biomass_to_FLIndex))

# ============================================================
# Alternative input set:
#   - drop Chile_AcousCS, DEPM, Peru_Acoustic entirely
#   - retain Chile_AcousN only from 2000 onward
# Apply before q-splitting so downstream objects inherit the change.
# ============================================================

drop_indices <- c("Chile_AcousCS", "DEPM", "Peru_Acoustic")
idx_h1 <- subset_indices(idx_h1, setdiff(names(idx_h1), drop_indices))
idx_h1[["Chile_AcousN"]] <- window(idx_h1[["Chile_AcousN"]], start = 2000)

# ============================================================
# Split indices with time-varying q (break years from JJM h1_1.14.ctl)
#   Chile_AcousCS : break 2002
#   Chile_CPUE    : break 2000
#   Offshore_CPUE : break 2021
# Split pairs share obs.var but receive separate catchability parameters.
# biomassTreat is preserved via explicit assignment below.
# ============================================================

idx_h1 <- split_indices_by_q_break(idx_h1)

stk_h1 <- fix_spwn(stk_h1)

# ============================================================
# FLSAM control
# ============================================================

ctrl_h1 <- build_ctrl_alt(stk_h1, idx_h1, residuals = TRUE)

# ============================================================
# Fit FLSAM model
# ============================================================

sam_h1 <- FLSAM(stk_h1, idx_h1, ctrl_h1)

# ============================================================
# Retrospective analysis (sequential loop)
# retro() cannot handle sub-series that become empty mid-peel.
# For each peel: window stock + indices, detect sub-series that
# would have < 2 years of data and drop them, then fit.
# ============================================================
ctrl_h1@residuals <- TRUE     # enable residuals for diagnostic plots
max_yr  <- range(stk_h1)["maxyear"]
n_retro <- 10

retro_h1 <- setNames(vector("list", n_retro),
                     as.character(max_yr - seq_len(n_retro)))

for (peel in seq_len(n_retro)) {
  yr_end   <- max_yr - peel
  stk_peel <- window(stk_h1, end = yr_end)

  in_range <- sapply(idx_h1, function(fl) range(fl)["minyear"] <= yr_end - 1)
  dropped  <- names(idx_h1)[!in_range]

  if (length(dropped) > 0) {
    cat(sprintf("Retro peel %i (->%i): dropping %s\n",
                peel, yr_end, paste(dropped, collapse = ", ")))
  } else {
    cat(sprintf("Retro peel %i (->%i)\n", peel, yr_end))
  }

  keep_peel           <- names(idx_h1)[in_range]
  idx_peel            <- FLIndices(setNames(lapply(keep_peel, function(nm) window(idx_h1[[nm]], end = yr_end)), keep_peel))
  ctrl_peel           <- build_ctrl_alt(stk_peel, idx_peel)
  ctrl_peel@residuals <- TRUE

  retro_h1[[as.character(yr_end)]] <- tryCatch(
    FLSAM(stk_peel, idx_peel, ctrl_peel, newtonsteps = 0),
    error = function(e) { message("  FAILED: ", e$message); NULL }
  )
}
retro_h1 <- Filter(Negate(is.null), retro_h1)
retro_h1[[as.character(max_yr)]] <- sam_h1   # append full model

# ============================================================
# Leave-one-out (LOO) analysis
# Drop one survey group at a time; split pairs dropped together.
# ============================================================

loo_groups <- list(
  Chile_AcousCS = c("Chile_AcousCS_early", "Chile_AcousCS_late"),
  Chile_AcousN  = "Chile_AcousN",
  Chile_CPUE    = c("Chile_CPUE_early",    "Chile_CPUE_late"),
  DEPM          = "DEPM",
  Peru_Acoustic = "Peru_Acoustic",
  Peru_CPUE     = "Peru_CPUE",
  Offshore_CPUE = c("Offshore_CPUE_early", "Offshore_CPUE_late")
)
loo_groups <- loo_groups[vapply(loo_groups, function(x) any(x %in% names(idx_h1)), logical(1))]

loo_fits <- lapply(names(loo_groups), function(grp) {
  drop     <- loo_groups[[grp]]
  idx_loo  <- subset_indices(idx_h1, names(idx_h1)[!names(idx_h1) %in% drop])
  ctrl_loo <- build_ctrl(stk_h1, idx_loo)
  cat(sprintf("LOO: dropping %s ...\n", grp))
  tryCatch(
    FLSAM(stk_h1, idx_loo, ctrl_loo, newtonsteps = 0),
    error = function(e) { message("  FAILED: ", e$message); NULL }
  )
})
names(loo_fits) <- names(loo_groups)

# ============================================================
# Save model outputs and input data
# ============================================================
save(stk_h1,      # input FLStock
     idx_h1,      # input FLIndices (with early/late splits)
     ctrl_h1,     # FLSAM control
     sam_h1,      # full assessment fit
     retro_h1,    # retrospective fits keyed by final year (incl. full model)
     loo_fits,    # leave-one-out fits (one per survey group)
     file = "output_h1.RData")
cat("Saved output_h1.RData\n")

# ============================================================
# Shared setup for all diagnostic scripts
# ============================================================

dir.create("diagnostics", showWarnings = FALSE)

inject_context(prepare_diagnostics_context(stk_h1, idx_h1, sam_h1))
if (!is.null(jjm_ssb))
  cat(sprintf("JJM SSB loaded: %i years (%i\u2013%i)\n",
              nrow(jjm_ssb), min(jjm_ssb$year), max(jjm_ssb$year)))

# ============================================================
# Source diagnostic scripts
# ============================================================

source("02_standard_diagnostics.R")
source("03_retro_diagnostics.R")
source("04_catchability_diagnostics.R")
source("05_productivity.R")
source("06_selectivity_blocks.R")

# ============================================================
# Render Quarto report
# ============================================================

render_quarto_report("sam_diagnostics.qmd")
