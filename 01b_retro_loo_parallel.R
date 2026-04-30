# ============================================================
# 01b_retro_loo_parallel.R
# Parallel versions of the retrospective and LOO analyses.
#
# Run AFTER 01_run_assessment.R has produced output_h1.RData
# (which contains stk_h1, idx_h1, ctrl_h1, sam_h1).
# Results overwrite retro_h1 and loo_fits in output_h1.RData
# so that downstream diagnostic scripts (02_, 03_, …) work
# without modification.
#
# Uses mclapply() (fork-based, Linux/macOS only).
# Workers inherit the parent environment — no cluster setup,
# no clusterExport, no package reloading needed.
# ============================================================

library(FLCore)
library(FLSAM)
library(parallel)
library(R.utils)   # withTimeout()
source("R/sam_helpers.R")

skip_flsam_osa_residuals()

# ---- load main model outputs --------------------------------
load("output_h1.RData")   # stk_h1, idx_h1, ctrl_h1, sam_h1

# ---- LOO groups (must match 01_run_assessment.R) ------------
loo_groups <- list(
  Chile_AcousCS = c("Chile_AcousCS_early", "Chile_AcousCS_late"),
  Chile_AcousN  = "Chile_AcousN",
  Chile_CPUE    = c("Chile_CPUE_early",    "Chile_CPUE_late"),
  DEPM          = "DEPM",
  Peru_Acoustic = "Peru_Acoustic",
  Peru_CPUE     = "Peru_CPUE",
  Offshore_CPUE = c("Offshore_CPUE_early", "Offshore_CPUE_late")
)

# ---- timeouts -----------------------------------------------
# Per-job wall-clock limit in seconds.  Adjust to your machine.
retro_timeout_sec <- 2 * 3600   # 2 hours per retro peel
loo_timeout_sec   <- 1 * 3600   # 1 hour  per LOO fit

# ---- worker count -------------------------------------------
# Cap at n_retro (10) for retro; LOO only has 7 groups.
# Leave one core free for the OS.
n_workers <- max(1, detectCores() - 3)
cat(sprintf("Available workers: %i (using up to %i)\n",
            n_workers, n_workers))

# ============================================================
# Parallel retrospective
# ============================================================
max_yr  <- range(stk_h1)["maxyear"]
n_retro <- 10
peels   <- seq_len(n_retro)

cat(sprintf("Running %i retro peels in parallel (mc.cores = %i) ...\n",
            n_retro, min(n_workers, n_retro)))
t_retro <- system.time({
  retro_list <- mclapply(peels, function(peel) {
    yr_end   <- max_yr - peel
    stk_peel <- window(stk_h1, end = yr_end)

    in_range <- sapply(idx_h1, function(fl) range(fl)["minyear"] <= yr_end - 1)
    dropped  <- names(idx_h1)[!in_range]
    if (length(dropped) > 0)
      message(sprintf("Peel %i (->%i): dropping %s", peel, yr_end,
                      paste(dropped, collapse = ", ")))
    else
      message(sprintf("Peel %i (->%i)", peel, yr_end))

    idx_peel            <- FLIndices(lapply(idx_h1[in_range], window, end = yr_end))
    ctrl_peel           <- build_ctrl_alt(stk_peel, idx_peel)
    ctrl_peel@residuals <- TRUE

    tryCatch(
      withTimeout(FLSAM(stk_peel, idx_peel, ctrl_peel, newtonsteps = 0),
                  timeout = retro_timeout_sec, onTimeout = "error"),
      error = function(e) { message("FAILED peel ", peel, ": ", e$message); NULL }
    )
  }, mc.cores = min(n_workers, n_retro))
})
cat(sprintf("Retro done in %.1f min\n", t_retro["elapsed"] / 60))

retro_h1 <- setNames(retro_list, as.character(max_yr - peels))
retro_h1 <- Filter(Negate(is.null), retro_h1)
retro_h1[[as.character(max_yr)]] <- sam_h1   # append full model

# ============================================================
# Parallel LOO
# ============================================================
cat(sprintf("Running %i LOO fits in parallel (mc.cores = %i) ...\n",
            length(loo_groups), min(n_workers, length(loo_groups))))
t_loo <- system.time({
  loo_list <- mclapply(names(loo_groups), function(grp) {
    drop     <- loo_groups[[grp]]
    idx_loo  <- FLIndices(idx_h1[!names(idx_h1) %in% drop])
    ctrl_loo <- build_ctrl_alt(stk_h1, idx_loo)
    ctrl_loo@residuals <- TRUE
    message(sprintf("LOO: dropping %s", grp))
    tryCatch(
      withTimeout(FLSAM(stk_h1, idx_loo, ctrl_loo, newtonsteps = 0),
                  timeout = loo_timeout_sec, onTimeout = "error"),
      error = function(e) { message("LOO FAILED ", grp, ": ", e$message); NULL }
    )
  }, mc.cores = min(n_workers, length(loo_groups)))
})
cat(sprintf("LOO done in %.1f min\n", t_loo["elapsed"] / 60))

loo_fits <- setNames(loo_list, names(loo_groups))

# ============================================================
# Save — overwrites retro_h1 and loo_fits in output_h1.RData
# ============================================================
save(stk_h1, idx_h1, ctrl_h1, sam_h1, retro_h1, loo_fits,
     file = "output_h1.RData")
cat("Saved output_h1.RData\n")
