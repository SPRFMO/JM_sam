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

# ---- load main model outputs --------------------------------
load("output_h1.RData")   # stk_h1, idx_h1, ctrl_h1, sam_h1

# ---- build_ctrl (must be identical to 01_run_assessment.R) --
build_ctrl <- function(stk, idx_sub) {
  ctrl <- FLSAM.control(stk, idx_sub)
  ctrl@plus.group[] <- 1

  ctrl@states["catch N_Chile",]          <- c(1:7, rep(8,5))
  ctrl@states["catch SC_Chile_PS",]      <- c(1:8, rep(9,4))  + 101
  ctrl@states["catch FarNorth",]         <- c(1:5, rep(6,7))  + 201
  ctrl@states["catch Offshore_Trawl",]   <- c(1:7, rep(8,5))  + 301

  ctrl@f.vars["catch N_Chile",]          <- c(0,0,0,1,1,rep(2,7))
  ctrl@f.vars["catch SC_Chile_PS",]      <- c(0,0,0,rep(2,9))    + 101
  ctrl@f.vars["catch FarNorth",]         <- c(0,0,1,2,rep(3,8))  + 201
  ctrl@f.vars["catch Offshore_Trawl",]   <- c(rep(0,5),rep(1,7)) + 301

  ctrl@obs.vars["catch N_Chile",]        <- c(1, rep(2,11))
  ctrl@obs.vars["catch SC_Chile_PS",]    <- c(rep(1,12))       + 101
  ctrl@obs.vars["catch FarNorth",]       <- c(1,1,1,rep(2,9))  + 201
  ctrl@obs.vars["catch Offshore_Trawl",] <- c(rep(1,12))       + 301

  obs_map <- c(Chile_AcousCS_early=401, Chile_AcousCS_late=401,
               Chile_AcousN=402,
               Chile_CPUE_early=403,   Chile_CPUE_late=403,
               DEPM=404, Peru_Acoustic=405, Peru_CPUE=406,
               Offshore_CPUE_early=407, Offshore_CPUE_late=407)
  bio_map <- c(Chile_AcousCS_early=5,  Chile_AcousCS_late=5,
               Chile_AcousN=5,
               Chile_CPUE_early=2,     Chile_CPUE_late=2,
               DEPM=0, Peru_Acoustic=5, Peru_CPUE=2,
               Offshore_CPUE_early=2,  Offshore_CPUE_late=2)

  for (nm in names(idx_sub)) {
    ctrl@obs.vars[nm, 1]                                    <- obs_map[nm]
    ctrl@biomassTreat[which(names(ctrl@fleets) == nm) - 4] <- bio_map[nm]
  }

  ctrl <- update(ctrl)
  ctrl@residuals <- FALSE
  ctrl
}

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

# ---- worker count -------------------------------------------
# Cap at n_retro (10) for retro; LOO only has 7 groups.
# Leave one core free for the OS.
n_workers <- max(1, detectCores() - 1)
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
    ctrl_peel           <- build_ctrl(stk_peel, idx_peel)
    ctrl_peel@residuals <- TRUE

    tryCatch(
      FLSAM(stk_peel, idx_peel, ctrl_peel),
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
    ctrl_loo <- build_ctrl(stk_h1, idx_loo)
    message(sprintf("LOO: dropping %s", grp))
    tryCatch(
      FLSAM(stk_h1, idx_loo, ctrl_loo),
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
