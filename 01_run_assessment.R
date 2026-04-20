library(FLCore)
library(FLSAM)
library(ggplot2)
library(patchwork)
library(RColorBrewer)

# ============================================================
# Load FLStock and FLIndices objects from RDA files
# ============================================================

# H1 - Single stock hypothesis
load("h1_1.14_2025.rda")   # loads: stk (FLStock), idx (FLIndices)
stk_h1 <- stk
idx_h1 <- idx
rm(stk, idx)

# ============================================================
# Convert FLIndexBiomass -> FLIndex (type = "biomass")
# FLSAM.control requires FLIndex; FLIndexBiomass lacks the
# 'type' slot so must be coerced before building controls.
# ============================================================

biomass_to_FLIndex <- function(x) {
  empty <- x@index; empty[] <- NA

  r               <- x@range
  r["min"]        <- NA
  r["max"]        <- NA
  r["plusgroup"]  <- NA

  out <- new("FLIndex",
    distribution = x@distribution,
    index        = x@index,
    index.var    = empty,
    catch.n      = empty,
    catch.wt     = empty,
    effort       = empty,
    sel.pattern  = empty,
    index.q      = empty,
    name         = x@name,
    desc         = x@desc,
    range        = r
  )
  type(out) <- "biomass"
  out
}

idx_h1 <- FLIndices(lapply(idx_h1, biomass_to_FLIndex))

# ============================================================
# Split indices with time-varying q (break years from JJM h1_1.14.ctl)
#   Chile_AcousCS : break 2002
#   Chile_CPUE    : break 2000
#   Offshore_CPUE : break 2021
# Split pairs share obs.var but receive separate catchability parameters.
# biomassTreat is preserved via explicit assignment below.
# ============================================================

q_breaks <- list(
  Chile_AcousCS = 2002,
  Chile_CPUE    = 2000,
  Offshore_CPUE = 2021
)

split_idx <- function(fl, break_yr) {
  early      <- window(fl, end   = break_yr - 1)
  late       <- window(fl, start = break_yr)
  early@name <- paste0(fl@name, "_early")
  late@name  <- paste0(fl@name, "_late")
  list(early = early, late = late)
}

idx_split <- list()
for (nm in names(idx_h1)) {
  if (nm %in% names(q_breaks)) {
    sp <- split_idx(idx_h1[[nm]], q_breaks[[nm]])
    idx_split[[paste0(nm, "_early")]] <- sp$early
    idx_split[[paste0(nm, "_late")]]  <- sp$late
  } else {
    idx_split[[nm]] <- idx_h1[[nm]]
  }
}
idx_h1 <- FLIndices(idx_split)

# ============================================================
# Expand harvest.spwn / m.spwn to match fleet (area) count
# ============================================================

expand_to_fleets <- function(flq, n) {
  if (dim(flq)[5] == n) return(flq)
  new_dim          <- dim(flq);      new_dim[5]          <- n
  new_dn           <- dimnames(flq); new_dn[["area"]]    <- seq_len(n)
  out <- FLQuant(NA_real_, dim = new_dim, dimnames = new_dn, units = units(flq))
  for (a in seq_len(n)) out[,,,,a] <- flq[,,,,1]
  out
}

fix_spwn <- function(stk) {
  n <- dims(stk)$area
  stk@harvest.spwn <- expand_to_fleets(stk@harvest.spwn, n)
  stk@m.spwn       <- expand_to_fleets(stk@m.spwn, n)
  stk
}

stk_h1 <- fix_spwn(stk_h1)

# ============================================================
# FLSAM control
# ============================================================

ctrl_h1 <- FLSAM.control(stk_h1, idx_h1)

ctrl_h1@plus.group[] <- 1

ctrl_h1@states["catch N_Chile",]          <- c(1:7, rep(8,5))
ctrl_h1@states["catch SC_Chile_PS",]      <- c(1:8, rep(9,4))   + 101
ctrl_h1@states["catch FarNorth",]         <- c(1:5, rep(6,7))   + 201
ctrl_h1@states["catch Offshore_Trawl",]   <- c(1:7, rep(8,5))   + 301

ctrl_h1@f.vars["catch N_Chile",]          <- c(0,0,0,1,1,rep(2,7))
ctrl_h1@f.vars["catch SC_Chile_PS",]      <- c(0,0,0,rep(2,9)) + 101
ctrl_h1@f.vars["catch FarNorth",]         <- c(0,0,1,2,rep(3,8))   + 201
ctrl_h1@f.vars["catch Offshore_Trawl",]   <- c(rep(0,5),rep(1,7))+ 301

ctrl_h1@obs.vars["catch N_Chile",]        <- c(1, rep(2,11))
ctrl_h1@obs.vars["catch SC_Chile_PS",]    <- c(rep(1,12))       + 101
ctrl_h1@obs.vars["catch FarNorth",]       <- c(1,1,1,rep(2,9))  + 201
ctrl_h1@obs.vars["catch Offshore_Trawl",] <- c(rep(1,12))       + 301

# Split pairs share the same obs.var index; each gets its own catchability.
# Index order after splitting (10 surveys):
#   Chile_AcousCS_early, Chile_AcousCS_late, Chile_AcousN,
#   Chile_CPUE_early, Chile_CPUE_late, DEPM,
#   Peru_Acoustic, Peru_CPUE, Offshore_CPUE_early, Offshore_CPUE_late
ctrl_h1@obs.vars["Chile_AcousCS_early", 1] <- 401
ctrl_h1@obs.vars["Chile_AcousCS_late",  1] <- 401   # shared with early
ctrl_h1@obs.vars["Chile_AcousN",        1] <- 402
ctrl_h1@obs.vars["Chile_CPUE_early",    1] <- 403
ctrl_h1@obs.vars["Chile_CPUE_late",     1] <- 403   # shared with early
ctrl_h1@obs.vars["DEPM",               1]  <- 404
ctrl_h1@obs.vars["Peru_Acoustic",       1] <- 405
ctrl_h1@obs.vars["Peru_CPUE",          1]  <- 406
ctrl_h1@obs.vars["Offshore_CPUE_early", 1] <- 407
ctrl_h1@obs.vars["Offshore_CPUE_late",  1] <- 407   # shared with early

# biomassTreat: 0=SSB (DEPM), 5=total biomass (acoustics + CPUEs)
ctrl_h1@biomassTreat[which(names(ctrl_h1@fleets) == "Chile_AcousCS_early")] <- 5
ctrl_h1@biomassTreat[which(names(ctrl_h1@fleets) == "Chile_AcousCS_late")]  <- 5
ctrl_h1@biomassTreat[which(names(ctrl_h1@fleets) == "Chile_AcousN")]        <- 5
ctrl_h1@biomassTreat[which(names(ctrl_h1@fleets) == "Chile_CPUE_early")]    <- 5
ctrl_h1@biomassTreat[which(names(ctrl_h1@fleets) == "Chile_CPUE_late")]     <- 5
ctrl_h1@biomassTreat[which(names(ctrl_h1@fleets) == "DEPM")]                <- 0
ctrl_h1@biomassTreat[which(names(ctrl_h1@fleets) == "Peru_Acoustic")]       <- 5
ctrl_h1@biomassTreat[which(names(ctrl_h1@fleets) == "Peru_CPUE")]           <- 5
ctrl_h1@biomassTreat[which(names(ctrl_h1@fleets) == "Offshore_CPUE_early")] <- 5
ctrl_h1@biomassTreat[which(names(ctrl_h1@fleets) == "Offshore_CPUE_late")]  <- 5

ctrl_h1 <- update(ctrl_h1)
ctrl_h1@residuals <- TRUE

# ============================================================
# Helper: rebuild ctrl for any subset of indices
# Used by both the retrospective loop and leave-one-out fits.
# MUST be defined before first use (retro loop and LOO).
# ============================================================
build_ctrl <- function(stk, idx_sub) {
  ctrl <- FLSAM.control(stk, idx_sub)
  ctrl@plus.group[] <- 1

  ctrl@states["catch N_Chile",]          <- c(1:7, rep(8,5))
  ctrl@states["catch SC_Chile_PS",]      <- c(1:8, rep(9,4))  + 101
  ctrl@states["catch FarNorth",]         <- c(1:5, rep(6,7))  + 201
  ctrl@states["catch Offshore_Trawl",]   <- c(1:7, rep(8,5))  + 301

  ctrl@f.vars["catch N_Chile",]          <- c(0,0,0,1,1,rep(2,7))
  ctrl@f.vars["catch SC_Chile_PS",]      <- c(0,0,0,rep(2,9))     + 101
  ctrl@f.vars["catch FarNorth",]         <- c(0,0,1,2,rep(3,8))   + 201
  ctrl@f.vars["catch Offshore_Trawl",]   <- c(rep(0,5),rep(1,7))  + 301

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
               Chile_CPUE_early=5,     Chile_CPUE_late=5,
               DEPM=0, Peru_Acoustic=5, Peru_CPUE=5,
               Offshore_CPUE_early=5,  Offshore_CPUE_late=5)

  for (nm in names(idx_sub)) {
    ctrl@obs.vars[nm, 1]                      <- obs_map[nm]
    ctrl@biomassTreat[which(names(ctrl@fleets) == nm)] <- bio_map[nm]
  }

  ctrl <- update(ctrl)
  ctrl@residuals <- FALSE
  ctrl
}

# ============================================================
# Alternative helper: FLSAM defaults for catch f.vars / obs.vars
# Catch fleets get one shared f.var and one shared obs.var per
# fleet (all ages same code) — the FLSAM.control() defaults.
# States, survey obs.vars and biomassTreat are identical to
# build_ctrl.  Use for model comparisons / sensitivity runs.
# ============================================================
build_ctrl_alt <- function(stk, idx_sub) {
  ctrl <- FLSAM.control(stk, idx_sub)
  ctrl@plus.group[] <- 1

  # States: same age-grouping as main model
  ctrl@states["catch N_Chile",]          <- c(1:7, rep(8,5))
  ctrl@states["catch SC_Chile_PS",]      <- c(1:8, rep(9,4))  + 101
  ctrl@states["catch FarNorth",]         <- c(1:5, rep(6,7))  + 201
  ctrl@states["catch Offshore_Trawl",]   <- c(1:7, rep(8,5))  + 301

  # f.vars and obs.vars for catch fleets: FLSAM defaults retained
  # (FLSAM.control assigns one shared parameter per fleet across
  # all ages — no override needed here)

  # Survey obs.vars: shared pairs for split series (same as build_ctrl)
  obs_map <- c(Chile_AcousCS_early=401, Chile_AcousCS_late=401,
               Chile_AcousN=402,
               Chile_CPUE_early=403,   Chile_CPUE_late=403,
               DEPM=404, Peru_Acoustic=405, Peru_CPUE=406,
               Offshore_CPUE_early=407, Offshore_CPUE_late=407)
  bio_map <- c(Chile_AcousCS_early=5,  Chile_AcousCS_late=5,
               Chile_AcousN=5,
               Chile_CPUE_early=5,     Chile_CPUE_late=5,
               DEPM=0, Peru_Acoustic=5, Peru_CPUE=5,
               Offshore_CPUE_early=5,  Offshore_CPUE_late=5)

  for (nm in names(idx_sub)) {
    ctrl@obs.vars[nm, 1]                                    <- obs_map[nm]
    ctrl@biomassTreat[which(names(ctrl@fleets) == nm)] <- bio_map[nm]
  }

  ctrl <- update(ctrl)
  ctrl@residuals <- FALSE
  ctrl
}

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

  idx_peel            <- FLIndices(lapply(idx_h1[in_range], window, end = yr_end))
  ctrl_peel           <- build_ctrl_alt(stk_peel, idx_peel)
  ctrl_peel@residuals <- TRUE

  retro_h1[[as.character(yr_end)]] <- tryCatch(
    FLSAM(stk_peel, idx_peel, ctrl_peel),
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

loo_fits <- lapply(names(loo_groups), function(grp) {
  drop     <- loo_groups[[grp]]
  idx_loo  <- FLIndices(idx_h1[!names(idx_h1) %in% drop])
  ctrl_loo <- build_ctrl(stk_h1, idx_loo)
  cat(sprintf("LOO: dropping %s ...\n", grp))
  tryCatch(
    FLSAM(stk_h1, idx_loo, ctrl_loo),
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

catch_fleets  <- c("catch N_Chile", "catch SC_Chile_PS",
                   "catch FarNorth", "catch Offshore_Trawl")
survey_fleets <- c("Chile_AcousCS_early", "Chile_AcousCS_late",
                   "Chile_AcousN",
                   "Chile_CPUE_early", "Chile_CPUE_late",
                   "DEPM", "Peru_Acoustic", "Peru_CPUE",
                   "Offshore_CPUE_early", "Offshore_CPUE_late")

# Strip _early/_late suffixes for display (split pairs merge into one panel)
base_fleet         <- function(x) sub("_(early|late)$", "", x)
survey_fleets_base <- unique(base_fleet(survey_fleets))   # 7 base names

base_yr_chr <- as.character(max_yr)

# Common model outputs used across multiple diagnostic scripts
stk_h1_fit  <- stk_h1 + sam_h1
res_all     <- residuals(sam_h1)
ssb_df      <- ssb(sam_h1)

# Selectivity pattern: needed in 02 and 06
sel.pat        <- merge(f(sam_h1), fbar(sam_h1), by = "year", suffixes = c(".f", ".fbar"))
sel.pat$sel    <- sel.pat$value.f / sel.pat$value.fbar
sel.pat$age    <- as.numeric(as.character(sel.pat$age))
sel.pat$pentad <- sprintf("%i's", floor(sel.pat$year / 5) * 5)

# ============================================================
# Parse JJM h1_1.14 SSB for model comparison overlay
# jjm_ssb  : data.frame(year, value, lbnd, ubnd)  — or NULL if file absent
# ============================================================
.parse_jjm_ssb <- function(rep_file, max_year) {
  if (!file.exists(rep_file)) {
    warning("JJM rep file not found: ", rep_file)
    return(NULL)
  }
  lines     <- readLines(rep_file, warn = FALSE)
  sec_idx   <- grep("^\\$[A-Za-z_]+$", lines)
  sec_names <- sub("^\\$", "", lines[sec_idx])

  i_ssb <- sec_idx[sec_names == "SSB"]
  if (length(i_ssb) == 0) { warning("$SSB not found in rep file"); return(NULL) }

  pos  <- which(sec_idx == i_ssb)
  end  <- if (pos < length(sec_idx)) sec_idx[pos + 1] - 1 else length(lines)
  dlines <- trimws(lines[(i_ssb + 1):end])
  dlines <- dlines[nchar(dlines) > 0]
  mat    <- do.call(rbind, lapply(dlines,
              function(l) as.numeric(strsplit(l, "\\s+")[[1]])))

  if (ncol(mat) < 5) { warning("Unexpected $SSB format"); return(NULL) }
  # columns: year | estimate | SE | lower_95 | upper_95
  mat <- mat[mat[, 1] <= max_year, , drop = FALSE]
  data.frame(year  = as.integer(mat[, 1]),
             value = mat[, 2],
             lbnd  = mat[, 4],
             ubnd  = mat[, 5])
}

jjm_ssb <- .parse_jjm_ssb(
  rep_file  = "../assessment/results/h1_1.14_1_R.rep",
  max_year  = as.integer(range(stk_h1)["maxyear"])
)
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

Sys.setenv(QUARTO_PATH = "C:/Program Files/Quarto/bin/quarto.exe")
quarto::quarto_render("sam_diagnostics.qmd")
