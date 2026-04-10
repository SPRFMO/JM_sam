library(FLCore)
library(FLSAM)

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
  # FLIndexBiomass has catch.n/catch.wt with age 1:12, but index with age
  # "all". dims() picks the largest age dim, so it returns min=1 instead of
  # NA. FLSAM2SAM checks is.na(dims$min) to detect biomass indices, so all
  # slots must share the same (1-row, age="all") dimensions as index.
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
# These slots are stock-level (1 area) but FLSAM2SAM loops over
# all areas; replicate the single layer to each fleet.
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

stk_h1           <- fix_spwn(stk_h1)

# ============================================================
# FLSAM control objects
# ============================================================

# --- H1 ---
ctrl_h1 <- FLSAM.control(stk_h1, idx_h1)

ctrl_h1@plus.group[] <- 1

ctrl_h1@states["catch N_Chile",] 		        <- c(1:7,rep(8,5))
ctrl_h1@states["catch SC_Chile_PS",] 	      <- c(1:8,rep(9,4)) 		+ 101
ctrl_h1@states["catch FarNorth",] 		      <- c(1:5,rep(6,7)) 		+ 201
ctrl_h1@states["catch Offshore_Trawl",]     <- c(1:7,rep(8,5)) 		+ 301

ctrl_h1@f.vars["catch N_Chile",] 		        <- c(0,0,1,1,rep(2,8))
ctrl_h1@f.vars["catch SC_Chile_PS",]	      <- c(0,1,1,1,rep(2,8)) 	+ 101
ctrl_h1@f.vars["catch FarNorth",] 		      <- c(0,1,2,rep(3,9)) 	+ 201
ctrl_h1@f.vars["catch Offshore_Trawl",]     <- c(0) 				+ 301

ctrl_h1@obs.vars["catch N_Chile",] 			    <- c(1,rep(2,11))
ctrl_h1@obs.vars["catch SC_Chile_PS",] 		  <- c(rep(1,12)) 		+ 101
ctrl_h1@obs.vars["catch FarNorth",] 		    <- c(1,1,1,rep(2,9)) 		+ 201
ctrl_h1@obs.vars["catch Offshore_Trawl",] 	<- c(rep(1,12)) 		+ 301

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
ctrl_h1@obs.vars["DEPM",                1] <- 404
ctrl_h1@obs.vars["Peru_Acoustic",       1] <- 405
ctrl_h1@obs.vars["Peru_CPUE",           1] <- 406
ctrl_h1@obs.vars["Offshore_CPUE_early", 1] <- 407
ctrl_h1@obs.vars["Offshore_CPUE_late",  1] <- 407   # shared with early

# biomassTreat: 0=SSB, 2=exploitable biomass, 5=total biomass
ctrl_h1@biomassTreat["Chile_AcousCS_early"] <- 5
ctrl_h1@biomassTreat["Chile_AcousCS_late"]  <- 5
ctrl_h1@biomassTreat["Chile_AcousN"]        <- 5
ctrl_h1@biomassTreat["Chile_CPUE_early"]    <- 2
ctrl_h1@biomassTreat["Chile_CPUE_late"]     <- 2
ctrl_h1@biomassTreat["DEPM"]               <- 0
ctrl_h1@biomassTreat["Peru_Acoustic"]       <- 5
ctrl_h1@biomassTreat["Peru_CPUE"]          <- 2
ctrl_h1@biomassTreat["Offshore_CPUE_early"] <- 2
ctrl_h1@biomassTreat["Offshore_CPUE_late"]  <- 2

ctrl_h1 <- update(ctrl_h1)
ctrl_h1@residuals <- TRUE

# ============================================================
# Helper: rebuild ctrl for any subset of indices
# Used by both the retrospective loop and leave-one-out fits.
# ============================================================
build_ctrl <- function(stk, idx_sub) {
  ctrl <- FLSAM.control(stk, idx_sub)
  ctrl@plus.group[] <- 1

  ctrl@states["catch N_Chile",]          <- c(1:7, rep(8,5))
  ctrl@states["catch SC_Chile_PS",]      <- c(1:8, rep(9,4))  + 101
  ctrl@states["catch FarNorth",]         <- c(1:5, rep(6,7))  + 201
  ctrl@states["catch Offshore_Trawl",]   <- c(1:7, rep(8,5))  + 301

  ctrl@f.vars["catch N_Chile",]          <- c(0,0,1,1,rep(2,8))
  ctrl@f.vars["catch SC_Chile_PS",]      <- c(0,1,1,1,rep(2,8)) + 101
  ctrl@f.vars["catch FarNorth",]         <- c(0,1,2,rep(3,9))   + 201
  ctrl@f.vars["catch Offshore_Trawl",]   <- c(0)                + 301

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
    ctrl@obs.vars[nm, 1]  <- obs_map[nm]
    ctrl@biomassTreat[nm] <- bio_map[nm]
  }

  ctrl <- update(ctrl)
  ctrl@residuals <- FALSE
  ctrl
}

# ============================================================
# Fit FLSAM models
# ============================================================

sam_h1   <- FLSAM(stk_h1,             idx_h1, ctrl_h1)

# ============================================================
# Retrospective analysis (manual loop)
# retro() cannot handle sub-series that become empty mid-peel (q-break splits).
# For each peel: window stock + indices, detect empty sub-series, drop them
# from ctrl with drop.from.control(), then fit. The named list retro_h1 is
# compatible with the downstream ssb()/fbar() extraction code.
# ============================================================
ctrl_h1@residuals <- FALSE
max_yr  <- range(stk_h1)["maxyear"]
n_retro <- 10

retro_h1 <- setNames(vector("list", n_retro),
                     as.character(max_yr - seq_len(n_retro)))

for (peel in seq_len(n_retro)) {
  yr_end   <- max_yr - peel
  stk_peel <- window(stk_h1, end = yr_end)

  # Only include indices with at least 2 years of data in [minyear, yr_end].
  # A single-year sub-series (e.g. Offshore_CPUE_late at peel->2021) causes
  # SAM's internal data.frame construction to fail with a row-count mismatch.
  in_range <- sapply(idx_h1, function(fl) range(fl)["minyear"] <= yr_end - 1)
  dropped  <- names(idx_h1)[!in_range]

  if (length(dropped) > 0) {
    cat(sprintf("Retro peel %i (->%i): dropping out-of-range: %s\n",
                peel, yr_end, paste(dropped, collapse = ", ")))
  } else {
    cat(sprintf("Retro peel %i (->%i)\n", peel, yr_end))
  }

  idx_peel  <- FLIndices(lapply(idx_h1[in_range], window, end = yr_end))
  # Rebuild ctrl from scratch for this subset — avoids drop.from.control() issues
  ctrl_peel <- build_ctrl(stk_peel, idx_peel)

  retro_h1[[as.character(yr_end)]] <- tryCatch(
    FLSAM(stk_peel, idx_peel, ctrl_peel),
    error = function(e) { message("  FAILED: ", e$message); NULL }
  )
}
retro_h1 <- Filter(Negate(is.null), retro_h1)

# Add full assessment keyed by its final year
retro_h1[[as.character(max_yr)]] <- sam_h1

# ============================================================
# Leave-one-out (LOO) analysis
# Drop one survey group at a time; split pairs dropped together.
# ============================================================

# Groups to drop — split pairs removed together
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
  fit <- tryCatch(
    FLSAM(stk_h1, idx_loo, ctrl_loo),
    error = function(e) { message("  FAILED: ", e$message); NULL }
  )
  fit
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
# Diagnostics  (H1)
# ============================================================

library(ggplot2)
dir.create("diagnostics", showWarnings = FALSE)

catch_fleets  <- c("catch N_Chile", "catch SC_Chile_PS",
                   "catch FarNorth", "catch Offshore_Trawl")
survey_fleets <- c("Chile_AcousCS_early", "Chile_AcousCS_late",
                   "Chile_AcousN",
                   "Chile_CPUE_early", "Chile_CPUE_late",
                   "DEPM", "Peru_Acoustic", "Peru_CPUE",
                   "Offshore_CPUE_early", "Offshore_CPUE_late")

# Helper: strip _early/_late suffixes for display — split pairs merge into one panel
base_fleet          <- function(x) sub("_(early|late)$", "", x)
survey_fleets_base  <- unique(base_fleet(survey_fleets))   # 7 base names

##################################################
## Catchability shift for split indices
##################################################
split_base <- names(q_breaks)   # Chile_AcousCS, Chile_CPUE, Offshore_CPUE

q_df <- catchabilities(sam_h1)

q_split <- do.call(rbind, lapply(split_base, function(b) {
  early <- subset(q_df, fleet == paste0(b, "_early"))
  late  <- subset(q_df, fleet == paste0(b, "_late"))
  if (nrow(early) == 0 || nrow(late) == 0) return(NULL)
  data.frame(
    index      = b,
    period     = factor(c("Early", "Late"), levels = c("Early", "Late")),
    log_q      = c(early$value, late$value),
    break_year = q_breaks[[b]]
  )
}))
q_split$q <- exp(q_split$log_q)

q_split$label <- sprintf("log(q) = %.3f\nq = %.4f", q_split$log_q, q_split$q)

png("diagnostics/h1_catchability_change.png", width = 1400, height = 600, res = 150)
print(
  ggplot(q_split, aes(x = period, y = log_q, colour = period, group = index)) +
    geom_line(colour = "grey60", linewidth = 0.9) +
    geom_point(size = 4) +
    geom_text(aes(label = label), vjust = -0.8, size = 3, lineheight = 0.9) +
    geom_vline(xintercept = 1.5, linetype = "dashed", colour = "grey40") +
    facet_wrap(~index, scales = "free_y") +
    scale_colour_manual(values = c("Early" = "steelblue", "Late" = "tomato"),
                        guide = "none") +
    labs(title = "Catchability shift at q-break year (JJM h1_1.14 breaks)",
         subtitle = paste(sapply(split_base, function(b)
           sprintf("%s: break %i", b, q_breaks[[b]])), collapse = "   |   "),
         x = NULL, y = "log(q)") +
    theme_bw() +
    theme(plot.subtitle = element_text(size = 9))
)
dev.off()

##################################################
## Full residual diagnostics (SAM default)
##################################################
pdf("diagnostics/h1_fit.pdf")
residual.diagnostics(sam_h1)
dev.off()

##################################################
## Stock trajectory
##################################################
png("diagnostics/h1_stock_trajectory.png", width = 1600, height = 1200, res = 150)
print(plot(sam_h1))
dev.off()

##################################################
## Stock trajectory zoom (SSB / Fbar / Rec)
##################################################
df.ssb        <- ssb(sam_h1);  df.ssb$quant  <- "SSB"
df.fbar       <- fbar(sam_h1); df.fbar$quant <- "Fbar"
df.rec        <- rec(sam_h1);  df.rec$quant  <- "Recruitment"
df.traj       <- rbind(df.ssb, df.fbar, df.rec)

png("diagnostics/h1_stock_trajectory_zoom.png", width = 1200, height = 1600, res = 150)
print(
  ggplot(subset(df.traj, year > 2002), aes(x = year, y = value)) +
    geom_ribbon(aes(ymin = lbnd, ymax = ubnd), alpha = 0.3) +
    geom_line() +
    ylim(0, NA) +
    facet_wrap(~quant, scales = "free", ncol = 1) +
    theme_bw()
)
dev.off()

##################################################
## F vars
##################################################

fv <- f.var(sam_h1)

png("diagnostics/h1_fvars.png", width = 1600, height = 900, res = 150)
print(
  ggplot(fv, aes(x = age, y = value)) +
    geom_ribbon(aes(ymin = lbnd, ymax = ubnd), alpha = 0.2) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.8) +
    facet_wrap(~ fleet) +
    labs(title = "F variance by fleet and age", x = "Age", y = "Value") +
    theme_bw()
)
dev.off()

##################################################
## Observation variance by data source
##################################################
obv     <- obs.var(sam_h1)
obv$str <- paste(obv$fleet, ifelse(is.na(obv$age), "", obv$age))
obv     <- obv[order(obv$value), ]

png("diagnostics/h1_observation_var.png", width = 1600, height = 900, res = 150)
bp <- barplot(obv$value, ylab = "Observation Variance",
              main = "Observation variances by data source",
              col  = as.integer(factor(obv$fleet)))
axis(1, at = bp, labels = obv$str, las = 3, lty = 0, mgp = c(0, 0, 0))
legend("topleft", levels(factor(obv$fleet)), pch = 15,
       col = seq_len(nlevels(factor(obv$fleet))), pt.cex = 1.5)
dev.off()

##################################################
## Observation variance vs CV
##################################################
png("diagnostics/h1_variance_vs_cv.png", width = 1200, height = 900, res = 150)
plot(obv$value, obv$CV, xlab = "Observation variance", ylab = "CV of estimate",
     log = "x", pch = 16, col = as.integer(factor(obv$fleet)),
     main = "Observation variance vs uncertainty")
text(obv$value, obv$CV, obv$str, pos = 4, cex = 0.75, xpd = NA)
dev.off()

##################################################
## Selectivity (F / Fbar at age) - one panel per fleet x pentad
## f() has fleet + age columns; fbar() is a scalar per year,
## so merge on year only, then facet by fleet.
##################################################
sel.pat       <- merge(f(sam_h1), fbar(sam_h1), by = "year", suffixes = c(".f", ".fbar"))
sel.pat$sel   <- sel.pat$value.f / sel.pat$value.fbar
sel.pat$age   <- as.numeric(as.character(sel.pat$age))
sel.pat$pentad <- sprintf("%i's", floor(sel.pat$year / 5) * 5)

png("diagnostics/h1_selectivity.png", width = 2000, height = 1400, res = 150)
print(
  ggplot(sel.pat, aes(x = age, y = sel, group = year, colour = year)) +
    geom_line(alpha = 0.6) +
    geom_point(size = 1.2) +
    facet_grid(fleet ~ pentad, scales = "free_y") +
    scale_colour_viridis_c() +
    scale_x_continuous(breaks = seq(min(sel.pat$age), max(sel.pat$age), by = 1)) +
    labs(title = "Selectivity (F/Fbar) by fleet and pentad", x = "Age", y = "F/Fbar") +
    theme_bw() +
    theme(axis.text.x = element_text(size = 7))
)
dev.off()

##################################################
## Correlation matrix of model parameters
##################################################
png("diagnostics/h1_cor_params.png", width = 1400, height = 1400, res = 150)
cor.plot(sam_h1)
dev.off()

##################################################
## Catch residuals (bubble plot, one panel per fleet)
##################################################
res_all <- residuals(sam_h1)

png("diagnostics/h1_catch_residuals.png", width = 1800, height = 1200, res = 150)
print(
  ggplot(subset(res_all, fleet %in% catch_fleets),
         aes(x = year, y = age, size = 10 * abs(std.res), fill = std.res > 0)) +
    geom_point(shape = 21) +
    scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "red")) +
    facet_wrap(~fleet, ncol = 2) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(title = "Standardised residuals - Catch fleets")
)
dev.off()

##################################################
## Survey residuals (biomass indices: year only)
##################################################
res_survey_disp <- subset(res_all, fleet %in% survey_fleets)
res_survey_disp$fleet <- base_fleet(res_survey_disp$fleet)

png("diagnostics/h1_survey_residuals.png", width = 1800, height = 1400, res = 150)
print(
  ggplot(res_survey_disp, aes(x = year, y = std.res)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(aes(fill = std.res > 0), shape = 21, size = 3) +
    geom_segment(aes(xend = year, yend = 0)) +
    scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "red")) +
    facet_wrap(~fleet, scales = "free_y", ncol = 2) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(title = "Standardised residuals - Survey indices", y = "Std. residual")
)
dev.off()

##################################################
## Process error - N / M
## procerr.plot merges slots by area; after stk + sam_h1,
## harvest has 4 areas (one per fleet) but m/stock.n/stock.wt
## have 1 area -> non-conformable arrays.
## Fix: collapse harvest to total F (sum over fleets) so all
## slots share a single area before passing to procerr.plot.
##################################################
stk_h1_fit     <- stk_h1 + sam_h1
stk_h1_procerr <- stk_h1_fit
h_dn           <- dimnames(harvest(stk_h1_fit)); h_dn[["area"]] <- "unique"
harvest(stk_h1_procerr) <- FLQuant(
  apply(harvest(stk_h1_fit), c(1, 2, 3, 4, 6), sum, na.rm = TRUE),
  dimnames = h_dn
)

png("diagnostics/h1_procerr_N.png", width = 1600, height = 900, res = 150)
procerr.plot(stk_h1_procerr, weight = "stock.wt", type = "n", rel = TRUE)
dev.off()

png("diagnostics/h1_procerr_M.png", width = 1600, height = 900, res = 150)
procerr.plot(stk_h1_procerr, weight = "stock.wt", type = "mort", rel = TRUE)
dev.off()

##################################################
## F at age - one panel per catch fleet (area)
## harvest() has 4 areas matching the catch fleets;
## the area dimname carries the fleet label.
##################################################
stk_h1_fit <- stk_h1 + sam_h1

png("diagnostics/h1_fatage.png", width = 1800, height = 1200, res = 150)
print(
  ggplot(as.data.frame(window(harvest(stk_h1_fit), start = 2000)),
         aes(x = year, y = data, colour = factor(age), group = age)) +
    geom_line() +
    facet_wrap(~area, ncol = 2, scales = "free_y") +
    labs(title = "F at age by fleet", x = "Year", y = "F", colour = "Age") +
    theme_bw()
)
dev.off()

##################################################
## Stock-recruit relationship
##################################################
df.sr       <- merge(
  setNames(ssb(sam_h1)[, c("year", "value")],  c("year", "ssb")),
  setNames(rec(sam_h1)[, c("year", "value")],  c("year", "rec")),
  by = "year"
)

png("diagnostics/h1_stock_recruit.png", width = 1200, height = 900, res = 150)
print(
  ggplot(df.sr, aes(x = ssb, y = rec, colour = year)) +
    geom_path(colour = "grey70") +
    geom_point(size = 3) +
    scale_colour_viridis_c() +
    labs(title = "Stock-recruit relationship", x = "SSB", y = "Recruitment age-1") +
    theme_bw()
)
dev.off()

##################################################
## Retrospective - SSB and Fbar
## retro_h1 is an FLSAMs list (n=10 peels)
##################################################
retro_ssb  <- do.call(rbind, lapply(names(retro_h1), function(nm) {
  d <- ssb(retro_h1[[nm]]); d$peel <- nm; d
}))
retro_fbar <- do.call(rbind, lapply(names(retro_h1), function(nm) {
  d <- fbar(retro_h1[[nm]]); d$peel <- nm; d
}))
retro_ssb$peel  <- factor(retro_ssb$peel,  levels = names(retro_h1))
retro_fbar$peel <- factor(retro_fbar$peel, levels = names(retro_h1))

base_ssb  <- ssb(sam_h1);  base_ssb$peel  <- as.character(max_yr)
base_fbar <- fbar(sam_h1); base_fbar$peel <- as.character(max_yr)

retro_ssb_all  <- rbind(retro_ssb,  base_ssb)
retro_fbar_all <- rbind(retro_fbar, base_fbar)

base_yr_chr <- as.character(max_yr)

png("diagnostics/h1_retro_ssb.png", width = 1600, height = 900, res = 150)
print(
  ggplot(retro_ssb_all, aes(x = year, y = value, colour = peel, group = peel)) +
    geom_line(data = subset(retro_ssb_all, peel != base_yr_chr), alpha = 0.7) +
    geom_line(data = subset(retro_ssb_all, peel == base_yr_chr), linewidth = 1.2, colour = "black") +
    geom_ribbon(data = subset(retro_ssb_all, peel == base_yr_chr),
                aes(ymin = lbnd, ymax = ubnd), alpha = 0.2, colour = NA, fill = "black") +
    coord_cartesian(xlim = c(max(retro_ssb_all$year) - 15, max(retro_ssb_all$year))) +
    ylim(0, NA) +
    labs(title = "Retrospective - SSB", x = "Year", y = "SSB", colour = "Final year") +
    theme_bw()
)
dev.off()

png("diagnostics/h1_retro_fbar.png", width = 1600, height = 900, res = 150)
print(
  ggplot(retro_fbar_all, aes(x = year, y = value, colour = peel, group = peel)) +
    geom_line(data = subset(retro_fbar_all, peel != base_yr_chr), alpha = 0.7) +
    geom_line(data = subset(retro_fbar_all, peel == base_yr_chr), linewidth = 1.2, colour = "black") +
    geom_ribbon(data = subset(retro_fbar_all, peel == base_yr_chr),
                aes(ymin = lbnd, ymax = ubnd), alpha = 0.2, colour = NA, fill = "black") +
    coord_cartesian(xlim = c(max(retro_fbar_all$year) - 15, max(retro_fbar_all$year))) +
    ylim(0, NA) +
    labs(title = "Retrospective - Fbar", x = "Year", y = "Fbar", colour = "Final year") +
    theme_bw()
)
dev.off()

##################################################
## Mohn's rho
##################################################
mohn_rho <- function(base_df, retro_list_df) {
  peels <- setdiff(unique(retro_list_df$peel), as.character(max_yr))
  rhos  <- sapply(peels, function(p) {
    sub  <- subset(retro_list_df, peel == p)
    yr   <- max(sub$year)
    base_val <- base_df$value[base_df$year == yr]
    retro_val <- sub$value[sub$year == yr]
    if (length(base_val) == 0 || length(retro_val) == 0) return(NA)
    (retro_val - base_val) / base_val
  })
  mean(rhos, na.rm = TRUE)
}

rho_ssb  <- mohn_rho(base_ssb,  retro_ssb)
rho_fbar <- mohn_rho(base_fbar, retro_fbar)
cat(sprintf("Mohn's rho - SSB: %.4f | Fbar: %.4f\n", rho_ssb, rho_fbar))

sink("diagnostics/h1_mohns_rho.txt")
cat(sprintf("Mohn's rho\n  SSB : %.4f\n  Fbar: %.4f\n", rho_ssb, rho_fbar))
sink()

##################################################
## Leave-one-out plots
##################################################
extract_loo <- function(fun, fits, base_fit) {
  base_df         <- fun(base_fit)
  base_df$dropped <- "Full model"
  loo_df <- do.call(rbind, lapply(names(fits), function(grp) {
    if (is.null(fits[[grp]])) return(NULL)
    df         <- fun(fits[[grp]])
    df$dropped <- grp
    df
  }))
  rbind(base_df, loo_df)
}

loo_ssb_all <- extract_loo(ssb,  loo_fits, sam_h1)
loo_fbar_all <- extract_loo(fbar, loo_fits, sam_h1)
loo_rec_all  <- extract_loo(rec,  loo_fits, sam_h1)

plot_loo <- function(dat, ylab, title) {
  base  <- subset(dat, dropped == "Full model")
  loodf <- subset(dat, dropped != "Full model")
  ggplot() +
    geom_ribbon(data = base,
                aes(x = year, ymin = lbnd, ymax = ubnd),
                fill = "grey70", alpha = 0.5) +
    geom_line(data = base,
              aes(x = year, y = value),
              colour = "black", linewidth = 1) +
    geom_ribbon(data = loodf,
                aes(x = year, ymin = lbnd, ymax = ubnd),
                fill = "tomato", alpha = 0.25) +
    geom_line(data = loodf,
              aes(x = year, y = value),
              colour = "tomato", linewidth = 0.9) +
    facet_wrap(~dropped, ncol = 2) +
    ylim(0, NA) +
    labs(title = title,
         subtitle = "Black/grey = full model; red = model with survey group dropped",
         x = "Year", y = ylab) +
    theme_bw()
}

png("diagnostics/h1_loo_ssb.png",  width = 1800, height = 1400, res = 150)
print(plot_loo(loo_ssb_all,  "SSB",         "Leave-one-out — SSB"))
dev.off()

png("diagnostics/h1_loo_fbar.png", width = 1800, height = 1400, res = 150)
print(plot_loo(loo_fbar_all, "Mean F",      "Leave-one-out — Fbar"))
dev.off()

png("diagnostics/h1_loo_rec.png",  width = 1800, height = 1400, res = 150)
print(plot_loo(loo_rec_all,  "Recruitment", "Leave-one-out — Recruitment"))
dev.off()

##################################################
## Index fits - observed vs model-predicted
##################################################
fit_idx <- residuals(sam_h1)

# For biomass surveys age is NA; use year on x-axis
survey_fit <- subset(fit_idx, fleet %in% survey_fleets)
survey_fit$fleet <- base_fleet(survey_fit$fleet)

png("diagnostics/h1_index_fits.png", width = 1800, height = 1400, res = 150)
print(
  ggplot(survey_fit, aes(x = year)) +
    geom_point(aes(y = exp(log.obs)), colour = "black", size = 1.5) +
    geom_line(aes(y = exp(log.mdl)), colour = "red",   linewidth = 0.8) +
    facet_wrap(~fleet, scales = "free_y", ncol = 2) +
    labs(title = "Survey index fits (obs = points, fitted = red line)",
         x = "Year", y = "Index") +
    theme_bw()
)
dev.off()

##################################################
## Input: catch proportions at age by fleet
##################################################
catch_prop <- as.data.frame(stk_h1@catch.n)
catch_prop$data[catch_prop$data == -1] <- NA

png("diagnostics/h1_catch_proportions.png", width = 1800, height = 1200, res = 150)
print(
  ggplot(subset(catch_prop, year >= 2000 & !is.na(data)),
         aes(x = year, y = data, fill = factor(age))) +
    geom_bar(stat = "identity", position = "fill") +
    facet_wrap(~area, ncol = 2) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(title = "Catch proportions at age by fleet",
         x = "Year", y = "Proportion", fill = "Age") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
)
dev.off()

##################################################
## Input: stock weight at age timeseries
##################################################
png("diagnostics/h1_stock_wt.png", width = 1600, height = 900, res = 150)
print(
  ggplot(as.data.frame(window(stk_h1@stock.wt, start = 2000)),
         aes(x = year, y = data, colour = factor(age), group = age)) +
    geom_line() +
    labs(title = "Stock weight at age", x = "Year", y = "Weight (kg)", colour = "Age") +
    theme_bw()
)
dev.off()

##################################################
## Input: catch weight at age by fleet
##################################################
png("diagnostics/h1_catch_wt.png", width = 1800, height = 1200, res = 150)
print(
  ggplot(subset(as.data.frame(window(stk_h1@catch.wt, start = 2000)),
                data > 0 & !is.na(data)),
         aes(x = year, y = data, colour = factor(age), group = age)) +
    geom_line() +
    facet_wrap(~area, ncol = 2, scales = "free_y") +
    labs(title = "Catch weight at age by fleet", x = "Year", y = "Weight (kg)", colour = "Age") +
    theme_bw()
)
dev.off()

##################################################
## Input: maturity at age
##################################################
png("diagnostics/h1_maturity.png", width = 1600, height = 900, res = 150)
print(
  ggplot(as.data.frame(window(stk_h1@mat, start = 2000)),
         aes(x = year, y = data, colour = factor(age), group = age)) +
    geom_line() +
    ylim(0, 1) +
    labs(title = "Maturity at age", x = "Year", y = "Proportion mature", colour = "Age") +
    theme_bw()
)
dev.off()

##################################################
## Input: natural mortality at age
##################################################
png("diagnostics/h1_natmort.png", width = 1600, height = 900, res = 150)
print(
  ggplot(as.data.frame(stk_h1@m),
         aes(x = year, y = data, colour = factor(age), group = age)) +
    geom_line() +
    labs(title = "Natural mortality (M) at age", x = "Year", y = "M", colour = "Age") +
    theme_bw()
)
dev.off()

##################################################
## Survey raw index timeseries (z-normalised overlay)
##################################################
survey_raw <- do.call(rbind, lapply(names(idx_h1), function(nm) {
  d        <- as.data.frame(idx_h1[[nm]]@index)
  d$fleet  <- base_fleet(nm)   # merge early/late into one series
  d$data[d$data <= 0] <- NA
  d
}))
survey_raw <- subset(survey_raw, !is.na(data))
# z-normalise per base fleet (early+late treated as one continuous series)
survey_raw <- do.call(rbind, lapply(split(survey_raw, survey_raw$fleet), function(df) {
  df$z <- (df$data - mean(df$data, na.rm = TRUE)) / sd(df$data, na.rm = TRUE)
  df
}))

png("diagnostics/h1_survey_timeseries.png", width = 1800, height = 1200, res = 150)
print(
  ggplot(survey_raw, aes(x = year, y = z, colour = fleet, group = fleet)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.5) +
    labs(title = "Survey indices (z-normalised)", x = "Year", y = "Standardised index",
         colour = "Survey") +
    theme_bw()
)
dev.off()

##################################################
## N at age - stock numbers timeseries by fleet
##################################################
png("diagnostics/h1_natage.png", width = 1600, height = 900, res = 150)
print(
  ggplot(as.data.frame(window(stk_h1_fit@stock.n, start = 2000)),
         aes(x = year, y = data, colour = factor(age), group = age)) +
    geom_line() +
    scale_y_log10() +
    labs(title = "Stock numbers at age (log scale)", x = "Year",
         y = "Numbers (log10)", colour = "Age") +
    theme_bw()
)
dev.off()

################################################################################
## CATCHABILITY & FLEET CONFLICT DIAGNOSTICS
################################################################################

# ------------------------------------------------------------------
# Shared prep: get residuals and SSB once
# ------------------------------------------------------------------
res_all  <- residuals(sam_h1)           # fleet, year, age, observed, fitted, std.res
ssb_df   <- ssb(sam_h1)                 # year, value (model SSB)

# ------------------------------------------------------------------
# 1. IMPLIED CATCHABILITY OVER TIME
# ------------------------------------------------------------------
# In SAM log-space: fitted = log(q) + log(B_pred)
# So residual = log(obs) - log(q*B_pred)
# A rolling mean of residuals reveals a drift in effective q.
# Cumulative sum of residuals = cumulative log-q change from model q.
# ------------------------------------------------------------------
survey_res <- subset(res_all, fleet %in% survey_fleets)
survey_res$fleet <- base_fleet(survey_res$fleet)

# Rolling mean of std residuals (window = 5 years)
roll_mean <- function(x, k = 5) {
  n <- length(x); out <- rep(NA, n)
  for (i in seq_len(n)) {
    idx <- max(1, i - floor(k/2)):min(n, i + floor(k/2))
    out[i] <- mean(x[idx], na.rm = TRUE)
  }
  out
}

survey_res_roll <- do.call(rbind, lapply(split(survey_res, survey_res$fleet), function(df) {
  df <- df[order(df$year), ]
  df$roll_res <- roll_mean(df$std.res, k = 5)
  df$cum_res  <- cumsum(ifelse(is.na(df$std.res), 0, df$std.res))
  df
}))

png("diagnostics/h1_q_rolling_residuals.png", width = 1800, height = 1400, res = 150)
print(
  ggplot(survey_res_roll, aes(x = year)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_col(aes(y = std.res, fill = std.res > 0), alpha = 0.4, width = 0.8) +
    geom_line(aes(y = roll_res), colour = "black", linewidth = 1) +
    scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "tomato"), guide = "none") +
    facet_wrap(~fleet, scales = "free_y", ncol = 2) +
    labs(title = "Survey residuals with 5-yr rolling mean (trend = q drift)",
         x = "Year", y = "Std. residual") +
    theme_bw()
)
dev.off()

# Cumulative residuals - sustained positive = q declining, negative = q increasing
png("diagnostics/h1_q_cumulative_residuals.png", width = 1800, height = 1200, res = 150)
print(
  ggplot(survey_res_roll, aes(x = year, y = cum_res)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_line(linewidth = 1, colour = "steelblue") +
    geom_point(size = 1.5) +
    facet_wrap(~fleet, scales = "free_y", ncol = 2) +
    labs(title = "Cumulative standardised residuals per survey (slope = q trend)",
         x = "Year", y = "Cumulative std. residual") +
    theme_bw()
)
dev.off()

# ------------------------------------------------------------------
# 2. RUNS TEST (SS3-style)
# Tests whether signs of residuals are random; a non-random
# pattern (too few runs) indicates a systematic trend in q or fit.
# ------------------------------------------------------------------
runs_test <- function(x) {
  x   <- x[!is.na(x)]
  sgn <- sign(x)
  n   <- length(sgn)
  if (n < 4) return(data.frame(n = n, n_runs = NA, p_value = NA, result = "too few obs"))
  runs <- sum(sgn[-1] != sgn[-n]) + 1
  n1   <- sum(sgn > 0); n2 <- sum(sgn < 0)
  if (n1 == 0 || n2 == 0) return(data.frame(n = n, n_runs = runs, p_value = NA, result = "all same sign"))
  mu   <- (2 * n1 * n2 / n) + 1
  sig2 <- (2 * n1 * n2 * (2 * n1 * n2 - n)) / (n^2 * (n - 1))
  z    <- (runs - mu) / sqrt(sig2)
  p    <- 2 * pnorm(-abs(z))
  data.frame(n = n, n_runs = runs, expected_runs = round(mu, 1),
             z = round(z, 3), p_value = round(p, 4),
             result = ifelse(p < 0.05, "FAIL (trend)", "pass"))
}

runs_results <- do.call(rbind, lapply(c(catch_fleets, survey_fleets_base), function(fl) {
  # Match on base fleet name so early+late are combined into one series
  sub <- subset(res_all, base_fleet(fleet) == fl)
  if (fl %in% catch_fleets) {
    sub <- aggregate(std.res ~ year, data = sub, FUN = mean, na.rm = TRUE)
    rt  <- runs_test(sub$std.res)
  } else {
    sub <- sub[order(sub$year), ]
    rt  <- runs_test(sub$std.res)
  }
  cbind(fleet = fl, rt)
}))

write.csv(runs_results, "diagnostics/h1_runs_test.csv", row.names = FALSE)

png("diagnostics/h1_runs_test.png", width = 1400, height = 700, res = 150)
print(
  ggplot(runs_results, aes(x = fleet, y = p_value, fill = result)) +
    geom_col() +
    geom_hline(yintercept = 0.05, linetype = "dashed", colour = "red") +
    scale_fill_manual(values = c("pass" = "steelblue", "FAIL (trend)" = "tomato",
                                 "too few obs" = "grey70", "all same sign" = "grey50")) +
    coord_flip() +
    labs(title = "Runs test p-value by fleet (p < 0.05 = non-random residuals)",
         x = NULL, y = "p-value", fill = NULL) +
    theme_bw()
)
dev.off()

# ------------------------------------------------------------------
# 3. SURVEY Q CHANGE: index / model-predicted SSB over time
# A slope in this ratio reveals a trend in effective catchability.
# Fitted = log(q) + log(B), so obs/fitted * q_est = implied q_t.
# We approximate: implied_q_t ~ exp(observed) / SSB_model_t
# (only valid if the survey covers total SSB, approximate otherwise)
# ------------------------------------------------------------------
implied_q <- merge(survey_res, ssb_df[, c("year", "value")], by = "year")
implied_q$implied_log_q <- implied_q$log.obs - log(implied_q$value)

# Normalise within fleet so curves start at 0 (relative change in q)
implied_q <- do.call(rbind, lapply(split(implied_q, implied_q$fleet), function(df) {
  df <- df[order(df$year), ]
  df$rel_log_q <- df$implied_log_q - mean(df$implied_log_q, na.rm = TRUE)
  df
}))

png("diagnostics/h1_q_trend.png", width = 1800, height = 1200, res = 150)
print(
  ggplot(implied_q, aes(x = year, y = rel_log_q)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(aes(colour = rel_log_q > 0), size = 2) +
    geom_smooth(method = "lm", se = TRUE, colour = "black", linewidth = 0.8) +
    scale_colour_manual(values = c("TRUE" = "steelblue", "FALSE" = "tomato"), guide = "none") +
    facet_wrap(~fleet, scales = "free_y", ncol = 2) +
    labs(title = "Implied catchability trend: log(obs/SSB) centred (slope ≠ 0 = q drift)",
         x = "Year", y = "Relative log(q)") +
    theme_bw()
)
dev.off()

# ------------------------------------------------------------------
# 4. RESIDUAL CROSS-CORRELATION MATRIX (fleet conflict / compensation)
# Per-year residuals (mean over ages for catch fleets) are correlated
# across all fleets.  A strong NEGATIVE correlation between two fleets
# means one is consistently over-predicting when the other under-predicts
# -> they are pulling the model in opposite directions (conflict).
# ------------------------------------------------------------------
annual_res <- do.call(rbind, lapply(c(catch_fleets, survey_fleets), function(fl) {
  sub <- subset(res_all, fleet == fl)
  ag  <- aggregate(std.res ~ year, data = sub, FUN = mean, na.rm = TRUE)
  ag$fleet <- fl
  ag
}))

res_wide <- reshape(annual_res, idvar = "year", timevar = "fleet", direction = "wide")
colnames(res_wide) <- gsub("std.res\\.", "", colnames(res_wide))

res_mat  <- res_wide[, -1]                 # drop year
cor_mat  <- cor(res_mat, use = "pairwise.complete.obs")

cor_df   <- as.data.frame(as.table(cor_mat))
colnames(cor_df) <- c("fleet1", "fleet2", "correlation")

png("diagnostics/h1_residual_crosscor.png", width = 1600, height = 1400, res = 150)
print(
  ggplot(cor_df, aes(x = fleet1, y = fleet2, fill = correlation)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = round(correlation, 2)), size = 3) +
    scale_fill_gradient2(low = "tomato", mid = "white", high = "steelblue",
                         midpoint = 0, limits = c(-1, 1)) +
    labs(title = "Residual cross-correlation matrix\n(negative = fleets conflict; positive = co-vary)",
         x = NULL, y = NULL) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
)
dev.off()

# ------------------------------------------------------------------
# 5. RETROSPECTIVE RESIDUALS PER SURVEY (SS3-style hindcast)
# Using retro_h1: in each peeled fit, what was the one-step-ahead
# prediction for the terminal year of that peel?
# Large deviations across peels = poor predictive skill for that fleet.
# ------------------------------------------------------------------
# Full hindcast residual extraction for all surveys
hindcast_res <- do.call(rbind, lapply(names(retro_h1), function(nm) {
  r <- residuals(retro_h1[[nm]])
  r <- subset(r, fleet %in% survey_fleets)
  r$peel      <- nm
  r$terminal  <- max(r$year)
  r$is_terminal <- r$year == max(r$year)
  r
}))

# Extract only terminal-year predictions from each peel
terminal_preds <- subset(hindcast_res, is_terminal)

# Compare to base run residual for same year
base_res_survey <- subset(res_all, fleet %in% survey_fleets)[, c("fleet", "year", "std.res")]
colnames(base_res_survey)[3] <- "base_std_res"
terminal_preds <- merge(terminal_preds, base_res_survey, by = c("fleet", "year"), all.x = TRUE)

png("diagnostics/h1_hindcast_residuals.png", width = 1800, height = 1200, res = 150)
print(
  ggplot(terminal_preds, aes(x = year, y = std.res)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_segment(aes(xend = year, yend = 0), colour = "grey70") +
    geom_point(aes(colour = abs(std.res) > 2), size = 3) +
    geom_point(aes(y = base_std_res), shape = 4, size = 3, colour = "black") +
    scale_colour_manual(values = c("FALSE" = "steelblue", "TRUE" = "tomato"),
                        name = "|residual| > 2") +
    facet_wrap(~fleet, scales = "free_y", ncol = 2) +
    labs(title = "Hindcast terminal-year residuals per survey\n(dot = peel prediction, cross = base model)",
         x = "Year", y = "Std. residual") +
    theme_bw()
)
dev.off()

# ------------------------------------------------------------------
# 6. INTER-SURVEY CONSISTENCY: pairwise index scatter
# Plot pairs of survey indices over shared years.
# Divergence or negative correlation = surveys conflict on stock trend.
# ------------------------------------------------------------------
survey_pairs <- merge(
  reshape(survey_res[, c("fleet", "year", "log.obs")],
          idvar = "year", timevar = "fleet", direction = "wide"),
  data.frame(year = ssb_df$year), by = "year"
)
colnames(survey_pairs) <- gsub("log\\.obs\\.", "", colnames(survey_pairs))

# survey_res$fleet already has base names; use survey_fleets_base for pairs
pair_combos <- combn(survey_fleets_base, 2, simplify = FALSE)
pair_plots  <- lapply(pair_combos, function(fl) {
  tmp <- survey_pairs[, c("year", fl[1], fl[2])]
  colnames(tmp) <- c("year", "x", "y")
  tmp <- tmp[complete.cases(tmp), ]
  if (nrow(tmp) < 4) return(NULL)
  r   <- round(cor(tmp$x, tmp$y), 2)
  ggplot(tmp, aes(x = x, y = y, colour = year)) +
    geom_point(size = 2) +
    scale_colour_viridis_c() +
    geom_smooth(method = "lm", se = FALSE, colour = "black", linewidth = 0.6) +
    labs(title = sprintf("%s vs %s\nr = %s", fl[1], fl[2], r),
         x = fl[1], y = fl[2]) +
    theme_bw() +
    theme(legend.position = "none", plot.title = element_text(size = 8))
})
pair_plots <- Filter(Negate(is.null), pair_plots)

library(patchwork)
n_cols <- 3
png("diagnostics/h1_survey_pairwise.png",
    width = n_cols * 600, height = ceiling(length(pair_plots) / n_cols) * 600, res = 150)
print(wrap_plots(pair_plots, ncol = n_cols))
dev.off()

# ------------------------------------------------------------------
# 7. COHORT TRACKING THROUGH CATCH-AT-AGE FLEETS
# For a cohort c, plot catch[age, c+age] across all fleets.
# Diverging trends between fleets for the same cohort reveal conflicts
# in the age-composition data used to estimate selectivity / F.
# ------------------------------------------------------------------
catch_n_df <- as.data.frame(stk_h1@catch.n)
catch_n_df$cohort <- catch_n_df$year - catch_n_df$age
catch_n_df$data[catch_n_df$data <= 0] <- NA
catch_n_df <- subset(catch_n_df, !is.na(data) & year >= 2000)

# Normalise within fleet x cohort so relative trajectory is visible
catch_n_df <- do.call(rbind, lapply(
  split(catch_n_df, list(catch_n_df$area, catch_n_df$cohort)),
  function(df) {
    df$rel <- df$data / mean(df$data, na.rm = TRUE)
    df
  }
))

recent_cohorts <- sort(unique(catch_n_df$cohort))
recent_cohorts <- tail(recent_cohorts[recent_cohorts >= max(catch_n_df$cohort) - 10], 8)

png("diagnostics/h1_cohort_tracking.png", width = 2000, height = 1400, res = 150)
print(
  ggplot(subset(catch_n_df, cohort %in% recent_cohorts),
         aes(x = age, y = rel, colour = area, group = area)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.5) +
    facet_wrap(~cohort, ncol = 4, scales = "free_y") +
    labs(title = "Cohort tracking across catch fleets (normalised catch at age)\nDivergence between fleets = selectivity / q conflict",
         x = "Age", y = "Relative catch (fleet mean = 1)", colour = "Fleet") +
    theme_bw()
)
dev.off()

# ------------------------------------------------------------------
# 8. FLEET-SPECIFIC MEAN AGE OF CATCH: observed vs fitted
# A persistent offset between observed and model-predicted mean age
# signals that the model is misallocating F across ages for that fleet.
# ------------------------------------------------------------------
mean_age_obs  <- aggregate(age * abs(std.res) ~ year + fleet,
                           data = subset(res_all, fleet %in% catch_fleets & !is.na(age)),
                           FUN = sum)
mean_age_wt   <- aggregate(abs(std.res) ~ year + fleet,
                           data = subset(res_all, fleet %in% catch_fleets & !is.na(age)),
                           FUN = sum)
colnames(mean_age_obs)[3] <- "wsum"
colnames(mean_age_wt)[3]  <- "wt"
mean_age_df <- merge(mean_age_obs, mean_age_wt, by = c("year", "fleet"))
mean_age_df$wmean_age_resid <- mean_age_df$wsum / mean_age_df$wt

png("diagnostics/h1_mean_age_residual.png", width = 1800, height = 900, res = 150)
print(
  ggplot(mean_age_df, aes(x = year, y = wmean_age_resid)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_line(colour = "steelblue", linewidth = 0.8) +
    geom_point(aes(fill = wmean_age_resid > 0), shape = 21, size = 2.5) +
    scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "tomato"), guide = "none") +
    facet_wrap(~fleet, ncol = 2, scales = "free_y") +
    labs(title = "Residual-weighted mean age discrepancy by catch fleet\n(positive = model underfits older fish)",
         x = "Year", y = "Residual-weighted mean age") +
    theme_bw()
)
dev.off()

cat("Diagnostics written to diagnostics/\n")

################################################################################
## PRODUCTIVITY ANALYSIS: BEVERTON-HOLT S-R WITH CHANGE-POINT DETECTION
##
## Uncertainty in both R and SSB estimates from SAM is carried through:
##   - Weights for NLS fitting = 1 / Var(log R), approximated from 95% CI
##   - Error bars on both axes in all plots
## Change points are found by minimising total weighted residual SS across
## candidate break years; significance is tested with an F-test (extra
## parameters = 2 per additional B-H curve). Up to 3 panels are produced:
##   Panel 1 – full time series, single B-H, steepness annotated
##   Panel 2 – if 1 break is significant (p<0.05): 2 B-H curves
##   Panel 3 – if a 2nd break is also significant: 3 B-H curves
################################################################################

# ------------------------------------------------------------------
# Collect S-R pairs with uncertainty from SAM
# ------------------------------------------------------------------
rec_sr <- rec(sam_h1)
ssb_sr <- ssb(sam_h1)

sr <- merge(
  data.frame(year  = rec_sr$year, R    = rec_sr$value,
             R_lo  = rec_sr$lbnd, R_hi = rec_sr$ubnd),
  data.frame(year    = ssb_sr$year, SSB    = ssb_sr$value,
             SSB_lo  = ssb_sr$lbnd, SSB_hi = ssb_sr$ubnd),
  by = "year"
)
sr <- sr[complete.cases(sr) & sr$R > 0 & sr$SSB > 0, ]
sr <- sr[order(sr$year), ]

# Observation weight = 1 / Var(log R), approximated from 95% CI
# Var(log R) ≈ ((R_hi - R_lo) / (4 * R))^2
sr$logR_se <- pmax((sr$R_hi - sr$R_lo) / (4 * sr$R), 0.01)
sr$wt      <- 1 / sr$logR_se^2
sr$wt      <- sr$wt / max(sr$wt)   # scale to [0,1]

# Proxy for unfished SSB: use max observed SSB as B0 reference for steepness
B0_ref <- max(sr$SSB)
n_sr   <- nrow(sr)

# ------------------------------------------------------------------
# Beverton-Holt NLS fit (weighted)
# Returns list: converged, a, b, h, wrss, fit
# B-H: R = a*S/(b+S)
# Steepness: h = 0.2*(b+B0)/(b+0.2*B0)  (bounded [0.21, 0.99])
# ------------------------------------------------------------------
fit_bh <- function(df, B0 = B0_ref) {
  if (nrow(df) < 6) return(list(converged = FALSE))
  starts <- list(a = quantile(df$R, 0.9) * 2,
                 b = median(df$SSB))
  tryCatch({
    fit <- nls(R ~ a * SSB / (b + SSB), data = df, weights = wt,
               start = starts,
               control = nls.control(maxiter = 500, warnOnly = FALSE))
    cf <- coef(fit)
    a  <- unname(cf["a"]); b <- unname(cf["b"])
    if (a <= 0 || b <= 0) return(list(converged = FALSE))
    h  <- pmin(pmax(0.2 * (b + B0) / (b + 0.2 * B0), 0.21), 0.99)
    list(converged = TRUE, a = a, b = b, h = round(h, 3),
         wrss = sum(df$wt * residuals(fit)^2), fit = fit)
  }, error = function(e) list(converged = FALSE))
}

# Predict B-H curve over a SSB vector
bh_pred <- function(a, b, ssb_vec) a * ssb_vec / (b + ssb_vec)

# Smooth SSB grid for curves (full observed range)
ssb_grid <- seq(0, max(sr$SSB) * 1.05, length.out = 300)

# ------------------------------------------------------------------
# 1. Full series fit
# ------------------------------------------------------------------
fit_all <- fit_bh(sr)
wrss_full <- if (fit_all$converged) fit_all$wrss else Inf

# ------------------------------------------------------------------
# 2. Single break-point search
# ------------------------------------------------------------------
min_seg_sr <- 10   # min years per segment

break1_res <- do.call(rbind, lapply(min_seg_sr:(n_sr - min_seg_sr), function(i) {
  f1 <- fit_bh(sr[seq_len(i), ])
  f2 <- fit_bh(sr[(i + 1):n_sr, ])
  if (!f1$converged || !f2$converged) return(NULL)
  data.frame(
    idx        = i,
    break_year = sr$year[i],
    wrss       = f1$wrss + f2$wrss,
    a1=f1$a, b1=f1$b, h1=f1$h,
    a2=f2$a, b2=f2$b, h2=f2$h,
    p1 = sprintf("%d\u2013%d", sr$year[1],   sr$year[i]),
    p2 = sprintf("%d\u2013%d", sr$year[i+1], sr$year[n_sr])
  )
}))

do_2panel <- FALSE; p1_val <- NA; best1 <- NULL
if (!is.null(break1_res) && nrow(break1_res) > 0 && fit_all$converged) {
  best1   <- break1_res[which.min(break1_res$wrss), ]
  f1stat  <- ((wrss_full - best1$wrss) / 2) / (best1$wrss / (n_sr - 4))
  p1_val  <- pf(f1stat, 2, n_sr - 4, lower.tail = FALSE)
  do_2panel <- isTRUE(p1_val < 0.05)
  cat(sprintf("Single break: year %d | F=%.2f | p=%.4f\n",
              best1$break_year, f1stat, p1_val))
}

# ------------------------------------------------------------------
# 3. Two break-point search (greedy: fix first break, find best second)
# ------------------------------------------------------------------
do_3panel <- FALSE; p2_val <- NA; best2 <- NULL
if (do_2panel) {
  i1 <- best1$idx

  # Second break in segment 1 (before first break)
  seg2a <- if (i1 >= 2 * min_seg_sr) {
    do.call(rbind, lapply(min_seg_sr:(i1 - min_seg_sr), function(j) {
      fa <- fit_bh(sr[seq_len(j), ])
      fb <- fit_bh(sr[(j + 1):i1, ])
      fc <- fit_bh(sr[(i1 + 1):n_sr, ])
      if (!fa$converged || !fb$converged || !fc$converged) return(NULL)
      data.frame(
        break1=sr$year[j], break2=sr$year[i1],
        wrss = fa$wrss + fb$wrss + fc$wrss,
        a1=fa$a,b1=fa$b,h1=fa$h,
        a2=fb$a,b2=fb$b,h2=fb$h,
        a3=fc$a,b3=fc$b,h3=fc$h,
        p1=sprintf("%d\u2013%d",sr$year[1],sr$year[j]),
        p2=sprintf("%d\u2013%d",sr$year[j+1],sr$year[i1]),
        p3=sprintf("%d\u2013%d",sr$year[i1+1],sr$year[n_sr])
      )
    }))
  } else NULL

  # Second break in segment 2 (after first break)
  seg2b <- if ((n_sr - i1) >= 2 * min_seg_sr) {
    do.call(rbind, lapply((i1 + min_seg_sr):(n_sr - min_seg_sr), function(j) {
      fa <- fit_bh(sr[seq_len(i1), ])
      fb <- fit_bh(sr[(i1 + 1):j, ])
      fc <- fit_bh(sr[(j + 1):n_sr, ])
      if (!fa$converged || !fb$converged || !fc$converged) return(NULL)
      data.frame(
        break1=sr$year[i1], break2=sr$year[j],
        wrss = fa$wrss + fb$wrss + fc$wrss,
        a1=fa$a,b1=fa$b,h1=fa$h,
        a2=fb$a,b2=fb$b,h2=fb$h,
        a3=fc$a,b3=fc$b,h3=fc$h,
        p1=sprintf("%d\u2013%d",sr$year[1],sr$year[i1]),
        p2=sprintf("%d\u2013%d",sr$year[i1+1],sr$year[j]),
        p3=sprintf("%d\u2013%d",sr$year[j+1],sr$year[n_sr])
      )
    }))
  } else NULL

  all2 <- rbind(seg2a, seg2b)
  if (!is.null(all2) && nrow(all2) > 0) {
    best2   <- all2[which.min(all2$wrss), ]
    f2stat  <- ((best1$wrss - best2$wrss) / 2) / (best2$wrss / (n_sr - 6))
    p2_val  <- pf(f2stat, 2, n_sr - 6, lower.tail = FALSE)
    do_3panel <- isTRUE(p2_val < 0.05)
    cat(sprintf("Second break: years %d & %d | F=%.2f | p=%.4f\n",
                best2$break1, best2$break2, f2stat, p2_val))
  }
}

# ------------------------------------------------------------------
# 4. Build plots
# ------------------------------------------------------------------
pal2 <- c("steelblue", "tomato")
pal3 <- c("steelblue", "tomato", "forestgreen")

# --- Panel 1: full series ---
crv1 <- data.frame(SSB = ssb_grid, R = bh_pred(fit_all$a, fit_all$b, ssb_grid))

p_sr1 <- ggplot(sr, aes(x = SSB, y = R)) +
  geom_errorbar(aes(ymin = R_lo, ymax = R_hi),
                colour = "grey60", linewidth = 0.35, width = 0) +
  geom_errorbarh(aes(xmin = SSB_lo, xmax = SSB_hi),
                 colour = "grey60", linewidth = 0.35, height = 0) +
  geom_point(aes(colour = year), size = 2.5) +
  geom_line(data = crv1, aes(x = SSB, y = R),
            colour = "black", linewidth = 1.1, inherit.aes = FALSE) +
  scale_colour_viridis_c(name = "Year") +
  annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.5,
           label = sprintf("h = %.3f  (%d\u2013%d)",
                           fit_all$h, min(sr$year), max(sr$year)),
           size = 4, fontface = "bold") +
  labs(title = "Full time series", x = "SSB", y = "Recruitment (age 1)") +
  theme_bw()

# --- Panel 2: one break ---
p_sr2 <- NULL
if (do_2panel) {
  sr$period2 <- factor(
    ifelse(sr$year <= best1$break_year, best1$p1, best1$p2),
    levels = c(best1$p1, best1$p2))

  leg2 <- setNames(
    c(sprintf("%s  (h=%.3f)", best1$p1, best1$h1),
      sprintf("%s  (h=%.3f)", best1$p2, best1$h2)),
    c(best1$p1, best1$p2))

  crv2 <- rbind(
    data.frame(SSB=ssb_grid, R=bh_pred(best1$a1,best1$b1,ssb_grid), period=best1$p1),
    data.frame(SSB=ssb_grid, R=bh_pred(best1$a2,best1$b2,ssb_grid), period=best1$p2))
  crv2$period <- factor(crv2$period, levels = c(best1$p1, best1$p2))

  p_sr2 <- ggplot(sr, aes(x = SSB, y = R)) +
    geom_errorbar(aes(ymin = R_lo, ymax = R_hi),
                  colour = "grey70", linewidth = 0.35, width = 0) +
    geom_errorbarh(aes(xmin = SSB_lo, xmax = SSB_hi),
                   colour = "grey70", linewidth = 0.35, height = 0) +
    geom_point(aes(colour = period2), size = 2.5) +
    geom_line(data = crv2, aes(x = SSB, y = R, colour = period),
              linewidth = 1.1, inherit.aes = FALSE) +
    scale_colour_manual(values = setNames(pal2, c(best1$p1, best1$p2)),
                        labels = leg2, name = "Period") +
    annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.5,
             label = sprintf("Break: %d  (p=%.3f)", best1$break_year + 1, p1_val),
             size = 3.5) +
    labs(title = "Single productivity break", x = "SSB", y = "Recruitment (age 1)") +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 8))
}

# --- Panel 3: two breaks ---
p_sr3 <- NULL
if (do_3panel) {
  sr$period3 <- factor(
    ifelse(sr$year <= best2$break1, best2$p1,
           ifelse(sr$year <= best2$break2, best2$p2, best2$p3)),
    levels = c(best2$p1, best2$p2, best2$p3))

  leg3 <- setNames(
    c(sprintf("%s  (h=%.3f)", best2$p1, best2$h1),
      sprintf("%s  (h=%.3f)", best2$p2, best2$h2),
      sprintf("%s  (h=%.3f)", best2$p3, best2$h3)),
    c(best2$p1, best2$p2, best2$p3))

  crv3 <- rbind(
    data.frame(SSB=ssb_grid, R=bh_pred(best2$a1,best2$b1,ssb_grid), period=best2$p1),
    data.frame(SSB=ssb_grid, R=bh_pred(best2$a2,best2$b2,ssb_grid), period=best2$p2),
    data.frame(SSB=ssb_grid, R=bh_pred(best2$a3,best2$b3,ssb_grid), period=best2$p3))
  crv3$period <- factor(crv3$period, levels=c(best2$p1, best2$p2, best2$p3))

  p_sr3 <- ggplot(sr, aes(x = SSB, y = R)) +
    geom_errorbar(aes(ymin = R_lo, ymax = R_hi),
                  colour = "grey70", linewidth = 0.35, width = 0) +
    geom_errorbarh(aes(xmin = SSB_lo, xmax = SSB_hi),
                   colour = "grey70", linewidth = 0.35, height = 0) +
    geom_point(aes(colour = period3), size = 2.5) +
    geom_line(data = crv3, aes(x = SSB, y = R, colour = period),
              linewidth = 1.1, inherit.aes = FALSE) +
    scale_colour_manual(values = setNames(pal3, c(best2$p1,best2$p2,best2$p3)),
                        labels = leg3, name = "Period") +
    annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.5,
             label = sprintf("Breaks: %d & %d  (p=%.3f)",
                             best2$break1+1, best2$break2+1, p2_val),
             size = 3.5) +
    labs(title = "Two productivity breaks", x = "SSB", y = "Recruitment (age 1)") +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 8))
}

# ------------------------------------------------------------------
# 5. Productivity over time: R/SSB with CI and smoother
# CI propagated via delta method on log scale:
#   log(R/SSB) SE ≈ sqrt(SE_logR^2 + SE_logSSB^2)
# ------------------------------------------------------------------
sr$logSSB_se <- pmax((sr$SSB_hi - sr$SSB_lo) / (4 * sr$SSB), 0.01)
sr$prod      <- sr$R / sr$SSB
sr$log_prod_se <- sqrt(sr$logR_se^2 + sr$logSSB_se^2)
sr$prod_lo   <- sr$prod * exp(-1.96 * sr$log_prod_se)
sr$prod_hi   <- sr$prod * exp( 1.96 * sr$log_prod_se)

png("diagnostics/h1_productivity_time.png", width = 1400, height = 800, res = 150)
print(
  ggplot(sr, aes(x = year, y = prod)) +
    geom_errorbar(aes(ymin = prod_lo, ymax = prod_hi),
                  colour = "grey65", linewidth = 0.4, width = 0) +
    geom_point(aes(colour = prod), size = 2.5) +
    geom_smooth(method = "loess", span = 0.4, se = FALSE,
                colour = "black", linewidth = 1.1) +
    scale_colour_viridis_c(option = "plasma", name = "R/SSB") +
    scale_y_continuous(limits = c(0, NA)) +
    labs(title = "Productivity over time (R / SSB)",
         subtitle = "Points = annual estimate with 95% CI; line = LOESS smoother",
         x = "Year", y = "Recruitment per Spawner") +
    theme_bw() +
    theme(legend.position = "right")
)
dev.off()

# ------------------------------------------------------------------
# 6. Combine B-H panels with patchwork and save
# ------------------------------------------------------------------
panel_plots <- Filter(Negate(is.null), list(p_sr1, p_sr2, p_sr3))
n_panels    <- length(panel_plots)

png("diagnostics/h1_productivity_sr.png",
    width = n_panels * 900, height = 900, res = 150)
print(wrap_plots(panel_plots, nrow = 1))
dev.off()

# Save text summary
sink("diagnostics/h1_productivity_summary.txt")
cat(sprintf("Beverton-Holt productivity analysis\n"))
cat(sprintf("B0 reference (max observed SSB): %.0f\n\n", B0_ref))
cat(sprintf("Full series: h = %.3f (%d-%d)\n",
            fit_all$h, min(sr$year), max(sr$year)))
if (!is.na(p1_val))
  cat(sprintf("1-break test: break year %d, p = %.4f %s\n",
              best1$break_year+1, p1_val, ifelse(do_2panel,"[significant]","[not significant]")))
if (do_2panel)
  cat(sprintf("  Period 1 (%s): h = %.3f\n  Period 2 (%s): h = %.3f\n",
              best1$p1, best1$h1, best1$p2, best1$h2))
if (!is.na(p2_val))
  cat(sprintf("2-break test: breaks %d & %d, p = %.4f %s\n",
              best2$break1+1, best2$break2+1, p2_val,
              ifelse(do_3panel,"[significant]","[not significant]")))
if (do_3panel)
  cat(sprintf("  Period 1 (%s): h = %.3f\n  Period 2 (%s): h = %.3f\n  Period 3 (%s): h = %.3f\n",
              best2$p1, best2$h1, best2$p2, best2$h2, best2$p3, best2$h3))
sink()

cat("Productivity analysis written to diagnostics/\n")

################################################################################
## SELECTIVITY BLOCKING ANALYSIS FOR JJM
##
## Rationale:
##   JJM assumes selectivity is constant within user-defined year blocks.
##   SAM estimates selectivity annually (as F/Fbar at age) using a random walk,
##   so the SAM trajectory reveals when the selectivity *shape* changed
##   substantially enough to warrant a new block in JJM.
##
## Method:
##   1. For each catch fleet, extract the normalised selectivity-at-age vector
##      per year (shape only, peak = 1) and the total catch (weight) per year.
##   2. Reduce to one dimension via PCA (PC1 captures the dominant mode of
##      selectivity change over time).
##   3. Apply CATCH-WEIGHTED recursive binary segmentation on PC1: the RSS at
##      each candidate split uses catch as observation weights so that years
##      with negligible catch have little influence on where breaks are placed.
##      An F-test (p < 0.05, min 5 years per segment) gates each split.
##   4. Rank breakpoints by importance = F-statistic × fraction of total fleet
##      catch in the split segment. Breaks in high-catch periods rank highest;
##      label 1 is most important and should be retained last when reducing
##      the number of blocks.
##
## Output:
##   diagnostics/h1_selectivity_blocks.csv      — block suggestion table
##   diagnostics/h1_selectivity_pc1.png         — PC1 + ranked breakpoints
##   diagnostics/h1_selectivity_byblock.png     — selectivity curves by block
##   diagnostics/h1_selectivity_catch_weight.png — catch weights used
################################################################################

# ------------------------------------------------------------------
# Step 1 – Total catch per fleet per year (used as weights)
# ------------------------------------------------------------------
catch_wt_df <- as.data.frame(stk_h1@catch.n * stk_h1@catch.wt)
catch_wt_df <- aggregate(data ~ year + area, data = catch_wt_df, FUN = sum, na.rm = TRUE)
catch_wt_df$data[catch_wt_df$data <= 0] <- NA

# Map area dimname to fleet name positionally
area_names <- dimnames(stk_h1@catch.n)$area
catch_wt_df$fleet <- catch_fleets[match(as.character(catch_wt_df$area), area_names)]

# Normalise within fleet: weight = catch / max(catch) so max weight = 1
catch_wt_df <- do.call(rbind, lapply(split(catch_wt_df, catch_wt_df$fleet), function(df) {
  mx        <- max(df$data, na.rm = TRUE)
  df$weight <- ifelse(is.na(df$data), 0, df$data / mx)
  df
}))

# ------------------------------------------------------------------
# Step 2 – Weighted binary segmentation helpers
# ------------------------------------------------------------------
min_seg <- 5
alpha   <- 0.05

find_cp_w <- function(x, w) {
  n  <- length(x)
  w  <- pmax(w, 1e-6)           # avoid zero-weight years collapsing segment
  w  <- w / sum(w)
  if (n < 2 * min_seg) return(list(cp = NA, p = 1, f = 0))

  wmean <- function(xi, wi) sum(wi * xi) / sum(wi)
  wrss  <- function(xi, wi) sum(wi * (xi - wmean(xi, wi))^2)

  rss0    <- wrss(x, w)
  cands   <- min_seg:(n - min_seg)
  rss_vec <- sapply(cands, function(i)
    wrss(x[seq_len(i)], w[seq_len(i)]) +
    wrss(x[(i + 1):n],  w[(i + 1):n]))

  best_pos <- which.min(rss_vec)
  best_i   <- cands[best_pos]
  best_rss <- rss_vec[best_pos]
  f_stat   <- (rss0 - best_rss) / (best_rss / (n - 2))
  p_val    <- pf(f_stat, 1, n - 2, lower.tail = FALSE)
  list(cp = best_i, p = p_val, f = f_stat)
}

# Returns data frame: idx_pos (index into yrs), f_stat, catch_sum (raw catch
# in split segment — used later for importance weighting)
bin_seg_w <- function(x, w, idx, catch_raw, depth = 0) {
  empty <- data.frame(idx_pos = integer(0), f_stat = numeric(0),
                      catch_sum = numeric(0))
  if (depth > 4 || length(x) < 2 * min_seg) return(empty)
  res <- find_cp_w(x, w)
  if (is.na(res$cp) || res$p > alpha)        return(empty)
  cp_global <- idx[res$cp]
  rbind(
    data.frame(idx_pos   = cp_global,
               f_stat    = res$f,
               catch_sum = sum(catch_raw, na.rm = TRUE)),
    bin_seg_w(x[seq_len(res$cp)],        w[seq_len(res$cp)],
              idx[seq_len(res$cp)],      catch_raw[seq_len(res$cp)],      depth + 1),
    bin_seg_w(x[(res$cp + 1):length(x)], w[(res$cp + 1):length(x)],
              idx[(res$cp + 1):length(x)], catch_raw[(res$cp + 1):length(x)], depth + 1)
  )
}

# ------------------------------------------------------------------
# Step 3 – Run per fleet
# ------------------------------------------------------------------
sel_block_list <- lapply(catch_fleets, function(fl) {

  df <- subset(sel.pat, fleet == fl & is.finite(sel) & sel >= 0)
  df <- df[order(df$year, df$age), ]

  sel_wide <- reshape(df[, c("year", "age", "sel")],
                      idvar = "year", timevar = "age", direction = "wide")
  yrs     <- sel_wide$year
  sel_mat <- as.matrix(sel_wide[, -1])
  rownames(sel_mat) <- yrs

  # Normalise shape: peak = 1 per year
  sel_norm <- t(apply(sel_mat, 1, function(x) {
    x[!is.finite(x)] <- NA
    mx <- max(x, na.rm = TRUE)
    if (is.finite(mx) && mx > 0) x / mx else x
  }))
  ok_col   <- colSums(!is.na(sel_norm)) > 1
  sel_norm <- sel_norm[, ok_col, drop = FALSE]

  n <- nrow(sel_norm)
  if (n < 2 * min_seg) {
    message("Fleet ", fl, ": too few years for blocking analysis")
    return(NULL)
  }

  # PCA
  pca  <- prcomp(sel_norm, center = TRUE, scale. = FALSE)
  pc1  <- pca$x[, 1]
  var1 <- round(summary(pca)$importance[2, 1] * 100, 1)

  # Match catch weights to years in sel (some years may be missing)
  wt_fleet   <- subset(catch_wt_df, fleet == fl)
  w_vec      <- wt_fleet$weight[match(yrs, wt_fleet$year)]
  w_vec[is.na(w_vec)] <- 0
  catch_raw  <- wt_fleet$data[match(yrs, wt_fleet$year)]
  catch_raw[is.na(catch_raw)] <- 0

  # Total fleet catch — used to compute importance fractions
  total_catch <- sum(catch_raw, na.rm = TRUE)

  # Weighted binary segmentation
  cps_df <- bin_seg_w(pc1, w_vec, seq_along(yrs), catch_raw)

  if (nrow(cps_df) == 0) {
    break_yrs  <- integer(0)
    importance <- numeric(0)
  } else {
    # Importance = F-stat × (catch in split segment / total fleet catch)
    # Higher = more important; rank 1 = keep last when reducing blocks
    cps_df$importance <- cps_df$f_stat * (cps_df$catch_sum / total_catch)
    cps_df <- cps_df[order(cps_df$idx_pos), ]
    # Global rank across all breaks for this fleet (1 = most important)
    cps_df$rank <- rank(-cps_df$importance, ties.method = "first")
    break_yrs  <- yrs[cps_df$idx_pos]
    cps_df$break_year <- break_yrs
    importance <- cps_df$importance
  }

  # Build block table
  starts <- c(min(yrs), break_yrs + 1L)
  ends   <- c(break_yrs, max(yrs))

  block_df <- do.call(rbind, lapply(seq_along(starts), function(b) {
    yr_b    <- yrs[yrs >= starts[b] & yrs <= ends[b]]
    mat_b   <- sel_norm[rownames(sel_norm) %in% yr_b, , drop = FALSE]
    mn      <- colMeans(mat_b, na.rm = TRUE)
    wb_var  <- mean(apply(mat_b, 2, var, na.rm = TRUE), na.rm = TRUE)
    pk_age  <- as.integer(gsub("sel\\.", "", names(which.max(mn))))
    ct_b    <- sum(catch_raw[yrs %in% yr_b], na.rm = TRUE)
    ct_frac <- round(ct_b / total_catch * 100, 1)

    # Breakpoint rank that *ends* this block (NA for last block)
    bp_rank <- if (b < length(starts)) {
      cps_df$rank[cps_df$break_year == ends[b]]
    } else NA_integer_

    prev_pk <- if (b > 1) {
      yr_p  <- yrs[yrs >= starts[b - 1] & yrs <= ends[b - 1]]
      mat_p <- sel_norm[rownames(sel_norm) %in% yr_p, , drop = FALSE]
      mn_p  <- colMeans(mat_p, na.rm = TRUE)
      as.integer(gsub("sel\\.", "", names(which.max(mn_p))))
    } else NA_integer_

    rationale <- if (b == 1) {
      "Initial block"
    } else if (!is.na(prev_pk) && abs(pk_age - prev_pk) >= 2) {
      sprintf("Peak age shifted %i → %i (%.1f%% of catch)", prev_pk, pk_age, ct_frac)
    } else {
      sprintf("Selectivity shape change (%.1f%% of fleet catch in block)", ct_frac)
    }

    data.frame(fleet         = fl,
               block         = b,
               start_year    = starts[b],
               end_year      = ends[b],
               n_years       = length(yr_b),
               catch_pct     = ct_frac,
               peak_age      = pk_age,
               within_var    = round(wb_var, 5),
               pc1_pct_var   = var1,
               boundary_rank = bp_rank,
               rationale     = rationale,
               stringsAsFactors = FALSE)
  }))

  list(blocks    = block_df,
       pc1_df    = data.frame(year = yrs, pc1 = pc1, fleet = fl,
                              catch_weight = w_vec),
       cps_df    = if (nrow(cps_df) > 0) cps_df else NULL,
       break_yrs = break_yrs,
       sel_norm  = data.frame(year = yrs, sel_norm),
       fl        = fl)
})
names(sel_block_list) <- catch_fleets
sel_block_list        <- Filter(Negate(is.null), sel_block_list)

# ------------------------------------------------------------------
# Step 4 – Combined block table
# ------------------------------------------------------------------
sel_blocks_all <- do.call(rbind, lapply(sel_block_list, `[[`, "blocks"))
rownames(sel_blocks_all) <- NULL
write.csv(sel_blocks_all, "diagnostics/h1_selectivity_blocks.csv", row.names = FALSE)

cat("\n--- Suggested JJM selectivity blocks ---\n")
print(sel_blocks_all[, c("fleet", "block", "start_year", "end_year",
                          "n_years", "catch_pct", "peak_age",
                          "within_var", "boundary_rank", "rationale")])

# ------------------------------------------------------------------
# Step 5 – Plot: PC1 with ranked breakpoints and catch weight ribbon
# ------------------------------------------------------------------
pc1_all <- do.call(rbind, lapply(sel_block_list, `[[`, "pc1_df"))

bp_all <- do.call(rbind, lapply(sel_block_list, function(x) {
  if (is.null(x$cps_df) || nrow(x$cps_df) == 0) return(NULL)
  data.frame(fleet      = x$fl,
             break_year = x$cps_df$break_year,
             rank       = x$cps_df$rank,
             importance = x$cps_df$importance)
}))

# y-range per fleet for placing rank labels
pc1_range <- do.call(rbind, lapply(split(pc1_all, pc1_all$fleet), function(df) {
  data.frame(fleet = df$fleet[1],
             ymin  = min(df$pc1, na.rm = TRUE),
             ymax  = max(df$pc1, na.rm = TRUE))
}))

if (!is.null(bp_all) && nrow(bp_all) > 0) {
  bp_all <- merge(bp_all, pc1_range, by = "fleet")
  bp_all$label_y <- bp_all$ymax   # place rank label at top of panel
}

# Pre-compute per-fleet ribbon bounds: catch weight ribbon fills the bottom
# 35% of each panel's y-range, so it is visible regardless of free_y scaling.
pc1_all <- do.call(rbind, lapply(split(pc1_all, pc1_all$fleet), function(df) {
  y_min   <- min(df$pc1, na.rm = TRUE)
  y_range <- max(df$pc1, na.rm = TRUE) - y_min
  df$ribbon_lo <- y_min
  df$ribbon_hi <- y_min + df$catch_weight * y_range * 0.35
  df
}))

png("diagnostics/h1_selectivity_pc1.png", width = 1800, height = 1200, res = 150)
p_pc1 <- ggplot(pc1_all, aes(x = year, y = pc1)) +
  # Catch weight ribbon anchored to the bottom of each panel's y-range
  geom_ribbon(aes(ymin = ribbon_lo, ymax = ribbon_hi),
              fill = "grey75", alpha = 0.6) +
  geom_line(colour = "steelblue", linewidth = 0.9) +
  geom_point(size = 1.2, colour = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  facet_wrap(~fleet, scales = "free_y", ncol = 2) +
  labs(title = "PC1 of annual selectivity shape by fleet",
       subtitle = paste("Vertical lines = breakpoints; number = importance rank",
                        "(1 = most important, keep last when reducing blocks).",
                        "Grey shading = relative catch weight."),
       x = "Year", y = "PC1 score") +
  theme_bw()

if (!is.null(bp_all) && nrow(bp_all) > 0) {
  p_pc1 <- p_pc1 +
    geom_vline(data = bp_all,
               aes(xintercept = break_year + 0.5),
               colour = "tomato", linewidth = 0.9) +
    geom_label(data = bp_all,
               aes(x = break_year + 0.5, y = label_y, label = rank),
               colour = "tomato", fill = "white", size = 3.5,
               label.padding = unit(0.15, "lines"), inherit.aes = FALSE)
}
print(p_pc1)
dev.off()

# ------------------------------------------------------------------
# Step 6 – Plot: selectivity curves coloured by block
# ------------------------------------------------------------------
sel_with_block <- do.call(rbind, lapply(sel_block_list, function(x) {
  bd <- x$blocks
  df <- subset(sel.pat, fleet == x$fl & is.finite(sel))
  df$block <- NA_integer_
  for (b in seq_len(nrow(bd))) {
    df$block[df$year >= bd$start_year[b] & df$year <= bd$end_year[b]] <- b
  }
  df$block_label <- sprintf("Block %i (%i–%i)", df$block,
                             bd$start_year[df$block], bd$end_year[df$block])
  df
}))

png("diagnostics/h1_selectivity_byblock.png", width = 2000, height = 1400, res = 150)
print(
  ggplot(sel_with_block, aes(x = age, y = sel,
                              group = year, colour = factor(block))) +
    geom_line(alpha = 0.5, linewidth = 0.55) +
    stat_summary(aes(group = factor(block)), fun = mean,
                 geom = "line", linewidth = 1.6, linetype = "dashed") +
    scale_colour_brewer(palette = "Set1", name = "Block") +
    facet_wrap(~fleet, scales = "free_y", ncol = 2) +
    scale_x_continuous(breaks = seq(1, 12, by = 1)) +
    labs(title = "Annual selectivity curves by fleet, coloured by block",
         subtitle = "Thin lines = individual years; dashed = block mean",
         x = "Age", y = "F / Fbar") +
    theme_bw()
)
dev.off()

# ------------------------------------------------------------------
# Step 7 – Plot: catch weight timeseries per fleet
# ------------------------------------------------------------------
png("diagnostics/h1_selectivity_catch_weight.png", width = 1800, height = 900, res = 150)
print(
  ggplot(pc1_all, aes(x = year, y = catch_weight)) +
    geom_area(fill = "steelblue", alpha = 0.5) +
    geom_line(colour = "steelblue", linewidth = 0.7) +
    facet_wrap(~fleet, scales = "free_y", ncol = 2) +
    labs(title = "Relative catch weight used in selectivity block analysis",
         subtitle = "Weight = annual catch / fleet maximum catch. Low-catch years have little influence on block boundaries.",
         x = "Year", y = "Relative catch weight") +
    theme_bw()
)
dev.off()

cat("Selectivity blocking output written to diagnostics/\n")

Sys.setenv(QUARTO_PATH = "C:/Program Files/Quarto/bin/quarto.exe")  # adjust to whichever path was found
quarto::quarto_render("sam_diagnostics.qmd")
