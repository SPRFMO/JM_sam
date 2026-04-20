# ============================================================
# 03_retro_diagnostics.R
# Retrospective analysis, Mohn's rho, leave-one-out plots,
# and in-depth retrospective parameter / residual diagnostics.
# Sourced from 01_run_assessment.R; requires:
#   sam_h1, retro_h1, loo_fits, res_all, ssb_df
#   catch_fleets, survey_fleets, survey_fleets_base, base_fleet
#   max_yr, base_yr_chr
# NOTE: retro peels are run with ctrl@residuals = TRUE so that
#       residual-based stability diagnostics are available.
# ============================================================

# If run standalone (not via 01_run_assessment.R), parse jjm_ssb here.
# When sourced from 01_run_assessment.R the object already exists.
if (!exists("jjm_ssb")) {
  .parse_jjm_ssb <- function(rep_file, max_year) {
    if (!file.exists(rep_file)) return(NULL)
    lines     <- readLines(rep_file, warn = FALSE)
    sec_idx   <- grep("^\\$[A-Za-z_]+$", lines)
    sec_names <- sub("^\\$", "", lines[sec_idx])
    i_ssb <- sec_idx[sec_names == "SSB"]
    if (length(i_ssb) == 0) return(NULL)
    pos  <- which(sec_idx == i_ssb)
    end  <- if (pos < length(sec_idx)) sec_idx[pos + 1] - 1 else length(lines)
    dlines <- trimws(lines[(i_ssb + 1):end])
    dlines <- dlines[nchar(dlines) > 0]
    mat    <- do.call(rbind, lapply(dlines,
                function(l) as.numeric(strsplit(l, "\\s+")[[1]])))
    if (ncol(mat) < 5) return(NULL)
    mat <- mat[mat[, 1] <= max_year, , drop = FALSE]
    data.frame(year = as.integer(mat[, 1]), value = mat[, 2],
               lbnd = mat[, 4], ubnd = mat[, 5])
  }
  jjm_ssb <- .parse_jjm_ssb(
    rep_file = "../assessment/results/h1_1.14_1_R.rep",
    max_year = as.integer(range(stk_h1)["maxyear"])
  )
}

##################################################
## Retrospective - SSB and Fbar
##################################################
retro_ssb  <- do.call(rbind, lapply(names(retro_h1), function(nm) {
  fit <- retro_h1[[nm]]
  if (is.null(fit) || inherits(fit, "try-error")) return(NULL)
  d <- ssb(fit); d$peel <- nm; d
}))
retro_fbar <- do.call(rbind, lapply(names(retro_h1), function(nm) {
  fit <- retro_h1[[nm]]
  if (is.null(fit) || inherits(fit, "try-error")) return(NULL)
  d <- fbar(fit); d$peel <- nm; d
}))
retro_ssb$peel  <- factor(retro_ssb$peel,  levels = names(retro_h1))
retro_fbar$peel <- factor(retro_fbar$peel, levels = names(retro_h1))

base_ssb  <- ssb(sam_h1);  base_ssb$peel  <- base_yr_chr
base_fbar <- fbar(sam_h1); base_fbar$peel <- base_yr_chr

retro_ssb_all  <- rbind(retro_ssb,  base_ssb)
retro_fbar_all <- rbind(retro_fbar, base_fbar)

p_retro_ssb <-
  ggplot(retro_ssb_all, aes(x = year, y = value, colour = peel, group = peel)) +
    geom_line(data = subset(retro_ssb_all, peel != base_yr_chr), alpha = 0.7) +
    geom_line(data = subset(retro_ssb_all, peel == base_yr_chr),
              linewidth = 1.2, colour = "black") +
    geom_ribbon(data = subset(retro_ssb_all, peel == base_yr_chr),
                aes(ymin = lbnd, ymax = ubnd), alpha = 0.2, colour = NA, fill = "black") +
    coord_cartesian(xlim = c(max(retro_ssb_all$year) - 15, max(retro_ssb_all$year))) +
    ylim(0, NA) +
    labs(title = "Retrospective \u2014 SSB",
         x = "Year", y = "SSB", colour = "Final year") +
    theme_bw()

# Overlay JJM h1_1.14 SSB as a dashed green reference line
if (!is.null(jjm_ssb)) {
  jjm_retro_sub <- subset(jjm_ssb,
    year >= max(retro_ssb_all$year) - 15 & year <= max(retro_ssb_all$year))
  p_retro_ssb <- p_retro_ssb +
    geom_ribbon(data = jjm_retro_sub,
                aes(x = year, ymin = lbnd, ymax = ubnd),
                fill = "forestgreen", alpha = 0.12, colour = NA,
                inherit.aes = FALSE) +
    geom_line(data = jjm_retro_sub,
              aes(x = year, y = value),
              colour = "forestgreen", linetype = "dashed", linewidth = 1.0,
              inherit.aes = FALSE) +
    labs(caption = "Dashed green = JJM h1_1.14 SSB \u00b1 95% CI")
}

png("diagnostics/h1_retro_ssb.png", width = 1600, height = 900, res = 150)
print(p_retro_ssb)
dev.off()

png("diagnostics/h1_retro_fbar.png", width = 1600, height = 900, res = 150)
print(
  ggplot(retro_fbar_all, aes(x = year, y = value, colour = peel, group = peel)) +
    geom_line(data = subset(retro_fbar_all, peel != base_yr_chr), alpha = 0.7) +
    geom_line(data = subset(retro_fbar_all, peel == base_yr_chr),
              linewidth = 1.2, colour = "black") +
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
mohn_rho <- function(base_df, retro_df) {
  peels <- setdiff(unique(retro_df$peel), base_yr_chr)
  rhos  <- sapply(peels, function(p) {
    sub       <- subset(retro_df, peel == p)
    yr        <- max(sub$year)
    base_val  <- base_df$value[base_df$year == yr]
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
    fit <- fits[[grp]]
    if (is.null(fit) || inherits(fit, "try-error")) return(NULL)
    df         <- fun(fit)
    df$dropped <- grp
    df
  }))
  rbind(base_df, loo_df)
}

loo_ssb_all  <- extract_loo(ssb,  loo_fits, sam_h1)
loo_fbar_all <- extract_loo(fbar, loo_fits, sam_h1)
loo_rec_all  <- extract_loo(rec,  loo_fits, sam_h1)

plot_loo <- function(dat, ylab, title, ref_df = NULL) {
  # All LOO lines overlaid; no CI ribbons; full model thicker black line.
  # Colour is mapped to `dropped` so each LOO run gets a distinct hue.
  base  <- subset(dat, dropped == "Full model")
  loodf <- subset(dat, dropped != "Full model")

  # Okabe-Ito palette: 7 colorblind-safe, high-contrast hues.
  # Black is reserved for the full-model line so is excluded here.
  okabe_ito <- c("#E69F00", "#56B4E9", "#009E73",
                 "#D55E00", "#CC79A7", "#0072B2", "#F0E442")
  grps <- unique(loodf$dropped)
  pal  <- setNames(okabe_ito[seq_along(grps)], grps)

  p <- ggplot() +
    # LOO lines — coloured by dropped survey
    geom_line(data = loodf,
              aes(x = year, y = value, colour = dropped),
              linewidth = 0.9, alpha = 0.9) +
    # Full model on top — thicker black
    geom_line(data = base,
              aes(x = year, y = value),
              colour = "black", linewidth = 1.5) +
    scale_colour_manual(values = pal, name = "Dropped survey") +
    ylim(0, NA) +
    labs(title = title,
         subtitle = "Black = full model; coloured lines = model with that survey group dropped",
         x = "Year", y = ylab) +
    theme_bw() +
    theme(legend.position = "right")

  # Optional JJM reference line
  if (!is.null(ref_df)) {
    p <- p +
      geom_line(data = ref_df,
                aes(x = year, y = value),
                colour = "forestgreen", linetype = "dashed", linewidth = 1.0,
                inherit.aes = FALSE) +
      labs(caption = "Dashed green = JJM h1_1.14 SSB")
  }
  p
}

png("diagnostics/h1_loo_ssb.png",  width = 1600, height = 900, res = 150)
print(plot_loo(loo_ssb_all, "SSB", "Leave-one-out \u2014 SSB",
               ref_df = jjm_ssb))
dev.off()

png("diagnostics/h1_loo_fbar.png", width = 1600, height = 900, res = 150)
print(plot_loo(loo_fbar_all, "Mean F",      "Leave-one-out \u2014 Fbar"))
dev.off()

png("diagnostics/h1_loo_rec.png",  width = 1600, height = 900, res = 150)
print(plot_loo(loo_rec_all,  "Recruitment", "Leave-one-out \u2014 Recruitment"))
dev.off()

################################################################################
## IN-DEPTH RETROSPECTIVE DIAGNOSTICS
##
##  1. retroParams — catchability (log-q) vs peel year
##     Catchability is re-estimated in each peel. If log-q shifts substantially
##     as years are added, that survey's catchability is poorly identified or
##     driven by the most recent observations.
##
##  2. retroParams — observation variance (sigma) vs peel year
##     Observation variance parameters should be stable across peels if the
##     data are consistent over time. Drift signals changing noise levels.
##
##  3. Residual spaghetti across peels (per survey)
##     For each survey, overlay the standardised residuals from every peel that
##     covered that year. Lines should be near-identical in the middle of the
##     series and may diverge only near each peel's terminal year. Systematic
##     divergence throughout indicates that the survey is in conflict with the
##     stock trajectory implied by other fleets.
##
##  4. Residual stability heatmap
##     For each survey-year combination, the SD of std.residuals across all
##     retro peels that included that year. High SD = the model re-interprets
##     that observation substantially as more data are added, i.e. that year
##     is a high-influence point.
################################################################################

# ------------------------------------------------------------------
# 1 & 2. retroParams: estimated catchability and obs.var per peel
# ------------------------------------------------------------------

# -- Catchabilities --
retro_q <- do.call(rbind, lapply(names(retro_h1), function(nm) {
  q <- tryCatch(catchabilities(retro_h1[[nm]]), error = function(e) NULL)
  if (is.null(q) || nrow(q) == 0) return(NULL)
  q$peel    <- nm
  q$peel_yr <- as.integer(nm)
  q$base    <- base_fleet(q$fleet)
  q
}))

if (!is.null(retro_q) && nrow(retro_q) > 0) {
  retro_q <- retro_q[order(retro_q$peel_yr), ]

  # Colour by period (early / late / single) — at most 2 lines per facet,
  # so 3 colours suffice regardless of how many fleets there are in total.
  retro_q$period <- ifelse(grepl("_early$", retro_q$fleet), "Early",
                    ifelse(grepl("_late$",  retro_q$fleet), "Late", "Single"))

  png("diagnostics/h1_retro_q_params.png", width = 1800, height = 1200, res = 150)
  print(
    ggplot(retro_q, aes(x = peel_yr, y = value,
                        colour = period, group = fleet)) +
      geom_line(linewidth = 0.9) +
      geom_point(size = 2.5) +
      scale_colour_manual(
        values = c("Early" = "steelblue", "Late" = "tomato", "Single" = "grey30"),
        name   = "Period") +
      scale_x_continuous(breaks = unique(retro_q$peel_yr)) +
      facet_wrap(~base, scales = "free_y", ncol = 2) +
      labs(title = "Catchability parameter estimates across retrospective peels",
           subtitle = "Early/late split series shown separately within each panel",
           x = "Terminal year of peel", y = "log(q)") +
      theme_bw() +
      theme(axis.text.x    = element_text(angle = 45, hjust = 1),
            legend.position = "right",
            strip.text      = element_text(face = "bold"))
  )
  dev.off()
}

# -- Observation variances --
retro_obv <- do.call(rbind, lapply(names(retro_h1), function(nm) {
  ov <- tryCatch(obs.var(retro_h1[[nm]]), error = function(e) NULL)
  if (is.null(ov) || nrow(ov) == 0) return(NULL)
  ov$peel    <- nm
  ov$peel_yr <- as.integer(nm)
  ov$base    <- base_fleet(ov$fleet)
  ov$group   <- ifelse(ov$fleet %in% catch_fleets, "Catch fleets", "Survey indices")
  ov$label   <- paste0(sub("catch ", "", ov$base),
                       ifelse(is.na(ov$age), "", paste0(" a", ov$age)))
  ov
}))

if (!is.null(retro_obv) && nrow(retro_obv) > 0) {
  retro_obv <- retro_obv[order(retro_obv$peel_yr), ]

  # Deduplicate: ages sharing the same estimated parameter have identical
  # values within a fleet × peel. Keep one row per unique value and label
  # it with the age range it represents.
  retro_obv_uniq <- do.call(rbind, lapply(
    split(retro_obv, list(retro_obv$fleet, retro_obv$peel_yr), drop = TRUE),
    function(df) {
      seen <- c()
      rows <- list()
      for (i in order(df$age, na.last = FALSE)) {
        v <- round(df$value[i], 10)
        if (v %in% seen) next
        seen <- c(seen, v)
        match_ages <- sort(df$age[round(df$value, 10) == v])
        match_ages <- match_ages[!is.na(match_ages)]
        if (length(match_ages) == 0) {
          age_lbl <- "idx"
        } else if (length(match_ages) == 1) {
          age_lbl <- paste0("a", match_ages)
        } else if (all(diff(match_ages) == 1)) {
          age_lbl <- paste0("a", min(match_ages), "-", max(match_ages))
        } else {
          age_lbl <- paste0("a", paste(match_ages, collapse = ","))
        }
        row              <- df[i, , drop = FALSE]
        row$age_label    <- age_lbl
        row$param_label  <- paste0(sub("catch ", "", row$base),
                                   " (", age_lbl, ")")
        rows[[length(rows) + 1]] <- row
      }
      do.call(rbind, rows)
    }
  ))

  png("diagnostics/h1_retro_obsvar_params.png", width = 1800, height = 1200, res = 150)
  print(
    ggplot(retro_obv_uniq,
           aes(x = peel_yr, y = value,
               colour = param_label, group = param_label)) +
      geom_line(linewidth = 0.9) +
      geom_point(size = 2) +
      scale_x_continuous(breaks = unique(retro_obv_uniq$peel_yr)) +
      scale_colour_discrete(name = NULL) +
      facet_wrap(~group, scales = "free_y", ncol = 1) +
      labs(title = "Observation variance estimates across retrospective peels",
           subtitle = "One line per unique estimated parameter; legend shows which ages share that parameter",
           x = "Terminal year of peel", y = "Observation variance") +
      theme_bw() +
      theme(axis.text.x    = element_text(angle = 45, hjust = 1),
            legend.position = "right",
            strip.text      = element_text(face = "bold"))
  )
  dev.off()
}

# ------------------------------------------------------------------
# 3 & 4. Residual spaghetti and stability
# Requires retro peels run with ctrl@residuals = TRUE
# ------------------------------------------------------------------

retro_res_all <- do.call(rbind, lapply(names(retro_h1), function(nm) {
  r <- tryCatch(residuals(retro_h1[[nm]]), error = function(e) NULL)
  if (is.null(r) || nrow(r) == 0) return(NULL)
  r$peel    <- nm
  r$peel_yr <- as.integer(nm)
  r
}))

if (!is.null(retro_res_all) && nrow(retro_res_all) > 0) {

  # --- 3. Survey residual spaghetti ---
  svy_retro <- subset(retro_res_all, fleet %in% survey_fleets)
  svy_retro$fleet <- base_fleet(svy_retro$fleet)

  png("diagnostics/h1_retro_residual_spaghetti.png", width = 1800, height = 1400, res = 150)
  print(
    ggplot() +
      # Peel lines (all except base), coloured by terminal year
      geom_line(data   = subset(svy_retro, peel != base_yr_chr),
                aes(x = year, y = std.res, group = peel, colour = peel_yr),
                linewidth = 0.55, alpha = 0.7) +
      # Base model on top in black
      geom_line(data   = subset(svy_retro, peel == base_yr_chr),
                aes(x = year, y = std.res),
                colour = "black", linewidth = 1.3) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
      scale_colour_viridis_c(name = "Peel terminal\nyear") +
      facet_wrap(~fleet, scales = "free_y", ncol = 2) +
      labs(title = "Survey residuals across retrospective peels",
           subtitle = "Each line = one retro peel; black = full model. Divergence from black = influence of terminal data.",
           x = "Year", y = "Std. residual") +
      theme_bw()
  )
  dev.off()

  # --- 4. Residual stability heatmap ---
  # For each survey x year: SD of std.res across all peels that included that year.
  # Peels are ordered from most to least data; years near the base terminal year
  # are only present in few peels, so limit to years appearing in >= 2 peels.
  stab_df <- do.call(rbind, lapply(
    split(svy_retro, list(svy_retro$fleet, svy_retro$year)),
    function(df) {
      valid <- df[!is.na(df$std.res), ]
      if (nrow(valid) < 2) return(NULL)
      data.frame(fleet    = df$fleet[1],
                 year     = df$year[1],
                 n_peels  = nrow(valid),
                 mean_res = mean(valid$std.res),
                 sd_res   = sd(valid$std.res))
    }
  ))

  if (!is.null(stab_df) && nrow(stab_df) > 0) {
    stab_df$fleet <- factor(stab_df$fleet, levels = rev(survey_fleets_base))

    png("diagnostics/h1_retro_residual_stability.png", width = 1800, height = 900, res = 150)
    print(
      ggplot(stab_df, aes(x = year, y = fleet, fill = sd_res)) +
        geom_tile(colour = "white", linewidth = 0.3) +
        scale_fill_viridis_c(option = "plasma", name = "SD(std.res)\nacross peels",
                             limits = c(0, NA)) +
        labs(title = "Residual stability across retrospective peels",
             subtitle = "SD of standardised residuals across all peels containing that year.\nHigh values = observation strongly re-interpreted as more data are added.",
             x = "Year", y = NULL) +
        theme_bw() +
        theme(panel.grid = element_blank(),
              axis.text.y = element_text(size = 9))
    )
    dev.off()
  }

  # --- Terminal-year residuals: what is driving the retrospective pattern ---
  # For each retro peel, extract residuals at its terminal year only.
  # These are the observations the model is struggling to reconcile when
  # it has just that much data — the direct cause of the retrospective pattern.
  # For catch fleets (age-structured): aggregate to obs/pred ratio using
  # the geometric mean across ages.  For surveys: one obs per year directly.
  terminal_res <- do.call(rbind, lapply(
    setdiff(names(retro_h1), base_yr_chr),   # exclude full model
    function(nm) {
      peel_yr <- as.integer(nm)
      r <- tryCatch(residuals(retro_h1[[nm]]), error = function(e) NULL)
      if (is.null(r) || nrow(r) == 0) return(NULL)
      term <- r[r$year == peel_yr & !is.na(r$std.res), ]
      if (nrow(term) == 0) return(NULL)
      do.call(rbind, lapply(split(term, term$fleet), function(df) {
        data.frame(
          fleet       = df$fleet[1],
          peel_yr     = peel_yr,
          # obs/pred ratio: exp(mean(log.obs)) / exp(mean(log.mdl))
          # = geometric mean ratio, works for both surveys and catch-at-age
          bio_obs     = exp(mean(df$log.obs, na.rm = TRUE)),
          bio_pred    = exp(mean(df$log.mdl, na.rm = TRUE)),
          mean_stdres = mean(df$std.res,    na.rm = TRUE)
        )
      }))
    }
  ))

  if (!is.null(terminal_res) && nrow(terminal_res) > 0) {
    # Use full fleet name for faceting — base_fleet() would merge early/late
    # into one panel causing their points to overlap or stack incorrectly.
    terminal_res$fleet_label <- sub("^catch ", "", terminal_res$fleet)
    terminal_res$type        <- ifelse(terminal_res$fleet %in% catch_fleets,
                                       "Catch", "Survey")
    # obs/pred ratio: > 1 = model underpredicts, < 1 = overpredicts
    terminal_res$ratio       <- terminal_res$bio_obs / terminal_res$bio_pred

    png("diagnostics/h1_retro_residual_topchange.png", width = 1800, height = 1400, res = 150)
    print(
      # Lollipop anchored at ratio = 1 (perfect fit).
      # Segment runs from 1 to the ratio; point sits at the ratio.
      # Colour encodes direction; no stacking ambiguity.
      ggplot(terminal_res,
             aes(x = peel_yr, colour = ratio > 1)) +
        geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
        geom_segment(aes(xend = peel_yr, y = 1, yend = ratio),
                     linewidth = 1.0) +
        geom_point(aes(y = ratio), size = 3) +
        facet_wrap(~fleet_label, ncol = 3, scales = "free_y") +
        scale_colour_manual(values = c("TRUE"  = "steelblue",
                                       "FALSE" = "tomato"),
                            labels = c("TRUE"  = "Obs > Pred (underpredicted)",
                                       "FALSE" = "Obs < Pred (overpredicted)"),
                            name = NULL) +
        scale_x_continuous(breaks = unique(terminal_res$peel_yr)) +
        labs(title = "Observed / predicted biomass ratio at the terminal year of each retrospective peel",
             subtitle = paste0("Each point = one peel; dashed line = perfect fit (ratio = 1).\n",
                               "Blue = model underpredicts; red = overpredicts.\n",
                               "Consistent direction across peels = systematic mis-fit driving the retrospective."),
             x = "Peel terminal year", y = "Obs / Predicted ratio") +
        theme_bw() +
        theme(axis.text.x    = element_text(angle = 45, hjust = 1, size = 7),
              strip.text      = element_text(size = 8),
              legend.position = "top")
    )
    dev.off()
  }

  # --- Catch-fleet residual spaghetti ---
  # Catch-at-age fleets: summarise across ages (mean std.res per year) then overlay peels.
  catch_retro <- subset(retro_res_all, fleet %in% catch_fleets)
  catch_retro_agg <- aggregate(std.res ~ fleet + year + peel + peel_yr,
                               data = catch_retro, FUN = mean, na.rm = TRUE)

  png("diagnostics/h1_retro_catch_residual_spaghetti.png", width = 1800, height = 1000, res = 150)
  print(
    ggplot() +
      geom_line(data   = subset(catch_retro_agg, peel != base_yr_chr),
                aes(x = year, y = std.res, group = peel, colour = peel_yr),
                linewidth = 0.55, alpha = 0.7) +
      geom_line(data   = subset(catch_retro_agg, peel == base_yr_chr),
                aes(x = year, y = std.res),
                colour = "black", linewidth = 1.3) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
      scale_colour_viridis_c(name = "Peel terminal\nyear") +
      facet_wrap(~fleet, scales = "free_y", ncol = 2) +
      labs(title = "Catch-fleet mean residuals across retrospective peels (mean across ages)",
           subtitle = "Black = full model; colour = peel terminal year.",
           x = "Year", y = "Mean std. residual (across ages)") +
      theme_bw()
  )
  dev.off()

} else {
  message("No residuals found in retro peels — run with ctrl@residuals = TRUE to enable spaghetti diagnostics.")
}

cat("Retrospective diagnostics written to diagnostics/\n")
