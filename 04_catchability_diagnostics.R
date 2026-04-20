# ============================================================
# 04_catchability_diagnostics.R
# Fleet conflict and catchability diagnostics.
# Sourced from 01_run_assessment.R; requires:
#   sam_h1, stk_h1, retro_h1, res_all, ssb_df
#   catch_fleets, survey_fleets, survey_fleets_base, base_fleet
# ============================================================

library(patchwork)

# ------------------------------------------------------------------
# Shared prep: survey residuals with base-name merging
# ------------------------------------------------------------------
survey_res        <- subset(res_all, fleet %in% survey_fleets)
survey_res$fleet  <- base_fleet(survey_res$fleet)

# ------------------------------------------------------------------
# 1. ROLLING MEAN AND CUMULATIVE RESIDUALS
# Rolling mean of std residuals (window = 5 years) reveals drift in
# effective catchability. Cumulative sum: sustained positive = q
# declining over time; sustained negative = q increasing.
# ------------------------------------------------------------------
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
# Tests whether signs of residuals are random; too few runs indicate
# a systematic trend in catchability or model fit.
# ------------------------------------------------------------------
runs_test <- function(x) {
  x   <- x[!is.na(x)]
  sgn <- sign(x)
  n   <- length(sgn)
  if (n < 4) return(data.frame(n = n, n_runs = NA, p_value = NA, result = "too few obs"))
  runs <- sum(sgn[-1] != sgn[-n]) + 1
  n1   <- sum(sgn > 0); n2 <- sum(sgn < 0)
  if (n1 == 0 || n2 == 0)
    return(data.frame(n = n, n_runs = runs, p_value = NA, result = "all same sign"))
  mu   <- (2 * n1 * n2 / n) + 1
  sig2 <- (2 * n1 * n2 * (2 * n1 * n2 - n)) / (n^2 * (n - 1))
  z    <- (runs - mu) / sqrt(sig2)
  p    <- 2 * pnorm(-abs(z))
  data.frame(n = n, n_runs = runs, expected_runs = round(mu, 1),
             z = round(z, 3), p_value = round(p, 4),
             result = ifelse(p < 0.05, "FAIL (trend)", "pass"))
}

runs_results <- do.call(rbind, lapply(c(catch_fleets, survey_fleets_base), function(fl) {
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

# NA p-values (all same sign / too few obs) are invisible with geom_col.
# Give them a sentinel bar height of 1.0 (top of axis) so they always render;
# a text label inside the bar explains the actual situation.
runs_plot    <- runs_results
runs_plot$p_display <- ifelse(is.na(runs_plot$p_value), 1.0, runs_plot$p_value)
runs_plot$na_label  <- ifelse(is.na(runs_plot$p_value), runs_plot$result, NA_character_)

png("diagnostics/h1_runs_test.png", width = 1400, height = 700, res = 150)
print(
  ggplot(runs_plot, aes(x = fleet, y = p_display, fill = result)) +
    geom_col() +
    geom_text(data = subset(runs_plot, !is.na(na_label)),
              aes(label = na_label, y = 0.5),
              colour = "white", size = 3, fontface = "italic",
              hjust = 0.5, inherit.aes = FALSE) +
    geom_hline(yintercept = 0.05, linetype = "dashed", colour = "red") +
    scale_fill_manual(values = c("pass"          = "steelblue",
                                 "FAIL (trend)"  = "tomato",
                                 "too few obs"   = "grey70",
                                 "all same sign" = "grey50")) +
    scale_y_continuous(limits = c(0, 1),
                       breaks = c(0, 0.05, 0.25, 0.5, 0.75, 1)) +
    coord_flip() +
    labs(title = "Runs test p-value by fleet (p < 0.05 = non-random residuals)",
         subtitle = "Grey bar = test not applicable; label shows reason",
         x = NULL, y = "p-value", fill = NULL) +
    theme_bw()
)
dev.off()

# ------------------------------------------------------------------
# 3. IMPLIED CATCHABILITY TREND OVER TIME
# log(obs/SSB) centred within fleet. A slope reveals whether effective
# catchability is increasing or decreasing relative to model assumption.
# ------------------------------------------------------------------
implied_q <- merge(survey_res, ssb_df[, c("year", "value")], by = "year")
implied_q$implied_log_q <- implied_q$log.obs - log(implied_q$value)

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
    labs(title = "Implied catchability trend: log(obs/SSB) centred (slope \u2260 0 = q drift)",
         x = "Year", y = "Relative log(q)") +
    theme_bw()
)
dev.off()

# ------------------------------------------------------------------
# 4. RESIDUAL CROSS-CORRELATION MATRIX
# Strong NEGATIVE correlation between two fleets = they pull the model
# in opposite directions (conflict). Positive = co-vary.
# ------------------------------------------------------------------
annual_res <- do.call(rbind, lapply(c(catch_fleets, survey_fleets), function(fl) {
  sub <- subset(res_all, fleet == fl)
  ag  <- aggregate(std.res ~ year, data = sub, FUN = mean, na.rm = TRUE)
  ag$fleet <- fl
  ag
}))

res_wide <- reshape(annual_res, idvar = "year", timevar = "fleet", direction = "wide")
colnames(res_wide) <- gsub("std.res\\.", "", colnames(res_wide))
res_mat  <- res_wide[, -1]
cor_mat  <- cor(res_mat, use = "pairwise.complete.obs")
cor_df   <- as.data.frame(as.table(cor_mat))
colnames(cor_df) <- c("fleet1", "fleet2", "correlation")
# Use base names for display
cor_df$fleet1 <- base_fleet(cor_df$fleet1)
cor_df$fleet2 <- base_fleet(cor_df$fleet2)

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
# 5. HINDCAST TERMINAL-YEAR RESIDUALS (SS3-style)
# From each peeled retrospective fit: what was the one-step-ahead
# prediction residual in the terminal year of that peel?
# Large deviations across peels = poor predictive skill for that fleet.
# ------------------------------------------------------------------
hindcast_res <- do.call(rbind, lapply(names(retro_h1), function(nm) {
  r <- tryCatch(residuals(retro_h1[[nm]]), error = function(e) NULL)
  if (is.null(r)) return(NULL)
  r <- subset(r, fleet %in% survey_fleets)
  r$peel        <- nm
  r$terminal    <- max(r$year)
  r$is_terminal <- r$year == max(r$year)
  r
}))

terminal_preds        <- subset(hindcast_res, is_terminal)
terminal_preds$fleet  <- base_fleet(terminal_preds$fleet)

base_res_survey        <- subset(res_all, fleet %in% survey_fleets)[, c("fleet", "year", "std.res")]
base_res_survey$fleet  <- base_fleet(base_res_survey$fleet)
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
# Divergence or negative correlation = surveys conflict on stock trend.
# ------------------------------------------------------------------
survey_pairs <- merge(
  reshape(survey_res[, c("fleet", "year", "log.obs")],
          idvar = "year", timevar = "fleet", direction = "wide"),
  data.frame(year = ssb_df$year), by = "year"
)
colnames(survey_pairs) <- gsub("log\\.obs\\.", "", colnames(survey_pairs))

pair_combos <- combn(survey_fleets_base, 2, simplify = FALSE)
pair_plots  <- lapply(pair_combos, function(fl) {
  if (!all(fl %in% colnames(survey_pairs))) return(NULL)
  tmp <- survey_pairs[, c("year", fl[1], fl[2])]
  colnames(tmp) <- c("year", "x", "y")
  tmp <- tmp[complete.cases(tmp), ]
  if (nrow(tmp) < 4) return(NULL)
  r <- round(cor(tmp$x, tmp$y), 2)
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

n_cols <- 3
png("diagnostics/h1_survey_pairwise.png",
    width  = n_cols * 600,
    height = ceiling(length(pair_plots) / n_cols) * 600, res = 150)
print(wrap_plots(pair_plots, ncol = n_cols))
dev.off()

# ------------------------------------------------------------------
# 7. COHORT TRACKING THROUGH CATCH-AT-AGE FLEETS
# For each cohort: catch[age, cohort+age] across all fleets.
# Diverging trends between fleets = conflicts in age-composition data.
# ------------------------------------------------------------------
catch_n_df <- as.data.frame(stk_h1@catch.n)
catch_n_df$cohort <- catch_n_df$year - catch_n_df$age
catch_n_df$data[catch_n_df$data <= 0] <- NA
catch_n_df <- subset(catch_n_df, !is.na(data) & year >= 2000)

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
# 8. FLEET-SPECIFIC MEAN AGE OF CATCH: residual-weighted
# Persistent offset between observed and model-predicted mean age
# signals misallocation of F across ages for that fleet.
# ------------------------------------------------------------------
mean_age_obs <- aggregate(age * abs(std.res) ~ year + fleet,
                          data = subset(res_all, fleet %in% catch_fleets & !is.na(age)),
                          FUN = sum)
mean_age_wt  <- aggregate(abs(std.res) ~ year + fleet,
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

cat("Catchability diagnostics written to diagnostics/\n")
