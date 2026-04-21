# ============================================================
# 02_standard_diagnostics.R
# Standard SAM diagnostics: residuals, trajectories, inputs.
# Sourced from 01_run_assessment.R; requires:
#   sam_h1, stk_h1, stk_h1_fit, idx_h1, res_all, sel.pat
#   catch_fleets, survey_fleets, survey_fleets_base, base_fleet
#   q_breaks
# ============================================================

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
q_split$q     <- exp(q_split$log_q)
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
df.ssb  <- ssb(sam_h1);  df.ssb$quant  <- "SSB"
df.fbar <- fbar(sam_h1); df.fbar$quant <- "Fbar"
df.rec  <- rec(sam_h1);  df.rec$quant  <- "Recruitment"
df.traj_all <- rbind(df.ssb, df.fbar, df.rec)

png("diagnostics/h1_stock_trajectory.png", width = 1600, height = 1200, res = 150)
print(
  ggplot(df.traj_all, aes(x = year, y = value)) +
    geom_ribbon(aes(ymin = lbnd, ymax = ubnd), alpha = 0.25, fill = "grey70") +
    geom_line(linewidth = 0.9, colour = "black") +
    facet_wrap(~quant, scales = "free_y", ncol = 1) +
    labs(title = "Model-estimated stock trajectory", x = "Year", y = NULL) +
    theme_bw()
)
dev.off()

##################################################
## Stock trajectory zoom (SSB / Fbar / Rec)
##################################################
df.traj <- df.traj_all

png("diagnostics/h1_stock_trajectory_zoom.png", width = 1200, height = 1600, res = 150)
p_traj <- ggplot(subset(df.traj, year > 2002), aes(x = year, y = value)) +
  geom_ribbon(aes(ymin = lbnd, ymax = ubnd), alpha = 0.3) +
  geom_line() +
  ylim(0, NA) +
  facet_wrap(~quant, scales = "free", ncol = 1) +
  labs(x = "Year", y = NULL) +
  theme_bw()

# Overlay JJM h1_1.14 SSB on the SSB panel (dashed tomato line)
if (!is.null(jjm_ssb)) {
  jjm_ssb_sub        <- subset(jjm_ssb, year > 2002)
  jjm_ssb_sub$quant  <- "SSB"   # must match facet label
  p_traj <- p_traj +
    geom_ribbon(data = jjm_ssb_sub,
                aes(x = year, ymin = lbnd, ymax = ubnd),
                fill = "tomato", alpha = 0.15, colour = NA,
                inherit.aes = FALSE) +
    geom_line(data = jjm_ssb_sub,
              aes(x = year, y = value),
              colour = "tomato", linetype = "dashed", linewidth = 0.9,
              inherit.aes = FALSE) +
    labs(caption = "SAM: solid black; JJM h1_1.14: dashed red (SSB panel only)")
}
print(p_traj)
dev.off()

##################################################
## F variance by fleet and age
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
## Observation variance by fleet and age
## One facet per fleet, fixed y-axis across all panels.
##################################################
obv             <- obs.var(sam_h1)
obv$fleet_clean <- sub("^catch ", "", obv$fleet)
obv$label       <- paste0(obv$fleet_clean,
                           ifelse(is.na(obv$age), "", paste0(" a", obv$age)))
obv$age_label   <- ifelse(is.na(obv$age), "idx", as.character(obv$age))
obv$age_label   <- factor(obv$age_label,
                           levels = c(as.character(sort(unique(na.omit(obv$age)))), "idx"))

y_max <- max(obv$ubnd, na.rm = TRUE)

png("diagnostics/h1_observation_var.png", width = 2400, height = 1600, res = 150)
print(
  ggplot(obv, aes(x = age_label, y = value)) +
    geom_col(aes(fill = fleet_clean), width = 0.7, show.legend = FALSE) +
    geom_errorbar(aes(ymin = lbnd, ymax = ubnd), width = 0.25) +
    facet_wrap(~fleet_clean, scales = "fixed") +
    scale_y_continuous(limits = c(0, y_max * 1.05)) +
    labs(title = "Observation variances by fleet and age",
         x = "Age", y = "Observation variance") +
    theme_bw() +
    theme(strip.text  = element_text(size = 8),
          axis.text.x = element_text(size = 7))
)
dev.off()

##################################################
## Observation variance vs CV
##################################################
png("diagnostics/h1_variance_vs_cv.png", width = 1200, height = 900, res = 150)
print(
  ggplot(subset(obv, is.finite(value) & is.finite(CV) & value > 0),
         aes(x = value, y = CV, colour = fleet_clean, label = label)) +
    geom_point(size = 2.5) +
    geom_text(hjust = 0, nudge_x = 0.02, size = 2.8, show.legend = FALSE) +
    scale_x_log10() +
    labs(title = "Observation variance vs uncertainty",
         x = "Observation variance",
         y = "CV of estimate",
         colour = "Fleet") +
    theme_bw() +
    theme(legend.position = "bottom")
)
dev.off()

##################################################
## Selectivity (F/Fbar at age) by fleet and pentad
## sel.pat pre-computed in 01_run_assessment.R
##################################################
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
try(dev.off())
##################################################
## Correlation matrix of model parameters
##################################################
param_vcov <- tryCatch(vcov(sam_h1), error = function(e) NULL)
if (!is.null(param_vcov)) {
  param_cor <- cov2cor(param_vcov)
  param_names <- make.unique(colnames(param_cor))
  colnames(param_cor) <- param_names
  rownames(param_cor) <- make.unique(rownames(param_cor))
  cor_df <- as.data.frame(as.table(param_cor))
  names(cor_df) <- c("param_x", "param_y", "correlation")
  cor_df$param_x <- factor(cor_df$param_x, levels = param_names)
  cor_df$param_y <- factor(cor_df$param_y, levels = rev(param_names))

  png("diagnostics/h1_cor_params.png", width = 1800, height = 1800, res = 120)
  print(
    ggplot(cor_df, aes(x = param_x, y = param_y, fill = correlation)) +
      geom_tile(colour = "white", linewidth = 0.15) +
      scale_fill_gradient2(low = "#b2182b", mid = "white", high = "#2166ac",
                           midpoint = 0, limits = c(-1, 1)) +
      labs(title = "Pairwise correlation matrix of estimated parameters",
           x = NULL, y = NULL, fill = "Correlation") +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6),
        panel.grid = element_blank()
      )
  )
  dev.off()
}


##################################################
## Catch residuals (bubble plot, one panel per fleet)
##################################################
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
res_survey_disp        <- subset(res_all, fleet %in% survey_fleets)
res_survey_disp$fleet  <- base_fleet(res_survey_disp$fleet)

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
##################################################
stock_n_df  <- as.data.frame(stk_h1_fit@stock.n)
stock_wt_df <- as.data.frame(stk_h1_fit@stock.wt)
m_df        <- as.data.frame(stk_h1_fit@m)
f_df        <- as.data.frame(harvest(stk_h1_fit))

total_f_df <- aggregate(data ~ age + year, data = f_df, FUN = sum, na.rm = TRUE)
names(total_f_df)[names(total_f_df) == "data"] <- "total_f"

n_curr <- stock_n_df[, c("age", "year", "data")]
names(n_curr)[3] <- "stock_n"
n_prev <- n_curr
n_prev$age  <- n_prev$age + 1
n_prev$year <- n_prev$year + 1
names(n_prev)[3] <- "stock_n_prev"

z_prev <- merge(
  total_f_df,
  m_df[, c("age", "year", "data")],
  by = c("age", "year"),
  all = FALSE
)
names(z_prev)[names(z_prev) == "data"] <- "m"
z_prev$z <- z_prev$total_f + z_prev$m
z_prev$age  <- z_prev$age + 1
z_prev$year <- z_prev$year + 1
z_prev <- z_prev[, c("age", "year", "z")]

proc_n_df <- merge(n_curr, n_prev, by = c("age", "year"))
proc_n_df <- merge(proc_n_df, z_prev, by = c("age", "year"))
proc_n_df <- merge(proc_n_df, stock_wt_df[, c("age", "year", "data")], by = c("age", "year"))
names(proc_n_df)[names(proc_n_df) == "data"] <- "stock_wt"
proc_n_df <- subset(proc_n_df,
                    stock_n > 0 & stock_n_prev > 0 & is.finite(z) & is.finite(stock_wt))
proc_n_df$proc_err <- log(proc_n_df$stock_n) - (log(proc_n_df$stock_n_prev) - proc_n_df$z)
proc_n_df$weight_scale <- proc_n_df$stock_wt / max(proc_n_df$stock_wt, na.rm = TRUE)
proc_n_df$proc_err_scaled <- proc_n_df$proc_err * proc_n_df$weight_scale

png("diagnostics/h1_procerr_N.png", width = 1600, height = 900, res = 150)
print(
  ggplot(proc_n_df, aes(x = year, y = age, fill = proc_err_scaled)) +
    geom_tile() +
    scale_fill_gradient2(low = "#b2182b", mid = "white", high = "#2166ac", midpoint = 0) +
    labs(title = "Process errors in log-N", x = "Year", y = "Age", fill = "Scaled\nerror") +
    theme_bw()
)
dev.off()

proc_m_df <- merge(total_f_df, total_f_df, by = "age", suffixes = c("", "_prev"))
proc_m_df <- subset(proc_m_df, year == year_prev + 1 & total_f > 0 & total_f_prev > 0)
proc_m_df$proc_err <- log(proc_m_df$total_f) - log(proc_m_df$total_f_prev)
proc_m_df$year <- proc_m_df$year

png("diagnostics/h1_procerr_M.png", width = 1600, height = 900, res = 150)
print(
  ggplot(proc_m_df, aes(x = year, y = age, fill = proc_err)) +
    geom_tile() +
    scale_fill_gradient2(low = "#b2182b", mid = "white", high = "#2166ac", midpoint = 0) +
    labs(title = "Process errors in log-mortality", x = "Year", y = "Age", fill = "Delta\nlog(F)") +
    theme_bw()
)
dev.off()

##################################################
## F at age - one panel per catch fleet
##################################################
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
df.sr <- merge(
  setNames(ssb(sam_h1)[, c("year", "value")], c("year", "ssb")),
  setNames(rec(sam_h1)[, c("year", "value")], c("year", "rec")),
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
## Index fits - observed vs model-predicted
##################################################
fit_idx     <- residuals(sam_h1)
survey_fit  <- subset(fit_idx, fleet %in% survey_fleets)
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
  d$fleet  <- base_fleet(nm)
  d$data[d$data <= 0] <- NA
  d
}))
survey_raw <- subset(survey_raw, !is.na(data))
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
## N at age (log scale)
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

cat("Standard diagnostics written to diagnostics/\n")
