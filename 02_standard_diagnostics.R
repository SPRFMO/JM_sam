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
png("diagnostics/h1_stock_trajectory.png", width = 1600, height = 1200, res = 150)
print(plot(sam_h1))
dev.off()

##################################################
## Stock trajectory zoom (SSB / Fbar / Rec)
##################################################
df.ssb  <- ssb(sam_h1);  df.ssb$quant  <- "SSB"
df.fbar <- fbar(sam_h1); df.fbar$quant <- "Fbar"
df.rec  <- rec(sam_h1);  df.rec$quant  <- "Recruitment"
df.traj <- rbind(df.ssb, df.fbar, df.rec)

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
## Observation variance by data source
## Two-panel figure (catch / survey), colour-coded by dataset.
## Each unique obs.var parameter appears once, labelled with the
## fleet/survey name and (for catch) the age or age-range it covers.
##################################################
obv            <- obs.var(sam_h1)
obv$base       <- base_fleet(obv$fleet)
obv$group      <- ifelse(obv$fleet %in% catch_fleets, "Catch fleets", "Survey indices")
# Short label: drop "catch " prefix; tag age for catch fleets
obv$label      <- paste0(
  sub("catch ", "", obv$base),
  ifelse(is.na(obv$age), "", paste0(" (age ", obv$age, ")"))
)

# Colour palette: one colour per base dataset (4 catch + 7 survey = 11 groups)
all_bases   <- unique(obv$base)
n_bases     <- length(all_bases)
base_colors <- setNames(
  c(RColorBrewer::brewer.pal(4, "Dark2"),           # catch fleets
    RColorBrewer::brewer.pal(min(7, n_bases - 4), "Set1")),   # surveys
  all_bases
)

# Sort within each panel by value (ascending)
obv_catch  <- obv[obv$group == "Catch fleets",   ]
obv_survey <- obv[obv$group == "Survey indices",  ]
obv_catch  <- obv_catch[order(obv_catch$value),   ]
obv_survey <- obv_survey[order(obv_survey$value), ]
obv_catch$label  <- factor(obv_catch$label,  levels = obv_catch$label)
obv_survey$label <- factor(obv_survey$label, levels = obv_survey$label)

p_obv_catch <- ggplot(obv_catch, aes(x = label, y = value, fill = base)) +
  geom_col(width = 0.75) +
  coord_flip() +
  scale_fill_manual(values = base_colors, name = "Fleet") +
  labs(title = "Catch fleets", x = NULL, y = "Observation variance") +
  theme_bw() +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 10))

p_obv_survey <- ggplot(obv_survey, aes(x = label, y = value, fill = base)) +
  geom_col(width = 0.75) +
  coord_flip() +
  scale_fill_manual(values = base_colors, name = "Survey") +
  labs(title = "Survey indices", x = NULL, y = "Observation variance") +
  theme_bw() +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 10))

# Height proportional to number of rows in each panel
n_catch_rows  <- nrow(obv_catch)
n_survey_rows <- nrow(obv_survey)

png("diagnostics/h1_observation_var.png",
    width  = 1600,
    height = 200 + (n_catch_rows + n_survey_rows) * 55,
    res    = 150)
print(
  p_obv_catch / p_obv_survey +
    plot_layout(heights = c(n_catch_rows, n_survey_rows)) +
    plot_annotation(title = "Observation variances by data source",
                    theme = theme(plot.title = element_text(face = "bold", size = 13)))
)
dev.off()

##################################################
## Observation variance vs CV
##################################################
png("diagnostics/h1_variance_vs_cv.png", width = 1200, height = 900, res = 150)
plot(obv$value, obv$CV, xlab = "Observation variance", ylab = "CV of estimate",
     log = "x", pch = 16, col = as.integer(factor(obv$fleet)),
     main = "Observation variance vs uncertainty")
text(obv$value, obv$CV, obv$label, pos = 4, cex = 0.75, xpd = NA)
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

##################################################
## Correlation matrix of model parameters
## Large canvas + lower res so parameter labels have room to breathe.
##################################################
png("diagnostics/h1_cor_params.png", width = 3000, height = 3000, res = 200)
par(mar = c(8, 8, 4, 2))   # generous margins for axis labels
cor.plot(sam_h1)
dev.off()

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
## Collapse harvest to total F so all slots share one area.
##################################################
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
