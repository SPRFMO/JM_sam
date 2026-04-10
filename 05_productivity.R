# ============================================================
# 05_productivity.R
# Beverton-Holt stock-recruit analysis with change-point detection.
# Sourced from 01_run_assessment.R; requires:
#   sam_h1
# ============================================================

library(patchwork)

################################################################################
## PRODUCTIVITY ANALYSIS: BEVERTON-HOLT S-R WITH CHANGE-POINT DETECTION
##
## Uncertainty in both R and SSB estimates from SAM is carried through:
##   - Weights for NLS fitting = 1 / Var(log R), approximated from 95% CI
##   - Error bars on both axes in all plots
## Change points are found by minimising total weighted residual SS across
## candidate break years; significance is tested with an F-test (extra
## parameters = 2 per additional B-H curve). Up to 3 panels are produced:
##   Panel 1 - full time series, single B-H, steepness annotated
##   Panel 2 - if 1 break is significant (p<0.05): 2 B-H curves
##   Panel 3 - if a 2nd break is also significant: 3 B-H curves
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
sr$logR_se <- pmax((sr$R_hi - sr$R_lo) / (4 * sr$R), 0.01)
sr$wt      <- 1 / sr$logR_se^2
sr$wt      <- sr$wt / max(sr$wt)   # scale to [0,1]

B0_ref <- max(sr$SSB)
n_sr   <- nrow(sr)

# ------------------------------------------------------------------
# Beverton-Holt NLS fit (weighted)
# B-H: R = a*S/(b+S)
# Steepness: h = 0.2*(b+B0)/(b+0.2*B0), bounded [0.21, 0.99]
# ------------------------------------------------------------------
fit_bh <- function(df, B0 = B0_ref) {
  if (nrow(df) < 6) return(list(converged = FALSE))
  starts <- list(a = quantile(df$R, 0.9) * 2, b = median(df$SSB))
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

bh_pred  <- function(a, b, ssb_vec) a * ssb_vec / (b + ssb_vec)
ssb_grid <- seq(0, max(sr$SSB) * 1.05, length.out = 300)

# ------------------------------------------------------------------
# 1. Full series fit
# ------------------------------------------------------------------
fit_all   <- fit_bh(sr)
wrss_full <- if (fit_all$converged) fit_all$wrss else Inf

# ------------------------------------------------------------------
# 2. Single break-point search
# ------------------------------------------------------------------
min_seg_sr <- 10

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
  best1     <- break1_res[which.min(break1_res$wrss), ]
  f1stat    <- ((wrss_full - best1$wrss) / 2) / (best1$wrss / (n_sr - 4))
  p1_val    <- pf(f1stat, 2, n_sr - 4, lower.tail = FALSE)
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
    best2     <- all2[which.min(all2$wrss), ]
    f2stat    <- ((best1$wrss - best2$wrss) / 2) / (best2$wrss / (n_sr - 6))
    p2_val    <- pf(f2stat, 2, n_sr - 6, lower.tail = FALSE)
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

# Panel 1: full series
crv1  <- data.frame(SSB = ssb_grid, R = bh_pred(fit_all$a, fit_all$b, ssb_grid))
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

# Panel 2: one break
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
    theme(legend.position = "bottom", legend.text = element_text(size = 8))
}

# Panel 3: two breaks
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
  crv3$period <- factor(crv3$period, levels = c(best2$p1, best2$p2, best2$p3))

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
    theme(legend.position = "bottom", legend.text = element_text(size = 8))
}

# ------------------------------------------------------------------
# 5. Productivity over time: R/SSB with CI (delta method)
# ------------------------------------------------------------------
sr$logSSB_se   <- pmax((sr$SSB_hi - sr$SSB_lo) / (4 * sr$SSB), 0.01)
sr$prod        <- sr$R / sr$SSB
sr$log_prod_se <- sqrt(sr$logR_se^2 + sr$logSSB_se^2)
sr$prod_lo     <- sr$prod * exp(-1.96 * sr$log_prod_se)
sr$prod_hi     <- sr$prod * exp( 1.96 * sr$log_prod_se)

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
# 6. Combine B-H panels and save
# ------------------------------------------------------------------
panel_plots <- Filter(Negate(is.null), list(p_sr1, p_sr2, p_sr3))
n_panels    <- length(panel_plots)

png("diagnostics/h1_productivity_sr.png",
    width = n_panels * 900, height = 900, res = 150)
print(wrap_plots(panel_plots, nrow = 1))
dev.off()

# Save text summary
sink("diagnostics/h1_productivity_summary.txt")
cat("Beverton-Holt productivity analysis\n")
cat(sprintf("B0 reference (max observed SSB): %.0f\n\n", B0_ref))
cat(sprintf("Full series: h = %.3f (%d-%d)\n", fit_all$h, min(sr$year), max(sr$year)))
if (!is.na(p1_val))
  cat(sprintf("1-break test: break year %d, p = %.4f %s\n",
              best1$break_year+1, p1_val,
              ifelse(do_2panel, "[significant]", "[not significant]")))
if (do_2panel)
  cat(sprintf("  Period 1 (%s): h = %.3f\n  Period 2 (%s): h = %.3f\n",
              best1$p1, best1$h1, best1$p2, best1$h2))
if (!is.na(p2_val))
  cat(sprintf("2-break test: breaks %d & %d, p = %.4f %s\n",
              best2$break1+1, best2$break2+1, p2_val,
              ifelse(do_3panel, "[significant]", "[not significant]")))
if (do_3panel)
  cat(sprintf("  Period 1 (%s): h = %.3f\n  Period 2 (%s): h = %.3f\n  Period 3 (%s): h = %.3f\n",
              best2$p1, best2$h1, best2$p2, best2$h2, best2$p3, best2$h3))
sink()

cat("Productivity analysis written to diagnostics/\n")
