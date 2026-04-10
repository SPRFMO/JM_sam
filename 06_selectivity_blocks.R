# ============================================================
# 06_selectivity_blocks.R
# Selectivity blocking analysis for JJM.
# Sourced from 01_run_assessment.R; requires:
#   stk_h1, sel.pat, catch_fleets
# ============================================================

library(patchwork)

################################################################################
## SELECTIVITY BLOCKING ANALYSIS FOR JJM
##
## Rationale:
##   JJM assumes selectivity is constant within user-defined year blocks.
##   SAM estimates selectivity annually (as F/Fbar at age) using a random walk,
##   so the SAM trajectory reveals when the selectivity shape changed
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
##   4. Rank breakpoints by importance = F-statistic x fraction of total fleet
##      catch in the split segment. Breaks in high-catch periods rank highest;
##      label 1 is most important and should be retained last when reducing
##      the number of blocks.
##
## Output:
##   diagnostics/h1_selectivity_blocks.csv       - block suggestion table
##   diagnostics/h1_selectivity_pc1.png          - PC1 + ranked breakpoints
##   diagnostics/h1_selectivity_byblock.png      - selectivity curves by block
##   diagnostics/h1_selectivity_catch_weight.png - catch weights used
################################################################################

# ------------------------------------------------------------------
# Step 1 - Total catch per fleet per year (used as weights)
# ------------------------------------------------------------------
catch_wt_df <- as.data.frame(stk_h1@catch.n * stk_h1@catch.wt)
catch_wt_df <- aggregate(data ~ year + area, data = catch_wt_df, FUN = sum, na.rm = TRUE)
catch_wt_df$data[catch_wt_df$data <= 0] <- NA

# Map area dimname to fleet name positionally
area_names  <- dimnames(stk_h1@catch.n)$area
catch_wt_df$fleet <- catch_fleets[match(as.character(catch_wt_df$area), area_names)]

# Normalise within fleet: weight = catch / max(catch) so max weight = 1
catch_wt_df <- do.call(rbind, lapply(split(catch_wt_df, catch_wt_df$fleet), function(df) {
  mx        <- max(df$data, na.rm = TRUE)
  df$weight <- ifelse(is.na(df$data), 0, df$data / mx)
  df
}))

# ------------------------------------------------------------------
# Step 2 - Weighted binary segmentation helpers
# ------------------------------------------------------------------
min_seg <- 5
alpha   <- 0.05

find_cp_w <- function(x, w) {
  n <- length(x)
  w <- pmax(w, 1e-6)
  w <- w / sum(w)
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
    bin_seg_w(x[seq_len(res$cp)],         w[seq_len(res$cp)],
              idx[seq_len(res$cp)],        catch_raw[seq_len(res$cp)],      depth + 1),
    bin_seg_w(x[(res$cp + 1):length(x)],  w[(res$cp + 1):length(x)],
              idx[(res$cp + 1):length(x)], catch_raw[(res$cp + 1):length(x)], depth + 1)
  )
}

# ------------------------------------------------------------------
# Step 3 - Run per fleet
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

  # Match catch weights to years in sel
  wt_fleet  <- subset(catch_wt_df, fleet == fl)
  w_vec     <- wt_fleet$weight[match(yrs, wt_fleet$year)]
  w_vec[is.na(w_vec)] <- 0
  catch_raw <- wt_fleet$data[match(yrs, wt_fleet$year)]
  catch_raw[is.na(catch_raw)] <- 0

  total_catch <- sum(catch_raw, na.rm = TRUE)

  # Weighted binary segmentation
  cps_df <- bin_seg_w(pc1, w_vec, seq_along(yrs), catch_raw)

  if (nrow(cps_df) == 0) {
    break_yrs  <- integer(0)
    importance <- numeric(0)
  } else {
    cps_df$importance <- cps_df$f_stat * (cps_df$catch_sum / total_catch)
    cps_df <- cps_df[order(cps_df$idx_pos), ]
    cps_df$rank       <- rank(-cps_df$importance, ties.method = "first")
    break_yrs         <- yrs[cps_df$idx_pos]
    cps_df$break_year <- break_yrs
    importance        <- cps_df$importance
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
      sprintf("Peak age shifted %i \u2192 %i (%.1f%% of catch)", prev_pk, pk_age, ct_frac)
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
# Step 4 - Combined block table
# ------------------------------------------------------------------
sel_blocks_all <- do.call(rbind, lapply(sel_block_list, `[[`, "blocks"))
rownames(sel_blocks_all) <- NULL
write.csv(sel_blocks_all, "diagnostics/h1_selectivity_blocks.csv", row.names = FALSE)

cat("\n--- Suggested JJM selectivity blocks ---\n")
print(sel_blocks_all[, c("fleet", "block", "start_year", "end_year",
                          "n_years", "catch_pct", "peak_age",
                          "within_var", "boundary_rank", "rationale")])

# ------------------------------------------------------------------
# Step 5 - Plot: PC1 with ranked breakpoints and catch weight ribbon
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
  bp_all$label_y <- bp_all$ymax
}

# Pre-compute per-fleet ribbon bounds: catch weight fills bottom 35% of y-range
pc1_all <- do.call(rbind, lapply(split(pc1_all, pc1_all$fleet), function(df) {
  y_min   <- min(df$pc1, na.rm = TRUE)
  y_range <- max(df$pc1, na.rm = TRUE) - y_min
  df$ribbon_lo <- y_min
  df$ribbon_hi <- y_min + df$catch_weight * y_range * 0.35
  df
}))

png("diagnostics/h1_selectivity_pc1.png", width = 1800, height = 1200, res = 150)
p_pc1 <- ggplot(pc1_all, aes(x = year, y = pc1)) +
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
# Step 6 - Plot: selectivity curves coloured by block
# ------------------------------------------------------------------
sel_with_block <- do.call(rbind, lapply(sel_block_list, function(x) {
  bd <- x$blocks
  df <- subset(sel.pat, fleet == x$fl & is.finite(sel))
  df$block <- NA_integer_
  for (b in seq_len(nrow(bd))) {
    df$block[df$year >= bd$start_year[b] & df$year <= bd$end_year[b]] <- b
  }
  df$block_label <- sprintf("Block %i (%i\u2013%i)", df$block,
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
# Step 7 - Plot: catch weight timeseries per fleet
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
