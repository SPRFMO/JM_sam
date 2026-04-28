library(FLCore)
library(FLSAM)
library(ggplot2)
library(patchwork)
library(RColorBrewer)

load("output_h1.RData")

dir.create("diagnostics", showWarnings = FALSE)

catch_fleets <- c(
  "catch N_Chile", "catch SC_Chile_PS",
  "catch FarNorth", "catch Offshore_Trawl"
)
survey_fleets <- names(idx_h1)

q_breaks <- list(
  Chile_AcousCS = 2002,
  Chile_CPUE = 2000,
  Offshore_CPUE = 2021
)

base_fleet <- function(x) sub("_(early|late)$", "", x)
survey_fleets_base <- unique(base_fleet(survey_fleets))

max_yr <- range(stk_h1)["maxyear"]
base_yr_chr <- as.character(max_yr)

stk_h1_fit <- stk_h1 + sam_h1
res_all <- residuals(sam_h1)
ssb_df <- ssb(sam_h1)

sel.pat <- merge(f(sam_h1), fbar(sam_h1), by = "year", suffixes = c(".f", ".fbar"))
sel.pat$sel <- sel.pat$value.f / sel.pat$value.fbar
sel.pat$age <- as.numeric(as.character(sel.pat$age))
sel.pat$pentad <- sprintf("%i's", floor(sel.pat$year / 5) * 5)

.parse_jjm_ssb <- function(rep_file, max_year) {
  if (!file.exists(rep_file)) {
    warning("JJM rep file not found: ", rep_file)
    return(NULL)
  }

  lines <- readLines(rep_file, warn = FALSE)
  sec_idx <- grep("^\\$[A-Za-z_]+$", lines)
  sec_names <- sub("^\\$", "", lines[sec_idx])

  i_ssb <- sec_idx[sec_names == "SSB"]
  if (length(i_ssb) == 0) {
    warning("$SSB not found in rep file")
    return(NULL)
  }

  pos <- which(sec_idx == i_ssb)
  end <- if (pos < length(sec_idx)) sec_idx[pos + 1] - 1 else length(lines)
  dlines <- trimws(lines[(i_ssb + 1):end])
  dlines <- dlines[nchar(dlines) > 0]
  mat <- do.call(rbind, lapply(dlines, function(l) as.numeric(strsplit(l, "\\s+")[[1]])))

  if (ncol(mat) < 5) {
    warning("Unexpected $SSB format")
    return(NULL)
  }

  mat <- mat[mat[, 1] <= max_year, , drop = FALSE]
  data.frame(
    year = as.integer(mat[, 1]),
    value = mat[, 2],
    lbnd = mat[, 4],
    ubnd = mat[, 5]
  )
}

jjm_ssb <- .parse_jjm_ssb(
  rep_file = "../assessment/results/h1_1.14_1_R.rep",
  max_year = as.integer(range(stk_h1)["maxyear"])
)
if (!is.null(jjm_ssb)) {
  cat(sprintf(
    "JJM SSB loaded: %i years (%i-%i)\n",
    nrow(jjm_ssb), min(jjm_ssb$year), max(jjm_ssb$year)
  ))
}

source("02_standard_diagnostics.R")
source("03_retro_diagnostics.R")
source("04_catchability_diagnostics.R")

rec_sr <- rec(sam_h1)
ssb_sr <- ssb(sam_h1)
sr_check <- merge(
  data.frame(year = rec_sr$year, R = rec_sr$value, R_lo = rec_sr$lbnd, R_hi = rec_sr$ubnd),
  data.frame(year = ssb_sr$year, SSB = ssb_sr$value, SSB_lo = ssb_sr$lbnd, SSB_hi = ssb_sr$ubnd),
  by = "year"
)
sr_check <- sr_check[complete.cases(sr_check) & sr_check$R > 0 & sr_check$SSB > 0, ]

write_productivity_placeholders <- function(msg) {
  placeholder <- ggplot(data.frame(x = 0, y = 0), aes(x, y)) +
    annotate("text", x = 0, y = 0,
             label = paste("Productivity analysis unavailable", msg, sep = "\n"),
             size = 5) +
    theme_void()

  png("diagnostics/h1_productivity_time.png", width = 1400, height = 800, res = 150)
  print(placeholder + labs(title = "Productivity over time"))
  dev.off()

  png("diagnostics/h1_productivity_sr.png", width = 900, height = 900, res = 150)
  print(placeholder + labs(title = "Beverton-Holt productivity analysis"))
  dev.off()

  writeLines(
    c(
      "Beverton-Holt productivity analysis",
      paste("Unavailable:", msg)
    ),
    "diagnostics/h1_productivity_summary.txt"
  )
  cat("Productivity analysis skipped:", msg, "\n")
}

if (nrow(sr_check) >= 6) {
  prod_ok <- tryCatch({
    source("05_productivity.R")
    TRUE
  }, error = function(e) {
    write_productivity_placeholders(conditionMessage(e))
    FALSE
  })
} else {
  write_productivity_placeholders("No finite stock-recruit pairs in this run.")
}
source("06_selectivity_blocks.R")

status <- system2("quarto", c("render", "sam_diagnostics.qmd"))
if (!identical(status, 0L)) {
  stop("quarto render failed with status ", status)
}
