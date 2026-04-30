library(FLCore)
library(FLSAM)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
source("R/sam_helpers.R")

load("output_h1.RData")

dir.create("diagnostics", showWarnings = FALSE)

inject_context(prepare_diagnostics_context(stk_h1, idx_h1, sam_h1))
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

render_quarto_report("sam_diagnostics.qmd")
