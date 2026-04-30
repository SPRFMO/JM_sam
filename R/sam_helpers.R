survey_obs_groups <- list(
  Chile_AcousCS = c("Chile_AcousCS_early", "Chile_AcousCS_late", "Chile_AcousCS"),
  Chile_AcousN = "Chile_AcousN",
  Chile_CPUE = c("Chile_CPUE_early", "Chile_CPUE_late", "Chile_CPUE"),
  DEPM = "DEPM",
  Peru_Acoustic = "Peru_Acoustic",
  Peru_CPUE = "Peru_CPUE",
  Offshore_CPUE = c("Offshore_CPUE_early", "Offshore_CPUE_late", "Offshore_CPUE")
)

survey_biomass_treat <- c(
  Chile_AcousCS_early = 5, Chile_AcousCS_late = 5, Chile_AcousCS = 5,
  Chile_AcousN = 5,
  Chile_CPUE_early = 2, Chile_CPUE_late = 2, Chile_CPUE = 2,
  DEPM = 0,
  Peru_Acoustic = 5,
  Peru_CPUE = 2,
  Offshore_CPUE_early = 2, Offshore_CPUE_late = 2, Offshore_CPUE = 2
)

q_breaks <- list(
  Chile_AcousCS = 2002,
  Chile_CPUE = 2000,
  Offshore_CPUE = 2021
)

catch_fleets <- c(
  "catch N_Chile",
  "catch SC_Chile_PS",
  "catch FarNorth",
  "catch Offshore_Trawl"
)

base_fleet <- function(x) sub("_(early|late)$", "", x)

skip_flsam_osa_residuals <- function() {
  sam2flr <- get("SAM2FLR", envir = asNamespace("FLSAM"))

  replace_osa_call <- function(expr) {
    if (!is.call(expr)) return(expr)
    if (identical(expr[[1]], as.name("oneStepPredict"))) {
      return(quote(data.frame(residual = fit$data$logobs - fit$rep$predObs)))
    }
    as.call(lapply(as.list(expr), replace_osa_call))
  }

  body(sam2flr) <- replace_osa_call(body(sam2flr))
  assignInNamespace("SAM2FLR", sam2flr, ns = "FLSAM")
}

biomass_to_FLIndex <- function(x) {
  empty <- x@index
  empty[] <- NA

  r <- x@range
  r["min"] <- NA
  r["max"] <- NA
  r["plusgroup"] <- NA

  out <- new("FLIndex",
    distribution = x@distribution,
    index = x@index,
    index.var = empty,
    catch.n = empty,
    catch.wt = empty,
    effort = empty,
    sel.pattern = empty,
    index.q = empty,
    name = x@name,
    desc = x@desc,
    range = r
  )
  type(out) <- "biomass"
  out
}

subset_indices <- function(x, nms) {
  FLIndices(setNames(lapply(nms, function(nm) x[[nm]]), nms))
}

build_obs_map <- function(idx_names, start_code = 401) {
  out <- numeric(0)
  code <- start_code
  for (grp in survey_obs_groups) {
    present <- grp[grp %in% idx_names]
    if (length(present) > 0) {
      out[present] <- code
      code <- code + 1
    }
  }
  out
}

split_index_at_q_break <- function(fl, break_yr) {
  early <- window(fl, end = break_yr - 1)
  late <- window(fl, start = break_yr)
  early@name <- paste0(fl@name, "_early")
  late@name <- paste0(fl@name, "_late")
  list(early = early, late = late)
}

split_indices_by_q_break <- function(idx, breaks = q_breaks) {
  idx_split <- list()
  for (nm in names(idx)) {
    if (nm %in% names(breaks) &&
        range(idx[[nm]])["minyear"] < breaks[[nm]] &&
        range(idx[[nm]])["maxyear"] >= breaks[[nm]]) {
      sp <- split_index_at_q_break(idx[[nm]], breaks[[nm]])
      idx_split[[paste0(nm, "_early")]] <- sp$early
      idx_split[[paste0(nm, "_late")]] <- sp$late
    } else {
      idx_split[[nm]] <- idx[[nm]]
    }
  }
  FLIndices(idx_split)
}

expand_to_fleets <- function(flq, n) {
  if (dim(flq)[5] == n) return(flq)
  new_dim <- dim(flq)
  new_dim[5] <- n
  new_dn <- dimnames(flq)
  new_dn[["area"]] <- seq_len(n)
  out <- FLQuant(NA_real_, dim = new_dim, dimnames = new_dn, units = units(flq))
  for (a in seq_len(n)) out[,,,,a] <- flq[,,,,1]
  out
}

fix_spwn <- function(stk) {
  n <- dims(stk)$area
  stk@harvest.spwn <- expand_to_fleets(stk@harvest.spwn, n)
  stk@m.spwn <- expand_to_fleets(stk@m.spwn, n)
  stk
}

apply_catch_state_groups <- function(ctrl) {
  ctrl@plus.group[] <- TRUE

  ctrl@states["catch N_Chile",] <- c(1:7, rep(8, 5))
  ctrl@states["catch SC_Chile_PS",] <- c(1:8, rep(9, 4)) + 101
  ctrl@states["catch FarNorth",] <- c(1:5, rep(6, 7)) + 201
  ctrl@states["catch Offshore_Trawl",] <- c(1:7, rep(8, 5)) + 301

  ctrl
}

apply_catch_variance_groups <- function(ctrl) {
  ctrl@f.vars["catch N_Chile",] <- c(0, 0, 0, 1, 1, rep(2, 7))
  ctrl@f.vars["catch SC_Chile_PS",] <- c(0, 0, 0, rep(2, 9)) + 101
  ctrl@f.vars["catch FarNorth",] <- c(0, 0, 1, 2, rep(3, 8)) + 201
  ctrl@f.vars["catch Offshore_Trawl",] <- c(rep(0, 5), rep(1, 7)) + 301

  ctrl@obs.vars["catch N_Chile",] <- c(1, rep(2, 11))
  ctrl@obs.vars["catch SC_Chile_PS",] <- c(rep(1, 12)) + 101
  ctrl@obs.vars["catch FarNorth",] <- c(1, 1, 1, rep(2, 9)) + 201
  ctrl@obs.vars["catch Offshore_Trawl",] <- c(rep(1, 12)) + 301

  ctrl
}

apply_survey_groups <- function(ctrl, idx_sub) {
  obs_map <- build_obs_map(names(idx_sub))
  bio_map <- survey_biomass_treat[names(idx_sub)]

  for (nm in names(idx_sub)) {
    ctrl@obs.vars[nm, 1] <- obs_map[nm]
    ctrl@biomassTreat[which(names(ctrl@fleets) == nm)] <- bio_map[nm]
  }

  ctrl
}

build_ctrl <- function(stk, idx_sub, residuals = FALSE,
                       use_catch_var_groups = TRUE) {
  ctrl <- FLSAM.control(stk, idx_sub)
  ctrl <- apply_catch_state_groups(ctrl)
  if (use_catch_var_groups) {
    ctrl <- apply_catch_variance_groups(ctrl)
  }
  ctrl <- apply_survey_groups(ctrl, idx_sub)
  ctrl <- update(ctrl)
  ctrl@residuals <- residuals
  ctrl
}

build_ctrl_alt <- function(stk, idx_sub, residuals = FALSE) {
  build_ctrl(stk, idx_sub,
    residuals = residuals,
    use_catch_var_groups = FALSE
  )
}

parse_jjm_ssb <- function(rep_file, max_year) {
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
  mat <- do.call(rbind, lapply(dlines, function(line) {
    as.numeric(strsplit(line, "\\s+")[[1]])
  }))

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

prepare_diagnostics_context <- function(stk_h1, idx_h1, sam_h1,
                                        jjm_rep_file = "../assessment/results/h1_1.14_1_R.rep") {
  survey_fleets <- names(idx_h1)
  survey_fleets_base <- unique(base_fleet(survey_fleets))
  max_yr <- range(stk_h1)["maxyear"]

  stk_h1_fit <- stk_h1 + sam_h1
  res_all <- residuals(sam_h1)
  ssb_df <- ssb(sam_h1)

  sel.pat <- merge(f(sam_h1), fbar(sam_h1), by = "year", suffixes = c(".f", ".fbar"))
  sel.pat$sel <- sel.pat$value.f / sel.pat$value.fbar
  sel.pat$age <- as.numeric(as.character(sel.pat$age))
  sel.pat$pentad <- sprintf("%i's", floor(sel.pat$year / 5) * 5)

  jjm_ssb <- parse_jjm_ssb(jjm_rep_file, as.integer(max_yr))

  list(
    catch_fleets = catch_fleets,
    survey_fleets = survey_fleets,
    survey_fleets_base = survey_fleets_base,
    max_yr = max_yr,
    base_yr_chr = as.character(max_yr),
    stk_h1_fit = stk_h1_fit,
    res_all = res_all,
    ssb_df = ssb_df,
    sel.pat = sel.pat,
    jjm_ssb = jjm_ssb
  )
}

inject_context <- function(ctx, envir = parent.frame()) {
  list2env(ctx, envir = envir)
  invisible(ctx)
}

render_quarto_report <- function(input = "sam_diagnostics.qmd") {
  if (requireNamespace("quarto", quietly = TRUE)) {
    return(quarto::quarto_render(input))
  }

  status <- system2("quarto", c("render", input))
  if (!identical(status, 0L)) {
    stop("quarto render failed with status ", status)
  }
  invisible(status)
}
