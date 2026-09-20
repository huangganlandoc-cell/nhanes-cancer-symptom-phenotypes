## =============================================================================
## compare_outputs.R  --  compare every regenerated CSV with the copy in supporting/
##
## Reads : supporting/<target>.csv                      (the original outputs; read-only)
##         review/reproduction/output/<target>.csv      (written by the other scripts here)
## Writes: review/reproduction/output/comparison.csv   one row per target
##         review/reproduction/output/comparison_by_column.csv   one row per target x column
##
## Rules
##   * Both files are read with read.csv(check.names = FALSE); a nameless first column (row
##     names written by R's write.csv or a pandas index) is compared like any text column.
##   * Numeric columns: maximum absolute difference and maximum relative difference
##     |a - b| / max(|a|, |b|) over all cells (NA in both files counts as equal; NA in only one
##     file counts as a mismatch).
##   * Text / logical columns: exact match, cell by cell.
##   * Status: "完全一致" (identical)       all text identical and max relative difference <= 1e-10
##                                          (i.e. only floating-point noise in the 11th+ digit)
##             "差异可忽略" (negligible)    text identical, numeric max relative difference <= 1e-4
##                                          (cannot change any figure reported to 3 significant digits)
##             "不一致" (different)         anything larger, a text mismatch, or different shape
##             "无法运行" (not produced)    the regenerated file does not exist
## Run   : /opt/homebrew/bin/Rscript code/reconstructed/compare_outputs.R  (from the project root)
## =============================================================================
if (!dir.exists("supporting")) stop("run from the project root")
invisible(suppressWarnings(Sys.setlocale("LC_CTYPE", "en_US.UTF-8")))   # so the Chinese status labels print readably
OUT <- "review/reproduction/output"

targets <- c(
  "RR4_inclusive_profiles.csv", "RR4_inclusive_pseudoclass.csv", "RR4_localdep_cox.csv",
  "RR4_measurement_fit.csv", "RR4_threshold_cox.csv", "RR_attenuation_comparison.csv",
  "RR_contrast_pseudoclass.csv", "RR_crosstab_vs_prior.csv", "RR_dimensional.csv",
  "RR_dimensional_affective.csv", "RR_dimensional_multiplicity.csv", "RR_efa_loadings.csv",
  "RR_fullsample_lca_profiles.csv", "RR_fullsample_lca_pseudoclass.csv", "RR_incremental_joint.csv",
  "RR_model4.csv", "RR_model4_contrast.csv", "RR_model4_evalues.csv", "RR_model4_multiplicity.csv",
  "RR_model4_pseudoclass.csv", "RR_power_confirmatory.csv", "RR_vs_prior_classification.csv",
  "lca_fit_all_adults.csv", "lca_fit_indices_cancer.csv", "ph_period_specific.csv",
  "ph_stratified_baseline.csv", "ph_tests.csv", "sens_evalues.csv", "sens_pseudoclass.csv",
  "sens_pseudoclass_M100.csv")
script_of <- c(
  RR4_inclusive_profiles.csv = "measurement_sensitivity_RR4.R", RR4_inclusive_pseudoclass.csv = "measurement_sensitivity_RR4.R",
  RR4_localdep_cox.csv = "measurement_sensitivity_RR4.R", RR4_measurement_fit.csv = "measurement_sensitivity_RR4.R",
  RR4_threshold_cox.csv = "measurement_sensitivity_RR4.R", RR_attenuation_comparison.csv = "model4.R",
  RR_contrast_pseudoclass.csv = "model4.R", RR_crosstab_vs_prior.csv = "prior_classification.R",
  RR_dimensional.csv = "dimensional.R", RR_dimensional_affective.csv = "dimensional.R",
  RR_dimensional_multiplicity.csv = "dimensional.R", RR_efa_loadings.csv = "dimensional.R",
  RR_fullsample_lca_profiles.csv = "all_adult_lca.R", RR_fullsample_lca_pseudoclass.csv = "all_adult_lca.R",
  RR_incremental_joint.csv = "incremental_value.R", RR_model4.csv = "model4.R",
  RR_model4_contrast.csv = "model4.R", RR_model4_evalues.csv = "model4.R",
  RR_model4_multiplicity.csv = "model4.R", RR_model4_pseudoclass.csv = "model4.R",
  RR_power_confirmatory.csv = "model4.R", RR_vs_prior_classification.csv = "prior_classification.R",
  lca_fit_all_adults.csv = "all_adult_lca.R", lca_fit_indices_cancer.csv = "fit_indices.R",
  ph_period_specific.csv = "proportional_hazards.R", ph_stratified_baseline.csv = "proportional_hazards.R",
  ph_tests.csv = "proportional_hazards.R", sens_evalues.csv = "pseudoclass.R",
  sens_pseudoclass.csv = "pseudoclass.R", sens_pseudoclass_M100.csv = "pseudoclass.R")

rd <- function(f) read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
relmax <- function(a, b) {
  d <- abs(a - b); s <- pmax(abs(a), abs(b)); r <- ifelse(d == 0, 0, d / s); max(c(0, r), na.rm = TRUE)
}
cmp_one <- function(t) {
  fo <- file.path("supporting", t); fn <- file.path(OUT, t)
  base <- data.frame(target = t, script = script_of[[t]], rows = NA, cols = NA, numeric_cols = NA,
                     text_cols = NA, text_cells_mismatched = NA, na_pattern_mismatches = NA,
                     max_abs_diff = NA, max_rel_diff = NA, byte_identical = NA, status = NA, stringsAsFactors = FALSE)
  if (!file.exists(fn)) { base$status <- "无法运行"; return(list(summary = base, cols = NULL)) }
  base$byte_identical <- identical(readBin(fo, "raw", file.size(fo)), readBin(fn, "raw", file.size(fn)))
  A <- rd(fo); B <- rd(fn)
  base$rows <- nrow(A); base$cols <- ncol(A)
  if (!identical(dim(A), dim(B)) || !identical(names(A), names(B))) {
    base$status <- "不一致"; base$text_cells_mismatched <- NA
    return(list(summary = base, cols = data.frame(target = t, column = "(shape or header differs)",
      type = NA, max_abs_diff = NA, max_rel_diff = NA, mismatches = NA))) }
  cr <- list(); nm_bad <- 0; na_bad <- 0; mabs <- 0; mrel <- 0; ncol_num <- 0; ncol_txt <- 0
  for (j in seq_along(A)) {
    a <- A[[j]]; b <- B[[j]]
    if (is.numeric(a) && is.numeric(b)) {
      ncol_num <- ncol_num + 1
      nab <- sum(is.na(a) != is.na(b)); ok <- !is.na(a) & !is.na(b)
      ma <- if (any(ok)) max(abs(a[ok] - b[ok])) else 0; mr <- if (any(ok)) relmax(a[ok], b[ok]) else 0
      na_bad <- na_bad + nab; mabs <- max(mabs, ma); mrel <- max(mrel, mr)
      cr[[j]] <- data.frame(target = t, column = names(A)[j], type = "numeric", max_abs_diff = ma,
                            max_rel_diff = mr, mismatches = nab)
    } else {
      ncol_txt <- ncol_txt + 1
      mm <- sum(!(as.character(a) == as.character(b) | (is.na(a) & is.na(b))), na.rm = TRUE) +
            sum(xor(is.na(a), is.na(b)))
      nm_bad <- nm_bad + mm
      cr[[j]] <- data.frame(target = t, column = names(A)[j], type = "text", max_abs_diff = NA,
                            max_rel_diff = NA, mismatches = mm)
    }
  }
  base$numeric_cols <- ncol_num; base$text_cols <- ncol_txt
  base$text_cells_mismatched <- nm_bad; base$na_pattern_mismatches <- na_bad
  base$max_abs_diff <- mabs; base$max_rel_diff <- mrel
  base$status <- if (nm_bad > 0 || na_bad > 0) "不一致" else if (mrel <= 1e-10) "完全一致" else
                 if (mrel <= 1e-4) "差异可忽略" else "不一致"
  list(summary = base, cols = do.call(rbind, cr))
}
res <- lapply(targets, cmp_one)
S <- do.call(rbind, lapply(res, `[[`, "summary")); C <- do.call(rbind, lapply(res, `[[`, "cols"))
S$status_en <- c("完全一致" = "identical", "差异可忽略" = "negligible difference",
                 "不一致" = "different", "无法运行" = "not produced")[S$status]
write.csv(S, file.path(OUT, "comparison.csv"), row.names = FALSE)
write.csv(C, file.path(OUT, "comparison_by_column.csv"), row.names = FALSE)
op <- options(width = 200)
print(transform(S, max_abs_diff = signif(max_abs_diff, 3), max_rel_diff = signif(max_rel_diff, 3))[,
      c("target", "script", "rows", "text_cells_mismatched", "na_pattern_mismatches", "max_abs_diff",
        "max_rel_diff", "byte_identical", "status")], row.names = FALSE)
options(op)
cat("\nstatus counts:\n"); print(table(S$status))
