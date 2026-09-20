## =============================================================================
## compare_harmonised.R  --  old (supporting/) versus harmonised (review/reproduction/harmonised/)
##
## Writes review/reproduction/harmonised/value_changes.csv: one row per numeric cell, with the
## old and new value, the absolute difference, and both values rounded the way such numbers are
## displayed in the manuscript and tables (hazard ratios, limits, E-values: 2 decimals; P values:
## 3 and 4 decimals; fraction of missing information: 2 and 3 decimals). The flag columns say
## whether the rounded values differ. Section B pairs: RR3_* files (2,557 frame) versus their
## harmonised counterparts (RR4_* files, already on the 2,564 frame, and diagnostics_*.csv).
## Run   : /opt/homebrew/bin/Rscript code/reconstructed/harmonised/compare_harmonised.R (project root)
## =============================================================================
if (!dir.exists("supporting")) stop("run from the project root")
HOUT <- "review/reproduction/harmonised"
rd <- function(f) read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
pairs <- list(
  list("supporting/RR_model4.csv",                     "RR_model4.csv",                     c("model","term")),
  list("supporting/RR_model4_contrast.csv",            "RR_model4_contrast.csv",            "term"),
  list("supporting/RR_model4_pseudoclass.csv",         "RR_model4_pseudoclass.csv",         "term"),
  list("supporting/RR_model4_multiplicity.csv",        "RR_model4_multiplicity.csv",        c("outcome","term")),
  list("supporting/RR_model4_evalues.csv",             "RR_model4_evalues.csv",             "estimate"),
  list("supporting/RR_contrast_pseudoclass.csv",       "RR_contrast_pseudoclass.csv",       c("cov","term")),
  list("supporting/RR_attenuation_comparison.csv",     "RR_attenuation_comparison.csv",     c("comparison","model")),
  list("supporting/RR_power_confirmatory.csv",         "RR_power_confirmatory.csv",         ".row"),
  list("supporting/RR_dimensional.csv",                "RR_dimensional.csv",                c("model","term")),
  list("supporting/RR_dimensional_affective.csv",      "RR_dimensional_affective.csv",      c("block","model","term")),
  list("supporting/RR_dimensional_multiplicity.csv",   "RR_dimensional_multiplicity.csv",   c("outcome","term")),
  list("supporting/RR_fullsample_lca_pseudoclass.csv", "RR_fullsample_lca_pseudoclass.csv", "term"),
  list("supporting/ph_period_specific.csv",            "ph_period_specific.csv",            c("period","term")),
  list("supporting/ph_stratified_baseline.csv",        "ph_stratified_baseline.csv",        "term"),
  list("supporting/ph_tests.csv",                      "ph_tests.csv",                      "test"),
  list("supporting/RR_incremental_joint.csv",          "RR_incremental_joint.csv",          "model"),
  list("supporting/RR_incremental_value.csv",          "RR_incremental_value.csv",          c("model","term")),
  ## section B: RR3_* (2,557 frame) versus the same analyses on the 2,564 frame (the RR4_* files)
  list("supporting/RR3_localdep_cox.csv",              "supporting/RR4_localdep_cox.csv",   "term"),
  list("supporting/RR3_threshold_cox.csv",             "supporting/RR4_threshold_cox.csv",  "term"),
  list("supporting/RR3_inclusive_pseudoclass.csv",     "supporting/RR4_inclusive_pseudoclass.csv", c("model","term")))
kind <- function(col) if (col %in% c("HR","lo","hi","approx_RR","E_value_point","E_value_CI","CI_lower")) "hr" else
  if (col %in% c("p","p_bh","p_bonf","joint_p","joint_p_phenotype","joint_p_11_indicators")) "p" else
  if (col == "fmi") "fmi" else "other"
out <- list()
for (pp in pairs) {
  O <- rd(pp[[1]]); N <- rd(if (startsWith(pp[[2]], "supporting/")) pp[[2]] else file.path(HOUT, pp[[2]]))
  key <- pp[[3]]
  if (identical(key, ".row")) { O$.row <- seq_len(nrow(O)); N$.row <- seq_len(nrow(N)) }
  ko <- do.call(paste, c(O[key], sep = " | ")); kn <- do.call(paste, c(N[key], sep = " | "))
  for (i in seq_len(nrow(N))) {
    j <- match(kn[i], ko)
    for (col in setdiff(names(N), key)) {
      nv <- N[[col]][i]; ov <- if (is.na(j) || !col %in% names(O)) NA else O[[col]][j]
      if (!is.numeric(nv)) { if (!identical(as.character(nv), as.character(ov)))
        out[[length(out)+1]] <- data.frame(old_file = pp[[1]], new_file = pp[[2]], row = kn[i], column = col,
          kind = "text", old = as.character(ov), new = as.character(nv), abs_diff = NA, old_disp = as.character(ov),
          new_disp = as.character(nv), changed_2dp = NA, changed_3dp = NA, changed_4dp = NA); next }
      k <- kind(col)
      r2 <- function(x) sprintf("%.2f", x); r3 <- function(x) sprintf("%.3f", x); r4 <- function(x) sprintf("%.4f", x)
      out[[length(out)+1]] <- data.frame(old_file = pp[[1]], new_file = pp[[2]], row = kn[i], column = col, kind = k,
        old = if (is.na(ov)) NA_character_ else format(ov, digits = 15), new = format(nv, digits = 15),
        abs_diff = if (is.na(ov)) NA else abs(nv - ov),
        old_disp = if (is.na(ov)) NA else if (k == "hr") r2(ov) else if (k == "p") r4(ov) else if (k == "fmi") r3(ov) else format(ov),
        new_disp = if (k == "hr") r2(nv) else if (k == "p") r4(nv) else if (k == "fmi") r3(nv) else format(nv),
        changed_2dp = if (is.na(ov)) NA else r2(ov) != r2(nv),
        changed_3dp = if (is.na(ov)) NA else r3(ov) != r3(nv),
        changed_4dp = if (is.na(ov)) NA else r4(ov) != r4(nv))
    }
  }
}
V <- do.call(rbind, out)
write.csv(V, file.path(HOUT, "value_changes.csv"), row.names = FALSE)
cat("cells compared:", nrow(V), "| changed at 2 dp:", sum(V$changed_2dp, na.rm = TRUE),
    "| at 3 dp:", sum(V$changed_3dp, na.rm = TRUE), "| at 4 dp:", sum(V$changed_4dp, na.rm = TRUE), "\n")
op <- options(width = 250)
print(V[(V$kind == "hr" & V$changed_2dp %in% TRUE) | (V$kind == "p" & V$changed_3dp %in% TRUE) |
        (V$kind %in% c("fmi","other") & V$changed_3dp %in% TRUE) | is.na(V$old),
        c("new_file","row","column","old","new","old_disp","new_disp")], row.names = FALSE)
options(op)
