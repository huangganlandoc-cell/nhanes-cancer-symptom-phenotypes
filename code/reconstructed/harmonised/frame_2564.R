## =============================================================================
## frame_2564.R  --  the ONE harmonised Model 4 analysis frame, sourced by every script in
##                   code/reconstructed/harmonised/ (not meant to be run on its own)
##
## Construction (identical to the RR4_* frame, steps 406-409, and to
## code/analysis_three_step.R):
##   * data/derived/analysis_frame_primary.rds (all 70,190 NHANES participants; educ, married and
##     smoke already modally imputed for the survivors; phenotype = 4-class modal assignment)
##   * + the Model 4 covariates of data/derived/nhanes_model4_design_frame.csv.gz
##     (plus the seven single-condition indicators, needed only for Model 4b)
##   * bio_m = pmax(LBXSAL_m, egfr_m)  (one missingness indicator for the biochemistry panel)
##   * analysed = primary cohort AND complete cases on time, event and every Model 4 (C4) covariate
##     -> 2,564 survivors, 675 deaths (the 5 excluded have missing hypertension or diabetes status)
##   * survey design built on all 70,190 participants, then subset(design, analysed == 1)
## Every model in the harmonised scripts is fitted on this same 2,564-person domain, including
## Model 2 and Model 3 fits, so that all estimates in a table share one analytic sample.
## Outputs go to review/reproduction/harmonised/ only.
## =============================================================================
suppressPackageStartupMessages({library(survey); library(survival)})
if (!file.exists("data/derived/analysis_frame_primary.rds")) stop("run from the project root")
options(survey.lonely.psu = "adjust")
DER  <- "data/derived"
HOUT <- "review/reproduction/harmonised"
dir.create(file.path(HOUT, "logs"), recursive = TRUE, showWarnings = FALSE)

A  <- readRDS(file.path(DER, "analysis_frame_primary.rds"))
G4 <- read.csv(file.path(DER, "nhanes_model4_design_frame.csv.gz"))
M4KEEP <- c("SEQN", "comorb_n_i", "LBXHGB_i", "LBXHGB_m", "LBXSAL_i", "LBXSAL_m", "egfr_i", "egfr_m",
            "func_lim_i", "n_rx_i", "n_rx_m", "antidep_i",
            "hf_i", "stroke_i", "emphysema_i", "bronchitis_i", "liver_i", "arthritis_i", "ckd_i")
A <- merge(A, G4[, M4KEEP], by = "SEQN", all.x = TRUE)
A$cvd   <- as.numeric(A$cvd)
A$bio_m <- pmax(A$LBXSAL_m, A$egfr_m)
## phenotype labels as in the original output files (only the renamed third class is changed)
levels(A$lca)[levels(A$lca) == "Hypersomnia-somatic"] <- "Somatic-depressive"

C2 <- "age+sex+race4+educ+married+pir_i+pir_m"
C3 <- paste(C2, "bmi_i+bmi_m+smoke+pa3+htn+dm+cvd+ydx_i+ydx_m+multi_primary+catype", sep = "+")
B_com  <- "comorb_n_i+LBXHGB_i+LBXHGB_m+LBXSAL_i+egfr_i+bio_m+func_lim_i+n_rx_i+n_rx_m"
B_full <- paste("hf_i+stroke_i+emphysema_i+bronchitis_i+liver_i+arthritis_i+ckd_i",
                "LBXHGB_i+LBXHGB_m+LBXSAL_i+egfr_i+bio_m+func_lim_i+n_rx_i+n_rx_m", sep = "+")
C4 <- paste(C3, B_com, "antidep_i", sep = "+")

cc <- complete.cases(A[, c("time", "event", all.vars(as.formula(paste("~", C4))))])
A$analysed <- as.numeric(A$primary == 1 & cc)
stopifnot(sum(A$analysed) == 2564, sum(A$event[A$analysed == 1]) == 675)

mk  <- function(d) svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt, nest = TRUE, data = d)
des <- mk(A)
P   <- subset(des, analysed == 1)

tid <- function(m, pat) { s <- summary(m)$coef; k <- grep(pat, rownames(s)); ci <- confint(m)[k, , drop = FALSE]
  data.frame(term = gsub(pat, "", rownames(s)[k]), HR = exp(s[k, "coef"]), lo = exp(ci[, 1]), hi = exp(ci[, 2]),
             p = s[k, ncol(s)]) }

## analysed n / deaths of every fitted model, collected per script
.NLOG <- list()
log_fit <- function(script, output, model, m) {
  .NLOG[[length(.NLOG) + 1]] <<- data.frame(script = script, output = output, model = model,
                                            n = as.character(m$n), deaths = as.character(m$nevent))
  invisible(m)
}
log_draws <- function(script, output, model, n, d)      # pseudo-class draws: n / deaths of every draw
  .NLOG[[length(.NLOG) + 1]] <<- data.frame(script = script, output = output, model = model,
                                            n = paste(unique(n), collapse = "/"), deaths = paste(unique(d), collapse = "/"))
write_nlog <- function(script) {
  x <- unique(do.call(rbind, .NLOG)); f <- file.path(HOUT, sprintf("analysed_n_%s.csv", script))
  write.csv(x, f, row.names = FALSE); cat("\nanalysed n / deaths per model ->", f, "\n"); print(x, row.names = FALSE)
}

## the way numbers were typed into the original Python cells: HR / CI as printed with "%.2f",
## P values as printed by R with digits = 3 within the printed block
typed_hr <- function(x) round(x, 2)
typed_p  <- function(p) as.numeric(format(p, digits = 3))
py_round <- function(x, d) round(x * 10^d) / 10^d      # numpy / pandas round (half to even on x*10^d)
py_num <- function(x) vapply(x, function(v) {          # Python repr() of a float, as pandas writes it
  if (is.na(v)) return("")
  for (d in 1:17) { s <- formatC(v, digits = d, format = "g"); if (as.numeric(s) == v) break }
  s <- trimws(s); if (!grepl("[.e]", s)) s <- paste0(s, ".0"); s }, "")
write_py_csv <- function(df, file) {                    # pandas DataFrame.to_csv(index=False)
  cols <- lapply(df, function(col) {
    if (is.logical(col)) ifelse(is.na(col), "", ifelse(col, "True", "False"))
    else if (is.integer(col)) ifelse(is.na(col), "", as.character(col))
    else if (is.numeric(col)) py_num(col)
    else { s <- as.character(col); s[is.na(s)] <- ""; q <- grepl('[,"\n]', s)
           s[q] <- paste0('"', gsub('"', '""', s[q]), '"'); s } })
  writeLines(c(paste(names(df), collapse = ","), do.call(paste, c(cols, sep = ","))), file)
}
bh_statsmodels <- function(p) { n <- length(p); o <- order(p); raw <- p[o] / (seq_len(n) / n)
  adj <- pmin(rev(cummin(rev(raw))), 1); out <- numeric(n); out[o] <- adj; out }
cat(sprintf("harmonised frame: %d analysed survivors, %d deaths (primary cohort %d / %d)\n",
            sum(A$analysed), sum(A$event[A$analysed == 1]), sum(A$primary == 1), sum(A$event[A$primary == 1])))
