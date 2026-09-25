## Somatic-depressive versus sleep-fatigue phenotype at a fixed PHQ-9 total score (post hoc, 2026-09-25).
## The two intermediate phenotypes differ in severity as well as in pattern; this is the 1-df contrast
## between them with the summed score in the model (linear per 5 points, and a 4-knot spline), under
## Model 3 and Model 4a, modal assignment, primary cohort.
## Run from the project root: Rscript code/analysis_fixed_score_contrast.R
## Output: supporting/RR7_fixed_score_contrast.csv
suppressPackageStartupMessages({library(survey); library(survival); library(rms)})
options(survey.lonely.psu = "adjust")
DER <- "data/derived"; if (!dir.exists(DER)) stop("run from the project root")
OCSV <- Sys.getenv("OUT", "supporting"); dir.create(OCSV, recursive = TRUE, showWarnings = FALSE)

A  <- readRDS(file.path(DER, "analysis_frame_primary.rds"))
G4 <- read.csv(file.path(DER, "nhanes_model4_design_frame.csv.gz"))
A  <- merge(A, G4[, c("SEQN", "comorb_n_i", "egfr_i", "egfr_m")], by = "SEQN", all.x = TRUE)
A$cvd <- as.numeric(A$cvd)
C3  <- "age+sex+race4+educ+married+pir_i+pir_m+bmi_i+bmi_m+smoke+pa3+htn+dm+cvd+ydx_i+ydx_m+multi_primary+catype"
C4a <- paste(C3, "comorb_n_i+egfr_i+egfr_m", sep = "+")
P <- subset(svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt, nest = TRUE, data = A), primary == 1)

SD <- "lcaHypersomnia-somatic"; SF <- "lcaInsomnia-fatigue"
rows <- list()
for (cov in c("C3", "C4a")) for (sev in c("I(phq9_score/5)", "rcs(phq9_score,4)")) {
  m <- svycoxph(as.formula(paste("Surv(time, event) ~ lca +", sev, "+", get(cov))), design = P)
  stopifnot(all(c(SD, SF) %in% names(coef(m))))
  for (k in list(c(SD, "", "somatic-depressive vs low symptom burden"),
                 c(SD, SF, "somatic-depressive vs sleep-fatigue"))) {
    w <- setNames(rep(0, length(coef(m))), names(coef(m))); w[k[1]] <- 1; if (nzchar(k[2])) w[k[2]] <- -1
    ct <- svycontrast(m, w); b <- as.numeric(coef(ct)); se <- as.numeric(SE(ct))
    rows[[length(rows) + 1]] <- data.frame(
      model = c(C3 = "Model 3", C4a = "Model 4a")[cov],
      severity = c("I(phq9_score/5)" = "PHQ-9 score, linear per 5 points", "rcs(phq9_score,4)" = "PHQ-9 score, 4-knot spline")[sev],
      contrast = k[3], HR = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se), p = 2 * pnorm(-abs(b / se)),
      n = m$n, deaths = m$nevent, row.names = NULL)
  }
}
out <- do.call(rbind, rows)
print(transform(out, HR95 = sprintf("%.2f (%.2f-%.2f)", HR, lo, hi))[, c("model", "severity", "contrast", "HR95", "p")],
      row.names = FALSE, digits = 3)
write.csv(out, file.path(OCSV, "RR7_fixed_score_contrast.csv"), row.names = FALSE)
