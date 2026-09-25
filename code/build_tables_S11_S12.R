## Supplementary Tables S11 (correction for misclassification of class membership) and
## S12 (bivariate residuals), assembled from saved outputs; no model is fitted here.
## Run from the project root after analysis_three_step.R (fit, sim, rep, rep3only, boot),
## analysis_inclusive_draws_score.R and analysis_bvr.R:  Rscript code/build_tables_S11_S12.R
OUT <- "supporting"
f3 <- function(x) sprintf("%.2f (%.2f-%.2f)", x$HR, x$lo, x$hi)
fp <- function(p) ifelse(p < 0.001, formatC(p, format = "e", digits = 1), sprintf("%.3f", p))
LAB <- c(Low = "Low symptom burden", SleepFatigue = "Sleep-fatigue",
         SomaticDepressive = "Somatic-depressive", High = "High symptom burden")
rows <- list()
add <- function(block, row, estimate, detail = "", p = "") rows[[length(rows) + 1]] <<-
  data.frame(Block = block, Row = row, Estimate = estimate, Detail = detail, P = p)

## A. classification-error matrix
D <- readRDS(file.path(OUT, "RR5_error_matrices.rds"))$Dp
for (t in rownames(D)) add("A. Classification-error matrix, primary cohort (probability of assigned class given true class)",
                           paste("True", LAB[t]),
                           paste(sprintf("%s %.3f", LAB[colnames(D)], D[t, ]), collapse = " | "))

## B-C. bias-adjusted three-step estimates
E <- read.csv(file.path(OUT, "RR5_three_step_estimates.csv"))
blk <- c(C3 = "B. Bias-adjusted three-step estimator, Model 3", C4 = "B. Bias-adjusted three-step estimator, Model 4",
         C3_fullD = "C. Sensitivity: error matrix from the full cohort, Model 3")
for (m in names(blk)) for (i in which(E$model == m))
  add(blk[m], E$term[i], f3(E[i, ]), sprintf("log HR %.3f, JKn SE %.3f", E$logHR[i], E$se[i]), fp(E$p[i]))

## D. implementation checks
ck <- read.csv(file.path(OUT, "RR5_fit_checks.csv"))
for (i in seq_len(nrow(ck))) add("D. Implementation checks", ck$check[i], ck$result[i], ck$detail[i])
MLAB <- c(C3 = "Model 3", C4 = "Model 4", C3_fullD = "Model 3, full-cohort error matrix")
for (v in list(c("", "Jackknife"), c("_bootstrap", "Rao-Wu bootstrap"))) {
  JD <- read.csv(file.path(OUT, paste0("RR5_jackknife_diagnostics", v[1], ".csv")))
  for (i in seq_len(nrow(JD))) add("D. Implementation checks",
    sprintf("%s (steps 1-3 repeated in each of %d replicates), %s: step-3 EM stopping", v[2], JD$replicates[i], MLAB[JD$model[i]]),
    sprintf("pseudo-log-likelihood %d | hazard ratios stable %d | not converged %d", JD$stopped_loglik[i],
            JD$stopped_coefficients[i], JD$not_converged[i]),
    sprintf("median iterations %g; largest membership coefficient %.1f", JD$median_iterations[i], JD$max_membership_coef[i]))
  JL <- read.csv(file.path(OUT, paste0("RR5_jackknife_lca_diagnostics", v[1], ".csv")))
  add("D. Implementation checks", sprintf("%s: measurement model re-fitted in every replicate", v[2]),
      sprintf("largest item-response difference from the main solution %.3f", JL$max_item_prob_difference),
      sprintf("assignments changed per replicate: median %g, maximum %g; somatic-depressive correctly assigned %s",
              JL$median_changed_assignments, JL$max_changed_assignments, JL$somatic_depressive_correctly_assigned_range))
}
E3 <- read.csv(file.path(OUT, "RR5_three_step_estimates_step3only.csv"))
EB <- read.csv(file.path(OUT, "RR5_three_step_estimates_bootstrap.csv"))
for (m in c("C3", "C4")) for (t in c("Somatic-depressive vs low", "Somatic-depressive vs sleep-fatigue")) {
  a <- E3[E3$model == m & E3$term == t, ]; b <- E[E$model == m & E$term == t, ]; k <- EB[EB$model == m & EB$term == t, ]
  add("D. Implementation checks", sprintf("Standard error of log HR, %s, %s", MLAB[m], t),
      sprintf("jackknife, step 3 only %.3f | jackknife, steps 1-3 %.3f | bootstrap, steps 1-3 %.3f", a$se, b$se, k$se),
      sprintf("95%% CI %.2f-%.2f | %.2f-%.2f | %.2f-%.2f", a$lo, a$hi, b$lo, b$hi, k$lo, k$hi))
}
## bootstrap SE without the unstable replicates: step-3 EM stopped at maxit, or a membership
## coefficient above 15 (the ridge trigger; Model 3 replicates are fitted without the penalty, as its
## full-sample fit). Variance as in analysis_three_step.R: sum((theta_r - estimate)^2) / (R - 1).
FB <- readRDS(file.path(OUT, "RR5_three_step_fits.rds")); BB <- readRDS(file.path(OUT, "RR5_jackknife_full_bootstrap.rds"))
for (m in c("C3", "C4")) {
  est <- FB[[m]]$fit$coef["clsSomaticDepressive"]
  x <- BB$theta[, paste0(m, ":clsSomaticDepressive")]
  nc <- BB$diag[, paste0(m, ":stop_rule")] == 3; dv <- BB$diag[, paste0(m, ":max_abs_mn")] > 15
  sev <- function(v) sqrt(sum((v - est)^2) / (length(v) - 1))
  ci <- function(s) sprintf("%.2f-%.2f", exp(est - 1.96 * s), exp(est + 1.96 * s))
  add("D. Implementation checks", sprintf("Bootstrap without unstable replicates, %s, Somatic-depressive vs low", MLAB[m]),
      sprintf("all 500: SE %.3f | without %d not converged: SE %.3f | without %d not converged or with a membership coefficient above 15: SE %.3f",
              sev(x), sum(nc), sev(x[!nc]), sum(nc | dv), sev(x[!(nc | dv)])),
      sprintf("95%% CI %s | %s | %s", ci(sev(x)), ci(sev(x[!nc])), ci(sev(x[!(nc | dv)]))))
}
se <- read.csv(file.path(OUT, "RR5_check2_se_comparison.csv"))
for (i in seq_len(nrow(se))) add("D. Implementation checks",
  sprintf("Modal model, %s, %s: SE of log HR", ifelse(se$model[i] == "C3", "Model 3", "Model 4"),
          LAB[sub("^lca", "", se$term[i])]), sprintf("Taylor %.3f | jackknife %.3f", se$se_taylor[i], se$se_jkn[i]))
sim <- read.csv(file.path(OUT, "RR5_simulation_check.csv"))
for (i in seq_len(nrow(sim))) add("D. Implementation checks",
  sprintf("Simulation (60 data sets), %s covariates, %s: true HR %.1f", sim$covariates[i],
          LAB[sim$phenotype[i]], sim$true_HR[i]),
  sprintf("modal %.2f | three-step %.2f", sim$modal_HR[i], sim$ml3_HR[i]),
  sprintf("three-step bias in log HR %.3f (Monte Carlo SE %.3f)", sim$ml3_bias_logHR[i], sim$mcse_ml3_bias[i]))

## E. inclusive draws
I <- read.csv(file.path(OUT, "RR5_inclusive_draws_score.csv"))
for (i in seq_len(nrow(I))) add(sprintf("E. Inclusive draws, %s", ifelse(I$model[i] == "C3", "Model 3", "Model 4")),
  I$term[i], f3(I[i, ]), sprintf("fraction of missing information %.2f", I$fmi[i]), fp(I$p[i]))
pr <- read.csv(file.path(OUT, "RR5_inclusive_profiles_somatic_depressive_score.csv"))
for (m in unique(pr$model)) add(sprintf("E. Inclusive draws, %s", ifelse(m == "C3", "Model 3", "Model 4")),
  "Somatic-depressive item profile: largest difference from the main solution",
  sprintf("%.3f", max(abs(pr$main[pr$model == m] - pr$inclusive[pr$model == m]))))

S11 <- do.call(rbind, rows)
write.csv(S11, "supplementary/TableS11_misclassification_correction.csv", row.names = FALSE)

## S12: bivariate residuals
b <- read.csv(file.path(OUT, "RR5_bvr.csv"))
S12 <- data.frame(`Indicator 1` = b$item1, `Indicator 2` = b$item2, df = b$df,
                  `Pearson X2` = sprintf("%.2f", b$X2), BVR = sprintf("%.2f", b$BVR),
                  Flag = ifelse(b$BVR > 3.84, "> 3.84", ""), check.names = FALSE)
write.csv(S12, "supplementary/TableS12_bivariate_residuals.csv", row.names = FALSE)
cat("Table S11:", nrow(S11), "rows | Table S12:", nrow(S12), "rows\n")
