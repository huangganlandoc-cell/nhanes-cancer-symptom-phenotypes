## Writes the classification-error matrices and the implementation checks of the ML three-step
## analysis to CSV (for Supplementary Table S11), recomputing each check rather than copying logs.
## Run from the project root after `Rscript code/analysis_three_step.R fit`:
##   Rscript code/export_three_step_checks.R
commandArgs <- function(...) "none"           # load data and helpers only
source("code/analysis_three_step.R")
rm(commandArgs)
res <- readRDS(file.path(OUT, "RR5_three_step_fits.rds"))
write.csv(round(Dp, 4), file.path(OUT, "RR5_error_matrix_primary.csv"))
write.csv(round(Df, 4), file.path(OUT, "RR5_error_matrix_full_cohort.csv"))

ck <- list()
add <- function(check, result, detail = "") ck[[length(ck) + 1]] <<- data.frame(check, result, detail)
e0 <- ml3step(dat, C3, diag(4), sw, LEVS)
br3 <- coxph(as.formula(paste("Surv(time, event) ~ lca +", C3)), data = dat, weights = sw, ties = "breslow")
add("Error matrix set to the identity reproduces modal assignment (Model 3)",
    sprintf("max |difference in log HR| %.1e", max(abs(e0$coef - coef(br3)[grep("^lca", names(coef(br3)))]))))
for (m in c("C3", "C4")) {
  e <- res[[m]]$fit; lab <- c(C3 = "Model 3", C4 = "Model 4")[m]
  add(sprintf("Full-sample EM, %s", lab), sprintf("converged in %d iterations", e$iter),
      sprintf("pseudo-log-likelihood %.3f; largest membership coefficient %.1f; penalty %s", e$loglik,
              e$max_abs_mn, ifelse(res[[m]]$decay > 0, res[[m]]$decay, "none")))
  alt <- sapply(1:5, function(r) {
    set.seed(100 + r); mix <- runif(1, 0.2, 0.8)
    pp <- as.matrix(dat[, LEVS]); pp <- mix * pp + (1 - mix) * matrix(rgamma(length(pp), 1), nrow(pp))
    ee <- ml3step(dat, get(m), Dp, sw, LEVS, init_post = pp, decay = res[[m]]$decay)
    c(ee$loglik, exp(ee$coef["clsSomaticDepressive"]))
  })
  add(sprintf("Five perturbed starting values, %s", lab),
      sprintf("somatic-depressive HR %.3f to %.3f (main fit %.3f)", min(alt[2, ]), max(alt[2, ]),
              exp(e$coef["clsSomaticDepressive"])),
      sprintf("pseudo-log-likelihood %.3f to %.3f", min(alt[1, ]), max(alt[1, ])))
}
write.csv(do.call(rbind, ck), file.path(OUT, "RR5_fit_checks.csv"), row.names = FALSE)
print(do.call(rbind, ck))
