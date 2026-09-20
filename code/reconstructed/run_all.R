## =============================================================================
## run_all.R  --  run every reconstruction script, then the comparison
##
## Each script runs in its own fresh R process (so no object can leak from one to the next),
## single-threaded, with its screen output saved to review/reproduction/output/logs/<script>.log.
## The all-adult LCA (about half an hour) runs last.
## Run   : /opt/homebrew/bin/Rscript code/reconstructed/run_all.R   (from the project root)
## =============================================================================
if (!dir.exists("code/reconstructed")) stop("run from the project root")
LOG <- "review/reproduction/output/logs"; dir.create(LOG, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
rscript <- file.path(R.home("bin"), "Rscript")
scripts <- c("fit_indices", "pseudoclass", "proportional_hazards", "incremental_value",
             "prior_classification", "derive_model4_covariates", "model4", "dimensional",
             "measurement_sensitivity_RR4", "all_adult_lca")
for (s in scripts) {
  t0 <- Sys.time(); cat(sprintf("%-30s ... ", s))
  st <- system2(rscript, file.path("code/reconstructed", paste0(s, ".R")),
                stdout = file.path(LOG, paste0(s, ".log")), stderr = file.path(LOG, paste0(s, ".log")))
  cat(if (st == 0) "ok" else paste("FAILED (exit", st, ") - see log"),
      sprintf("(%.1f min)\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
system2(rscript, "code/reconstructed/compare_outputs.R")
