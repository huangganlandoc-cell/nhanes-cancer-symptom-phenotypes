## =============================================================================
## run_all_harmonised.R  --  run every harmonised script, then the old-versus-new comparison
##
## Each script runs in its own R process, single-threaded; screen output goes to
## review/reproduction/harmonised/logs/<script>.log. Total about 2 minutes, provided
## review/reproduction/output/intermediate/lca_fits_all_adults.rds exists (written by
## code/reconstructed/all_adult_lca.R); otherwise add about 30 minutes for the all-adult refit.
## Run   : /opt/homebrew/bin/Rscript code/reconstructed/harmonised/run_all_harmonised.R  (project root)
## =============================================================================
if (!dir.exists("code/reconstructed/harmonised")) stop("run from the project root")
LOG <- "review/reproduction/harmonised/logs"; dir.create(LOG, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
rscript <- file.path(R.home("bin"), "Rscript")
scripts <- c("model4_harmonised", "dimensional_harmonised", "all_adult_lca_pseudoclass_harmonised",
             "proportional_hazards_fulldesign", "diagnostics_harmonised")
for (s in scripts) {
  t0 <- Sys.time(); cat(sprintf("%-40s ... ", s))
  lf <- file.path(LOG, paste0(s, ".log"))
  st <- system2(rscript, file.path("code/reconstructed/harmonised", paste0(s, ".R")), stdout = lf, stderr = lf)
  cat(if (st == 0) "ok" else paste("FAILED (exit", st, ") - see log"),
      sprintf("(%.1f min)\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
fs <- setdiff(list.files("review/reproduction/harmonised", pattern = "^analysed_n_.*[.]csv$", full.names = TRUE),
              "review/reproduction/harmonised/analysed_n_all.csv")
write.csv(do.call(rbind, lapply(fs, read.csv, colClasses = "character")),
          "review/reproduction/harmonised/analysed_n_all.csv", row.names = FALSE)
system2(rscript, "code/reconstructed/harmonised/compare_harmonised.R",
        stdout = file.path(LOG, "compare_harmonised.log"), stderr = file.path(LOG, "compare_harmonised.log"))
cat("comparison written to review/reproduction/harmonised/value_changes.csv\n")
