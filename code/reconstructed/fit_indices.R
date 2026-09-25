## =============================================================================
## fit_indices.R  --  latent class enumeration in the 3,268 cancer survivors
##
## Regenerates (into review/reproduction/output/):
##   lca_fit_indices_cancer.csv        (source of Supplementary Table S1)
## Extra check written to review/reproduction/output/input_checks/:
##   check_lca_fits_cancer_refit.csv   (refitted k = 1..5 vs data/derived/lca_fits_cancer.rds)
## Writes data/derived/lca_fits_cancer.rds ONLY if that file does not exist yet (a run from the
## raw data); an existing file is never overwritten, and the check then compares the refit with it.
##
## Reconstructed from the original analysis log:
##   step 54: poLCA k = 1..5, seed 20260907, 30 random starts,
##            builds the fit-index table `ic`, saves lca_fits_cancer.rds
##   step 70: write.csv(ic, "lca_fit_indices_cancer.csv")
##
## Input : data/derived/nhanes_cancer_design_frame.csv.gz (the original read the uncompressed
##         copy; the later version of this file only has two extra columns, only_nms and any_mel,
##         which do not enter the LCA)
## Run   : /opt/homebrew/bin/Rscript code/reconstructed/fit_indices.R   (from the project root)
## =============================================================================
suppressPackageStartupMessages({library(poLCA); library(survey); library(survival)})
if (!file.exists("data/derived/nhanes_cancer_design_frame.csv.gz")) stop("run from the project root")
OUT <- "review/reproduction/output"; dir.create(file.path(OUT, "input_checks"), recursive = TRUE, showWarnings = FALSE)

## ---- step 54 (verbatim apart from the file path) -------------------------------------
G <- read.csv("data/derived/nhanes_cancer_design_frame.csv.gz"); a <- subset(G, inAnalysis==1)
set.seed(20260907)
f <- cbind(phqi1,phqi2,phqi3,phqi4,phqi5,phqi6,phqi7,phqi8,phqi9,sleepcat,slq050b)~1
fits <- lapply(1:5, function(k) poLCA(f, a, nclass=k, maxiter=8000, nrep=30, verbose=FALSE, na.rm=TRUE))
## original: saveRDS(fits,"lca_fits_cancer.rds")  -- saved here only when the file is absent;
## otherwise the refit is compared with the saved object below.
if (!file.exists("data/derived/lca_fits_cancer.rds")) {
  saveRDS(fits, "data/derived/lca_fits_cancer.rds")
  cat("data/derived/lca_fits_cancer.rds did not exist: saved the refit (the check below is then trivial)\n")
}
ic <- data.frame(k=1:5, logLik=round(sapply(fits,`[[`,"llik"),1), AIC=round(sapply(fits,`[[`,"aic"),1),
  BIC=round(sapply(fits,`[[`,"bic"),1),
  cAIC=round(sapply(fits,function(x) -2*x$llik+x$npar*(log(x$Nobs)+1)),1),
  entropy=round(sapply(fits,function(x){p<-x$posterior; if(is.null(dim(p)))return(NA)
     1-sum(-p*log(pmax(p,1e-12)))/(nrow(p)*log(ncol(p)))}),3),
  min_class_pct=sapply(fits,function(x) round(min(x$P)*100,1)))
cat("n =",nrow(a)," deaths =",sum(a$event),"\n"); print(ic,row.names=FALSE)

## ---- step 70 -------------------------------------------------------------------------
write.csv(ic, file.path(OUT, "lca_fit_indices_cancer.csv"), row.names=FALSE)

## ---- check: does the refit reproduce the saved fits used by every later analysis? ----
saved <- readRDS("data/derived/lca_fits_cancer.rds")
chk <- do.call(rbind, lapply(1:5, function(k) data.frame(
  k = k,
  loglik_refit = fits[[k]]$llik, loglik_saved = saved[[k]]$llik,
  max_abs_diff_class_share = max(abs(fits[[k]]$P - saved[[k]]$P)),
  max_abs_diff_item_probs  = max(abs(unlist(fits[[k]]$probs) - unlist(saved[[k]]$probs))),
  max_abs_diff_posterior   = max(abs(fits[[k]]$posterior - saved[[k]]$posterior)))))
print(chk, row.names = FALSE, digits = 6)
write.csv(chk, file.path(OUT, "input_checks", "check_lca_fits_cancer_refit.csv"), row.names = FALSE)
cat("fit_indices.R done\n")
