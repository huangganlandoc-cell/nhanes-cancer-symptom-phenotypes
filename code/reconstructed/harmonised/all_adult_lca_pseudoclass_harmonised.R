## =============================================================================
## all_adult_lca_pseudoclass_harmonised.R  --  the Cox part of all_adult_lca.R on the harmonised frame
##
## Writes to review/reproduction/harmonised/:
##   RR_fullsample_lca_pseudoclass.csv, analysed_n_all_adult_lca.csv
## Not re-run (no Cox model, independent of the covariate frame): lca_fit_all_adults.csv and
##   RR_fullsample_lca_profiles.csv.
##
## The all-adult latent class fits are read from review/reproduction/output/intermediate/
## lca_fits_all_adults.rds (written by code/reconstructed/all_adult_lca.R, which reproduced the
## original fits exactly); if that file is missing they are refitted here with the original seed
## (about 30 minutes).
## Change relative to the original step 283: the Cox models use the
## harmonised 2,564-person domain (frame_2564.R) instead of the CSV-rebuilt frame (2,557 analysed).
## The draws themselves are unchanged: set.seed(20260907), one draw per primary-cohort member
## (2,569) per iteration, 50 iterations, reference = the most frequent class in each draw.
## Run   : /opt/homebrew/bin/Rscript code/reconstructed/harmonised/all_adult_lca_pseudoclass_harmonised.R
## =============================================================================
source("code/reconstructed/harmonised/frame_2564.R")
suppressPackageStartupMessages(library(poLCA))
S <- "all_adult_lca"
ALL <- read.csv(file.path(DER, "nhanes_all_adults_lca_input.csv.gz"))
ff <- "review/reproduction/output/intermediate/lca_fits_all_adults.rds"
if (file.exists(ff)) {
  fitsA <- readRDS(ff); cat("all-adult fits read from", ff, "\n")
} else {
  cat("refitting the all-adult latent class models (about 30 minutes)\n")
  set.seed(20260907)
  f <- cbind(phqi1,phqi2,phqi3,phqi4,phqi5,phqi6,phqi7,phqi8,phqi9,sleepcat,slq050b)~1
  fitsA <- lapply(1:7, function(k) poLCA(f, ALL, nclass=k, maxiter=10000, nrep=20, verbose=FALSE, na.rm=TRUE))
}
stopifnot(round(fitsA[[4]]$llik, 1) == -107837.6)       # the original k = 4 solution

## ---- step 283 (Cox part), harmonised domain ---------------------------------------------------
ftA <- fitsA[[4]]; postA <- ftA$posterior
idx <- match(A$SEQN[A$primary==1], ALL$SEQN)
pS <- postA[idx[!is.na(idx)],,drop=FALSE]
cat("survivors matched:", sum(!is.na(idx)), "/", sum(A$primary==1),
    "| mean max posterior (primary cohort):", round(mean(apply(pS,1,max)),3), "\n")
set.seed(20260907)
seqn_ok <- A$SEQN[A$primary==1][!is.na(idx)]
poolA <- function(M, rhs){
  est <- lapply(1:M, function(i){
    cl <- apply(pS,1,function(pp) sample.int(length(pp),1,prob=pp))
    dd <- data.frame(SEQN=seqn_ok, lcaA=factor(paste0("C",cl), levels=paste0("C",1:4)))
    Hh <- merge(A[,setdiff(names(A),"lcaA")], dd, by="SEQN", all.x=TRUE)
    Hh$lcaA <- relevel(Hh$lcaA, ref=names(which.max(table(dd$lcaA))))
    mm <- svycoxph(as.formula(paste0("Surv(time,event)~lcaA+",rhs)), design=subset(mk(Hh), analysed==1 & !is.na(lcaA)))
    k <- grep("^lcaA", names(coef(mm)))
    list(b=coef(mm)[k], v=diag(vcov(mm))[k], nm=gsub("^lcaA","",names(coef(mm))[k]), n=mm$n, d=mm$nevent)})
  log_draws(S, "RR_fullsample_lca_pseudoclass.csv", "Model 4 + antidepressant, 50 draws (each)",
            sapply(est,`[[`,"n"), sapply(est,`[[`,"d"))
  B<-do.call(rbind,lapply(est,`[[`,"b")); V<-do.call(rbind,lapply(est,`[[`,"v"))
  qb<-colMeans(B); ub<-colMeans(V); bv<-apply(B,2,var); tot<-ub+(1+1/M)*bv
  data.frame(term=est[[1]]$nm, HR=exp(qb), lo=exp(qb-1.96*sqrt(tot)), hi=exp(qb+1.96*sqrt(tot)),
             p=2*pnorm(-abs(qb/sqrt(tot))), fmi=round((1+1/M)*bv/tot,3))}
pcA <- poolA(50, C4)
print(transform(pcA,HR95=sprintf("%.2f (%.2f-%.2f)",HR,lo,hi))[,c("term","HR95","p","fmi")],row.names=FALSE,digits=3)
write.csv(pcA, file.path(HOUT, "RR_fullsample_lca_pseudoclass.csv"), row.names=FALSE)
write_nlog(S)
cat("all_adult_lca_pseudoclass_harmonised.R done\n")
