## =============================================================================
## model4_harmonised.R  --  model4.R re-run on the harmonised 2,564-person Model 4 frame
##
## Writes to review/reproduction/harmonised/ (same file names and columns as the originals):
##   RR_model4.csv, RR_model4_contrast.csv, RR_model4_pseudoclass.csv, RR_model4_multiplicity.csv,
##   RR_model4_evalues.csv, RR_contrast_pseudoclass.csv, RR_attenuation_comparison.csv,
##   RR_power_confirmatory.csv, model4_cause_specific.csv (cause-specific Model 4, printed only in
##   the original), model4_typed_vs_exact.csv, analysed_n_model4.csv
##
## What changes relative to code/reconstructed/model4.R (original steps 252-257, 278):
##   1. the analysis frame: frame_2564.R (analysis_frame_primary.rds + Model 4 covariates,
##      bio_m = pmax(LBXSAL_m, egfr_m), design on all 70,190, domain = 2,564 complete cases on C4)
##      instead of the CSV-rebuilt frame without modal imputation (2,557 analysed);
##   2. RR_model4_evalues.csv: the E-value for a confidence limit is set to 1 when the interval
##      includes 1 (original finding 4).
## Everything else is kept as in the original, including how the Python cells took their inputs:
## hazard ratios and limits as printed with 2 decimals and P values as printed by R (3 significant
## digits within the printed block) -- only now they are taken from the harmonised fits by code
## instead of being typed by hand. model4_typed_vs_exact.csv shows what unrounded inputs would give.
## Random numbers: same seeds and call order as the original (set.seed(20260907) before the 100
## Model 4 draws; again before the 2 x 100 contrast draws), so the draws are the original ones.
## Run   : /opt/homebrew/bin/Rscript code/reconstructed/harmonised/model4_harmonised.R  (project root)
## =============================================================================
source("code/reconstructed/harmonised/frame_2564.R")
suppressPackageStartupMessages(library(poLCA))
S <- "model4"

## ---- step 254: modal-assignment specifications ------------------------------------------------
specs <- list(
 "Model 3 (reference)"                              = C3,
 "Model 3 + antidepressant use"                     = paste(C3,"antidep_i",sep="+"),
 "Model 4: + comorbidity and frailty block"         = paste(C3,B_com,sep="+"),
 "Model 4 + antidepressant use"                     = paste(C3,B_com,"antidep_i",sep="+"),
 "Model 4b: individual conditions + antidepressant" = paste(C3,B_full,"antidep_i",sep="+"))
out <- list()
for(nm in names(specs)){
  m <- svycoxph(as.formula(paste0("Surv(time,event)~lca+",specs[[nm]])), design=P)
  log_fit(S, "RR_model4.csv", nm, m)
  jp <- tryCatch(regTermTest(m,~lca)$p, error=function(e) NA)
  out[[nm]] <- cbind(model=nm, tid(m,"^lca"), joint_p=as.numeric(jp), npar=length(coef(m)))
}
M <- do.call(rbind,out); rownames(M) <- NULL
print(transform(M, HR95=sprintf("%.2f (%.2f-%.2f)",HR,lo,hi))[,c("model","term","HR95","p","joint_p","npar")],
      row.names=FALSE, digits=3)
write.csv(M, file.path(HOUT, "RR_model4.csv"), row.names=FALSE)

## ---- step 255: 100 pseudo-class draws under Model 4, contrast, cause-specific -----------------
set.seed(20260907)
fits <- readRDS(file.path(DER, "lca_fits_cancer.rds")); post <- fits[[4]]$posterior
LAB <- c("High symptom burden","Low symptom burden","Insomnia-fatigue","Somatic-depressive")
LEV <- c("Low symptom burden","Insomnia-fatigue","Somatic-depressive","High symptom burden")
aid <- readRDS(file.path(DER, "analysis_frame_primary.rds")); aid <- aid[aid$inAnalysis==1,]
stopifnot(nrow(post)==nrow(aid))
pool_draws <- function(M, rhs, ref = NULL, output, label) {
  est <- lapply(1:M, function(i){
    cl <- apply(post,1,function(pp) sample.int(length(pp),1,prob=pp))
    f <- factor(LAB[cl], levels=LEV); if (!is.null(ref)) f <- relevel(f, ref=ref)
    dd <- data.frame(SEQN=aid$SEQN, lca_d=f)
    Hh <- merge(A[,setdiff(names(A),"lca_d")], dd, by="SEQN", all.x=TRUE)
    mm <- svycoxph(as.formula(paste0("Surv(time,event)~lca_d+",rhs)), design=subset(mk(Hh), analysed==1))
    k <- grep("^lca_d", names(coef(mm)))
    list(b=coef(mm)[k], v=diag(vcov(mm))[k], nm=gsub("^lca_d","",names(coef(mm))[k]), n=mm$n, d=mm$nevent)})
  log_draws(S, output, label, sapply(est,`[[`,"n"), sapply(est,`[[`,"d"))
  B<-do.call(rbind,lapply(est,`[[`,"b")); V<-do.call(rbind,lapply(est,`[[`,"v"))
  qb<-colMeans(B); ub<-colMeans(V); bv<-apply(B,2,var); tot<-ub+(1+1/M)*bv
  data.frame(term=est[[1]]$nm, HR=exp(qb), lo=exp(qb-1.96*sqrt(tot)), hi=exp(qb+1.96*sqrt(tot)),
             p=2*pnorm(-abs(qb/sqrt(tot))), fmi=round((1+1/M)*bv/tot,3))}
pc4 <- pool_draws(100, C4, output="RR_model4_pseudoclass.csv", label="Model 4 + antidepressant, 100 draws (each)")
cat("\n=== Model 4 + antidepressant + 100 pseudo-class draws ===\n")
print(transform(pc4, HR95=sprintf("%.2f (%.2f-%.2f)",HR,lo,hi))[,c("term","HR95","p","fmi")], row.names=FALSE, digits=3)
write.csv(pc4, file.path(HOUT, "RR_model4_pseudoclass.csv"), row.names=FALSE)
## direct contrast under Model 4 (and under Model 3, needed for the attenuation table)
A2 <- A; A2$lca_if <- relevel(A2$lca, ref="Insomnia-fatigue")
P_if <- subset(mk(A2), analysed==1)
mif4 <- log_fit(S, "RR_model4_contrast.csv", "Model 4 + antidepressant, reference insomnia-fatigue",
                svycoxph(as.formula(paste0("Surv(time,event)~lca_if+",C4)), design=P_if))
mif3 <- log_fit(S, "RR_attenuation_comparison.csv", "Model 3, reference insomnia-fatigue",
                svycoxph(as.formula(paste0("Surv(time,event)~lca_if+",C3)), design=P_if))
ct4 <- tid(mif4,"^lca_if"); ct3 <- tid(mif3,"^lca_if")
cat("\n=== Model 4: contrasts against insomnia-fatigue ===\n")
print(transform(ct4,HR95=sprintf("%.2f (%.2f-%.2f)",HR,lo,hi))[,c("term","HR95","p")], row.names=FALSE, digits=3)
write.csv(ct4, file.path(HOUT, "RR_model4_contrast.csv"), row.names=FALSE)
## cause-specific Model 4 (printed only in the original; typed into step 256)
cs4 <- list()
for(oc in c("ev_ca","ev_cvd")){
  m <- log_fit(S, "RR_model4_multiplicity.csv", paste("Model 4 + antidepressant,", oc),
               svycoxph(as.formula(sprintf("Surv(time,%s)~lca+%s",oc,C4)), design=P))
  cs4[[oc]] <- tid(m,"^lca")
  cat("\n",oc,": "); print(transform(cs4[[oc]],HR95=sprintf("%.2f (%.2f-%.2f)",HR,lo,hi))[,c("term","HR95","p")],
      row.names=FALSE,digits=3)}
write.csv(rbind(cbind(outcome="ev_ca", cs4$ev_ca), cbind(outcome="ev_cvd", cs4$ev_cvd)),
          file.path(HOUT, "model4_cause_specific.csv"), row.names=FALSE)

## ---- step 256 (Python -> R): inputs taken from the fits as the original typed them ------------
M4ad <- M[M$model=="Model 4 + antidepressant use",]
Mp_typed <- typed_p(M$p)                     # the 15-row RR_model4 print (3 significant digits)
allc <- data.frame(outcome="All-cause", term=M4ad$term, HR=typed_hr(M4ad$HR), lo=typed_hr(M4ad$lo),
                   hi=typed_hr(M4ad$hi), p=Mp_typed[M$model=="Model 4 + antidepressant use"])
csb <- function(x, lab) data.frame(outcome=lab, term=x$term, HR=typed_hr(x$HR), lo=typed_hr(x$lo),
                                   hi=typed_hr(x$hi), p=typed_p(x$p))
rows <- rbind(allc, csb(cs4$ev_ca, "Cancer"), csb(cs4$ev_cvd, "CVD/stroke"))
Fm <- rows; Fm$p_bh <- bh_statsmodels(Fm$p); Fm$p_bonf <- p.adjust(Fm$p, "bonferroni")
write_py_csv(Fm, file.path(HOUT, "RR_model4_multiplicity.csv"))
cat("\n=== Model 4 main family (9 tests) ===\n"); print(Fm, row.names=FALSE)
rr <- function(hr) (1-0.5^sqrt(hr))/(1-0.5^sqrt(1/hr))
ev <- function(r){ x <- max(r,1/r); x+sqrt(x*(x-1)) }
## E-value for the confidence limit closest to the null; 1 if the interval includes 1 (FIX 2)
ev_ci <- function(hr, lo, hi) if (hr >= 1) { if (lo <= 1) 1 else ev(rr(lo)) } else { if (hi >= 1) 1 else ev(rr(hi)) }
sdrow <- function(x) x[x$term=="Somatic-depressive",]
m3sd <- sdrow(M[M$model=="Model 3 (reference)",]); m4sd <- sdrow(M4ad)
m4bsd <- sdrow(M[M$model=="Model 4b: individual conditions + antidepressant",]); pcsd <- sdrow(pc4); ct4sd <- sdrow(ct4)
Ein <- data.frame(estimate=c("Model 3 (previous primary)","Model 4 + antidepressant (new primary)",
                             "Model 4b individual conditions","Model 4 + pseudo-class (most conservative)",
                             "Direct contrast under Model 4"),
                  HR=typed_hr(c(m3sd$HR, m4sd$HR, m4bsd$HR, pcsd$HR, ct4sd$HR)),
                  CI_lower=typed_hr(c(m3sd$lo, m4sd$lo, m4bsd$lo, pcsd$lo, ct4sd$lo)),
                  CI_upper=typed_hr(c(m3sd$hi, m4sd$hi, m4bsd$hi, pcsd$hi, ct4sd$hi)))
E <- data.frame(estimate=Ein$estimate, HR=Ein$HR, CI_lower=Ein$CI_lower,
                approx_RR=py_round(sapply(Ein$HR, rr),2),
                E_value_point=py_round(sapply(sapply(Ein$HR, rr), ev),2),
                E_value_CI=py_round(mapply(ev_ci, Ein$HR, Ein$CI_lower, Ein$CI_upper),2))
write_py_csv(E, file.path(HOUT, "RR_model4_evalues.csv"))
cat("\n=== E-values (CI limit set to 1 when the interval includes 1) ===\n"); print(E, row.names=FALSE)

## ---- step 257: contrasts with 2 x 100 draws; attenuation -----------------------------------------
set.seed(20260907)
cn3 <- pool_draws(100, C3, ref="Insomnia-fatigue", output="RR_contrast_pseudoclass.csv", label="Model 3, 100 draws (each)")
cn4 <- pool_draws(100, C4, ref="Insomnia-fatigue", output="RR_contrast_pseudoclass.csv", label="Model 4 + antidepressant, 100 draws (each)")
cat("\n=== contrasts (reference = insomnia-fatigue), 100 pseudo-class draws ===\n")
print(rbind(cbind(cov="Model 3",cn3), cbind(cov="Model 4 + antidepressant",cn4)), row.names=FALSE, digits=3)
write.csv(rbind(cbind(cov="Model 3",cn3),cbind(cov="Model 4 + antidepressant",cn4)),
          file.path(HOUT, "RR_contrast_pseudoclass.csv"), row.names=FALSE)
ct3sd <- sdrow(ct3)
cmp <- data.frame(comparison=c("vs low symptom burden","vs low symptom burden","vs insomnia-fatigue","vs insomnia-fatigue"),
                  model=c("Model 3","Model 4 + antidep","Model 3","Model 4 + antidep"),
                  HR=typed_hr(c(m3sd$HR, m4sd$HR, ct3sd$HR, ct4sd$HR)), lo=typed_hr(c(m3sd$lo, m4sd$lo, ct3sd$lo, ct4sd$lo)),
                  hi=typed_hr(c(m3sd$hi, m4sd$hi, ct3sd$hi, ct4sd$hi)))
cmp$attenuation_pct <- NA
cmp$attenuation_pct[2] <- round((cmp$HR[1]-cmp$HR[2])/(cmp$HR[1]-1)*100,1)
cmp$attenuation_pct[4] <- round((cmp$HR[3]-cmp$HR[4])/(cmp$HR[3]-1)*100,1)
print(cmp, row.names=FALSE)
write.csv(cmp, file.path(HOUT, "RR_attenuation_comparison.csv"), row.names=FALSE)

## ---- step 278, second half: sample size (inputs as printed, 2 decimals) -------------------------
## the original multiplied by the cohort size 2,569 / 676 deaths; kept (the analysed sample is 2,564 / 675)
power_tab <- function(hr, lo, hi, hr2, lo2, hi2) {
  b <- log(hr); se <- (log(hi)-log(lo))/(2*1.96); se_need <- abs(b)/(1.959964+0.8416212); ratio <- (se/se_need)^2
  b2 <- log(hr2); se2 <- (log(hi2)-log(lo2))/(2*1.96)
  data.frame(target=c(sprintf("assignment-corrected effect (HR %s)", format(hr)), sprintf("modal-assignment effect (HR %s)", format(hr2))),
             logHR=c(py_round(b,4), py_round(b2,4)), SE=c(py_round(se,4), py_round(se2,4)),
             z=c(py_round(b/se,2), py_round(b2/se2,2)),
             SE_for_80pct_power=c(py_round(se_need,4), py_round(abs(b2)/2.8016,4)),
             sample_multiplier=c(py_round(ratio,1), 1.0),
             participants_needed=c(as.integer(trunc(2569*ratio)), 2569L),
             deaths_needed=c(as.integer(trunc(676*ratio)), 676L)) }
PW <- power_tab(typed_hr(pcsd$HR), typed_hr(pcsd$lo), typed_hr(pcsd$hi), typed_hr(m4sd$HR), typed_hr(m4sd$lo), typed_hr(m4sd$hi))
write_py_csv(PW, file.path(HOUT, "RR_power_confirmatory.csv"))
print(PW, row.names=FALSE)

## ---- information only: the same quantities from unrounded inputs -----------------------------------
PWx <- power_tab(pcsd$HR, pcsd$lo, pcsd$hi, m4sd$HR, m4sd$lo, m4sd$hi)
TX <- data.frame(
  quantity=c("E-value point, Model 4 + antidepressant", "E-value CI, Model 4 + antidepressant",
             "E-value point, Model 4b", "E-value CI, Model 4b", "E-value point, Model 4 + pseudo-class",
             "E-value point, direct contrast Model 4", "E-value CI, direct contrast Model 4",
             "attenuation % vs low symptom burden", "attenuation % vs insomnia-fatigue",
             "sample multiplier for 80% power", "Bonferroni P, Model 4 all-cause somatic-depressive"),
  from_printed_inputs=c(E$E_value_point[2], E$E_value_CI[2], E$E_value_point[3], E$E_value_CI[3], E$E_value_point[4],
                        E$E_value_point[5], E$E_value_CI[5], cmp$attenuation_pct[2], cmp$attenuation_pct[4],
                        PW$sample_multiplier[1], Fm$p_bonf[2]),
  from_unrounded_inputs=c(round(ev(rr(m4sd$HR)),2), round(ev_ci(m4sd$HR,m4sd$lo,m4sd$hi),2),
                          round(ev(rr(m4bsd$HR)),2), round(ev_ci(m4bsd$HR,m4bsd$lo,m4bsd$hi),2), round(ev(rr(pcsd$HR)),2),
                          round(ev(rr(ct4sd$HR)),2), round(ev_ci(ct4sd$HR,ct4sd$lo,ct4sd$hi),2),
                          round((m3sd$HR-m4sd$HR)/(m3sd$HR-1)*100,1), round((ct3sd$HR-ct4sd$HR)/(ct3sd$HR-1)*100,1),
                          round(PWx$sample_multiplier[1],1), min(1, 9*m4sd$p)))
write.csv(TX, file.path(HOUT, "model4_typed_vs_exact.csv"), row.names=FALSE)
print(TX, row.names=FALSE)

## ---- analysed n / deaths ------------------------------------------------------------------------
write_nlog(S)
cat("model4_harmonised.R done\n")
