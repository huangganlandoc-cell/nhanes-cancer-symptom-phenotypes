## =============================================================================
## pseudoclass.R  --  pseudo-class draws from the symptom-only posteriors (Model 3),
##                    and the E-values that were computed from those results
##
## Regenerates (into review/reproduction/output/):
##   sens_pseudoclass.csv        20 draws, full cohort and primary cohort
##   sens_pseudoclass_M100.csv   100 draws, primary cohort, then 100 further draws, full cohort
##   sens_evalues.csv            E-values (typed-in estimates, see below)
##
## Reconstructed from the original analysis log:
##   step 96: builds G2 (analysis_frame.rds + only_nms/any_mel),
##            C2/C3 strings and mkdes(); its landmark and
##            non-melanoma-skin models wrote other files and are
##            not repeated here (they use no random numbers)
##   step 97: set.seed(20260907); M = 20 draws; sens_pseudoclass.csv
##   step 149: set.seed(20260907); poolM(100) on the primary cohort
##   step 151: poolM(100, primary_only = FALSE) WITHOUT re-seeding,
##             i.e. draws 101-200 of the same random stream;
##             writes the final sens_pseudoclass_M100.csv
##   step 150: Python/numpy E-values, translated to R
##
## Notes on faithfulness
##   * survey.lonely.psu was still at its default ("fail") when these cells ran; it is not set here.
##   * The draw order (all 20 draw sets first in step 97; draws interleaved with model fits in
##     steps 149/151) is kept exactly, so the draws are the same random numbers as in the original.
##   * sens_evalues.csv: the original Python cell typed the hazard ratios in by hand
##     (1.80/1.39, 1.45/1.05, 1.59/1.08, 2.11/1.55); the same numbers are used here, and the
##     script prints the corresponding model estimates next to them as a transcription check.
##     Python's round() (numpy, round-half-even on x*100) is emulated by py_round().
## Run   : /opt/homebrew/bin/Rscript code/reconstructed/pseudoclass.R   (from the project root)
## =============================================================================
suppressPackageStartupMessages({library(survey); library(survival); library(poLCA)})
if (!file.exists("data/derived/analysis_frame.rds")) stop("run from the project root")
OUT <- "review/reproduction/output"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
DER <- "data/derived"

## ---- step 96 (only the objects that step 97 uses) -------------------------------------
G2 <- readRDS(file.path(DER, "analysis_frame.rds"))
nf <- read.csv(file.path(DER, "nhanes_cancer_design_frame.csv.gz"))[,c("SEQN","only_nms","any_mel")]
G2 <- merge(G2, nf, by="SEQN", all.x=TRUE)
C2 <- "age+sex+race4+educ+married+pir_i+pir_m"
C3 <- paste(C2,"bmi_i+bmi_m+smoke+pa3+htn+dm+cvd+ydx_i+ydx_m+multi_primary+catype",sep="+")
mkdes <- function(d) svydesign(ids=~SDMVPSU,strata=~SDMVSTRA,weights=~wt,nest=TRUE,data=d)
hr <- function(m){s<-summary(m)$coef;k<-grep("^lca",rownames(s));ci<-confint(m)[k,,drop=FALSE]
  data.frame(term=gsub("^lca","",rownames(s)[k]),HR=exp(s[k,"coef"]),lo=exp(ci[,1]),hi=exp(ci[,2]),p=s[k,ncol(s)])}

## ---- step 97 (verbatim apart from file paths) ----------------------------------------
## 3) pseudo-class draws: correct for latent class assignment uncertainty
set.seed(20260907)
fits <- readRDS(file.path(DER, "lca_fits_cancer.rds")); post <- fits[[4]]$posterior
LAB <- c("High symptom burden","Low symptom burden","Insomnia-fatigue","Hypersomnia-somatic")
a_id <- readRDS(file.path(DER, "analysis_frame.rds")); a_id <- a_id[a_id$inAnalysis==1,]
stopifnot(nrow(post)==nrow(a_id))
M <- 20
draws <- lapply(1:M, function(m){
  cl <- apply(post,1,function(p) sample.int(length(p),1,prob=p))
  data.frame(SEQN=a_id$SEQN, lca_d=factor(LAB[cl],
    levels=c("Low symptom burden","Insomnia-fatigue","Hypersomnia-somatic","High symptom burden")))
})
pool <- function(keepexpr){
  est <- lapply(draws, function(dd){
    H <- merge(G2[,setdiff(names(G2),"lca_d")], dd, by="SEQN", all.x=TRUE)
    H$.k <- as.numeric(eval(keepexpr, H))
    s <- subset(mkdes(H), .k==1)
    m <- svycoxph(as.formula(paste0("Surv(time,event)~lca_d+",C3)), design=s)
    k <- grep("^lca_d", names(coef(m)))
    list(b=coef(m)[k], v=diag(vcov(m))[k], nm=gsub("^lca_d","",names(coef(m))[k]))
  })
  B  <- do.call(rbind, lapply(est,`[[`,"b")); V <- do.call(rbind, lapply(est,`[[`,"v"))
  qb <- colMeans(B); ub <- colMeans(V); bv <- apply(B,2,var)
  tot <- ub + (1+1/M)*bv
  data.frame(term=est[[1]]$nm, HR=exp(qb), lo=exp(qb-1.96*sqrt(tot)), hi=exp(qb+1.96*sqrt(tot)),
             p=2*pnorm(-abs(qb/sqrt(tot))), fmi=round((1+1/M)*bv/tot,3))
}
pc_full <- cbind(cohort="full", pool(quote(inAnalysis==1)))
pc_excl <- cbind(cohort="excl non-melanoma-skin-only",
                 pool(quote(inAnalysis==1 & !(only_nms %in% c(TRUE,"True")))))
pc <- rbind(pc_full, pc_excl); write.csv(pc, file.path(OUT, "sens_pseudoclass.csv"), row.names=FALSE)
print(transform(pc, HR95=sprintf("%.2f (%.2f-%.2f)",HR,lo,hi))[,c("cohort","term","HR95","p","fmi")],
      row.names=FALSE, digits=3)

## ---- step 149 ---------------------------------------------------------------------------
set.seed(20260907)
G3 <- readRDS(file.path(DER, "analysis_frame_primary.rds"))
fits <- readRDS(file.path(DER, "lca_fits_cancer.rds")); post <- fits[[4]]$posterior
LAB <- c("High symptom burden","Low symptom burden","Insomnia-fatigue","Somatic-depressive")
aid <- G3[G3$inAnalysis==1,]; stopifnot(nrow(post)==nrow(aid))
C2 <- "age+sex+race4+educ+married+pir_i+pir_m"
C3 <- paste(C2,"bmi_i+bmi_m+smoke+pa3+htn+dm+cvd+ydx_i+ydx_m+multi_primary+catype",sep="+")
mkdes <- function(d) svydesign(ids=~SDMVPSU,strata=~SDMVSTRA,weights=~wt,nest=TRUE,data=d)
poolM <- function(M, primary_only=TRUE){
  est <- lapply(1:M, function(m){
    cl <- apply(post,1,function(p) sample.int(length(p),1,prob=p))
    dd <- data.frame(SEQN=aid$SEQN, lca_d=factor(LAB[cl],
        levels=c("Low symptom burden","Insomnia-fatigue","Somatic-depressive","High symptom burden")))
    H <- merge(G3[,setdiff(names(G3),"lca_d")], dd, by="SEQN", all.x=TRUE)
    H$.k <- if(primary_only) as.numeric(H$primary==1) else as.numeric(H$inAnalysis==1)
    mm <- svycoxph(as.formula(paste0("Surv(time,event)~lca_d+",C3)), design=subset(mkdes(H), .k==1))
    k <- grep("^lca_d", names(coef(mm)))
    list(b=coef(mm)[k], v=diag(vcov(mm))[k], nm=gsub("^lca_d","",names(coef(mm))[k]))
  })
  B<-do.call(rbind,lapply(est,`[[`,"b")); V<-do.call(rbind,lapply(est,`[[`,"v"))
  qb<-colMeans(B); ub<-colMeans(V); bv<-apply(B,2,var); tot<-ub+(1+1/M)*bv
  data.frame(M=M, term=est[[1]]$nm, HR=exp(qb), lo=exp(qb-1.96*sqrt(tot)), hi=exp(qb+1.96*sqrt(tot)),
             p=2*pnorm(-abs(qb/sqrt(tot))), fmi=round((1+1/M)*bv/tot,3))
}
pc100 <- poolM(100)
## (step 149 wrote pc100 alone to sens_pseudoclass_M100.csv; step 151 overwrote that file)

## ---- step 151 (continues the random stream of step 149: no set.seed) --------------------
pc100_full <- poolM(100, primary_only=FALSE)
both <- rbind(cbind(cohort="full", pc100_full), cbind(cohort="primary", pc100))
write.csv(both, file.path(OUT, "sens_pseudoclass_M100.csv"), row.names=FALSE)
print(transform(both, HR95=sprintf("%.2f (%.2f-%.2f)",HR,lo,hi))[,c("cohort","term","HR95","p","fmi")],
      row.names=FALSE, digits=3)

## ---- step 150 (Python -> R) --------------------------------------------------------------
py_round <- function(x, d) round(x * 10^d) / 10^d      # numpy.round: x*10^d, round-half-even
py_num <- function(x) vapply(x, function(v) {           # Python repr() of a float, as pandas writes it
  if (is.na(v)) return("")
  for (d in 1:17) { s <- formatC(v, digits = d, format = "g"); if (as.numeric(s) == v) break }
  s <- trimws(s); if (!grepl("[.e]", s)) s <- paste0(s, ".0"); s }, "")
write_py_csv <- function(df, file) {                     # pandas DataFrame.to_csv(index=False)
  cols <- lapply(df, function(col) {
    if (is.logical(col)) ifelse(is.na(col), "", ifelse(col, "True", "False"))
    else if (is.integer(col)) ifelse(is.na(col), "", as.character(col))
    else if (is.numeric(col)) py_num(col)
    else { s <- as.character(col); s[is.na(s)] <- ""; q <- grepl('[,"\n]', s)
           s[q] <- paste0('"', gsub('"', '""', s[q]), '"'); s } })
  writeLines(c(paste(names(df), collapse = ","), do.call(paste, c(cols, sep = ","))), file)
}
rr_from_hr <- function(hr) (1-0.5^sqrt(hr))/(1-0.5^sqrt(1/hr))   # VanderWeele & Ding, common outcome
evalue <- function(rr){ rr <- max(rr,1/rr); rr+sqrt(rr*(rr-1)) }
typed <- data.frame(estimate=c("Somatic-depressive (primary, modal)","Somatic-depressive (pseudo-class M=100)",
                               "High symptom burden (primary, modal)",">=5 yr since dx, somatic-depressive"),
                    HR=c(1.80,1.45,1.59,2.11), CI_lower=c(1.39,1.05,1.08,1.55))
E <- do.call(rbind, lapply(seq_len(nrow(typed)), function(i){
  hr <- typed$HR[i]; lo <- typed$CI_lower[i]; rr <- rr_from_hr(hr); rrl <- rr_from_hr(lo)
  data.frame(estimate=typed$estimate[i], HR=hr, CI_lower=lo, approx_RR=py_round(rr,2),
             E_value_point=py_round(evalue(rr),2), E_value_CI=py_round(evalue(rrl),2))}))
write_py_csv(E, file.path(OUT, "sens_evalues.csv"))
print(E, row.names=FALSE)
cat(sprintf("\nevent rate %.1f%% (>15%%: common-outcome approximation used)\n", 676/2569*100))
## ---- transcription check: each typed HR / lower limit against the model it came from ------
## (run after all draws, so it cannot disturb the random-number stream)
sd100 <- pc100[pc100$term=="Somatic-depressive",]
m3p <- svycoxph(as.formula(paste0("Surv(time,event)~lca+",C3)), design=subset(mkdes(G3), primary==1))
r3 <- hr(m3p)
G3$.sg <- as.numeric(G3$primary==1 & G3$ydx_i>=5); G3$.sg[is.na(G3$.sg)] <- 0   # subgroup ">=5 yr since dx"
msg <- svycoxph(as.formula(paste0("Surv(time,event)~lca+",C2)), design=subset(mkdes(G3), .sg==1))  # Model 2, as all subgroups
rsg <- hr(msg)
comp <- rbind(unlist(r3[r3$term=="Hypersomnia-somatic", c("HR","lo")]), unlist(sd100[, c("HR","lo")]),
              unlist(r3[r3$term=="High symptom burden", c("HR","lo")]), unlist(rsg[rsg$term=="Hypersomnia-somatic", c("HR","lo")]))
TC <- data.frame(estimate=typed$estimate, typed_HR=typed$HR, computed_HR=comp[,1],
                 typed_lo=typed$CI_lower, computed_lo=comp[,2])
TC$agrees <- abs(round(TC$computed_HR,2)-TC$typed_HR) < 1e-9 & abs(round(TC$computed_lo,2)-TC$typed_lo) < 1e-9
cat("\n=== sens_evalues.csv transcription check ===\n"); print(TC, row.names=FALSE, digits=6)
dir.create(file.path(OUT, "input_checks"), showWarnings = FALSE)
write.csv(TC, file.path(OUT, "input_checks", "check_sens_evalues_transcription.csv"), row.names=FALSE)
cat("pseudoclass.R done\n")
