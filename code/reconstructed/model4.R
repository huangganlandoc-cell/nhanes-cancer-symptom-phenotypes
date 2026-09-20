## =============================================================================
## model4.R  --  Model 4 (comorbidity, frailty, medication burden, antidepressant use)
##
## Regenerates (into review/reproduction/output/):
##   RR_model4.csv                Model 3 / +antidepressant / Model 4 / Model 4+AD / Model 4b
##   RR_model4_contrast.csv       somatic-depressive vs insomnia(sleep)-fatigue under Model 4
##   RR_model4_pseudoclass.csv    100 pseudo-class draws under Model 4
##   RR_model4_multiplicity.csv   BH / Bonferroni over 3 phenotypes x 3 outcomes (typed-in HRs)
##   RR_model4_evalues.csv        E-values (typed-in HRs)
##   RR_contrast_pseudoclass.csv  phenotype contrasts vs insomnia-fatigue, 100 draws, Models 3 and 4
##   RR_attenuation_comparison.csv  attenuation of the excess hazard (typed-in HRs)
##   RR_power_confirmatory.csv    sample size for a confirmatory study (typed-in HRs)
##
## Reconstructed from the original analysis log:
##   step 252 (R): frame H from nhanes_model4_design_frame + phenotype;
##                 covariates re-derived on H (see NOTE 1)
##   step 253 (R): first attempt; Model 4 singular (LBXSAL_m and egfr_m identical);
##                 its output file was overwritten by step 254 -> not repeated
##   step 254 (R): bio_m replaces the two identical indicators -> RR_model4.csv
##   step 255 (R): set.seed(20260907); 100 draws -> RR_model4_pseudoclass.csv;
##                 contrast -> RR_model4_contrast.csv; cause-specific Model 4 (printed)
##   step 256 (Python/statsmodels): multiplicity and E-values -> translated to R
##   step 257 (R): set.seed(20260907); contrasts with draws -> RR_contrast_pseudoclass.csv;
##                 RR_attenuation_comparison.csv
##   step 278 (Python), second half: RR_power_confirmatory.csv -> translated to R
##                                   (placed here because its inputs are Model 4 results)
## Inputs: data/derived/nhanes_model4_design_frame.csv.gz, analysis_frame_primary.rds,
##                                   lca_fits_cancer.rds. The Model 4 covariates were derived from NHANES modules that are
##                                   not in data/raw_nhanes; see derive_model4_covariates.R (not run).
##
## NOTE 1 (reproduced, not fixed): H is rebuilt from the CSV design frame, and unlike
##   analysis_frame*.rds it does NOT get the modal imputation of educ/married/smoke; tumour type
##   blanks become "Unknown" with Breast as reference. Survivors with a missing educ, married,
##   smoke, htn or dm value are therefore dropped silently by svycoxph (the log prints how many).
## NOTE 2: RR_model4_multiplicity / _evalues / RR_attenuation_comparison / RR_power_confirmatory
##   use hazard ratios that were typed in by hand from printed output. They are kept as typed;
##   this script recomputes every typed number and prints a transcription check (also written to
##   output/input_checks/check_model4_transcription.csv).
## NOTE 3: In RR_model4_evalues.csv the row "Model 4 + pseudo-class" has a confidence interval that
##   includes 1 (0.90-1.73); the original code still converted the lower limit (E = 1.36), whereas
##   by the usual convention the E-value for such an interval is 1. Reproduced as is.
## Run   : /opt/homebrew/bin/Rscript code/reconstructed/model4.R   (from the project root)
## =============================================================================
suppressPackageStartupMessages({library(survey); library(survival); library(poLCA)})
if (!file.exists("data/derived/nhanes_model4_design_frame.csv.gz")) stop("run from the project root")
OUT <- "review/reproduction/output"; dir.create(file.path(OUT, "input_checks"), recursive = TRUE, showWarnings = FALSE)
DER <- "data/derived"

## helpers for the two translated Python cells ------------------------------------------
py_round <- function(x, d) round(x * 10^d) / 10^d   # numpy round(): x*10^d, round half to even
py_num <- function(x) vapply(x, function(v) {       # Python repr() of a float, as pandas writes it
  if (is.na(v)) return("")
  for (d in 1:17) { s <- formatC(v, digits = d, format = "g"); if (as.numeric(s) == v) break }
  s <- trimws(s); if (!grepl("[.e]", s)) s <- paste0(s, ".0"); s }, "")
write_py_csv <- function(df, file) {                 # pandas DataFrame.to_csv(index=False)
  cols <- lapply(df, function(col) {
    if (is.logical(col)) ifelse(is.na(col), "", ifelse(col, "True", "False"))
    else if (is.integer(col)) ifelse(is.na(col), "", as.character(col))
    else if (is.numeric(col)) py_num(col)
    else { s <- as.character(col); s[is.na(s)] <- ""; q <- grepl('[,"\n]', s)
           s[q] <- paste0('"', gsub('"', '""', s[q]), '"'); s } })
  writeLines(c(paste(names(df), collapse = ","), do.call(paste, c(cols, sep = ","))), file)
}

## ---- step 252 ----------------------------------------------------------------------------
options(survey.lonely.psu="adjust")
H <- read.csv(file.path(DER, "nhanes_model4_design_frame.csv.gz"))
lc <- readRDS(file.path(DER, "analysis_frame_primary.rds"))[,c("SEQN","primary","lca","time","event","ev_ca","ev_cvd")]
levels(lc$lca)[levels(lc$lca)=="Hypersomnia-somatic"] <- "Somatic-depressive"
H <- merge(H[,setdiff(names(H),c("time","event","ev_ca","ev_cvd"))], lc, by="SEQN", all.x=TRUE)
H$primary[is.na(H$primary)] <- 0
H$sex <- factor(H$female, levels=c(1,0), labels=c("Female","Male"))
H$race4 <- factor(ifelse(H$race=="NH White","NH White", ifelse(H$race=="NH Black","NH Black",
             ifelse(H$race %in% c("Mexican American","Other Hispanic"),"Hispanic","Other"))),
             levels=c("NH White","NH Black","Hispanic","Other"))
H$educ <- factor(H$educ, levels=c(">High school","High school","<High school"))
H$married <- factor(H$married, levels=c("Married/partnered","Not partnered"))
H$smoke <- factor(H$smoke, levels=c("Never","Former","Current"))
H$pa3 <- factor(ifelse(is.na(H$pa_active),"Unknown", ifelse(H$pa_active==1,"Active","Inactive")),
                levels=c("Active","Inactive","Unknown"))
H$pir_m <- as.numeric(is.na(H$INDFMPIR)); H$pir_i <- ifelse(is.na(H$INDFMPIR), median(H$INDFMPIR,na.rm=TRUE), H$INDFMPIR)
H$bmi_m <- as.numeric(is.na(H$BMXBMI));  H$bmi_i <- ifelse(is.na(H$BMXBMI),  median(H$BMXBMI,na.rm=TRUE),  H$BMXBMI)
H$ydx_m <- as.numeric(is.na(H$yrs_dx));  H$ydx_i <- ifelse(is.na(H$yrs_dx),  median(H$yrs_dx,na.rm=TRUE),  H$yrs_dx)
H$cvd <- as.numeric(H$cvd); H$catype <- factor(ifelse(is.na(H$ca_type)|H$ca_type=="","Unknown",as.character(H$ca_type)))
H$catype <- relevel(H$catype, ref="Breast")
C2 <- "age+sex+race4+educ+married+pir_i+pir_m"
C3 <- paste(C2,"bmi_i+bmi_m+smoke+pa3+htn+dm+cvd+ydx_i+ydx_m+multi_primary+catype",sep="+")
mk <- function(d) svydesign(ids=~SDMVPSU,strata=~SDMVSTRA,weights=~wt,nest=TRUE,data=d)
P <- subset(mk(H), primary==1)
tid <- function(m,pat){s<-summary(m)$coef;k<-grep(pat,rownames(s));ci<-confint(m)[k,,drop=FALSE]
  data.frame(term=gsub(pat,"",rownames(s)[k]),HR=exp(s[k,"coef"]),lo=exp(ci[,1]),hi=exp(ci[,2]),p=s[k,ncol(s)])}
m3 <- svycoxph(as.formula(paste0("Surv(time,event)~lca+",C3)), design=P)
cat("Model 3 check: ")
print(transform(tid(m3,"^lca"),HR95=sprintf("%.2f (%.2f-%.2f)",HR,lo,hi))[,c("term","HR95","p")],row.names=FALSE,digits=3)
cat("\nprimary n =", sum(H$primary), " deaths =", sum(H$event[H$primary==1]),
    "| rows used by Model 3 =", m3$n, " events used =", m3$nevent, "\n")
d0_ <- subset(H, primary==1)
cat("missing among primary: educ", sum(is.na(d0_$educ)), "married", sum(is.na(d0_$married)),
    "smoke", sum(is.na(d0_$smoke)), "htn", sum(is.na(d0_$htn)), "dm", sum(is.na(d0_$dm)), "\n")

## ---- step 254 ----------------------------------------------------------------------------
d0 <- subset(H, primary==1)
cat("LBXSAL_m identical to egfr_m:", identical(d0$LBXSAL_m, d0$egfr_m),
    " (missing", sum(d0$LBXSAL_m), "/", sum(d0$egfr_m), ")\n")
cat("overlap with LBXHGB_m:", sum(d0$LBXHGB_m==1 & d0$LBXSAL_m==1), "\n\n")
B_com  <- "comorb_n_i+LBXHGB_i+LBXHGB_m+LBXSAL_i+egfr_i+bio_m+func_lim_i+n_rx_i+n_rx_m"
B_full <- paste("hf_i+stroke_i+emphysema_i+bronchitis_i+liver_i+arthritis_i+ckd_i",
                "LBXHGB_i+LBXHGB_m+LBXSAL_i+egfr_i+bio_m+func_lim_i+n_rx_i+n_rx_m",sep="+")
H$bio_m <- H$LBXSAL_m
P <- subset(mk(H), primary==1)
specs <- list(
 "Model 3 (reference)"                              = C3,
 "Model 3 + antidepressant use"                     = paste(C3,"antidep_i",sep="+"),
 "Model 4: + comorbidity and frailty block"         = paste(C3,B_com,sep="+"),
 "Model 4 + antidepressant use"                     = paste(C3,B_com,"antidep_i",sep="+"),
 "Model 4b: individual conditions + antidepressant" = paste(C3,B_full,"antidep_i",sep="+"))
out <- list()
for(nm in names(specs)){
  m <- try(svycoxph(as.formula(paste0("Surv(time,event)~lca+",specs[[nm]])), design=P), silent=TRUE)
  if(inherits(m,"try-error")){ cat(nm,": FAILED -", attr(m,"condition")$message,"\n"); next }
  jp <- tryCatch(regTermTest(m,~lca)$p, error=function(e) NA)
  out[[nm]] <- cbind(model=nm, tid(m,"^lca"), joint_p=as.numeric(jp), npar=length(coef(m)))
}
M <- do.call(rbind,out)
print(transform(M, HR95=sprintf("%.2f (%.2f-%.2f)",HR,lo,hi))[,c("model","term","HR95","p","joint_p","npar")],
      row.names=FALSE, digits=3)
write.csv(M, file.path(OUT, "RR_model4.csv"), row.names=FALSE)

## ---- step 255 ----------------------------------------------------------------------------
set.seed(20260907)
fits <- readRDS(file.path(DER, "lca_fits_cancer.rds")); post <- fits[[4]]$posterior
LAB <- c("High symptom burden","Low symptom burden","Insomnia-fatigue","Somatic-depressive")
aid <- readRDS(file.path(DER, "analysis_frame_primary.rds")); aid <- aid[aid$inAnalysis==1,]
stopifnot(nrow(post)==nrow(aid))
C4 <- paste(C3,B_com,"antidep_i",sep="+")
poolM <- function(M, rhs, primary_only=TRUE){
  est <- lapply(1:M, function(i){
    cl <- apply(post,1,function(pp) sample.int(length(pp),1,prob=pp))
    dd <- data.frame(SEQN=aid$SEQN, lca_d=factor(LAB[cl],
        levels=c("Low symptom burden","Insomnia-fatigue","Somatic-depressive","High symptom burden")))
    Hh <- merge(H[,setdiff(names(H),"lca_d")], dd, by="SEQN", all.x=TRUE)
    Hh$.k <- as.numeric(if(primary_only) Hh$primary==1 else Hh$inAnalysis==1)
    mm <- svycoxph(as.formula(paste0("Surv(time,event)~lca_d+",rhs)), design=subset(mk(Hh), .k==1))
    k <- grep("^lca_d", names(coef(mm)))
    list(b=coef(mm)[k], v=diag(vcov(mm))[k], nm=gsub("^lca_d","",names(coef(mm))[k]))})
  B<-do.call(rbind,lapply(est,`[[`,"b")); V<-do.call(rbind,lapply(est,`[[`,"v"))
  qb<-colMeans(B); ub<-colMeans(V); bv<-apply(B,2,var); tot<-ub+(1+1/M)*bv
  data.frame(term=est[[1]]$nm, HR=exp(qb), lo=exp(qb-1.96*sqrt(tot)), hi=exp(qb+1.96*sqrt(tot)),
             p=2*pnorm(-abs(qb/sqrt(tot))), fmi=round((1+1/M)*bv/tot,3))}
pc4 <- poolM(100, C4)
cat("=== Model 4 + antidepressant + 100 pseudo-class draws ===\n")
print(transform(pc4, HR95=sprintf("%.2f (%.2f-%.2f)",HR,lo,hi))[,c("term","HR95","p","fmi")],
      row.names=FALSE, digits=3)
write.csv(pc4, file.path(OUT, "RR_model4_pseudoclass.csv"), row.names=FALSE)
## direct contrast under Model 4
H2 <- H; H2$lca_if <- relevel(H2$lca, ref="Insomnia-fatigue")
mif4 <- svycoxph(as.formula(paste0("Surv(time,event)~lca_if+",C4)), design=subset(mk(H2), primary==1))
cat("\n=== Model 4: somatic-depressive vs insomnia-fatigue ===\n")
print(transform(tid(mif4,"^lca_if"),HR95=sprintf("%.2f (%.2f-%.2f)",HR,lo,hi))[,c("term","HR95","p")],
      row.names=FALSE,digits=3)
write.csv(tid(mif4,"^lca_if"), file.path(OUT, "RR_model4_contrast.csv"), row.names=FALSE)
## cause-specific (printed only; typed into step 256)
cs4 <- list()
for(oc in c("ev_ca","ev_cvd")){
  m <- svycoxph(as.formula(sprintf("Surv(time,%s)~lca+%s",oc,C4)), design=P)
  cs4[[oc]] <- tid(m,"^lca")
  cat("\n",oc,": "); print(transform(cs4[[oc]],HR95=sprintf("%.2f (%.2f-%.2f)",HR,lo,hi))[,c("term","HR95","p")],
      row.names=FALSE,digits=3)}

## ---- step 256 (Python/statsmodels -> R) ----------------------------------------------------
## main family: Model 4 + antidepressant, 3 phenotypes x 3 outcomes -- numbers as typed in the cell
rows <- data.frame(
  outcome=rep(c("All-cause","Cancer","CVD/stroke"), each=3),
  term=rep(c("Insomnia-fatigue","Somatic-depressive","High symptom burden"), 3),
  HR=c(0.85,1.48,1.18, 0.69,0.95,0.85, 0.73,1.76,1.37),
  lo=c(0.62,1.13,0.80, 0.39,0.52,0.42, 0.40,0.95,0.49),
  hi=c(1.16,1.93,1.74, 1.22,1.75,1.71, 1.35,3.29,3.87),
  p =c(0.313,0.00460,0.407, 0.206,0.872,0.644, 0.3173,0.0731,0.5502))
## statsmodels multipletests(method="fdr_bh") divides p by rank/n (R's p.adjust multiplies by
## n/rank); both are the Benjamini-Hochberg adjustment and differ at most in the last binary
## digit. The statsmodels arithmetic is used so that the file is reproduced byte for byte.
bh_statsmodels <- function(p){ n <- length(p); o <- order(p); raw <- p[o]/(seq_len(n)/n)
  adj <- pmin(rev(cummin(rev(raw))), 1); out <- numeric(n); out[o] <- adj; out }
Fm <- rows
Fm$p_bh <- bh_statsmodels(Fm$p); Fm$p_bonf <- p.adjust(Fm$p, "bonferroni")
stopifnot(isTRUE(all.equal(Fm$p_bh, p.adjust(Fm$p, "BH"), tolerance = 1e-12)))
write_py_csv(Fm, file.path(OUT, "RR_model4_multiplicity.csv"))
cat("\n=== Model 4 main family (9 tests) ===\n"); print(Fm, row.names=FALSE)
rr <- function(hr) (1-0.5^sqrt(hr))/(1-0.5^sqrt(1/hr))
ev <- function(r){ x <- max(r,1/r); x+sqrt(x*(x-1)) }
Ein <- data.frame(estimate=c("Model 3 (previous primary)","Model 4 + antidepressant (new primary)",
                             "Model 4b individual conditions","Model 4 + pseudo-class (most conservative)",
                             "Direct contrast under Model 4"),
                  HR=c(1.80,1.48,1.52,1.25,1.73), CI_lower=c(1.39,1.13,1.14,0.90,1.15))
E <- transform(Ein, approx_RR=py_round(sapply(HR, rr),2),
               E_value_point=py_round(sapply(sapply(HR, rr), ev),2),
               E_value_CI=py_round(sapply(sapply(CI_lower, rr), ev),2))
write_py_csv(E, file.path(OUT, "RR_model4_evalues.csv"))
cat("\n=== E-values ===\n"); print(E, row.names=FALSE)

## ---- step 257 ----------------------------------------------------------------------------
poolC <- function(M, rhs){
  est <- lapply(1:M, function(i){
    cl <- apply(post,1,function(pp) sample.int(length(pp),1,prob=pp))
    dd <- data.frame(SEQN=aid$SEQN, lca_d=relevel(factor(LAB[cl],
        levels=c("Low symptom burden","Insomnia-fatigue","Somatic-depressive","High symptom burden")),
        ref="Insomnia-fatigue"))
    Hh <- merge(H[,setdiff(names(H),"lca_d")], dd, by="SEQN", all.x=TRUE)
    mm <- svycoxph(as.formula(paste0("Surv(time,event)~lca_d+",rhs)), design=subset(mk(Hh), primary==1))
    k <- grep("^lca_d", names(coef(mm)))
    list(b=coef(mm)[k], v=diag(vcov(mm))[k], nm=gsub("^lca_d","",names(coef(mm))[k]))})
  B<-do.call(rbind,lapply(est,`[[`,"b")); V<-do.call(rbind,lapply(est,`[[`,"v"))
  qb<-colMeans(B); ub<-colMeans(V); bv<-apply(B,2,var); tot<-ub+(1+1/M)*bv
  data.frame(term=est[[1]]$nm, HR=exp(qb), lo=exp(qb-1.96*sqrt(tot)), hi=exp(qb+1.96*sqrt(tot)),
             p=2*pnorm(-abs(qb/sqrt(tot))), fmi=round((1+1/M)*bv/tot,3))}
set.seed(20260907)
cn3 <- poolC(100, C3); cn4 <- poolC(100, C4)
cat("\n=== contrasts (reference = insomnia-fatigue), 100 pseudo-class draws ===\n-- Model 3 --\n")
print(transform(cn3,HR95=sprintf("%.2f (%.2f-%.2f)",HR,lo,hi))[,c("term","HR95","p","fmi")],row.names=FALSE,digits=3)
cat("-- Model 4 + antidepressant --\n")
print(transform(cn4,HR95=sprintf("%.2f (%.2f-%.2f)",HR,lo,hi))[,c("term","HR95","p","fmi")],row.names=FALSE,digits=3)
write.csv(rbind(cbind(cov="Model 3",cn3),cbind(cov="Model 4 + antidepressant",cn4)),
          file.path(OUT, "RR_contrast_pseudoclass.csv"), row.names=FALSE)
cmp <- rbind(
 data.frame(comparison="vs low symptom burden", model="Model 3", HR=1.80, lo=1.39, hi=2.34),
 data.frame(comparison="vs low symptom burden", model="Model 4 + antidep", HR=1.48, lo=1.13, hi=1.93),
 data.frame(comparison="vs insomnia-fatigue",  model="Model 3", HR=1.75, lo=1.19, hi=2.59),
 data.frame(comparison="vs insomnia-fatigue",  model="Model 4 + antidep", HR=1.73, lo=1.15, hi=2.60))
cmp$attenuation_pct <- NA
cmp$attenuation_pct[2] <- round((1.80-1.48)/(1.80-1)*100,1)
cmp$attenuation_pct[4] <- round((1.75-1.73)/(1.75-1)*100,1)
print(cmp, row.names=FALSE)
write.csv(cmp, file.path(OUT, "RR_attenuation_comparison.csv"), row.names=FALSE)

## ---- step 278, second half (Python -> R): sample size for a confirmatory study ---------------
lo <- 0.90; hi <- 1.73; hr <- 1.25
b <- log(hr); se <- (log(hi)-log(lo))/(2*1.96)
se_need <- abs(b)/(1.959964+0.8416212); ratio <- (se/se_need)^2
b2 <- log(1.48); lo2 <- 1.13; hi2 <- 1.93; se2 <- (log(hi2)-log(lo2))/(2*1.96)
cat(sprintf("\nassignment-corrected HR 1.25: logHR %.4f SE %.4f z %.2f; SE needed %.4f -> %.1f x sample (%s survivors / %s deaths)\n",
            b, se, b/se, se_need, ratio, format(2569*ratio, big.mark=","), format(676*ratio, big.mark=",")))
PW <- data.frame(target=c("assignment-corrected effect (HR 1.25)","modal-assignment effect (HR 1.48)"),
                 logHR=c(py_round(b,4), py_round(b2,4)), SE=c(py_round(se,4), py_round(se2,4)),
                 z=c(py_round(b/se,2), py_round(b2/se2,2)),
                 SE_for_80pct_power=c(py_round(se_need,4), py_round(abs(b2)/2.8016,4)),
                 sample_multiplier=c(py_round(ratio,1), 1.0),
                 participants_needed=c(as.integer(trunc(2569*ratio)), 2569L),
                 deaths_needed=c(as.integer(trunc(676*ratio)), 676L))
write_py_csv(PW, file.path(OUT, "RR_power_confirmatory.csv"))
print(PW, row.names=FALSE)

## ---- transcription check: every typed-in number against the value computed above ------------
## Model 3 contrast vs insomnia-fatigue (typed 1.75, 1.19-2.59) came from the design of
## analysis_reviewer_response.R (RR_contrast_vs_insomnia_fatigue.csv); recomputed here the same way.
G <- readRDS(file.path(DER, "analysis_frame_primary.rds")); G$cvd <- as.numeric(G$cvd)
levels(G$lca)[levels(G$lca)=="Hypersomnia-somatic"] <- "Somatic-depressive"
a2 <- subset(G, primary==1); a2$lca_if <- relevel(a2$lca, ref="Insomnia-fatigue")
P2 <- subset(mk(rbind(a2, transform(subset(G,primary==0), lca_if=NA))), primary==1)
ct3 <- tid(svycoxph(as.formula(paste0("Surv(time,event)~lca_if+",C3)), design=P2),"^lca_if")
Mrow <- function(model, term) M[M$model==model & M$term==term,]
sd <- "Somatic-depressive"
## digits: HR and confidence limits were typed with 2 decimals; P values with the decimals shown
ndec <- function(x) { s <- as.character(x); ifelse(grepl("\\.", s), nchar(sub("^[^.]*\\.", "", s)), 0) }
chk <- function(file, what, typed, computed, digits = 2) data.frame(file=file, quantity=what, typed=typed,
  computed=computed, digits=digits, agrees = abs(round(computed, digits) - typed) < 1e-9)
C <- rbind(
  do.call(rbind, lapply(1:3, function(i){ r <- Mrow("Model 4 + antidepressant use", rows$term[i])
    rbind(chk("RR_model4_multiplicity", paste("All-cause", rows$term[i], c("HR","lo","hi","p")),
              unlist(rows[i,c("HR","lo","hi","p")]), unlist(r[,c("HR","lo","hi","p")]),
              digits = c(2,2,2,ndec(rows$p[i])))) })),
  do.call(rbind, lapply(4:9, function(i){ r <- cs4[[ifelse(i<=6,"ev_ca","ev_cvd")]]; r <- r[r$term==rows$term[i],]
    chk("RR_model4_multiplicity", paste(rows$outcome[i], rows$term[i], c("HR","lo","hi","p")),
        unlist(rows[i,c("HR","lo","hi","p")]), unlist(r[,c("HR","lo","hi","p")]),
        digits = c(2,2,2,ndec(rows$p[i]))) })),
  chk("RR_model4_evalues / attenuation / sens_evalues", paste("Model 3, somatic-depressive", c("HR","lo","hi")),
      c(1.80,1.39,2.34), unlist(Mrow("Model 3 (reference)", sd)[,c("HR","lo","hi")])),
  chk("RR_model4_evalues / attenuation / power", paste("Model 4 + antidepressant, somatic-depressive", c("HR","lo","hi")),
      c(1.48,1.13,1.93), unlist(Mrow("Model 4 + antidepressant use", sd)[,c("HR","lo","hi")])),
  chk("RR_model4_evalues", paste("Model 4b, somatic-depressive", c("HR","lo")), c(1.52,1.14),
      unlist(Mrow("Model 4b: individual conditions + antidepressant", sd)[,c("HR","lo")])),
  chk("RR_model4_evalues / power", paste("Model 4 pseudo-class, somatic-depressive", c("HR","lo","hi")),
      c(1.25,0.90,1.73), unlist(pc4[pc4$term==sd, c("HR","lo","hi")])),
  chk("RR_model4_evalues / attenuation", paste("Model 4 contrast, somatic-depressive vs insomnia-fatigue", c("HR","lo","hi")),
      c(1.73,1.15,2.60), unlist(tid(mif4,"^lca_if")[2, c("HR","lo","hi")])),
  chk("RR_attenuation_comparison", paste("Model 3 contrast, somatic-depressive vs insomnia-fatigue", c("HR","lo","hi")),
      c(1.75,1.19,2.59), unlist(ct3[ct3$term==sd, c("HR","lo","hi")])))
rownames(C) <- NULL
cat("\n=== transcription check (typed value vs recomputed value rounded to the typed digits) ===\n")
print(C, row.names=FALSE, digits=6)
cat("all typed numbers agree:", all(C$agrees), "\n")
write.csv(C, file.path(OUT, "input_checks", "check_model4_transcription.csv"), row.names=FALSE)
cat("model4.R done\n")
