## =============================================================================
## measurement_sensitivity_RR4.R  --  class-assignment and measurement sensitivity ("RR4" round)
##
## Regenerates (into review/reproduction/output/):
##   RR4_inclusive_pseudoclass.csv  outcome-informed LCA (death + log follow-up in the membership
##                                  model), 50 draws each under Model 3 and Model 4
##   RR4_inclusive_profiles.csv     item profiles of that outcome-informed solution
##   RR4_threshold_cox.csv          items dichotomised at >= 1 ("several days"), Model 4
##   RR4_localdep_cox.csv           PHQ-9 sleep item removed, Model 4
##   RR4_measurement_fit.csv        entropy and Cramer's V of the two refitted solutions
##
## Reconstructed from the original analysis log:
##   step 406: H = analysis_frame_primary.rds (set.seed(20260908) here
##             is superseded: no random numbers are drawn before step 410 re-seeds)
##   step 407: G4 = nhanes_model4_design_frame.csv
##   step 408: H2 = H + Model 4 covariates of G4; phenotype labels renamed
##             (Sleep-fatigue, Somatic-depressive); its Model 4 with separate LBXSAL_m and egfr_m was
##             singular and the cell stopped there (not repeated)
##   step 409: bio_m; Model 4 formula C4 (43 terms) -> 1.48 (1.13-1.93)
##   step 410: set.seed(20260908); outcome-informed poLCA (10 starts); 2 x 50 draws
##   step 411: set.seed(20260908); threshold >= 1 poLCA, then (same random stream)
##             poLCA without the PHQ-9 sleep item
## These cells are the ones that produced the RR4_* files. code/analysis_class_assignment.R is a
## later tidy-up of the same analyses that sets the seed only once at the top, so it does NOT
## reproduce these exact files (it writes RR3_* names); this script keeps the per-cell seeds.
## Inputs: data/derived/analysis_frame_primary.rds, data/derived/nhanes_model4_design_frame.csv.gz,
##             supporting/nhanes_dpq_raw.csv (raw 0-3 PHQ-9 item scores, read as in the original)
## NOTE: this Model 4 frame starts from analysis_frame_primary.rds (educ/married/smoke already
##   modally imputed), so 2,564 survivors enter the Cox models, whereas model4.R (steps 252-257)
##   uses the CSV frame without that imputation and analyses 2,557. Both give 1.48 (1.13-1.93)
##   for the somatic-depressive Model 4 estimate.
## Run   : /opt/homebrew/bin/Rscript code/reconstructed/measurement_sensitivity_RR4.R (from the project root)
## =============================================================================
suppressPackageStartupMessages({library(poLCA); library(survey); library(survival)})
if (!file.exists("data/derived/analysis_frame_primary.rds")) stop("run from the project root")
OUT <- "review/reproduction/output"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
DER <- "data/derived"

## ---- step 406 ---------------------------------------------------------------------------------
options(survey.lonely.psu="adjust"); set.seed(20260908)
H  <- readRDS(file.path(DER, "analysis_frame_primary.rds"))
## ---- step 407 ---------------------------------------------------------------------------------
G4 <- read.csv(file.path(DER, "nhanes_model4_design_frame.csv.gz"))
## ---- step 408 (up to the failing model) ------------------------------------------------------
keep <- c("SEQN","comorb_n_i","LBXHGB_i","LBXHGB_m","LBXSAL_i","LBXSAL_m",
          "egfr_i","egfr_m","func_lim_i","n_rx_i","n_rx_m","antidep_i")
H2 <- merge(H, G4[,keep], by="SEQN", all.x=TRUE)
H2$cvd <- as.numeric(H2$cvd)
levels(H2$lca)[levels(H2$lca)=="Hypersomnia-somatic"] <- "Somatic-depressive"
levels(H2$lca)[levels(H2$lca)=="Insomnia-fatigue"]    <- "Sleep-fatigue"
cat("merged primary:", sum(H2$primary,na.rm=TRUE), "| deaths:", sum(H2$event[H2$primary==1],na.rm=TRUE), "\n")
mk <- function(d) svydesign(ids=~SDMVPSU,strata=~SDMVSTRA,weights=~wt,nest=TRUE,data=d)
tid <- function(m,pat){s<-summary(m)$coef;k<-grep(pat,rownames(s));ci<-confint(m)[k,,drop=FALSE]
  data.frame(term=gsub(pat,"",rownames(s)[k]),HR=exp(s[k,"coef"]),lo=exp(ci[,1]),hi=exp(ci[,2]),p=s[k,ncol(s)])}
C2 <- "age+sex+race4+educ+married+pir_i+pir_m"
C3 <- paste(C2,"bmi_i+bmi_m+smoke+pa3+htn+dm+cvd+ydx_i+ydx_m+multi_primary+catype",sep="+")

## ---- step 409 ---------------------------------------------------------------------------------
a <- subset(H2, primary==1)
cat("LBXSAL_m identical to egfr_m:", all(a$LBXSAL_m==a$egfr_m), "\n")
H2$bio_m <- pmax(H2$LBXSAL_m, H2$egfr_m)
B_com2 <- "comorb_n_i+LBXHGB_i+LBXHGB_m+LBXSAL_i+egfr_i+bio_m+func_lim_i+n_rx_i+n_rx_m"
C4 <- paste(C3, B_com2, "antidep_i", sep="+")
P  <- subset(mk(H2), primary==1)
m4 <- svycoxph(as.formula(paste0("Surv(time,event)~lca+",C4)), design=P)
cat("\n=== Model 4 check (manuscript: 1.48, 1.13-1.93) ===\n")
print(transform(tid(m4,"^lca"),HR95=sprintf("%.2f (%.2f-%.2f)",HR,lo,hi))[,c("term","HR95","p")],
      row.names=FALSE, digits=3)
cat("Model 4 terms:", length(coef(m4)), "| rows used:", m4$n, "| events used:", m4$nevent, "\n")
## (the original also saved analysis_frame_model4.rds and model4_formula.txt; not needed here)

## ---- step 410 ---------------------------------------------------------------------------------
RW <- read.csv("supporting/nhanes_dpq_raw.csv")
D  <- merge(subset(H2, inAnalysis==1), RW, by="SEQN", all.x=TRUE)
D$logt <- log(pmax(D$time,0.05))
ent <- function(p) 1-sum(-p*log(pmax(p,1e-12)))/(nrow(p)*log(ncol(p)))
draw_pool <- function(post, rhs, ids, M=50){
  est <- lapply(seq_len(M), function(i){
    cl <- apply(post,1,function(pp) sample.int(length(pp),1,prob=pp))
    dd <- data.frame(SEQN=ids, lz=factor(paste0("C",cl)))
    dd$lz <- relevel(dd$lz, ref=names(which.max(table(dd$lz))))
    Hh <- merge(H2[,setdiff(names(H2),"lz")], dd, by="SEQN", all.x=TRUE)
    mm <- svycoxph(as.formula(paste0("Surv(time,event)~lz+",rhs)),
                   design=subset(mk(Hh), primary==1 & !is.na(lz)))
    k <- grep("^lz", names(coef(mm)))
    list(b=coef(mm)[k], v=diag(vcov(mm))[k], nm=gsub("^lz","",names(coef(mm))[k]))})
  B <- do.call(rbind,lapply(est,`[[`,"b")); V <- do.call(rbind,lapply(est,`[[`,"v"))
  qb <- colMeans(B); ub <- colMeans(V); bv <- apply(B,2,var); tot <- ub+(1+1/M)*bv
  data.frame(term=est[[1]]$nm, HR=exp(qb), lo=exp(qb-1.96*sqrt(tot)), hi=exp(qb+1.96*sqrt(tot)),
             p=2*pnorm(-abs(qb/sqrt(tot))), fmi=round((1+1/M)*bv/tot,3))}
set.seed(20260908)
fI <- poLCA(cbind(phqi1,phqi2,phqi3,phqi4,phqi5,phqi6,phqi7,phqi8,phqi9,sleepcat,slq050b)~event+logt,
            D, nclass=4, maxiter=8000, nrep=10, verbose=FALSE, na.rm=TRUE)
cat("A. outcome-informed | entropy", round(ent(fI$posterior),3),
    "| class %", paste(round(fI$P*100,1),collapse="/"), "\n")
incl <- rbind(cbind(model="Model 3", draw_pool(fI$posterior, C3, D$SEQN)),
              cbind(model="Model 4", draw_pool(fI$posterior, C4, D$SEQN)))
print(transform(incl,HR95=sprintf("%.2f (%.2f-%.2f)",HR,lo,hi))[,c("model","term","HR95","p","fmi")],
      row.names=FALSE, digits=3)
write.csv(incl, file.path(OUT, "RR4_inclusive_pseudoclass.csv"), row.names=FALSE)
write.csv(round(rbind(do.call(rbind,lapply(1:9,function(j) fI$probs[[j]][,2]*100)),
  short=fI$probs[[10]][,1]*100, long=fI$probs[[10]][,3]*100, told=fI$probs[[11]][,1]*100),1),
  file.path(OUT, "RR4_inclusive_profiles.csv"))

## ---- step 411 ---------------------------------------------------------------------------------
set.seed(20260908)
## B. threshold >= 1
D2 <- D
for(i in 1:9) D2[[paste0("phqi",i)]] <- as.integer(D[[paste0("dpq",i,"_raw")]] >= 1) + 1L
f1 <- poLCA(cbind(phqi1,phqi2,phqi3,phqi4,phqi5,phqi6,phqi7,phqi8,phqi9,sleepcat,slq050b)~1,
            D2, nclass=4, maxiter=8000, nrep=10, verbose=FALSE, na.rm=TRUE)
D2$lt1 <- factor(paste0("T",apply(f1$posterior,1,which.max)))
ct <- table(original=D2$lca, threshold1=D2$lt1)
V1 <- sqrt(chisq.test(ct)$statistic/(sum(ct)*(min(dim(ct))-1)))
Ht <- merge(H2, D2[,c("SEQN","lt1")], by="SEQN", all.x=TRUE)
Ht$lt1 <- relevel(Ht$lt1, ref=names(which.max(table(D2$lt1))))
mt <- svycoxph(as.formula(paste0("Surv(time,event)~lt1+",C4)),
               design=subset(mk(Ht), primary==1 & !is.na(lt1)))
cat("=== B. threshold >= 1 (earlier: entropy 0.714, V 0.497, no class associated) ===\n")
cat("entropy", round(ent(f1$posterior),3), "| Cramer V", round(V1,3), "\n")
print(transform(tid(mt,"^lt1"),HR95=sprintf("%.2f (%.2f-%.2f)",HR,lo,hi))[,c("term","HR95","p")],
      row.names=FALSE, digits=3)
## C. PHQ-9 sleep item removed
fL <- poLCA(cbind(phqi1,phqi2,phqi4,phqi5,phqi6,phqi7,phqi8,phqi9,sleepcat,slq050b)~1,
            D, nclass=4, maxiter=8000, nrep=10, verbose=FALSE, na.rm=TRUE)
D$lns <- factor(paste0("L",apply(fL$posterior,1,which.max)))
ct2 <- table(original=D$lca, dropped=D$lns)
V2 <- sqrt(chisq.test(ct2)$statistic/(sum(ct2)*(min(dim(ct2))-1)))
Hn <- merge(H2, D[,c("SEQN","lns")], by="SEQN", all.x=TRUE)
Hn$lns <- relevel(Hn$lns, ref=names(which.max(table(D$lns))))
mL <- svycoxph(as.formula(paste0("Surv(time,event)~lns+",C4)),
               design=subset(mk(Hn), primary==1 & !is.na(lns)))
cat("\n=== C. PHQ-9 sleep item removed (earlier: V 0.736, HR 1.29, 0.98-1.69, p=0.065) ===\n")
cat("entropy", round(ent(fL$posterior),3), "| Cramer V", round(V2,3), "\n")
print(transform(tid(mL,"^lns"),HR95=sprintf("%.2f (%.2f-%.2f)",HR,lo,hi))[,c("term","HR95","p")],
      row.names=FALSE, digits=3)
cat("\ncorrespondence (original phenotype x item-removed class):\n"); print(ct2)
write.csv(tid(mt,"^lt1"), file.path(OUT, "RR4_threshold_cox.csv"), row.names=FALSE)
write.csv(tid(mL,"^lns"), file.path(OUT, "RR4_localdep_cox.csv"), row.names=FALSE)
write.csv(data.frame(analysis=c("threshold>=1","drop PHQ-9 sleep item"),
                     entropy=c(round(ent(f1$posterior),3),round(ent(fL$posterior),3)),
                     cramer_v=c(round(V1,3),round(V2,3))), file.path(OUT, "RR4_measurement_fit.csv"), row.names=FALSE)
cat("measurement_sensitivity_RR4.R done\n")
