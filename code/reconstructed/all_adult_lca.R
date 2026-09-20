## =============================================================================
## all_adult_lca.R  --  refit the measurement model in all adults (n = 34,022) and apply its
##                      4-class posteriors to the primary survivor cohort
##
## Regenerates (into review/reproduction/output/):
##   lca_fit_all_adults.csv             k = 1..7 fit indices in all adults
##   RR_fullsample_lca_profiles.csv     item profiles of the all-adult k = 4 solution
##   RR_fullsample_lca_pseudoclass.csv  50 pseudo-class draws from those posteriors, Model 4
## Also writes (review/reproduction/output/):
##   intermediate/lca_fits_all_adults.rds          the refitted models (the original file,
##                                                  lca_fits_all_adults.rds, was never copied
##                                                  into the project, so the fit must be redone)
##   input_checks/check_all_adults_lca_input.csv   re-derivation of the LCA input from raw NHANES
##
## Reconstructed from the original analysis log:
##   step 276 (Python): built nhanes_all_adults_lca_input.csv from the
##                      DEMO, DPQ and SLQ modules -> translated to R below, used ONLY as a check of the saved file
##   step 280 (R): set.seed(20260907); poLCA k = 1..7, 20 random starts,
##                 maxiter 10000 -> lca_fit_all_adults.csv
##   steps 252/254/255 (R): frame H, mk() and the Model 4 string C4 that
##                          step 283 relied on (re-created here exactly as in model4.R)
##   step 283 (R): profiles; set.seed(20260907); 50 draws, Model 4
## Runtime: the seven all-adult fits take roughly half an hour on one core.
## Run   : /opt/homebrew/bin/Rscript code/reconstructed/all_adult_lca.R   (from the project root)
## =============================================================================
suppressPackageStartupMessages({library(poLCA); library(survey); library(survival)})
if (!file.exists("data/derived/nhanes_all_adults_lca_input.csv.gz")) stop("run from the project root")
OUT <- "review/reproduction/output"
dir.create(file.path(OUT, "input_checks"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "intermediate"), recursive = TRUE, showWarnings = FALSE)
DER <- "data/derived"; RAWD <- "data/raw_nhanes"

## ---- check of the saved input: step 276 translated from pandas -----------------------------
CY <- c(`2005`="_D", `2007`="_E", `2009`="_F", `2011`="_G", `2013`="_H", `2015`="_I", `2017`="_J")
DPQ9 <- sprintf("DPQ0%d0", 1:9)
rawf <- function(mod, y) read.csv(file.path(RAWD, sprintf("%s%s_%s-%d.csv.gz", mod, CY[[y]], y, as.integer(y) + 1)))
rows <- lapply(names(CY), function(y) {
  dem <- rawf("DEMO", y); dpq <- rawf("DPQ", y); slq <- rawf("SLQ", y)
  d <- dem[, c("SEQN","RIDAGEYR","RIAGENDR","SDMVPSU","SDMVSTRA","WTMEC2YR")]
  q <- dpq[, c("SEQN", DPQ9)]; for (c in DPQ9) q[[c]][!(q[[c]] <= 3) | is.na(q[[c]])] <- NA
  h <- if ("SLD010H" %in% names(slq)) "SLD010H" else "SLD012"
  sl <- data.frame(SEQN = slq$SEQN, sleep_h = ifelse(slq[[h]] <= 24, slq[[h]], NA),
                   slq050 = ifelse(slq$SLQ050 <= 2, slq$SLQ050, NA))
  m <- merge(merge(d, q, by = "SEQN"), sl, by = "SEQN")        # pandas inner joins
  m$cycle <- sprintf("%s-%d", y, as.integer(y) + 1); m })
ALLAD <- do.call(rbind, rows)
ALLAD <- ALLAD[ALLAD$RIDAGEYR >= 20 & complete.cases(ALLAD[, DPQ9]) & !is.na(ALLAD$sleep_h) & !is.na(ALLAD$slq050), ]
ALLAD$phq9_score <- rowSums(ALLAD[, DPQ9])
for (i in 1:9) ALLAD[[paste0("phqi", i)]] <- as.integer(ALLAD[[DPQ9[i]]] >= 2) + 1L
ALLAD$sleepcat <- ifelse(ALLAD$sleep_h < 6, 1, ifelse(ALLAD$sleep_h >= 9, 3, 2))
ALLAD$slq050b <- ifelse(ALLAD$slq050 == 1, 1, 2)
ALLAD$wt <- ALLAD$WTMEC2YR / 7
surv <- subset(read.csv(file.path(DER, "nhanes_cancer_design_frame.csv.gz")), inAnalysis == 1)$SEQN
ALLAD$is_survivor <- as.integer(ALLAD$SEQN %in% surv)
SAVED <- read.csv(file.path(DER, "nhanes_all_adults_lca_input.csv.gz"))
ALLAD <- ALLAD[order(ALLAD$SEQN), ]; S2 <- SAVED[order(SAVED$SEQN), ]
chk <- data.frame(column = names(SAVED), n_saved = nrow(SAVED), n_rederived = nrow(ALLAD),
  same_ids = identical(as.numeric(S2$SEQN), as.numeric(ALLAD$SEQN)),
  max_abs_diff = sapply(names(SAVED), function(v) {
    a <- S2[[v]]; b <- ALLAD[[v]]
    if (is.null(b) || length(a) != length(b)) return(NA)
    if (is.numeric(a)) { if (!identical(is.na(a), is.na(b))) return(Inf); max(c(0, abs(a - b)), na.rm = TRUE) }
    else as.numeric(sum(a != b, na.rm = TRUE)) }))
cat("=== re-derivation of nhanes_all_adults_lca_input.csv from data/raw_nhanes ===\n")
print(chk, row.names = FALSE)
write.csv(chk, file.path(OUT, "input_checks", "check_all_adults_lca_input.csv"), row.names = FALSE)

## ---- step 280 (verbatim apart from file paths) ----------------------------------------------
A <- read.csv(file.path(DER, "nhanes_all_adults_lca_input.csv.gz"))
set.seed(20260907)
f <- cbind(phqi1,phqi2,phqi3,phqi4,phqi5,phqi6,phqi7,phqi8,phqi9,sleepcat,slq050b)~1
ent <- function(x){p<-x$posterior; if(is.null(dim(p))) return(NA)
  1-sum(-p*log(pmax(p,1e-12)))/(nrow(p)*log(ncol(p)))}
t0 <- Sys.time()
fitsA <- lapply(1:7, function(k) poLCA(f, A, nclass=k, maxiter=10000, nrep=20, verbose=FALSE, na.rm=TRUE))
cat("all-adult fits took", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
saveRDS(fitsA, file.path(OUT, "intermediate", "lca_fits_all_adults.rds"))
ic <- data.frame(k=1:7, logLik=round(sapply(fitsA,`[[`,"llik"),1),
  BIC=round(sapply(fitsA,`[[`,"bic"),1),
  cAIC=round(sapply(fitsA,function(x) -2*x$llik+x$npar*(log(x$Nobs)+1)),1),
  entropy=round(sapply(fitsA,ent),3),
  min_class_pct=sapply(fitsA,function(x) round(min(x$P)*100,1)),
  mean_maxpost=round(sapply(fitsA,function(x){p<-x$posterior
    if(is.null(dim(p))) NA else mean(apply(p,1,max))}),3))
write.csv(ic, file.path(OUT, "lca_fit_all_adults.csv"), row.names=FALSE)
cat("n =", nrow(A), "\n"); print(ic, row.names=FALSE)

## ---- steps 252 + 254 + 255: the Model 4 frame H, mk() and C4 used by step 283 -----------------
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
B_com  <- "comorb_n_i+LBXHGB_i+LBXHGB_m+LBXSAL_i+egfr_i+bio_m+func_lim_i+n_rx_i+n_rx_m"
H$bio_m <- H$LBXSAL_m
C4 <- paste(C3,B_com,"antidep_i",sep="+")

## ---- step 283 (verbatim apart from file paths) ----------------------------------------------
ftA <- fitsA[[4]]; postA <- ftA$posterior
items <- c("Anhedonia","Depressed mood","Sleep disturbance","Fatigue","Appetite change",
           "Worthlessness","Concentration","Psychomotor","Suicidal ideation")
tb <- t(sapply(1:9, function(j) round(ftA$probs[[j]][,2]*100,1))); rownames(tb) <- items
tb <- rbind(tb, "Short sleep"=round(ftA$probs[[10]][,1]*100,1),
                "Long sleep"=round(ftA$probs[[10]][,3]*100,1),
                "Insomnia complaint"=round(ftA$probs[[11]][,1]*100,1))
colnames(tb) <- sprintf("C%d (%.1f%%)",1:4,ftA$P*100)
cat("=== all-adult k=4 item-response probabilities ===\n"); print(tb)
idx <- match(H$SEQN[H$primary==1], A$SEQN)
cat("\nsurvivors matched to the all-adult fit:", sum(!is.na(idx)), "/", sum(H$primary==1), "\n")
pS <- postA[idx[!is.na(idx)],,drop=FALSE]
cat("mean max posterior in that subset:", round(mean(apply(pS,1,max)),3), " (0.886 within survivors)\n")
set.seed(20260907)
seqn_ok <- H$SEQN[H$primary==1][!is.na(idx)]
poolA <- function(M, rhs){
  est <- lapply(1:M, function(i){
    cl <- apply(pS,1,function(pp) sample.int(length(pp),1,prob=pp))
    dd <- data.frame(SEQN=seqn_ok, lcaA=factor(paste0("C",cl), levels=paste0("C",1:4)))
    Hh <- merge(H[,setdiff(names(H),"lcaA")], dd, by="SEQN", all.x=TRUE)
    Hh$lcaA <- relevel(Hh$lcaA, ref=names(which.max(table(dd$lcaA))))
    mm <- svycoxph(as.formula(paste0("Surv(time,event)~lcaA+",rhs)), design=subset(mk(Hh), primary==1 & !is.na(lcaA)))
    k <- grep("^lcaA", names(coef(mm))); list(b=coef(mm)[k], v=diag(vcov(mm))[k], nm=gsub("^lcaA","",names(coef(mm))[k]))})
  B<-do.call(rbind,lapply(est,`[[`,"b")); V<-do.call(rbind,lapply(est,`[[`,"v"))
  qb<-colMeans(B); ub<-colMeans(V); bv<-apply(B,2,var); tot<-ub+(1+1/M)*bv
  data.frame(term=est[[1]]$nm, HR=exp(qb), lo=exp(qb-1.96*sqrt(tot)), hi=exp(qb+1.96*sqrt(tot)),
             p=2*pnorm(-abs(qb/sqrt(tot))), fmi=round((1+1/M)*bv/tot,3))}
pcA <- poolA(50, C4)
cat("\n=== all-adult classes + Model 4 + pseudo-class draws ===\n")
print(transform(pcA,HR95=sprintf("%.2f (%.2f-%.2f)",HR,lo,hi))[,c("term","HR95","p","fmi")],row.names=FALSE,digits=3)
write.csv(pcA, file.path(OUT, "RR_fullsample_lca_pseudoclass.csv"), row.names=FALSE)
write.csv(as.data.frame(tb), file.path(OUT, "RR_fullsample_lca_profiles.csv"))
cat("all_adult_lca.R done\n")
