## =============================================================================
## proportional_hazards.R  --  proportional-hazards checks for the phenotype (Model 3)
##
## Regenerates (into review/reproduction/output/):
##   ph_period_specific.csv     phenotype HRs within 0-5 yr and >5 yr of follow-up
##   ph_stratified_baseline.csv baseline hazard stratified on age group x years since diagnosis
##   ph_tests.csv               Schoenfeld tests (unweighted) + design-based period interaction
##
## Reconstructed from the original analysis log:
##   step 162: G = analysis_frame_primary.rds, C2/C3
##             (its cox.zph on the WEIGHTED model was degenerate
##             and not written; skipped)
##   step 163: unweighted coxph + cox.zph; survSplit at 5 years
##   step 164: options(survey.lonely.psu = "adjust"); rebuilds sp
##             (its lca*period model was singular -> run stopped)
##   steps 165-168: failed attempts / diagnostics (not repeated);
##                  step 167 created sp$late, sp$L2, sp$L4 used below
##   step 169: phenotype-by-period terms L2/L3/L4, joint Wald test,
##             period-specific HRs -> ph_period_specific.csv
##   step 170: strata(agec, ydxc) model -> ph_stratified_baseline.csv;
##             ph_tests.csv written from numbers typed in by hand
##
## INCONSISTENCY KEPT ON PURPOSE (reproduced, not fixed):
##   Unlike every other analysis, the two design-based models here build the survey design on
##   the 2,569 primary-cohort rows only (svydesign(data = sp) / svydesign(data = a)) instead of on
##   all 70,190 participants followed by subset(). This is why survey.lonely.psu = "adjust" was
##   needed (stratum 89 has a single PSU inside the cohort). Variances therefore differ from the
##   full-design convention.
## ph_tests.csv: the original cell typed the statistics in (2.1781, 73.8367, 0.5363, 0.0001,
##   0.9683, df 3/33/3). They are kept as typed; the script recomputes each one and prints a
##   transcription check.
## Run   : /opt/homebrew/bin/Rscript code/reconstructed/proportional_hazards.R  (from the project root)
## =============================================================================
suppressPackageStartupMessages({library(survey); library(survival)})
if (!file.exists("data/derived/analysis_frame_primary.rds")) stop("run from the project root")
OUT <- "review/reproduction/output"; dir.create(file.path(OUT, "input_checks"), recursive = TRUE, showWarnings = FALSE)

## ---- step 162 (objects used later) --------------------------------------------------
G <- readRDS("data/derived/analysis_frame_primary.rds"); G$cvd <- as.numeric(G$cvd)
C2 <- "age+sex+race4+educ+married+pir_i+pir_m"
C3 <- paste(C2,"bmi_i+bmi_m+smoke+pa3+htn+dm+cvd+ydx_i+ydx_m+multi_primary+catype",sep="+")

## ---- step 163 ---------------------------------------------------------------------
a <- subset(G, primary==1)
## (1) Schoenfeld test on the UNWEIGHTED Cox model -- diagnostic only (ignores the design)
zu <- cox.zph(coxph(as.formula(paste0("Surv(time,event)~lca+",C3)), data=a))
tu <- as.data.frame(zu$table)
cat("=== unweighted cox.zph (diagnostic) ===\n")
print(round(tu[c(grep("^lca$",rownames(tu),value=TRUE),"GLOBAL"),,drop=FALSE],4))
cat("p<0.10:", paste(rownames(tu)[tu[,"p"]<0.10 & rownames(tu)!="GLOBAL"],collapse=", "), "\n")

## ---- step 164 (first lines; the rest of the cell failed) -----------------------------
options(survey.lonely.psu="adjust")
sp <- survSplit(Surv(time,event)~., data=a, cut=5, episode="period")
sp$period <- factor(sp$period, labels=c("0-5 yr",">5 yr"))

## ---- step 167 (the columns it created before failing) --------------------------------
sp$late <- as.numeric(sp$period==">5 yr")
sp$L2 <- as.numeric(sp$lca=="Insomnia-fatigue")   * sp$late
sp$L3 <- as.numeric(sp$lca=="Somatic-depressive") * sp$late   # all zero here: level not yet renamed
sp$L4 <- as.numeric(sp$lca=="High symptom burden")* sp$late

## ---- step 169 ----------------------------------------------------------------------
levels(sp$lca)[levels(sp$lca)=="Hypersomnia-somatic"] <- "Somatic-depressive"
sp$L3 <- as.numeric(sp$lca=="Somatic-depressive") * sp$late
d2 <- svydesign(ids=~SDMVPSU, strata=~SDMVSTRA, weights=~wt, nest=TRUE, data=sp)
mi <- svycoxph(as.formula(paste0("Surv(tstart,time,event)~lca+L2+L3+L4+",C3)), design=d2)
rt <- regTermTest(mi, ~L2+L3+L4)
cat("\n=== design-based: does the phenotype effect change with follow-up period? ===\njoint Wald p =",
    signif(rt$p,4), " (df =", rt$df, ")\n\n")
ss <- summary(mi)$coef
nm <- c(L2="Insomnia-fatigue", L3="Somatic-depressive", L4="High symptom burden")
for(v in names(nm)){ ci <- confint(mi)[v,]
  cat(sprintf("  %-20s >5yr vs 0-5yr ratio = %.2f (%.2f-%.2f), p = %.3f\n",
      nm[v], exp(ss[v,"coef"]), exp(ci[1]), exp(ci[2]), ss[v,ncol(ss)])) }
per <- do.call(rbind, lapply(levels(sp$period), function(pp){
  mm <- svycoxph(as.formula(paste0("Surv(tstart,time,event)~lca+",C3)), design=subset(d2, period==pp))
  s2 <- summary(mm)$coef; k <- grep("^lca", rownames(s2)); ci <- confint(mm)[k,,drop=FALSE]
  data.frame(period=pp, events=sum(sp$event[sp$period==pp]), term=gsub("^lca","",rownames(s2)[k]),
             HR=exp(s2[k,"coef"]), lo=exp(ci[,1]), hi=exp(ci[,2]), p=s2[k,ncol(s2)])}))
print(transform(per, HR95=sprintf("%.2f (%.2f-%.2f)",HR,lo,hi))[,c("period","events","term","HR95","p")],
      row.names=FALSE, digits=3)
write.csv(per, file.path(OUT, "ph_period_specific.csv"), row.names=FALSE)
## (step 169 also wrote ph_interaction_test.csv, which is not a target of this reconstruction)

## ---- step 170 ----------------------------------------------------------------------
a$lca <- factor(ifelse(as.character(a$lca)=="Hypersomnia-somatic","Somatic-depressive",as.character(a$lca)),
                levels=c("Low symptom burden","Insomnia-fatigue","Somatic-depressive","High symptom burden"))
a$agec <- cut(a$age,  c(0,55,65,75,120), labels=c("<55","55-64","65-74",">=75"))
a$ydxc <- cut(a$ydx_i, c(-1,2,5,10,100),  labels=c("0-2","3-5","6-10",">10"))
Ps <- svydesign(ids=~SDMVPSU, strata=~SDMVSTRA, weights=~wt, nest=TRUE, data=a)
C3s <- paste("sex+race4+educ+married+pir_i+pir_m+bmi_i+bmi_m+smoke+pa3+htn+dm+cvd+ydx_m+multi_primary+catype")
ms <- svycoxph(as.formula(paste0("Surv(time,event)~lca+strata(agec,ydxc)+",C3s)), design=Ps)
ss <- summary(ms)$coef; k <- grep("^lca",rownames(ss)); ci <- confint(ms)[k,,drop=FALSE]
strat <- data.frame(term=gsub("^lca","",rownames(ss)[k]), HR=exp(ss[k,"coef"]),
                    lo=exp(ci[,1]), hi=exp(ci[,2]), p=ss[k,ncol(ss)])
cat("\n=== baseline hazard stratified on age group x years since diagnosis (16 strata) ===\n")
print(transform(strat, HR95=sprintf("%.2f (%.2f-%.2f)",HR,lo,hi))[,c("term","HR95","p")],
      row.names=FALSE, digits=3)
write.csv(strat, file.path(OUT, "ph_stratified_baseline.csv"), row.names=FALSE)
## typed-in values, exactly as in the original cell
write.csv(data.frame(
  test=c("lca Schoenfeld residuals (unweighted diagnostic)","GLOBAL Schoenfeld (unweighted diagnostic)",
         "lca x follow-up period interaction (design-based joint Wald)"),
  statistic=c(2.1781,73.8367,NA), df=c(3,33,3), p=c(0.5363,0.0001,0.9683),
  note=c("exposure satisfies PH","driven by age and years since diagnosis","no evidence of time-varying phenotype effect")),
  file.path(OUT, "ph_tests.csv"), row.names=FALSE)

## ---- transcription check for ph_tests.csv -------------------------------------------
recomp <- data.frame(
  quantity = c("lca Schoenfeld chisq","lca Schoenfeld df","lca Schoenfeld p",
               "GLOBAL chisq","GLOBAL df","GLOBAL p","joint Wald p","joint Wald df"),
  typed    = c(2.1781, 3, 0.5363, 73.8367, 33, 0.0001, 0.9683, 3),
  recomputed = c(round(tu["lca","chisq"],4), tu["lca","df"], round(tu["lca","p"],4),
                 round(tu["GLOBAL","chisq"],4), tu["GLOBAL","df"], round(tu["GLOBAL","p"],4),
                 signif(rt$p,4), rt$df))
recomp$agrees <- abs(recomp$typed - recomp$recomputed) < 1e-9
cat("\n=== ph_tests.csv transcription check (typed value vs recomputed, same rounding) ===\n")
print(recomp, row.names = FALSE)
write.csv(recomp, file.path(OUT, "input_checks", "check_ph_tests_transcription.csv"), row.names = FALSE)
cat("proportional_hazards.R done\n")
