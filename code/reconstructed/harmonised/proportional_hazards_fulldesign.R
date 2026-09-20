## =============================================================================
## proportional_hazards_fulldesign.R  --  proportional_hazards.R with the survey design built on
##                                        all 70,190 participants and then subset (finding 1)
##
## Writes to review/reproduction/harmonised/:
##   ph_period_specific.csv, ph_stratified_baseline.csv, ph_tests.csv, analysed_n_proportional_hazards.csv
##
## Changes relative to code/reconstructed/proportional_hazards.R (original steps 163-170):
##   * the two design-based analyses (phenotype-by-period model and the stratified-baseline model)
##     now use svydesign() on all 70,190 participants and subset(design, analysed == 1), like every
##     other analysis, instead of svydesign() on the 2,569 cohort rows only. For the split-follow-up
##     model the cohort's split rows are stacked with the 67,621 other participants (one row each,
##     outside the domain) so that the design keeps every PSU and stratum;
##   * domain = the harmonised 2,564 (frame_2564.R). These Model 3 models already analysed exactly
##     these 2,564 people (the 5 with missing hypertension/diabetes were dropped as NA), so only the
##     design construction changes;
##   * ph_tests.csv is computed instead of typed (same rounding as the typed values).
## The unweighted Schoenfeld tests do not involve the design and are unchanged.
## The `events` column of ph_period_specific.csv keeps its original meaning (deaths in that period
## among all 2,569 cohort members); the deaths actually analysed are in analysed_n_proportional_hazards.csv.
## Run   : /opt/homebrew/bin/Rscript code/reconstructed/harmonised/proportional_hazards_fulldesign.R
## =============================================================================
source("code/reconstructed/harmonised/frame_2564.R")
S <- "proportional_hazards"

## ---- unweighted Schoenfeld test (step 163; unchanged, no design) ------------------------------------
a <- subset(A, primary==1)
cu <- log_fit(S, "ph_tests.csv", "unweighted coxph for cox.zph (Model 3)",
              coxph(as.formula(paste0("Surv(time,event)~lca+",C3)), data=a))
tu <- as.data.frame(cox.zph(cu)$table)
print(round(tu[c("lca","GLOBAL"),],4))

## ---- split follow-up at 5 years (steps 164, 167, 169) ----------------------------------------------
sp <- survSplit(Surv(time,event)~., data=a, cut=5, episode="period")
sp$period <- factor(sp$period, labels=c("0-5 yr",">5 yr"))
sp$late <- as.numeric(sp$period==">5 yr")
sp$L2 <- as.numeric(sp$lca=="Insomnia-fatigue")   * sp$late
sp$L3 <- as.numeric(sp$lca=="Somatic-depressive") * sp$late
sp$L4 <- as.numeric(sp$lca=="High symptom burden")* sp$late
rest <- subset(A, primary!=1)
for (v in setdiff(names(sp), names(rest))) rest[[v]] <- NA
full <- rbind(sp, rest[, names(sp)])
stopifnot(length(unique(full$SEQN)) == nrow(A))
d2 <- mk(full)                                               # design on all 70,190 participants
mi <- log_fit(S, "ph_tests.csv", "phenotype x period joint Wald (split data)",
              svycoxph(as.formula(paste0("Surv(tstart,time,event)~lca+L2+L3+L4+",C3)), design=subset(d2, analysed==1)))
rt <- regTermTest(mi, ~L2+L3+L4)
cat("\njoint Wald p =", signif(rt$p,4), " (df =", rt$df, ")\n")
per <- do.call(rbind, lapply(levels(sp$period), function(pp){
  mm <- log_fit(S, "ph_period_specific.csv", paste("period", pp),
                svycoxph(as.formula(paste0("Surv(tstart,time,event)~lca+",C3)), design=subset(d2, analysed==1 & period==pp)))
  s2 <- summary(mm)$coef; k <- grep("^lca", rownames(s2)); ci <- confint(mm)[k,,drop=FALSE]
  data.frame(period=pp, events=sum(sp$event[sp$period==pp]), term=gsub("^lca","",rownames(s2)[k]),
             HR=exp(s2[k,"coef"]), lo=exp(ci[,1]), hi=exp(ci[,2]), p=s2[k,ncol(s2)])}))
print(transform(per, HR95=sprintf("%.2f (%.2f-%.2f)",HR,lo,hi))[,c("period","events","term","HR95","p")], row.names=FALSE, digits=3)
write.csv(per, file.path(HOUT, "ph_period_specific.csv"), row.names=FALSE)

## ---- stratified baseline hazard (step 170) ---------------------------------------------------------
A$agec <- cut(A$age,  c(0,55,65,75,120), labels=c("<55","55-64","65-74",">=75"))
A$ydxc <- cut(A$ydx_i, c(-1,2,5,10,100),  labels=c("0-2","3-5","6-10",">10"))
C3s <- "sex+race4+educ+married+pir_i+pir_m+bmi_i+bmi_m+smoke+pa3+htn+dm+cvd+ydx_m+multi_primary+catype"
ms <- log_fit(S, "ph_stratified_baseline.csv", "strata(age group x years since diagnosis) + Model 3 without age, ydx_i",
              svycoxph(as.formula(paste0("Surv(time,event)~lca+strata(agec,ydxc)+",C3s)), design=subset(mk(A), analysed==1)))
ss <- summary(ms)$coef; k <- grep("^lca",rownames(ss)); ci <- confint(ms)[k,,drop=FALSE]
strat <- data.frame(term=gsub("^lca","",rownames(ss)[k]), HR=exp(ss[k,"coef"]), lo=exp(ci[,1]), hi=exp(ci[,2]), p=ss[k,ncol(ss)])
print(transform(strat, HR95=sprintf("%.2f (%.2f-%.2f)",HR,lo,hi))[,c("term","HR95","p")], row.names=FALSE, digits=3)
write.csv(strat, file.path(HOUT, "ph_stratified_baseline.csv"), row.names=FALSE)

## ---- ph_tests.csv, computed (rounded as the original typed values) --------------------------------
write.csv(data.frame(
  test=c("lca Schoenfeld residuals (unweighted diagnostic)","GLOBAL Schoenfeld (unweighted diagnostic)",
         "lca x follow-up period interaction (design-based joint Wald)"),
  statistic=c(round(tu["lca","chisq"],4), round(tu["GLOBAL","chisq"],4), NA),
  df=c(tu["lca","df"], tu["GLOBAL","df"], rt$df),
  p=c(round(tu["lca","p"],4), round(tu["GLOBAL","p"],4), signif(rt$p,4)),
  note=c("exposure satisfies PH","driven by age and years since diagnosis","no evidence of time-varying phenotype effect")),
  file.path(HOUT, "ph_tests.csv"), row.names=FALSE)
write_nlog(S)
cat("proportional_hazards_fulldesign.R done\n")
