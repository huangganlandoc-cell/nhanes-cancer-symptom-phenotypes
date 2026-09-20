## Additional analyses requested by adversarial review (nature-reviewer skill)
## R1-M1 incremental value; R1-M2 phenotype-vs-phenotype contrast; R1-M5 events per class;
## R1-M6 landmark detail; R1-M7 time-since-diagnosis interaction; R1-M13 survey cycle.
suppressPackageStartupMessages({library(survey);library(survival);library(rms);library(poLCA)})
options(survey.lonely.psu="adjust")

G <- readRDS("analysis_frame_primary.rds")
G$cvd <- as.numeric(G$cvd)
levels(G$lca)[levels(G$lca)=="Hypersomnia-somatic"] <- "Somatic-depressive"
C2 <- "age+sex+race4+educ+married+pir_i+pir_m"
C3 <- paste(C2,"bmi_i+bmi_m+smoke+pa3+htn+dm+cvd+ydx_i+ydx_m+multi_primary+catype",sep="+")
mk  <- function(d) svydesign(ids=~SDMVPSU, strata=~SDMVSTRA, weights=~wt, nest=TRUE, data=d)
des <- mk(G); P <- subset(des, primary==1)
tid <- function(m,pat){s<-summary(m)$coef;k<-grep(pat,rownames(s));ci<-confint(m)[k,,drop=FALSE]
  data.frame(term=gsub(pat,"",rownames(s)[k]),HR=exp(s[k,"coef"]),lo=exp(ci[,1]),hi=exp(ci[,2]),p=s[k,ncol(s)])}

cat("##### R1-M5: events, person-years and weighted rate per phenotype #####\n")
a <- subset(G, primary==1)
ev <- data.frame(phenotype=levels(a$lca),
  n = as.numeric(table(a$lca)),
  deaths = as.numeric(tapply(a$event, a$lca, sum)),
  person_years = round(as.numeric(tapply(a$time, a$lca, sum)),0))
ev$crude_rate_per_1000py <- round(ev$deaths/ev$person_years*1000, 1)
wr <- svyby(~event, ~lca, P, svymean, na.rm=TRUE)
ev$weighted_pct_dead <- round(wr[,2]*100,1); ev$se <- round(wr[,3]*100,2)
print(ev, row.names=FALSE)
write.csv(ev,"RR_events_per_phenotype.csv",row.names=FALSE)

cat("\n##### R1-M2: somatic-depressive vs insomnia-fatigue contrast (design-based) #####\n")
a2 <- a; a2$lca_if <- relevel(a2$lca, ref="Insomnia-fatigue")
P2 <- subset(mk(rbind(a2, transform(subset(G,primary==0), lca_if=NA))), primary==1)
m_if <- svycoxph(as.formula(paste0("Surv(time,event)~lca_if+",C3)), design=P2)
ct <- tid(m_if,"^lca_if")
print(transform(ct, HR95=sprintf("%.2f (%.2f-%.2f)",HR,lo,hi))[,c("term","HR95","p")],
      row.names=FALSE, digits=3)
write.csv(ct,"RR_contrast_vs_insomnia_fatigue.csv",row.names=FALSE)

cat("\n##### R1-M1: incremental value of phenotype over severity #####\n")
mods <- list(
  "phenotype only"                  = paste0("lca+",C3),
  "PHQ-9 continuous only"           = paste0("I(phq9_score/5)+",C3),
  "phenotype + PHQ-9 continuous"    = paste0("lca+I(phq9_score/5)+",C3),
  "phenotype + PHQ-9 spline"        = paste0("lca+rcs(phq9_score,4)+",C3),
  "all 11 indicators, no phenotype" = paste0("phqi1+phqi2+phqi3+phqi4+phqi5+phqi6+phqi7+phqi8+phqi9+sleep3+insomnia+",C3),
  "all 11 indicators + phenotype"   = paste0("lca+phqi1+phqi2+phqi3+phqi4+phqi5+phqi6+phqi7+phqi8+phqi9+sleep3+insomnia+",C3))
inc <- do.call(rbind, lapply(names(mods), function(nm){
  m <- try(svycoxph(as.formula(paste0("Surv(time,event)~",mods[[nm]])), design=P), silent=TRUE)
  if(inherits(m,"try-error")) return(data.frame(model=nm,term="FAILED",HR=NA,lo=NA,hi=NA,p=NA,joint_p=NA))
  jp <- if(grepl("phenotype", nm)) tryCatch(regTermTest(m, ~lca)$p, error=function(e) NA) else NA
  r <- if(grepl("phenotype", nm)) tid(m,"^lca") else tid(m,"phq9_score|^sleep3|^insomnia")
  cbind(model=nm, r, joint_p=as.numeric(jp))}))
print(transform(inc, HR95=ifelse(is.na(HR),"-",sprintf("%.2f (%.2f-%.2f)",HR,lo,hi)))[,c("model","term","HR95","p","joint_p")],
      row.names=FALSE, digits=3)
write.csv(inc,"RR_incremental_value.csv",row.names=FALSE)

cat("\n##### R1-M7: phenotype x time-since-diagnosis interaction (design-based) #####\n")
G$late5 <- as.numeric(G$ydx_i >= 5)
G$I2 <- as.numeric(G$lca=="Insomnia-fatigue")  * G$late5
G$I3 <- as.numeric(G$lca=="Somatic-depressive")* G$late5
G$I4 <- as.numeric(G$lca=="High symptom burden")*G$late5
Pi <- subset(mk(G), primary==1)
mi <- svycoxph(as.formula(paste0("Surv(time,event)~lca+late5+I2+I3+I4+",C3)), design=Pi)
rt <- regTermTest(mi, ~I2+I3+I4)
cat("joint interaction Wald p =", signif(rt$p,4), " (df =", rt$df, ")\n")
si <- summary(mi)$coef
ii <- do.call(rbind, lapply(c("I2","I3","I4"), function(v){ ci<-confint(mi)[v,]
  data.frame(term=c(I2="Insomnia-fatigue",I3="Somatic-depressive",I4="High symptom burden")[v],
             ratio_late_vs_early=exp(si[v,"coef"]), lo=exp(ci[1]), hi=exp(ci[2]), p=si[v,ncol(si)])}))
print(transform(ii, R95=sprintf("%.2f (%.2f-%.2f)",ratio_late_vs_early,lo,hi))[,c("term","R95","p")],
      row.names=FALSE, digits=3)
ii$joint_p <- rt$p
write.csv(ii,"RR_interaction_time_since_dx.csv",row.names=FALSE)
cat("events by stratum: <5yr =", sum(a$event[a$ydx_i<5]), " >=5yr =", sum(a$event[a$ydx_i>=5]),
    "| n:", sum(a$ydx_i<5), "/", sum(a$ydx_i>=5), "\n")

cat("\n##### R1-M6: landmark series with intervals and event counts #####\n")
land <- do.call(rbind, lapply(c(0,1,2,3,5), function(L){
  G$.k <- as.numeric(G$primary==1 & G$time > L)
  s <- subset(mk(G), .k==1)
  m <- svycoxph(as.formula(paste0("Surv(time,event)~lca+",C3)), design=s)
  cbind(exclude_yr=L, n=sum(G$.k), events=sum(G$event[G$.k==1],na.rm=TRUE),
        person_years=round(sum(G$time[G$.k==1]),0), tid(m,"^lca"))}))
print(transform(land, HR95=sprintf("%.2f (%.2f-%.2f)",HR,lo,hi))[,c("exclude_yr","n","events","term","HR95","p")],
      row.names=FALSE, digits=3)
write.csv(land,"RR_landmark_full.csv",row.names=FALSE)

cat("\n##### R1-M13: survey cycle #####\n")
G$cyc <- factor(G$cycle)
Pc <- subset(mk(G), primary==1)
m_cyc  <- svycoxph(as.formula(paste0("Surv(time,event)~lca+cyc+",C3)), design=Pc)
m_strc <- svycoxph(as.formula(paste0("Surv(time,event)~lca+strata(cyc)+",C3)), design=Pc)
G2 <- G[!(G$primary==1 & G$cycle=="2005-2006"),]
m_no05 <- svycoxph(as.formula(paste0("Surv(time,event)~lca+",C3)), design=subset(mk(G2), primary==1))
cyc <- rbind(cbind(spec="cycle as covariate", tid(m_cyc,"^lca")),
             cbind(spec="baseline hazard stratified on cycle", tid(m_strc,"^lca")),
             cbind(spec="excluding 2005-2006 cycle", tid(m_no05,"^lca")))
print(transform(cyc, HR95=sprintf("%.2f (%.2f-%.2f)",HR,lo,hi))[,c("spec","term","HR95","p")],
      row.names=FALSE, digits=3)
write.csv(cyc,"RR_cycle_sensitivity.csv",row.names=FALSE)
cat("\ndeaths by cycle:\n"); print(tapply(a$event, a$cycle, sum))
cat("\ndone\n")
