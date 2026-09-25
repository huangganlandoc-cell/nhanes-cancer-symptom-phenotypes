## Primary analysis on the restricted cohort (excluding non-melanoma-skin-only survivors)
## Symptom phenotypes are defined by the 4-class LCA fitted in the FULL cancer-survivor
## cohort (n=3,268); only the analytic population is restricted, so phenotype definitions
## do not shift with the analytic sample.
## Run    : /opt/homebrew/bin/Rscript code/analysis_primary_exclNMS.R   (from the project root)
## Inputs : data/derived/analysis_frame.rds, data/derived/nhanes_cancer_design_frame.csv.gz
## Outputs: P_*.csv -> supporting/, analysis_frame_primary.rds -> data/derived/
##          (all of them to <dir> instead when run with OUT=<dir>)
suppressPackageStartupMessages({library(survey);library(survival);library(rms)})
DER <- "data/derived"; if (!dir.exists(DER)) stop("run from the project root")
OUT <- Sys.getenv("OUT"); if (nzchar(OUT)) dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
OCSV <- if (nzchar(OUT)) OUT else "supporting"; ORDS <- if (nzchar(OUT)) OUT else DER
dir.create(OCSV, showWarnings=FALSE)

G <- readRDS(file.path(DER,"analysis_frame.rds"))
nf <- read.csv(file.path(DER,"nhanes_cancer_design_frame.csv.gz"))[,c("SEQN","only_nms","any_mel")]
G  <- merge(G, nf, by="SEQN", all.x=TRUE)
G$primary <- as.numeric(G$inAnalysis==1 & !(G$only_nms %in% c(TRUE,"True")))
G$cvd <- as.numeric(G$cvd)

C2 <- "age+sex+race4+educ+married+pir_i+pir_m"
C3 <- paste(C2,"bmi_i+bmi_m+smoke+pa3+htn+dm+cvd+ydx_i+ydx_m+multi_primary+catype",sep="+")
mkdes <- function(d) svydesign(ids=~SDMVPSU,strata=~SDMVSTRA,weights=~wt,nest=TRUE,data=d)
des <- mkdes(G); P <- subset(des, primary==1)
tidy <- function(m,pat){s<-summary(m)$coef;k<-grep(pat,rownames(s));ci<-confint(m)[k,,drop=FALSE]
  data.frame(term=rownames(s)[k],HR=exp(s[k,"coef"]),lo=exp(ci[,1]),hi=exp(ci[,2]),p=s[k,ncol(s)],row.names=NULL)}
prog <- function(exposure,pat,outcome="event",design=P){
  do.call(rbind,lapply(1:3,function(i){
    rhs<-c("",paste0("+",C2),paste0("+",C3))[i]
    m<-svycoxph(as.formula(sprintf("Surv(time,%s)~%s%s",outcome,exposure,rhs)),design=design)
    cbind(model=paste0("Model ",i),outcome=outcome,tidy(m,pat))}))}

## ---- Table 1: weighted baseline in the primary cohort ----------------------
cn <- c("age","phq9_score","sleep_h","bmi_i","ydx_i")
ct <- c("sex","race4","educ","married","smoke","pa3","catype")
bn <- c("htn","dm","cvd","multi_primary")
rows <- list()
for(v in cn){
  mn <- svyby(as.formula(paste0("~",v)),~lca,P,svymean,na.rm=TRUE)
  ov <- svymean(as.formula(paste0("~",v)),P,na.rm=TRUE)
  pv <- regTermTest(svyglm(as.formula(paste0(v,"~lca")),design=P),~lca)$p
  rows[[length(rows)+1]] <- data.frame(variable=v,level="mean (SE)",
    Overall=sprintf("%.1f (%.2f)",coef(ov),SE(ov)),
    t(setNames(sprintf("%.1f (%.2f)",mn[,2],mn[,3]),mn$lca)),p=signif(pv,3),check.names=FALSE)
}
for(v in c(ct,bn)){
  fm <- as.formula(paste0("~factor(",v,")"))
  mn <- svyby(fm,~lca,P,svymean,na.rm=TRUE); ov <- svymean(fm,P,na.rm=TRUE)
  pv <- svychisq(as.formula(paste0("~",v,"+lca")),P)$p.value
  lv <- gsub(paste0("factor\\(",v,"\\)"),"",names(coef(ov)))
  for(i in seq_along(lv)) rows[[length(rows)+1]] <- data.frame(
    variable=if(i==1) v else "",level=lv[i],Overall=sprintf("%.1f",coef(ov)[i]*100),
    t(setNames(sprintf("%.1f",mn[,1+i]*100),mn$lca)),p=if(i==1) signif(pv,3) else NA,check.names=FALSE)
}
write.csv(do.call(rbind,rows),file.path(OCSV,"P_table1_weighted_baseline.csv"),row.names=FALSE)

## ---- primary + cause-specific, with multiplicity ---------------------------
main <- do.call(rbind,lapply(c("event","ev_ca","ev_cvd"),function(o) prog("lca","^lca",o)))
m3 <- subset(main,model=="Model 3"); m3$p_bh <- p.adjust(m3$p,"BH"); m3$p_bonf <- p.adjust(m3$p,"bonferroni")
write.csv(main,file.path(OCSV,"P_cox_main.csv"),row.names=FALSE); write.csv(m3,file.path(OCSV,"P_multiplicity_main.csv"),row.names=FALSE)

## ---- secondary exposures ---------------------------------------------------
sec <- rbind(prog("dep_cat","^dep_cat"),prog("sleep3","^sleep3"),
             prog("insomnia","^insomnia"),prog("I(phq9_score/5)","phq9_score"))
s3 <- subset(sec,model=="Model 3"); s3$p_bh <- p.adjust(s3$p,"BH")
write.csv(sec,file.path(OCSV,"P_cox_secondary.csv"),row.names=FALSE); write.csv(s3,file.path(OCSV,"P_multiplicity_secondary.csv"),row.names=FALSE)

## ---- prespecified subgroups ------------------------------------------------
SG <- list(`Breast cancer`=quote(has_Breast %in% c(TRUE,"True")),
           `Prostate cancer`=quote(has_Prostate %in% c(TRUE,"True")),
           `Colorectal cancer`=quote(has_Colorectal %in% c(TRUE,"True")),
           `Melanoma`=quote(any_mel %in% c(TRUE,"True")),
           `Gynecologic cancer`=quote(has_Gynecologic %in% c(TRUE,"True")),
           `Female`=quote(female==1),`Male`=quote(female==0),
           `Age < 65 yr`=quote(age<65),`Age >= 65 yr`=quote(age>=65),
           `<5 yr since dx`=quote(ydx_i<5),`>=5 yr since dx`=quote(ydx_i>=5))
sexfree <- c("Prostate cancer","Gynecologic cancer","Female","Male")
sg <- do.call(rbind,lapply(names(SG),function(nm){
  keep <- with(G, primary==1 & eval(SG[[nm]])); keep[is.na(keep)] <- FALSE
  G$.sg <- as.numeric(keep); s2 <- subset(mkdes(G), .sg==1)
  cv <- if(nm %in% sexfree) sub("\\+sex","",C2) else C2
  m <- try(svycoxph(as.formula(paste0("Surv(time,event)~lca+",cv)),design=s2),silent=TRUE)
  if(inherits(m,"try-error")) return(NULL)
  cbind(subgroup=nm,n=sum(G$.sg),events=sum(G$event[G$.sg==1],na.rm=TRUE),
        model=if(nm %in% sexfree) "Model 2 (sex omitted)" else "Model 2", tidy(m,"^lca"))}))
sg$p_bh <- p.adjust(sg$p,"BH"); write.csv(sg,file.path(OCSV,"P_cox_subgroups.csv"),row.names=FALSE)

## ---- restricted cubic splines ---------------------------------------------
rcs_curve <- function(var,knots=c(.05,.35,.65,.95)){
  d0 <- G[G$primary==1,]
  kn <- quantile(d0[[var]],knots,na.rm=TRUE)
  m  <- svycoxph(as.formula(sprintf("Surv(time,event)~rcs(%s,c(%s))+%s",var,
        paste(round(kn,3),collapse=","),C3)),design=P)
  grid <- seq(min(d0[[var]],na.rm=TRUE),max(d0[[var]],na.rm=TRUE),length.out=120)
  nd <- d0[rep(1,length(grid)),]; num <- sapply(nd,is.numeric)
  for(v in names(nd)[num])  nd[[v]] <- median(d0[[v]],na.rm=TRUE)
  for(v in names(nd)[!num]) nd[[v]] <- d0[[v]][1]
  nd$lca <- factor("Low symptom burden",levels=levels(G$lca)); nd[[var]] <- grid
  ref <- nd[1,,drop=FALSE]; ref[[var]] <- median(d0[[var]],na.rm=TRUE)
  lp <- predict(m,newdata=nd,type="lp"); lpr <- as.numeric(predict(m,newdata=ref,type="lp"))
  ## Interval for the log hazard ratio against the reference (cohort median): only the spline terms
  ## differ between nd and ref, so the variance is that of the contrast, zero at the reference.
  ## (predict(se.fit = TRUE) gives the SE of the whole centred linear predictor, which is not this.)
  tt <- delete.response(terms(m))
  mm <- function(d) model.matrix(tt, model.frame(tt, d, xlev = m$xlevels))[, names(coef(m)), drop = FALSE]
  D <- mm(nd) - mm(ref)[rep(1, nrow(nd)), , drop = FALSE]
  est <- drop(D %*% coef(m)); se <- sqrt(pmax(rowSums((D %*% vcov(m)) * D), 0))
  stopifnot(max(abs(est - (lp - lpr))) < 1e-8)
  data.frame(var=var,x=grid,HR=exp(est),lo=exp(est-1.96*se),hi=exp(est+1.96*se))
}
write.csv(rbind(rcs_curve("phq9_score"),rcs_curve("sleep_h")),file.path(OCSV,"P_rcs_curves.csv"),row.names=FALSE)

## ---- weighted KM, prevalence, landmark, 3-class sensitivity ----------------
km <- svykm(Surv(time,event)~lca,design=P,se=FALSE)
write.csv(do.call(rbind,lapply(names(km),function(l)
  data.frame(lca=l,time=km[[l]]$time,surv=km[[l]]$surv))),file.path(OCSV,"P_km_curves.csv"),row.names=FALSE)
wp <- svymean(~lca,P,na.rm=TRUE); ci <- confint(wp)
write.csv(data.frame(class=gsub("^lca","",names(coef(wp))),weighted_pct=coef(wp)*100,
  lo=ci[,1]*100,hi=ci[,2]*100,n=as.numeric(table(G$lca[G$primary==1])),
  deaths=as.numeric(tapply(G$event[G$primary==1],G$lca[G$primary==1],sum)),
  mean_maxpost=as.numeric(tapply(G$maxpost[G$primary==1],G$lca[G$primary==1],mean))),
  file.path(OCSV,"P_lca_prevalence.csv"),row.names=FALSE)
land <- do.call(rbind,lapply(c(0,1,2,3,5),function(L){
  G$.k <- as.numeric(G$primary==1 & G$time>L); s <- subset(mkdes(G), .k==1)
  m <- svycoxph(as.formula(paste0("Surv(time,event)~lca+",C3)),design=s)
  cbind(exclude_yr=L,n=sum(G$.k),events=sum(G$event[G$.k==1],na.rm=TRUE),tidy(m,"^lca"))}))
m3c <- svycoxph(as.formula(paste0("Surv(time,event)~lca3+",C3)),design=P)
full <- svycoxph(as.formula(paste0("Surv(time,event)~lca+",C3)),design=subset(des,inAnalysis==1))
write.csv(land,file.path(OCSV,"P_sens_landmark.csv"),row.names=FALSE)
write.csv(rbind(cbind(analysis="3-class solution",tidy(m3c,"^lca3")),
                cbind(analysis="full cohort (incl. non-melanoma skin only)",tidy(full,"^lca"))),
          file.path(OCSV,"P_sens_other.csv"),row.names=FALSE)
saveRDS(G,file.path(ORDS,"analysis_frame_primary.rds"))
cat("primary n =",sum(G$primary),"deaths =",sum(G$event[G$primary==1]),
    "person-years =",round(sum(G$time[G$primary==1]),0),
    "median FU =",round(median(G$time[G$primary==1]),2),"\n")
