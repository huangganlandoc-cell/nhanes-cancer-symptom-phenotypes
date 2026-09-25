## Symptom phenotypes and mortality among US cancer survivors, NHANES 2005-2018
## Survey-weighted Cox proportional hazards; LCA-derived symptom phenotypes
## Run    : /opt/homebrew/bin/Rscript code/analysis_cancer_symptom_mortality.R   (from the project root)
## Inputs : data/derived/nhanes_cancer_design_frame.csv.gz, data/derived/lca_fits_cancer.rds
##          (its columns only_nms and any_mel were added after this script was first run; they
##          are dropped on reading, as analysis_primary_exclNMS.R merges them in itself)
## Outputs: cox_*.csv, rcs_curves.csv, km_curves.csv, lca_prevalence_cancer.csv,
##          lca_profiles_cancer.csv -> supporting/; analysis_frame.rds -> data/derived/
##          (all of them to <dir> instead when run with OUT=<dir>)
suppressPackageStartupMessages({library(poLCA);library(survey);library(survival);library(rms)})
DER <- "data/derived"; if (!dir.exists(DER)) stop("run from the project root")
OUT <- Sys.getenv("OUT"); if (nzchar(OUT)) dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
OCSV <- if (nzchar(OUT)) OUT else "supporting"; ORDS <- if (nzchar(OUT)) OUT else DER
dir.create(OCSV, showWarnings=FALSE)

G  <- read.csv(file.path(DER,"nhanes_cancer_design_frame.csv.gz"))
G  <- G[, setdiff(names(G), c("only_nms","any_mel"))]
fits <- readRDS(file.path(DER,"lca_fits_cancer.rds"))
K <- 4                                   # selected by BIC and cAIC minima
ft <- fits[[K]]
a  <- subset(G, inAnalysis==1)
cl <- apply(ft$posterior, 1, which.max)
LAB <- c("High symptom burden","Low symptom burden","Insomnia-fatigue","Hypersomnia-somatic")
a$lca     <- factor(LAB[cl], levels=c("Low symptom burden","Insomnia-fatigue",
                                      "Hypersomnia-somatic","High symptom burden"))
a$maxpost <- apply(ft$posterior, 1, max)
cl3 <- apply(fits[[3]]$posterior, 1, which.max)
a$lca3 <- factor(c("Moderate/sleep","High symptom burden","Low symptom burden")[cl3],
                 levels=c("Low symptom burden","Moderate/sleep","High symptom burden"))
G <- merge(G, a[,c("SEQN","lca","lca3","maxpost")], by="SEQN", all.x=TRUE)

## ---- covariate preparation -------------------------------------------------
G$race4 <- factor(ifelse(G$race=="NH White","NH White",
             ifelse(G$race=="NH Black","NH Black",
             ifelse(G$race %in% c("Mexican American","Other Hispanic"),"Hispanic","Other"))),
             levels=c("NH White","NH Black","Hispanic","Other"))
G$educ    <- factor(G$educ,    levels=c(">High school","High school","<High school"))
G$married <- factor(G$married, levels=c("Married/partnered","Not partnered"))
G$smoke   <- factor(G$smoke,   levels=c("Never","Former","Current"))
G$sex     <- factor(ifelse(G$female==1,"Female","Male"), levels=c("Female","Male"))
G$pa3     <- factor(ifelse(is.na(G$pa_active),"Unknown",ifelse(G$pa_active==1,"Active","Inactive")),
                    levels=c("Active","Inactive","Unknown"))
G$catype  <- factor(ifelse(is.na(G$ca_type),"Other",G$ca_type))
G$catype  <- relevel(G$catype, ref="Melanoma/skin")
mi <- function(x) ifelse(is.na(x), median(x, na.rm=TRUE), x)
G$pir_m <- as.numeric(!is.na(G$INDFMPIR)); G$pir_i <- mi(G$INDFMPIR)
G$bmi_m <- as.numeric(!is.na(G$BMXBMI));   G$bmi_i <- mi(G$BMXBMI)
G$ydx_m <- as.numeric(!is.na(G$yrs_dx));   G$ydx_i <- mi(G$yrs_dx)
for (v in c("educ","smoke","married")) {
  mo <- names(sort(table(G[[v]]), decreasing=TRUE))[1]
  G[[v]][is.na(G[[v]]) & G$inAnalysis==1] <- mo
}
G$dep_cat <- factor(G$dep_cat, levels=c("None/minimal (0-4)","Mild (5-9)","Moderate+ (>=10)"))
G$sleep3  <- factor(G$sleepcat, levels=c(2,1,3), labels=c("6-8.9 h","<6 h",">=9 h"))
G$insomnia<- factor(ifelse(G$slq050b==1,"Yes","No"), levels=c("No","Yes"))

des <- svydesign(ids=~SDMVPSU, strata=~SDMVSTRA, weights=~wt, nest=TRUE, data=G)
sub <- subset(des, inAnalysis==1)

C2 <- "age+sex+race4+educ+married+pir_i+pir_m"
C3 <- paste(C2, "bmi_i+bmi_m+smoke+pa3+htn+dm+cvd+ydx_i+ydx_m+multi_primary+catype", sep="+")

tidy_hr <- function(m, pat) {
  s <- summary(m)$coef; k <- grep(pat, rownames(s)); ci <- confint(m)[k,,drop=FALSE]
  data.frame(term=rownames(s)[k], HR=exp(s[k,"coef"]), lo=exp(ci[,1]), hi=exp(ci[,2]),
             p=s[k, ncol(s)], row.names=NULL)
}
fit_prog <- function(exposure, pat, outcome="event") {
  do.call(rbind, lapply(1:3, function(i) {
    rhs <- c("", paste0("+",C2), paste0("+",C3))[i]
    m <- svycoxph(as.formula(sprintf("Surv(time,%s)~%s%s", outcome, exposure, rhs)), design=sub)
    cbind(model=paste0("Model ",i), outcome=outcome, tidy_hr(m, pat))
  }))
}

## ---- primary and secondary outcomes ---------------------------------------
main <- do.call(rbind, lapply(c("event","ev_ca","ev_cvd"),
                              function(o) fit_prog("lca","^lca",o)))
cont <- rbind(fit_prog("dep_cat","^dep_cat"), fit_prog("sleep3","^sleep3"),
              fit_prog("insomnia","^insomnia"), fit_prog("I(phq9_score/5)","phq9_score"))
sens <- rbind(
  cbind(analysis="3-class solution",
        fit_prog("lca3","^lca3")[fit_prog("lca3","^lca3")$model=="Model 3",]),
  cbind(analysis="exclude first 2 yr",
        { s2 <- subset(des, inAnalysis==1 & time>2)
          m <- svycoxph(as.formula(paste0("Surv(time,event)~lca+",C3)), design=s2)
          cbind(model="Model 3", outcome="event", tidy_hr(m,"^lca")) }))
write.csv(main,file.path(OCSV,"cox_main.csv"),row.names=FALSE)
write.csv(cont,file.path(OCSV,"cox_secondary_exposures.csv"),row.names=FALSE)
write.csv(sens,file.path(OCSV,"cox_sensitivity.csv"),row.names=FALSE)

## ---- prespecified subgroups ------------------------------------------------
subgroups <- list(
  `Breast cancer`      = "has_Breast=='True'|has_Breast==TRUE",
  `Prostate cancer`    = "has_Prostate=='True'|has_Prostate==TRUE",
  `Colorectal cancer`  = "has_Colorectal=='True'|has_Colorectal==TRUE",
  `Melanoma/skin`      = "has_Melanoma.skin=='True'|has_Melanoma.skin==TRUE",
  `Gynecologic cancer` = "has_Gynecologic=='True'|has_Gynecologic==TRUE",
  `Female`             = "female==1", `Male` = "female==0",
  `Age < 65 yr`        = "age<65", `Age >= 65 yr` = "age>=65",
  `<5 yr since dx`     = "ydx_i<5", `>=5 yr since dx` = "ydx_i>=5")
sg <- do.call(rbind, lapply(names(subgroups), function(nm) {
  ex <- parse(text=subgroups[[nm]])
  keep <- with(G, inAnalysis==1 & eval(ex))
  G$.sg <- as.numeric(keep & !is.na(keep))
  d2 <- svydesign(ids=~SDMVPSU, strata=~SDMVSTRA, weights=~wt, nest=TRUE, data=G)
  s2 <- subset(d2, .sg==1)
  n  <- sum(G$.sg); ev <- sum(G$event[G$.sg==1], na.rm=TRUE)
  m <- try(svycoxph(as.formula(paste0("Surv(time,event)~lca+",C2)), design=s2), silent=TRUE)
  if (inherits(m,"try-error")) return(NULL)
  cbind(subgroup=nm, n=n, events=ev, tidy_hr(m,"^lca"))
}))
write.csv(sg,file.path(OCSV,"cox_subgroups.csv"),row.names=FALSE)

## ---- restricted cubic splines ---------------------------------------------
rcs_curve <- function(var, knots) {
  kn <- quantile(G[[var]][G$inAnalysis==1], knots, na.rm=TRUE)
  fm <- as.formula(sprintf("Surv(time,event)~rcs(%s,c(%s))+%s", var,
                           paste(round(kn,3),collapse=","), C3))
  m <- svycoxph(fm, design=sub)
  rng <- range(G[[var]][G$inAnalysis==1], na.rm=TRUE)
  grid <- seq(rng[1], rng[2], length.out=120)
  nd <- G[G$inAnalysis==1,][rep(1,length(grid)),]
  num <- sapply(nd, is.numeric)
  for (v in names(nd)[num]) nd[[v]] <- median(G[[v]][G$inAnalysis==1], na.rm=TRUE)
  for (v in names(nd)[!num]) nd[[v]] <- G[[v]][G$inAnalysis==1][1]
  nd$lca <- factor("Low symptom burden", levels=levels(G$lca)); nd[[var]] <- grid
  ref <- nd; ref[[var]] <- median(G[[var]][G$inAnalysis==1], na.rm=TRUE)
  ## interval for the log hazard ratio against the reference: variance of the contrast in the
  ## spline terms (predict(se.fit = TRUE) would give that of the whole centred linear predictor)
  tt <- delete.response(terms(m))
  mm <- function(d) model.matrix(tt, model.frame(tt, d, xlev = m$xlevels))[, names(coef(m)), drop = FALSE]
  D <- mm(nd) - mm(ref[rep(1, nrow(nd)), , drop = FALSE])
  est <- drop(D %*% coef(m)); se <- sqrt(pmax(rowSums((D %*% vcov(m)) * D), 0))
  data.frame(var=var, x=grid, HR=exp(est), lo=exp(est-1.96*se), hi=exp(est+1.96*se))
}
rcsdat <- rbind(try(rcs_curve("phq9_score", c(.05,.35,.65,.95)), silent=TRUE),
                try(rcs_curve("sleep_h",    c(.05,.35,.65,.95)), silent=TRUE))
write.csv(rcsdat,file.path(OCSV,"rcs_curves.csv"),row.names=FALSE)

## ---- weighted KM and prevalence -------------------------------------------
km <- svykm(Surv(time,event)~lca, design=sub, se=FALSE)
kmd <- do.call(rbind, lapply(names(km), function(l)
  data.frame(lca=l, time=km[[l]]$time, surv=km[[l]]$surv)))
write.csv(kmd,file.path(OCSV,"km_curves.csv"),row.names=FALSE)
wp <- svymean(~lca, sub, na.rm=TRUE); ci <- confint(wp)
prev <- data.frame(class=gsub("^lca","",names(coef(wp))), weighted_pct=coef(wp)*100,
                   lo=ci[,1]*100, hi=ci[,2]*100,
                   n=as.numeric(table(G$lca[G$inAnalysis==1])),
                   deaths=as.numeric(tapply(G$event[G$inAnalysis==1], G$lca[G$inAnalysis==1], sum)),
                   mean_maxpost=as.numeric(tapply(G$maxpost[G$inAnalysis==1], G$lca[G$inAnalysis==1], mean)))
write.csv(prev,file.path(OCSV,"lca_prevalence_cancer.csv"),row.names=FALSE)

## ---- item-response profiles ----------------------------------------------
items <- c("Anhedonia","Depressed mood","Sleep disturbance","Fatigue","Appetite change",
           "Worthlessness","Concentration","Psychomotor change","Suicidal ideation")
pf <- do.call(rbind, lapply(1:9, function(j)
        data.frame(item=items[j], class=LAB, prob=ft$probs[[j]][,2]*100)))
pf <- rbind(pf, data.frame(item="Short sleep (<6 h)",  class=LAB, prob=ft$probs[[10]][,1]*100),
                data.frame(item="Long sleep (>=9 h)",  class=LAB, prob=ft$probs[[10]][,3]*100),
                data.frame(item="Insomnia complaint",  class=LAB, prob=ft$probs[[11]][,1]*100))
write.csv(pf,file.path(OCSV,"lca_profiles_cancer.csv"),row.names=FALSE)
saveRDS(G,file.path(ORDS,"analysis_frame.rds"))
cat("done\n")
