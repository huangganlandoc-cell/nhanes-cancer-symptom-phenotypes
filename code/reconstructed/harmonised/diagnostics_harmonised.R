## =============================================================================
## diagnostics_harmonised.R  --  ADDITIONAL analyses found on the same 2,557-person frame
##
## While checking the frames I found that the numbers in supporting/RR3_diagnostics.csv and
## supporting/RR3_identifiability.csv (original steps 318, 319, 320 and 346) were computed on the frame of step 281 (CSV-rebuilt, no modal imputation,
## 2,557 analysed) AND with a survey design built on the cohort rows only
## (Pd <- subset(mk(a), primary == 1) with a = the 2,569 cohort rows). These analyses are the source
## of several manuscript numbers: the Model 4 incremental-value P values (0.061 / 0.062), the Model 3
## value (0.054), the eleven-indicator model (condition number 718, standard-error inflation, joint
## P 0.184 / 0.128), the unweighted survey-cycle model (HR 0.74-0.87, P 0.075-0.478, rank / condition
## number) and the Model 2 + phenotype sleep model (nonlinearity P, short / long sleep HRs).
## This script re-runs them on the harmonised frame (frame_2564.R: 2,564 analysed, design on all
## 70,190, then subset). No other change. It also re-runs the six incremental-value models of
## step 207 (RR_incremental_joint.csv) on the harmonised domain.
##
## Writes to review/reproduction/harmonised/:
##   RR_incremental_joint.csv, RR_incremental_value.csv   (step 207 models, harmonised domain)
##   diagnostics_incremental.csv     phenotype joint P with total score, Models 3 and 4 (step 318)
##   diagnostics_identifiability.csv phenotype + 11 indicators (steps 319, 346)
##   diagnostics_cycle.csv           survey cycle as a covariate (step 318)
##   diagnostics_sleep.csv           sleep spline and categories, Model 2 covariates + phenotype (steps 319-320)
##   analysed_n_diagnostics.csv
## Run   : /opt/homebrew/bin/Rscript code/reconstructed/harmonised/diagnostics_harmonised.R (project root)
## =============================================================================
source("code/reconstructed/harmonised/frame_2564.R")
suppressPackageStartupMessages({library(rms); library(Hmisc)})
S <- "diagnostics"
DIM <- read.csv(file.path(DER, "nhanes_dimension_scores.csv"))
H2 <- merge(A, DIM[,c("SEQN","som","aff","tot","som_resid")], by="SEQN", all.x=TRUE)
Pd <- subset(mk(H2), analysed==1)
a  <- subset(H2, analysed==1)

## ---- step 207 models on the harmonised domain -------------------------------------------------
grab <- function(m,pat){s<-summary(m)$coef;k<-grep(pat,rownames(s))
  if(!length(k)) return(NULL); ci<-confint(m)[k,,drop=FALSE]
  data.frame(term=rownames(s)[k],HR=exp(s[k,"coef"]),lo=exp(ci[,1]),hi=exp(ci[,2]),p=s[k,ncol(s)])}
specs <- list(
 "1 phenotype only"                  = list(f=paste0("lca+",C3),                          show="^lca"),
 "2 PHQ-9 continuous only"           = list(f=paste0("I(phq9_score/5)+",C3),              show="phq9_score"),
 "3 phenotype + PHQ-9 continuous"    = list(f=paste0("lca+I(phq9_score/5)+",C3),          show="^lca|phq9_score"),
 "4 phenotype + PHQ-9 spline"        = list(f=paste0("lca+rcs(phq9_score,4)+",C3),        show="^lca"),
 "5 all 11 indicators, no phenotype" = list(f=paste0("phqi1+phqi2+phqi3+phqi4+phqi5+phqi6+phqi7+phqi8+phqi9+sleep3+insomnia+",C3), show="ZZZ"),
 "6 all 11 indicators + phenotype"   = list(f=paste0("lca+phqi1+phqi2+phqi3+phqi4+phqi5+phqi6+phqi7+phqi8+phqi9+sleep3+insomnia+",C3), show="^lca"))
inc <- list(); jt <- list()
for(nm in names(specs)){
  m <- log_fit(S, "RR_incremental_joint.csv", nm, svycoxph(as.formula(paste0("Surv(time,event)~",specs[[nm]]$f)), design=P))
  jp_lca <- if(grepl("lca",specs[[nm]]$f)) tryCatch(regTermTest(m,~lca)$p, error=function(e) NA) else NA
  jp_sev <- tryCatch(regTermTest(m, ~phqi1+phqi2+phqi3+phqi4+phqi5+phqi6+phqi7+phqi8+phqi9+sleep3+insomnia)$p, error=function(e) NA)
  jt[[nm]] <- data.frame(model=nm, joint_p_phenotype=as.numeric(jp_lca), joint_p_11_indicators=as.numeric(jp_sev))
  g <- grab(m, specs[[nm]]$show); if(!is.null(g)) inc[[nm]] <- cbind(model=nm, g)
}
INC <- do.call(rbind,inc); JT <- do.call(rbind,jt); rownames(INC) <- rownames(JT) <- NULL
write.csv(INC, file.path(HOUT, "RR_incremental_value.csv"), row.names=FALSE)
write.csv(JT,  file.path(HOUT, "RR_incremental_joint.csv"), row.names=FALSE)
print(JT, row.names=FALSE, digits=4)

## ---- step 318 part 2: phenotype + total score, Models 3 and 4 ----------------------------------
incr <- do.call(rbind, lapply(c("Model 3","Model 4"), function(lab){
  cv <- if(lab=="Model 3") C3 else C4
  m_ph <- log_fit(S, "diagnostics_incremental.csv", paste(lab, "+ linear total"),
                  svycoxph(as.formula(paste0("Surv(time,event)~lca+I(tot/5)+",cv)), design=Pd))
  m_sp <- log_fit(S, "diagnostics_incremental.csv", paste(lab, "+ spline total"),
                  svycoxph(as.formula(paste0("Surv(time,event)~lca+rcs(tot,4)+",cv)), design=Pd))
  data.frame(covariates=lab, quantity=c("phenotype joint P, linear total","total score P, linear total",
                                        "phenotype joint P, spline total"),
             value=c(regTermTest(m_ph,~lca)$p, regTermTest(m_ph,~I(tot/5))$p, regTermTest(m_sp,~lca)$p))}))
print(incr, row.names=FALSE, digits=4)
write.csv(incr, file.path(HOUT, "diagnostics_incremental.csv"), row.names=FALSE)

## ---- steps 319 and 346: phenotype + all 11 manifest indicators --------------------------------
IT  <- paste0("phqi",1:9, collapse="+")
IND <- paste0(IT,"+factor(sleepcat)+factor(slq050b)")
Xi <- model.matrix(as.formula(paste0("~lca+",IND,"+",C3)), data=a)
m_both <- log_fit(S, "diagnostics_identifiability.csv", "phenotype + 11 indicators + Model 3",
                  svycoxph(as.formula(paste0("Surv(time,event)~lca+",IND,"+",C3)), design=Pd))
m_ind  <- log_fit(S, "diagnostics_identifiability.csv", "11 indicators + Model 3",
                  svycoxph(as.formula(paste0("Surv(time,event)~",IND,"+",C3)), design=Pd))
m_ph   <- log_fit(S, "diagnostics_identifiability.csv", "phenotype + Model 3",
                  svycoxph(as.formula(paste0("Surv(time,event)~lca+",C3)), design=Pd))
kl <- grep("^lca", names(coef(m_both)))
se_both <- sqrt(diag(vcov(m_both)))[kl]; se_alone <- sqrt(diag(vcov(m_ph)))[grep("^lca", names(coef(m_ph)))]
sdb <- summary(m_both)$coef; ks <- grep("^lcaSomatic", rownames(sdb)); cib <- confint(m_both)[ks,,drop=FALSE]
idf <- data.frame(quantity=c("design matrix columns","rank","condition number",
                             "somatic-depressive HR with all indicators",
                             "phenotype block joint P (with indicators)",
                             "PHQ-9 item block joint P (9 items, with phenotype)",
                             "eleven-indicator block joint P (with phenotype)",
                             "eleven-indicator block joint P (without phenotype)",
                             "phenotype SE alone (IF/SD/High)","phenotype SE with indicators (IF/SD/High)",
                             "SE ratio with/alone (IF/SD/High)"),
  value=c(ncol(Xi), qr(Xi)$rank, signif(kappa(Xi),4),
          sprintf("%.2f (%.2f-%.2f), P = %.3f", exp(sdb[ks,"coef"]), exp(cib[,1]), exp(cib[,2]), sdb[ks,ncol(sdb)]),
          sprintf("%.4f", regTermTest(m_both,~lca)$p),
          sprintf("%.4f", regTermTest(m_both, as.formula(paste0("~",IT)))$p),
          sprintf("%.4f", regTermTest(m_both, as.formula(paste0("~",IND)))$p),
          sprintf("%.4f", regTermTest(m_ind,  as.formula(paste0("~",IND)))$p),
          paste(sprintf("%.3f", se_alone), collapse="/"), paste(sprintf("%.3f", se_both), collapse="/"),
          paste(sprintf("%.2f", se_both/se_alone), collapse="/")))
print(idf, row.names=FALSE)
write.csv(idf, file.path(HOUT, "diagnostics_identifiability.csv"), row.names=FALSE)

## ---- step 318 part 1: survey cycle as a covariate -----------------------------------------------
H2$cyc <- factor(H2$cycle); a$cyc <- factor(a$cycle)
e1 <- try(svycoxph(as.formula(paste0("Surv(time,event)~lca+cyc+",C3)), design=subset(mk(H2), analysed==1)), silent=TRUE)
u1 <- log_fit(S, "diagnostics_cycle.csv", "unweighted coxph, Model 3 + cycle",
              coxph(as.formula(paste0("Surv(time,event)~lca+cyc+",C3)), data=a))
X <- model.matrix(as.formula(paste0("~lca+cyc+",C3)), data=a)
su <- summary(u1)$coef; kc <- grep("^cyc", rownames(su))
cyc <- data.frame(quantity=c("design-based model with cycle","unweighted model: cycle HR range","unweighted model: cycle P range",
                             "model matrix columns / rank","model matrix condition number",
                             paste("maximum follow-up (years),", levels(a$cyc))),
  value=c(if(inherits(e1,"try-error")) paste("FAILED:", sub("\n.*","",attr(e1,"condition")$message)) else "fits",
          paste(sprintf("%.3f", range(exp(su[kc,"coef"]))), collapse=" to "),
          paste(sprintf("%.3f", range(su[kc,ncol(su)])), collapse=" to "),
          paste(ncol(X), qr(X)$rank, sep=" / "), sprintf("%.3g", kappa(X)),
          sprintf("%.2f", tapply(a$time, a$cyc, max))))
print(cyc, row.names=FALSE)
print(round(cbind(HR=exp(su[kc,"coef"]), p=su[kc,ncol(su)]),3))
write.csv(cyc, file.path(HOUT, "diagnostics_cycle.csv"), row.names=FALSE)

## ---- steps 319-320: sleep duration, Model 2 covariates + phenotype -----------------------------
CVs <- "lca+age+sex+race4+educ+married+pir_i+pir_m"
kn <- rcspline.eval(a$sleep_h, nk=4, knots.only=TRUE)          # knots from the analysed survivors
B <- matrix(NA_real_, nrow(H2), 3, dimnames=list(NULL, paste0("sp",1:3)))
ok <- !is.na(H2$sleep_h); B[ok,] <- rcspline.eval(H2$sleep_h[ok], knots=kn, inclx=TRUE)
H3 <- cbind(H2, as.data.frame(B)); Pd3 <- subset(mk(H3), analysed==1)
ms <- log_fit(S, "diagnostics_sleep.csv", "sleep spline + Model 2 covariates + phenotype",
              svycoxph(as.formula(paste0("Surv(time,event)~sp1+sp2+sp3+",CVs)), design=Pd3))
H3$sleepc <- relevel(factor(H3$sleepcat), ref="2")            # reference 6 to < 9 h
mc <- log_fit(S, "diagnostics_sleep.csv", "sleep categories (ref 6-<9 h) + Model 2 covariates + phenotype",
              svycoxph(as.formula(paste0("Surv(time,event)~sleepc+",CVs)), design=subset(mk(H3), analysed==1)))
s2 <- summary(mc)$coef; kq <- grep("sleepc", rownames(s2)); ci <- confint(mc)[kq,,drop=FALSE]
slp <- data.frame(quantity=c("knots (h)","spline overall P","spline nonlinearity P",
                             "short sleep (<6 h) vs 6-<9 h","long sleep (>=9 h) vs 6-<9 h"),
  value=c(paste(round(kn,2), collapse=", "), sprintf("%.4f", regTermTest(ms,~sp1+sp2+sp3)$p),
          sprintf("%.4f", regTermTest(ms,~sp2+sp3)$p),
          sprintf("%.2f (%.2f-%.2f), P = %.4f", exp(s2[kq,"coef"]), exp(ci[,1]), exp(ci[,2]), s2[kq,ncol(s2)])))
print(slp, row.names=FALSE)
write.csv(slp, file.path(HOUT, "diagnostics_sleep.csv"), row.names=FALSE)
write_nlog(S)
cat("diagnostics_harmonised.R done\n")
