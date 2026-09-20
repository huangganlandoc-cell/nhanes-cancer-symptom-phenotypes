## =============================================================================
## incremental_value.R  --  does phenotype add to symptom severity? (joint Wald tests)
##
## Regenerates (into review/reproduction/output/):
##   RR_incremental_joint.csv   design-based joint Wald P for the phenotype block and for the
##                              block of 11 manifest indicators, six Model 3 specifications
##   (the companion RR_incremental_value.csv is written too, for convenience; it already has a
##    script in code/ and is not one of the targets)
##
## Reconstructed from the original analysis log:
##   step 205: created analysis_reviewer_response.R (first version)
##   step 206: source()d that script; it stopped with an error in its
##             own R1-M1 block, but had already defined G, C2, C3,
##             mk(), des and P in the global environment
##   step 207: the six models and joint tests -> RR_incremental_joint.csv
## Only the lines of the sourced script that define objects used by step 207 are repeated here.
## Run   : /opt/homebrew/bin/Rscript code/reconstructed/incremental_value.R  (from the project root)
## =============================================================================
suppressPackageStartupMessages({library(survey); library(survival); library(rms); library(poLCA)})
if (!file.exists("data/derived/analysis_frame_primary.rds")) stop("run from the project root")
OUT <- "review/reproduction/output"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

## ---- analysis_reviewer_response.R as sourced in step 206 (setup lines) ---------------
options(survey.lonely.psu="adjust")
G <- readRDS("data/derived/analysis_frame_primary.rds")
G$cvd <- as.numeric(G$cvd)
levels(G$lca)[levels(G$lca)=="Hypersomnia-somatic"] <- "Somatic-depressive"
C2 <- "age+sex+race4+educ+married+pir_i+pir_m"
C3 <- paste(C2,"bmi_i+bmi_m+smoke+pa3+htn+dm+cvd+ydx_i+ydx_m+multi_primary+catype",sep="+")
mk  <- function(d) svydesign(ids=~SDMVPSU, strata=~SDMVSTRA, weights=~wt, nest=TRUE, data=d)
des <- mk(G); P <- subset(des, primary==1)

## ---- step 207 (verbatim apart from file paths) ----------------------------------------
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
  m <- try(svycoxph(as.formula(paste0("Surv(time,event)~",specs[[nm]]$f)), design=P), silent=TRUE)
  if(inherits(m,"try-error")){ cat(nm,": FAILED\n"); next }
  jp_lca  <- if(grepl("lca",specs[[nm]]$f)) tryCatch(regTermTest(m,~lca)$p, error=function(e) NA) else NA
  jp_sev  <- tryCatch(regTermTest(m, ~phqi1+phqi2+phqi3+phqi4+phqi5+phqi6+phqi7+phqi8+phqi9+sleep3+insomnia)$p,
                      error=function(e) NA)
  jt[[nm]] <- data.frame(model=nm, joint_p_phenotype=as.numeric(jp_lca), joint_p_11_indicators=as.numeric(jp_sev))
  g <- grab(m, specs[[nm]]$show); if(!is.null(g)) inc[[nm]] <- cbind(model=nm, g)
}
INC <- do.call(rbind,inc); JT <- do.call(rbind,jt)
print(transform(INC, HR95=sprintf("%.2f (%.2f-%.2f)",HR,lo,hi))[,c("model","term","HR95","p")], row.names=FALSE, digits=3)
cat("\n--- joint tests ---\n"); print(JT, row.names=FALSE, digits=3)
write.csv(INC, file.path(OUT, "RR_incremental_value.csv"), row.names=FALSE)
write.csv(JT,  file.path(OUT, "RR_incremental_joint.csv"), row.names=FALSE)
cat("incremental_value.R done\n")
