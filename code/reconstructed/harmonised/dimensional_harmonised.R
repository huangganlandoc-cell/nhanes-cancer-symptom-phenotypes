## =============================================================================
## dimensional_harmonised.R  --  dimensional.R (Cox parts) on the harmonised 2,564-person frame
##
## Writes to review/reproduction/harmonised/:
##   RR_dimensional.csv, RR_dimensional_affective.csv, RR_dimensional_multiplicity.csv,
##   dimension_means_by_phenotype.csv, analysed_n_dimensional.csv
## Not re-run: RR_efa_loadings.csv (factor analysis of all adults; no Cox model, frame-independent).
##
## Changes relative to code/reconstructed/dimensional.R (original steps 281, 282, 284):
##   * frame: frame_2564.R instead of the CSV-rebuilt frame (2,557 analysed in Models 3 and 4);
##     every model, including the Model 2 row of the affective progression, uses the same 2,564
##     survivors (complete cases on the Model 4 covariates);
##   * specification "D total + predominance + sleep (Model 4)" failed in the original only because
##     the rebuilt frame had no variable `insomnia`; the harmonised frame has it, so D now fits and
##     is written as extra rows (no manuscript number uses it).
## RR_dimensional_multiplicity.csv takes its inputs from the harmonised cause-specific fits exactly
## as the original Python cell typed them (HR/limits with 2 decimals, P as printed by R, digits = 3).
## The weighted dimension means by phenotype are descriptive and are computed, as originally, on
## all 2,569 primary-cohort members (they do not depend on the covariate frame).
## Run   : /opt/homebrew/bin/Rscript code/reconstructed/harmonised/dimensional_harmonised.R (project root)
## =============================================================================
source("code/reconstructed/harmonised/frame_2564.R")
S <- "dimensional"

DIM <- read.csv(file.path(DER, "nhanes_dimension_scores.csv"))
H2 <- merge(A, DIM[,c("SEQN","som","aff","tot","som_resid")], by="SEQN", all.x=TRUE)
H2$sleep_dev <- abs(H2$sleep_h - 7)
stopifnot(all(!is.na(H2$tot[H2$analysed==1])))
P2 <- subset(mk(H2), analysed==1)
tt <- function(m,pat){s<-summary(m)$coef;k<-grep(pat,rownames(s));ci<-confint(m)[k,,drop=FALSE]
  data.frame(term=rownames(s)[k],HR=exp(s[k,"coef"]),lo=exp(ci[,1]),hi=exp(ci[,2]),p=s[k,ncol(s)])}

## ---- step 281 ----------------------------------------------------------------------------------
specs <- list(
 "A1 total score only (Model 3)"                    = list(e="I(tot/5)",                    cv=C3),
 "A2 total score only (Model 4)"                    = list(e="I(tot/5)",                    cv=C4),
 "B1 total + somatic predominance (Model 3)"        = list(e="I(tot/5)+som_resid",           cv=C3),
 "B2 total + somatic predominance (Model 4)"        = list(e="I(tot/5)+som_resid",           cv=C4),
 "C1 somatic and affective separately (Model 3)"    = list(e="I(som/3)+I(aff/3)",            cv=C3),
 "C2 somatic and affective separately (Model 4)"    = list(e="I(som/3)+I(aff/3)",            cv=C4),
 "D  total + predominance + sleep (Model 4)"        = list(e="I(tot/5)+som_resid+sleep_dev+insomnia", cv=C4))
res <- list()
for(nm in names(specs)){
  m <- try(svycoxph(as.formula(sprintf("Surv(time,event)~%s+%s",specs[[nm]]$e,specs[[nm]]$cv)), design=P2), silent=TRUE)
  if(inherits(m,"try-error")){cat(nm,": FAILED\n"); next}
  log_fit(S, "RR_dimensional.csv", nm, m)
  res[[nm]] <- cbind(model=nm, tt(m,"tot|som|aff|sleep_dev|insomnia"))
}
R <- do.call(rbind,res); rownames(R) <- NULL
print(transform(R, HR95=sprintf("%.2f (%.2f-%.2f)",HR,lo,hi))[,c("model","term","HR95","p")], row.names=FALSE, digits=3)
write.csv(R, file.path(HOUT, "RR_dimensional.csv"), row.names=FALSE)

## ---- step 282 ----------------------------------------------------------------------------------
means <- do.call(rbind, lapply(c("tot","som","aff"), function(v){
  m <- svyby(as.formula(paste0("~",v)), ~lca, subset(mk(H2),primary==1), svymean, na.rm=TRUE)
  data.frame(dimension=v, phenotype=as.character(m$lca), weighted_mean=m[,2], se=m[,3], sample="primary cohort, n = 2,569")}))
means2 <- do.call(rbind, lapply(c("tot","som","aff"), function(v){
  m <- svyby(as.formula(paste0("~",v)), ~lca, P2, svymean, na.rm=TRUE)
  data.frame(dimension=v, phenotype=as.character(m$lca), weighted_mean=m[,2], se=m[,3], sample="analysed, n = 2,564")}))
write.csv(rbind(means, means2), file.path(HOUT, "dimension_means_by_phenotype.csv"), row.names=FALSE)
prog <- do.call(rbind, lapply(list(c("Model 2",C2),c("Model 3",C3),c("Model 4",C4)), function(z){
  m <- log_fit(S, "RR_dimensional_affective.csv", paste(z[1], "(affective + somatic)"),
               svycoxph(as.formula(paste0("Surv(time,event)~I(aff/3)+I(som/3)+",z[2])), design=P2))
  cbind(model=z[1], tt(m,"aff|som"))}))
print(transform(prog, HR95=sprintf("%.2f (%.2f-%.2f)",HR,lo,hi))[,c("model","term","HR95","p")],row.names=FALSE,digits=3)
cs <- do.call(rbind, lapply(c("event","ev_ca","ev_cvd"), function(oc){
  m <- log_fit(S, "RR_dimensional_affective.csv", paste("Model 4 (affective + somatic),", oc),
               svycoxph(as.formula(sprintf("Surv(time,%s)~I(aff/3)+I(som/3)+%s",oc,C4)), design=P2))
  cbind(outcome=oc, tt(m,"aff|som"))}))
cat("\ncause-specific (Model 4):\n")
print(transform(cs, HR95=sprintf("%.2f (%.2f-%.2f)",HR,lo,hi))[,c("outcome","term","HR95","p")],row.names=FALSE,digits=3)
write.csv(rbind(cbind(block="model progression",prog), cbind(block="cause-specific",cs |> setNames(c("model",names(prog)[-1])))),
          file.path(HOUT, "RR_dimensional_affective.csv"), row.names=FALSE)

## ---- step 284 (Python -> R); inputs = harmonised cs as the original typed them --------------------
ord <- c("I(aff/3)","I(som/3)")
dimf <- do.call(rbind, lapply(list(c("event","All-cause"),c("ev_ca","Cancer"),c("ev_cvd","CVD/stroke")), function(z){
  x <- cs[cs$outcome==z[1],]; x <- x[match(ord, x$term),]
  data.frame(outcome=z[2], term=c("Affective dimension (per 3 pts)","Somatic dimension (per 3 pts)"),
             HR=typed_hr(x$HR), lo=typed_hr(x$lo), hi=typed_hr(x$hi), p=x$p)}))
dimf$p <- typed_p(dimf$p)                                  # the 6-row cause-specific print
dimf$p_bh <- bh_statsmodels(dimf$p); dimf$p_bonf <- p.adjust(dimf$p, "bonferroni")
write_py_csv(dimf, file.path(HOUT, "RR_dimensional_multiplicity.csv"))
print(dimf, row.names=FALSE)
rr <- function(h) (1-0.5^sqrt(h))/(1-0.5^sqrt(1/h)); ev <- function(x){ z <- max(x,1/x); z+sqrt(z*(z-1)) }
cat(sprintf("affective dimension E-value: %.2f (CI limit %.2f)\n", ev(rr(dimf$HR[1])), ev(rr(dimf$lo[1]))))
write_nlog(S)
cat("dimensional_harmonised.R done\n")
