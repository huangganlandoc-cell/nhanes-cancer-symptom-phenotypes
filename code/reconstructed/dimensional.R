## =============================================================================
## dimensional.R  --  continuous symptom dimensions instead of classes, and the exploratory
##                    two-factor check of the somatic / cognitive-affective split
##
## Regenerates (into review/reproduction/output/):
##   RR_dimensional.csv              total score, somatic predominance, somatic + affective (Models 3/4)
##   RR_dimensional_affective.csv    affective and somatic dimensions: Models 2-4 and cause-specific
##   RR_dimensional_multiplicity.csv BH / Bonferroni over the 6 dimension x outcome tests (typed-in HRs)
##   RR_efa_loadings.csv             2-factor ML factor analysis, varimax, all 34,022 adults
## Also writes input_checks/check_dimension_scores.csv (re-derivation of
##   data/derived/nhanes_dimension_scores.csv from the raw DPQ files) and
##   input_checks/check_dimensional_multiplicity_transcription.csv
##
## Reconstructed from the original analysis log:
##   step 278 (Python), first half: dimension scores and the somatic
##                                  residual -> translated to R and used ONLY to check the saved nhanes_dimension_scores.csv
##   steps 252/254 (R): frame H, mk(), C2/C3 and B_com (bio_m version) that
##                      steps 281-282 relied on (re-created exactly as in model4.R)
##   step 281 (R): RR_dimensional.csv (specification D failed in the
##                 original too: the variable `insomnia` does not exist in H; kept as a FAILED line)
##   step 282 (R): RR_dimensional_affective.csv
##   step 284 (Python/statsmodels): RR_dimensional_multiplicity.csv
##   step 288 (Python/scikit-learn): RR_efa_loadings.csv
##
## Translation of step 288 (scikit-learn FactorAnalysis(n_components=2, rotation="varimax",
## random_state=0)): R has no scikit-learn, so its algorithm is ported line by line below
## (maximum likelihood via the SVD fixed-point iteration of Barber 2012, stopping rule
## tol = 0.01 on the total log-likelihood; sign convention of sklearn's svd_flip; varimax
## without Kaiser normalisation, tol 1e-6, which is what sklearn's _ortho_rotation computes and
## equals stats::varimax(normalize = FALSE, eps = 1e-6)). With 9 variables and 2 + 10
## oversampled components, sklearn's randomised SVD spans the whole column space, so an exact
## SVD gives the same result. stats::factanal is printed alongside as an independent check.
## Run   : /opt/homebrew/bin/Rscript code/reconstructed/dimensional.R   (from the project root)
## =============================================================================
suppressPackageStartupMessages({library(survey); library(survival); library(poLCA); library(rms)})
if (!file.exists("data/derived/nhanes_dimension_scores.csv")) stop("run from the project root")
OUT <- "review/reproduction/output"; dir.create(file.path(OUT, "input_checks"), recursive = TRUE, showWarnings = FALSE)
DER <- "data/derived"; RAWD <- "data/raw_nhanes"

## helpers for the translated Python cells ----------------------------------------------------
py_round <- function(x, d) round(x * 10^d) / 10^d   # numpy / pandas round(): x*10^d, half to even
py_num <- function(x) vapply(x, function(v) {       # Python repr() of a float, as pandas writes it
  if (is.na(v)) return("")
  for (d in 1:17) { s <- formatC(v, digits = d, format = "g"); if (as.numeric(s) == v) break }
  s <- trimws(s); if (!grepl("[.e]", s)) s <- paste0(s, ".0"); s }, "")
write_py_csv <- function(df, file, index = NULL) {  # pandas DataFrame.to_csv(); index = row labels
  cols <- lapply(df, function(col) {
    if (is.logical(col)) ifelse(is.na(col), "", ifelse(col, "True", "False"))
    else if (is.integer(col)) ifelse(is.na(col), "", as.character(col))
    else if (is.numeric(col)) py_num(col)
    else { s <- as.character(col); s[is.na(s)] <- ""; q <- grepl('[,"\n]', s)
           s[q] <- paste0('"', gsub('"', '""', s[q]), '"'); s } })
  hdr <- names(df); if (!is.null(index)) { cols <- c(list(index), cols); hdr <- c("", hdr) }
  writeLines(c(paste(hdr, collapse = ","), do.call(paste, c(cols, sep = ","))), file)
}
bh_statsmodels <- function(p){ n <- length(p); o <- order(p); raw <- p[o]/(seq_len(n)/n)
  adj <- pmin(rev(cummin(rev(raw))), 1); out <- numeric(n); out[o] <- adj; out }

## ---- check of the saved input: step 278, first half (pandas/statsmodels -> R) ----------------
CY <- c(`2005`="_D", `2007`="_E", `2009`="_F", `2011`="_G", `2013`="_H", `2015`="_I", `2017`="_J")
DPQ9 <- sprintf("DPQ0%d0", 1:9)
SOM <- c("DPQ030","DPQ040","DPQ050","DPQ080")                 # sleep, fatigue, appetite, psychomotor
AFF <- c("DPQ010","DPQ020","DPQ060","DPQ070","DPQ090")        # anhedonia, mood, worthlessness, concentration, ideation
DQ <- do.call(rbind, lapply(names(CY), function(y) {
  d <- read.csv(file.path(RAWD, sprintf("DPQ%s_%s-%d.csv.gz", CY[[y]], y, as.integer(y) + 1)))[, c("SEQN", DPQ9)]
  for (c in DPQ9) d[[c]][!(d[[c]] <= 3) | is.na(d[[c]])] <- NA; d }))
DQ$som <- rowSums(DQ[, SOM]); DQ$aff <- rowSums(DQ[, AFF]); DQ$tot <- rowSums(DQ[, DPQ9])  # NA unless complete
gd <- read.csv(file.path(DER, "nhanes_model4_design_frame.csv.gz"))
sub <- gd[gd$inAnalysis == 1, "SEQN", drop = FALSE]
DIMr <- merge(sub, DQ[, c("SEQN","som","aff","tot")], by = "SEQN", all.x = TRUE)
DIMr <- DIMr[complete.cases(DIMr), ]
DIMr$som_resid <- residuals(lm(som ~ tot + I(tot^2), data = DIMr))
DIMs <- read.csv(file.path(DER, "nhanes_dimension_scores.csv"))
chkD <- data.frame(column = names(DIMs), n_saved = nrow(DIMs), n_rederived = nrow(DIMr),
                   same_ids = identical(as.numeric(DIMs$SEQN), as.numeric(DIMr$SEQN)),
                   max_abs_diff = sapply(names(DIMs), function(v) max(abs(DIMs[[v]] - DIMr[[v]]))))
cat("=== re-derivation of nhanes_dimension_scores.csv ===\n"); print(chkD, row.names = FALSE)
write.csv(chkD, file.path(OUT, "input_checks", "check_dimension_scores.csv"), row.names = FALSE)

## ---- steps 252 + 254: frame H, mk(), C2, C3, B_com --------------------------------------------
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

## ---- step 281 (verbatim apart from file paths) ------------------------------------------------
DIM <- read.csv(file.path(DER, "nhanes_dimension_scores.csv"))
H2 <- merge(H, DIM[,c("SEQN","som","aff","tot","som_resid")], by="SEQN", all.x=TRUE)
H2$sleep_dev <- abs(H2$sleep_h - 7)
P2 <- subset(mk(H2), primary==1)
tt <- function(m,pat){s<-summary(m)$coef;k<-grep(pat,rownames(s));ci<-confint(m)[k,,drop=FALSE]
  data.frame(term=rownames(s)[k],HR=exp(s[k,"coef"]),lo=exp(ci[,1]),hi=exp(ci[,2]),p=s[k,ncol(s)])}
C4 <- paste(C3,B_com,"antidep_i",sep="+")
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
  res[[nm]] <- cbind(model=nm, tt(m,"tot|som|aff|sleep_dev|insomnia"))
}
R <- do.call(rbind,res)
print(transform(R, HR95=sprintf("%.2f (%.2f-%.2f)",HR,lo,hi))[,c("model","term","HR95","p")],
      row.names=FALSE, digits=3)
write.csv(R, file.path(OUT, "RR_dimensional.csv"), row.names=FALSE)

## ---- step 282 (verbatim apart from file paths) ------------------------------------------------
a2 <- subset(H2, primary==1)
cat("=== weighted dimension means by phenotype ===\n")
for(v in c("tot","som","aff")){
  m <- svyby(as.formula(paste0("~",v)), ~lca, subset(mk(H2),primary==1), svymean, na.rm=TRUE)
  cat(sprintf("%-4s: %s\n", v, paste(sprintf("%s=%.2f", m$lca, m[,2]), collapse="  ")))
}
prog <- do.call(rbind, lapply(list(c("Model 2",C2),c("Model 3",C3),c("Model 4",C4)), function(z)
  cbind(model=z[1], tt(svycoxph(as.formula(paste0("Surv(time,event)~I(aff/3)+I(som/3)+",z[2])), design=P2),"aff|som"))))
print(transform(prog, HR95=sprintf("%.2f (%.2f-%.2f)",HR,lo,hi))[,c("model","term","HR95","p")],row.names=FALSE,digits=3)
cs <- do.call(rbind, lapply(c("event","ev_ca","ev_cvd"), function(oc)
  cbind(outcome=oc, tt(svycoxph(as.formula(sprintf("Surv(time,%s)~I(aff/3)+I(som/3)+%s",oc,C4)), design=P2),"aff|som"))))
cat("\ncause-specific (Model 4):\n")
print(transform(cs, HR95=sprintf("%.2f (%.2f-%.2f)",HR,lo,hi))[,c("outcome","term","HR95","p")],row.names=FALSE,digits=3)
write.csv(rbind(cbind(block="model progression",prog), cbind(block="cause-specific",cs |> setNames(c("model",names(prog)[-1])))),
          file.path(OUT, "RR_dimensional_affective.csv"), row.names=FALSE)

## ---- step 284 (Python/statsmodels -> R); numbers as typed in the cell ---------------------------
dimf <- data.frame(
  outcome=c("All-cause","All-cause","Cancer","Cancer","CVD/stroke","CVD/stroke"),
  term=rep(c("Affective dimension (per 3 pts)","Somatic dimension (per 3 pts)"), 3),
  HR=c(1.18,0.92,1.05,0.81,1.15,0.91), lo=c(1.02,0.78,0.82,0.63,0.81,0.64),
  hi=c(1.37,1.09,1.36,1.04,1.63,1.29), p=c(0.0257,0.3221,0.6795,0.1013,0.4427,0.5917))
dimf$p_bh <- bh_statsmodels(dimf$p); dimf$p_bonf <- p.adjust(dimf$p, "bonferroni")
stopifnot(isTRUE(all.equal(dimf$p_bh, p.adjust(dimf$p, "BH"), tolerance = 1e-12)))
write_py_csv(dimf, file.path(OUT, "RR_dimensional_multiplicity.csv"))
print(dimf, row.names=FALSE)
rr <- function(h) (1-0.5^sqrt(h))/(1-0.5^sqrt(1/h)); ev <- function(x){ z <- max(x,1/x); z+sqrt(z*(z-1)) }
cat(sprintf("affective dimension E-value: %.2f (CI limit %.2f)\n", ev(rr(1.18)), ev(rr(1.02))))
## transcription check against the cause-specific Model 4 fits of step 282
cmpd <- do.call(rbind, lapply(seq_len(nrow(dimf)), function(i){
  oc <- c("All-cause"="event","Cancer"="ev_ca","CVD/stroke"="ev_cvd")[[dimf$outcome[i]]]
  tm <- if (grepl("^Affective", dimf$term[i])) "I(aff/3)" else "I(som/3)"
  r <- cs[cs$outcome==oc & cs$term==tm,]
  data.frame(outcome=dimf$outcome[i], term=dimf$term[i], typed_HR=dimf$HR[i], HR=r$HR, typed_lo=dimf$lo[i], lo=r$lo,
             typed_hi=dimf$hi[i], hi=r$hi, typed_p=dimf$p[i], p=r$p,
             agrees = abs(round(r$HR,2)-dimf$HR[i])<1e-9 & abs(round(r$lo,2)-dimf$lo[i])<1e-9 &
                      abs(round(r$hi,2)-dimf$hi[i])<1e-9 & abs(round(r$p,4)-dimf$p[i])<1e-9) }))
cat("\n=== RR_dimensional_multiplicity transcription check ===\n"); print(cmpd, row.names=FALSE, digits=5)
write.csv(cmpd, file.path(OUT, "input_checks", "check_dimensional_multiplicity_transcription.csv"), row.names=FALSE)

## ---- step 288 (Python/scikit-learn -> R) ----------------------------------------------------
A2 <- read.csv(file.path(DER, "nhanes_all_adults_lca_input.csv.gz"))
X2 <- as.matrix(A2[, DPQ9])
sk_factor_analysis <- function(X, k, tol = 1e-2, max_iter = 1000) {   # port of sklearn FactorAnalysis.fit
  n <- nrow(X); p <- ncol(X); X <- sweep(X, 2, colMeans(X))
  nsqrt <- sqrt(n); llconst <- p * log(2 * pi) + k; var <- colMeans(X^2)
  psi <- rep(1, p); old_ll <- -Inf; SMALL <- 1e-12
  for (it in seq_len(max_iter)) {
    sqrt_psi <- sqrt(psi) + SMALL
    Y <- sweep(X, 2, sqrt_psi * nsqrt, "/")
    sv <- svd(Y, nu = k, nv = k)
    u <- sv$u; v <- sv$v; s <- sv$d[1:k]
    for (j in 1:k) { sg <- sign(u[which.max(abs(u[, j])), j]); v[, j] <- v[, j] * sg }  # svd_flip (u-based)
    unexp <- sum(Y^2) - sum(s^2)
    s <- s^2
    W <- t(v) * sqrt(pmax(s - 1, 0))                      # k x p
    W <- sweep(W, 2, sqrt_psi, "*")
    ll <- -n / 2 * (llconst + sum(log(s)) + unexp + sum(log(psi)))
    if ((ll - old_ll) < tol) break
    old_ll <- ll
    psi <- pmax(var - colSums(W^2), SMALL)
  }
  list(W = W, iter = it, loglik = ll)
}
sk_varimax <- function(comp, tol = 1e-6, max_iter = 100) {        # port of sklearn _ortho_rotation
  nrow <- nrow(comp); R <- diag(ncol(comp)); var <- 0
  for (i in seq_len(max_iter)) {
    cr <- comp %*% R
    tmp <- cr * matrix(colSums(cr^2) / nrow, nrow, ncol(comp), byrow = TRUE)
    sv <- svd(t(comp) %*% (cr^3 - tmp)); R <- sv$u %*% t(sv$v)
    var_new <- sum(sv$d); if (var != 0 && var_new < var * (1 + tol)) break; var <- var_new
  }
  comp %*% R
}
fa2 <- sk_factor_analysis(X2, 2)
Lr <- sk_varimax(t(fa2$W))                                          # items x factors
L2 <- data.frame(F1 = py_round(Lr[,1], 3), F2 = py_round(Lr[,2], 3))
rn <- c('Anhedonia','Depressed mood','Sleep disturbance','Fatigue','Appetite change',
        'Worthlessness','Concentration','Psychomotor','Suicidal ideation')
L2$dominant <- ifelse(abs(L2$F1) >= abs(L2$F2), "F1", "F2")
L2$literature <- c('AFF','AFF','SOM','SOM','SOM','AFF','AFF','SOM','AFF')
pd_mode <- function(x){ tb <- table(x); sort(names(tb)[tb == max(tb)])[1] }  # pandas Series.mode()[0]
grp <- sapply(c("F1","F2"), function(k) pd_mode(L2$literature[L2$dominant == k]))
lab <- setNames(c("F1","F2"), grp)                                  # literature label -> factor
L2$matches <- lab[L2$literature] == L2$dominant
L2$gap <- py_round(abs(L2$F1) - abs(L2$F2), 3)
cat(sprintf("\nEFA: n = %d, sklearn-style iterations = %d | F1<->%s  F2<->%s\n", nrow(A2), fa2$iter, grp[["F1"]], grp[["F2"]]))
print(cbind(item = rn, L2), row.names = FALSE)
cat("consistent with the literature split:", sum(L2$matches), "/ 9\n")
write_py_csv(L2, file.path(OUT, "RR_efa_loadings.csv"), index = rn)
## independent check with stats::factanal (ML on the correlation matrix, rescaled to raw units)
fa_r <- factanal(covmat = cov(X2), factors = 2, rotation = "none", n.obs = nrow(X2))
Lf <- varimax(unclass(fa_r$loadings) * sqrt(colMeans(sweep(X2, 2, colMeans(X2))^2)), normalize = FALSE, eps = 1e-6)$loadings
Lf <- unclass(Lf); Lf <- sweep(Lf, 2, sign(colSums(Lf)), "*")
cat("max |sklearn-port - factanal| loading difference:", signif(max(abs(Lf - sweep(Lr, 2, sign(colSums(Lr)), "*"))), 3), "\n")
cat("dimensional.R done\n")
