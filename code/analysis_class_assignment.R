## Class-assignment and measurement-sensitivity analyses (Table 4, Tables S7 to S8)
## Inputs : analysis_frame_primary.rds (or nhanes_model4_design_frame.csv + nhanes_dpq_raw.csv)
##          nhanes_model4_design_frame.csv (Model 4 covariate block)
##          nhanes_all_adults_lca_input.csv
## Verified: this Model 4 reproduces the manuscript's modal-assignment estimate
##           1.48 (1.13 to 1.93), the outcome-informed estimates 2.00 / 1.71,
##           the threshold analysis (entropy 0.714, Cramer V 0.497, no class
##           associated) and the sleep-item analysis (V 0.736, 1.29, 0.98 to 1.69).
## Outputs: RR3_inclusive_pseudoclass.csv  RR3_inclusive_profiles.csv
##          RR3_threshold_cox.csv  RR3_threshold_crosstab.csv  RR3_threshold_profiles.csv
##          RR3_localdep_cox.csv   RR3_localdep_crosstab.csv   RR3_localdep_profiles.csv
##          RR3_identifiability.csv  RR3_diagnostics.csv
suppressPackageStartupMessages({library(poLCA); library(survey); library(survival); library(rms)})
options(survey.lonely.psu = "adjust")
set.seed(20260908)

## ---------------------------------------------------------------- design helper
## Weights are divided by the number of pooled cycles; the design object is built
## on ALL participants and then subset to the survivor domain, as required for
## correct variance estimation under domain estimation.
mk <- function(d) svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt,
                            nest = TRUE, data = d)
tid <- function(m, pat) {
  s <- summary(m)$coef; k <- grep(pat, rownames(s)); ci <- confint(m)[k, , drop = FALSE]
  data.frame(term = gsub(pat, "", rownames(s)[k]), HR = exp(s[k, "coef"]),
             lo = exp(ci[, 1]), hi = exp(ci[, 2]), p = s[k, ncol(s)])
}
ent <- function(p) 1 - sum(-p * log(pmax(p, 1e-12))) / (nrow(p) * log(ncol(p)))

## analysis_frame_primary.rds carries the design, phenotype and primary flag but NOT the
## Model 4 covariate block; that block lives in nhanes_model4_design_frame.csv and is
## merged here. The merged object is saved as analysis_frame_model4.rds.
H  <- readRDS("analysis_frame_primary.rds")
G4 <- read.csv("nhanes_model4_design_frame.csv")
M4KEEP <- c("SEQN", "comorb_n_i", "LBXHGB_i", "LBXHGB_m", "LBXSAL_i", "LBXSAL_m",
            "egfr_i", "egfr_m", "func_lim_i", "n_rx_i", "n_rx_m", "antidep_i")
H2 <- merge(H, G4[, M4KEEP], by = "SEQN", all.x = TRUE)
H2$cvd <- as.numeric(H2$cvd)
## Phenotype labels in the saved frame predate the renaming; align them.
levels(H2$lca)[levels(H2$lca) == "Hypersomnia-somatic"] <- "Somatic-depressive"
levels(H2$lca)[levels(H2$lca) == "Insomnia-fatigue"]    <- "Sleep-fatigue"
H  <- H2
RW <- read.csv("nhanes_dpq_raw.csv")                 # raw 0-3 PHQ-9 item scores
D  <- merge(subset(H, inAnalysis == 1), RW, by = "SEQN", all.x = TRUE)
D$logt <- log(pmax(D$time, 0.05))
Pd <- subset(mk(H), primary == 1)

C2 <- "age+sex+race4+educ+married+pir_i+pir_m"
C3 <- paste(C2, "bmi_i+bmi_m+smoke+pa3+htn+dm+cvd+ydx_i+ydx_m+multi_primary+catype", sep = "+")
## Model 4 adds the comorbidity and frailty block plus antidepressant use. Continuous
## members are median-imputed with accompanying missingness indicators. Albumin and eGFR
## come from the same biochemistry panel, so their indicators are identical and are
## entered once as `bio_m`; entering them separately makes the information matrix exactly
## singular. `only_nms` is the primary-cohort DEFINITION, not a covariate.
H$bio_m <- pmax(H$LBXSAL_m, H$egfr_m)
H2 <- H
stopifnot(all(H2$LBXSAL_m == H2$egfr_m, na.rm = TRUE))
B_com <- paste("comorb_n_i", "LBXHGB_i", "LBXHGB_m", "LBXSAL_i", "egfr_i", "bio_m",
               "func_lim_i", "n_rx_i", "n_rx_m", sep = "+")
C4 <- paste(C3, B_com, "antidep_i", sep = "+")   # 43 terms with the phenotype

## ---- pseudo-class draws, pooled by Rubin's rules --------------------------------
## post : n x K posterior matrix, rows aligned to `ids`
## rhs  : covariate string
draw_pool <- function(post, rhs, ids, M = 50) {
  est <- lapply(seq_len(M), function(i) {
    cl <- apply(post, 1, function(pp) sample.int(length(pp), 1, prob = pp))
    dd <- data.frame(SEQN = ids, lz = factor(paste0("C", cl)))
    dd$lz <- relevel(dd$lz, ref = names(which.max(table(dd$lz))))
    Hh <- merge(H[, setdiff(names(H), "lz")], dd, by = "SEQN", all.x = TRUE)
    mm <- svycoxph(as.formula(paste0("Surv(time,event)~lz+", rhs)),
                   design = subset(mk(Hh), primary == 1 & !is.na(lz)))
    k <- grep("^lz", names(coef(mm)))
    list(b = coef(mm)[k], v = diag(vcov(mm))[k], nm = gsub("^lz", "", names(coef(mm))[k]))
  })
  B <- do.call(rbind, lapply(est, `[[`, "b")); V <- do.call(rbind, lapply(est, `[[`, "v"))
  qb <- colMeans(B); ub <- colMeans(V); bv <- apply(B, 2, var)
  tot <- ub + (1 + 1 / M) * bv                       # Rubin total variance
  data.frame(term = est[[1]]$nm, HR = exp(qb),
             lo = exp(qb - 1.96 * sqrt(tot)), hi = exp(qb + 1.96 * sqrt(tot)),
             p = 2 * pnorm(-abs(qb / sqrt(tot))),
             fmi = round((1 + 1 / M) * bv / tot, 3))
}

## ================================================================ A. outcome-informed
## The mortality indicator and log follow-up time are entered as predictors of class
## MEMBERSHIP (poLCA's multinomial membership model), so posteriors condition on the
## outcome. Class membership is then partly determined by mortality: this specification
## is NOT independent of the outcome and cannot be reproduced prospectively.
fI <- poLCA(cbind(phqi1, phqi2, phqi3, phqi4, phqi5, phqi6, phqi7, phqi8, phqi9,
                  sleepcat, slq050b) ~ event + logt,
            D, nclass = 4, maxiter = 8000, nrep = 10, verbose = FALSE, na.rm = TRUE)
cat("A. outcome-informed | entropy", round(ent(fI$posterior), 3),
    "| mean max posterior", round(mean(apply(fI$posterior, 1, max)), 3),
    "| class %", paste(round(fI$P * 100, 1), collapse = "/"), "\n")
incl <- rbind(cbind(model = "Model 3", draw_pool(fI$posterior, C3, D$SEQN)),
              cbind(model = "Model 4", draw_pool(fI$posterior, C4, D$SEQN)))
print(transform(incl, HR95 = sprintf("%.2f (%.2f-%.2f)", HR, lo, hi))[, c("model","term","HR95","p","fmi")],
      row.names = FALSE, digits = 3)
write.csv(incl, "RR3_inclusive_pseudoclass.csv", row.names = FALSE)
write.csv(round(rbind(do.call(rbind, lapply(1:9, function(j) fI$probs[[j]][, 2] * 100)),
                      short = fI$probs[[10]][, 1] * 100,
                      long  = fI$probs[[10]][, 3] * 100,
                      told  = fI$probs[[11]][, 1] * 100), 1),
          "RR3_inclusive_profiles.csv")

## ================================================================ B. item threshold
## Main analysis dichotomises at an item score >= 2 ("more than half the days").
## Here: >= 1 ("several days").
D2 <- D
for (i in 1:9) D2[[paste0("phqi", i)]] <- as.integer(D[[paste0("dpq", i, "_raw")]] >= 1) + 1L
f1 <- poLCA(cbind(phqi1, phqi2, phqi3, phqi4, phqi5, phqi6, phqi7, phqi8, phqi9,
                  sleepcat, slq050b) ~ 1,
            D2, nclass = 4, maxiter = 8000, nrep = 10, verbose = FALSE, na.rm = TRUE)
D2$lt1 <- factor(paste0("T", apply(f1$posterior, 1, which.max)))
ct <- table(original = D2$lca, threshold1 = D2$lt1)
cat("\nB. threshold >=1 | entropy", round(ent(f1$posterior), 3),
    "| Cramer V", round(sqrt(chisq.test(ct)$statistic / (sum(ct) * (min(dim(ct)) - 1))), 3), "\n")
Ht <- merge(H, D2[, c("SEQN", "lt1")], by = "SEQN", all.x = TRUE)
Ht$lt1 <- relevel(Ht$lt1, ref = names(which.max(table(D2$lt1))))
mt <- svycoxph(as.formula(paste0("Surv(time,event)~lt1+", C4)),
               design = subset(mk(Ht), primary == 1 & !is.na(lt1)))
print(transform(tid(mt, "^lt1"), HR95 = sprintf("%.2f (%.2f-%.2f)", HR, lo, hi))[, c("term","HR95","p")],
      row.names = FALSE, digits = 3)
write.csv(tid(mt, "^lt1"), "RR3_threshold_cox.csv", row.names = FALSE)
write.csv(as.data.frame.matrix(ct), "RR3_threshold_crosstab.csv")
write.csv(round(rbind(do.call(rbind, lapply(1:9, function(j) f1$probs[[j]][, 2] * 100)),
                      short = f1$probs[[10]][, 1] * 100, long = f1$probs[[10]][, 3] * 100,
                      told  = f1$probs[[11]][, 1] * 100), 1), "RR3_threshold_profiles.csv")

## ================================================================ C. local dependence
## The PHQ-9 sleep item (phqi3) overlaps in content with the two sleep indicators;
## refit without it.
fL <- poLCA(cbind(phqi1, phqi2, phqi4, phqi5, phqi6, phqi7, phqi8, phqi9,
                  sleepcat, slq050b) ~ 1,
            D, nclass = 4, maxiter = 8000, nrep = 10, verbose = FALSE, na.rm = TRUE)
D$lns <- factor(paste0("L", apply(fL$posterior, 1, which.max)))
ct2 <- table(original = D$lca, dropped = D$lns)
cat("\nC. PHQ-9 sleep item dropped | entropy", round(ent(fL$posterior), 3),
    "| Cramer V", round(sqrt(chisq.test(ct2)$statistic / (sum(ct2) * (min(dim(ct2)) - 1))), 3), "\n")
Hn <- merge(H, D[, c("SEQN", "lns")], by = "SEQN", all.x = TRUE)
Hn$lns <- relevel(Hn$lns, ref = names(which.max(table(D$lns))))
mL <- svycoxph(as.formula(paste0("Surv(time,event)~lns+", C4)),
               design = subset(mk(Hn), primary == 1 & !is.na(lns)))
print(transform(tid(mL, "^lns"), HR95 = sprintf("%.2f (%.2f-%.2f)", HR, lo, hi))[, c("term","HR95","p")],
      row.names = FALSE, digits = 3)
write.csv(tid(mL, "^lns"), "RR3_localdep_cox.csv", row.names = FALSE)
write.csv(as.data.frame.matrix(ct2), "RR3_localdep_crosstab.csv")
write.csv(round(rbind(do.call(rbind, lapply(1:8, function(j) fL$probs[[j]][, 2] * 100)),
                      short = fL$probs[[9]][, 1] * 100, long = fL$probs[[9]][, 3] * 100,
                      told  = fL$probs[[10]][, 1] * 100), 1), "RR3_localdep_profiles.csv")

## ================================================================ D. identifiability
## Does the phenotype + 11-indicator model actually fail to be identified?
IND <- paste0(paste0("phqi", 1:9, collapse = "+"), "+factor(sleepcat)+factor(slq050b)")
a  <- subset(H, primary == 1)
Xi <- model.matrix(as.formula(paste0("~lca+", IND, "+", C3)), data = a)
m_both <- svycoxph(as.formula(paste0("Surv(time,event)~lca+", IND, "+", C3)), design = Pd)
m_ind  <- svycoxph(as.formula(paste0("Surv(time,event)~", IND, "+", C3)),      design = Pd)
m_ph   <- svycoxph(as.formula(paste0("Surv(time,event)~lca+", C3)),            design = Pd)
cat(sprintf("\nD. design matrix: %d columns, rank %d, condition number %.3g\n",
            ncol(Xi), qr(Xi)$rank, kappa(Xi)))
cat(sprintf("   phenotype joint P = %.4f | indicator block P (with phenotype) = %.4f | without = %.4f\n",
            regTermTest(m_both, ~lca)$p,
            regTermTest(m_both, as.formula(paste0("~", IND)))$p,
            regTermTest(m_ind,  as.formula(paste0("~", IND)))$p))
se_alone <- sqrt(diag(vcov(m_ph)))[grep("^lca", names(coef(m_ph)))]
se_both  <- sqrt(diag(vcov(m_both)))[grep("^lca", names(coef(m_both)))]
cat("   phenotype SE alone:", paste(sprintf("%.3f", se_alone), collapse = "/"),
    "-> with indicators:",   paste(sprintf("%.3f", se_both),  collapse = "/"), "\n")
write.csv(data.frame(
  quantity = c("design matrix columns", "rank", "condition number",
               "phenotype joint P", "indicator block P (with phenotype)",
               "indicator block P (without phenotype)", "phenotype SE alone",
               "phenotype SE with indicators"),
  value = c(ncol(Xi), qr(Xi)$rank, signif(kappa(Xi), 4),
            signif(regTermTest(m_both, ~lca)$p, 4),
            signif(regTermTest(m_both, as.formula(paste0("~", IND)))$p, 4),
            signif(regTermTest(m_ind,  as.formula(paste0("~", IND)))$p, 4),
            paste(sprintf("%.3f", se_alone), collapse = "/"),
            paste(sprintf("%.3f", se_both),  collapse = "/"))),
  "RR3_identifiability.csv", row.names = FALSE)

## ================================================================ E. survey cycle
## Design-based model fails; the unweighted model converges. Maximum observable
## follow-up is a deterministic function of cycle (common administrative censoring).
a$cyc <- factor(a$cycle)
e_svy <- try(svycoxph(as.formula(paste0("Surv(time,event)~lca+cyc+", C3)),
                      design = subset(mk(a), primary == 1)), silent = TRUE)
u_unw <- try(coxph(as.formula(paste0("Surv(time,event)~lca+cyc+", C3)), data = a), silent = TRUE)
cat("\nE. cycle | design-based:", if (inherits(e_svy, "try-error")) "singular" else "fits",
    "| unweighted:", if (inherits(u_unw, "try-error")) "fails" else "fits", "\n")
cat("   maximum follow-up by cycle:\n"); print(round(tapply(a$time, a$cycle, max), 2))

## ================================================================ F. sleep spline
B <- rcspline.eval(a$sleep_h, nk = 4, inclx = TRUE); colnames(B) <- paste0("sp", 1:ncol(B))
a2 <- cbind(a, as.data.frame(B)); Pd2 <- subset(mk(a2), primary == 1)
CVs <- "lca+age+sex+race4+educ+married+pir_i+pir_m"
ms <- svycoxph(as.formula(paste0("Surv(time,event)~sp1+sp2+sp3+", CVs)), design = Pd2)
cat(sprintf("\nF. sleep spline: overall P = %.4f | nonlinearity P = %.4f\n",
            regTermTest(ms, ~sp1 + sp2 + sp3)$p, regTermTest(ms, ~sp2 + sp3)$p))
mc <- svycoxph(as.formula(paste0("Surv(time,event)~factor(sleepcat)+", CVs)), design = Pd2)
print(transform(tid(mc, "factor\\(sleepcat\\)"), HR95 = sprintf("%.2f (%.2f-%.2f)", HR, lo, hi))[, c("term","HR95","p")],
      row.names = FALSE, digits = 3)

cat("\nDONE\n")
