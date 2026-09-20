## NOT USED FOR THE MANUSCRIPT. With all covariates in the membership model this analysis could not
## be completed (one start ran > 45 min; two starts did not finish in 4 min; NaN variances), and a
## single start converged to a degenerate solution. The pre-specified fallback is
## code/analysis_inclusive_draws_score.R (review/step2_analysis_plan.md, deviations 8-9).
## Inclusive pseudo-class draws (Bray, Lanza and Tan 2015), specified as in
## review/step2_analysis_plan.md, section 1.2: the membership model of the latent class
## analysis includes the death indicator, the Nelson-Aalen cumulative hazard at each
## participant's follow-up time (White and Royston 2009) and all covariates of the
## analysis model (Model 3; Model 4). 100 draws, survey-weighted Cox, Rubin's rules.
## Run from the project root:  Rscript code/analysis_inclusive_draws.R
source("code/three_step_functions.R")
suppressPackageStartupMessages({library(poLCA); library(survey); library(parallel)})
options(survey.lonely.psu = "adjust")
DER <- "data/derived"; OUT <- "supporting"
NCORES <- as.integer(Sys.getenv("NCORES", max(1, parallel::detectCores() - 1)))

F0 <- readRDS(file.path(DER, "analysis_frame.rds"))
f4 <- readRDS(file.path(DER, "lca_fits_cancer.rds"))[[4]]
LEVS <- c("Low", "SleepFatigue", "SomaticDepressive", "High")
main_probs <- lapply(f4$probs, function(p) p[c(2, 3, 4, 1), , drop = FALSE])   # reorder to LEVS

A  <- readRDS(file.path(DER, "analysis_frame_primary.rds"))
G4 <- read.csv(file.path(DER, "nhanes_model4_design_frame.csv.gz"))
M4KEEP <- c("SEQN", "comorb_n_i", "LBXHGB_i", "LBXHGB_m", "LBXSAL_i", "LBXSAL_m",
            "egfr_i", "egfr_m", "func_lim_i", "n_rx_i", "n_rx_m", "antidep_i")
A <- merge(A, G4[, M4KEEP], by = "SEQN", all.x = TRUE)
A$cvd <- as.numeric(A$cvd); A$bio_m <- pmax(A$LBXSAL_m, A$egfr_m)
ITEMS <- c(paste0("phqi", 1:9), "sleepcat", "slq050b")
C3 <- "age+sex+race4+educ+married+pir_i+pir_m+bmi_i+bmi_m+smoke+pa3+htn+dm+cvd+ydx_i+ydx_m+multi_primary+catype"
C4 <- paste(C3, "comorb_n_i+LBXHGB_i+LBXHGB_m+LBXSAL_i+egfr_i+bio_m+func_lim_i+n_rx_i+n_rx_m+antidep_i", sep = "+")
CONT <- c("age", "pir_i", "bmi_i", "ydx_i", "comorb_n_i", "LBXHGB_i", "LBXSAL_i", "egfr_i", "n_rx_i")

## LCA sample: the 3,268 survivors of the measurement model, restricted to complete covariates
L <- subset(A, inAnalysis == 1)
L <- L[complete.cases(L[, all.vars(as.formula(paste("~", C4)))]), ]
km <- survfit(Surv(time, event) ~ 1, data = L)                     # unweighted Nelson-Aalen
L$H_NA <- stepfun(km$time, c(0, km$cumhaz))(L$time)
Ls <- L; for (v in CONT) Ls[[v]] <- as.numeric(scale(L[[v]]))
Ls$H_NA <- as.numeric(scale(L$H_NA))

des <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt, nest = TRUE, data = A)
cat("LCA sample:", nrow(L), "\n")

## match classes of a new fit to the main solution by item-response probabilities
match_classes <- function(fit) {
  perms <- as.matrix(expand.grid(1:4, 1:4, 1:4, 1:4)); perms <- perms[apply(perms, 1, function(r) length(unique(r)) == 4), ]
  dist <- function(p) sum(sapply(seq_along(main_probs), function(j) sum(abs(fit$probs[[j]][p, ] - main_probs[[j]]))))
  perms[which.min(apply(perms, 1, dist)), ]                 # new-class index for each of LEVS
}

rows <- list(); prof <- list()
for (m in c("C3", "C4")) {
  covs <- get(m)
  f <- as.formula(paste0("cbind(", paste(ITEMS, collapse = ","), ") ~ event + H_NA + ", covs))
  set.seed(20260918)
  t0 <- Sys.time()
  ## single start at the main solution's item-response probabilities: the inclusive model is
  ## meant to reproduce the main classes (a run with 10 random starts took hours per model)
  fi <- poLCA(f, data = Ls, nclass = 4, maxiter = 8000, nrep = 1, probs.start = f4$probs,
              verbose = FALSE, calc.se = FALSE)
  perm <- match_classes(fi)
  pst <- fi$posterior[, perm]; colnames(pst) <- LEVS
  maxdiff <- max(sapply(seq_along(main_probs), function(j) max(abs(fi$probs[[j]][perm, ] - main_probs[[j]]))))
  cat(sprintf("\n%s inclusive LCA: loglik %.1f, %.0f s, class sizes %s, max |item prob - main| %.3f\n", m,
              fi$llik, as.numeric(difftime(Sys.time(), t0, units = "secs")),
              paste(round(colMeans(pst), 3), collapse = "/"), maxdiff))
  prof[[m]] <- data.frame(model = m, item = rep(names(main_probs), sapply(main_probs, ncol)),
                          category = unlist(lapply(main_probs, function(p) seq_len(ncol(p)))),
                          do.call(rbind, lapply(seq_along(main_probs), function(j)
                            t(rbind(main = main_probs[[j]][3, ], inclusive = fi$probs[[j]][perm[3], ])))))
  ## 100 draws for the analysed primary-cohort members
  PD <- data.frame(SEQN = L$SEQN, pst)
  analysed <- A$primary == 1 & complete.cases(A[, all.vars(as.formula(paste("~", C4)))])
  draws <- mclapply(1:100, function(r) {
    set.seed(9000 + r)
    cl <- apply(as.matrix(PD[, LEVS]), 1, function(p) sample.int(4, 1, prob = p))
    A2 <- A; A2$lz <- factor(LEVS[cl][match(A$SEQN, PD$SEQN)], levels = LEVS)
    d2 <- subset(svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt, nest = TRUE, data = A2),
                 analysed & !is.na(lz))
    mm <- svycoxph(as.formula(paste("Surv(time, event) ~ lz +", covs)), design = d2)
    k <- paste0("lz", LEVS[-1]); b <- coef(mm)[k]; V <- vcov(mm)[k, k]
    cvec <- setNames(c(-1, 1, 0), k)                                  # SD - SF contrast
    c(b, SDvsSF = sum(cvec * b), vSF = V[1, 1], vSD = V[2, 2], vHigh = V[3, 3],
      vC = as.numeric(t(cvec) %*% V %*% cvec))
  }, mc.cores = NCORES)
  Dr <- do.call(rbind, draws)
  B <- Dr[, c(paste0("lz", LEVS[-1]), "SDvsSF")]; Vv <- Dr[, c("vSF", "vSD", "vHigh", "vC")]
  colnames(B) <- colnames(Vv) <- c("Sleep-fatigue vs low", "Somatic-depressive vs low",
                                   "High symptom burden vs low", "Somatic-depressive vs sleep-fatigue")
  rr <- rubin(B, Vv)
  rr <- transform(rr, model = m, HR = exp(logHR), lo = exp(logHR - 1.96 * se), hi = exp(logHR + 1.96 * se),
                  p = 2 * pnorm(-abs(logHR / se)))
  print(transform(rr, HR95 = sprintf("%.2f (%.2f-%.2f)", HR, lo, hi))[, c("term", "HR95", "p", "fmi")], row.names = FALSE)
  rows[[m]] <- rr
}
write.csv(do.call(rbind, rows), file.path(OUT, "RR5_inclusive_draws.csv"), row.names = FALSE)
write.csv(do.call(rbind, prof), file.path(OUT, "RR5_inclusive_profiles_somatic_depressive.csv"), row.names = FALSE)
