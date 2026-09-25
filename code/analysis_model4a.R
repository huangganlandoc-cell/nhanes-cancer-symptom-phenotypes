## Model 4a: Model 3 plus the diagnosed-condition count and eGFR only (post hoc, 2026-09-24).
## Model 4 also adjusts for haemoglobin, serum albumin, functional limitation, medication count and
## antidepressant use, some of which may lie on the causal pathway or restate the symptoms; Model 4a
## keeps only the covariates that are plausibly confounders. Same data, class assignment and
## three-step machinery as code/analysis_three_step.R, which is left untouched.
##
## Run from the project root:
##   Rscript code/analysis_model4a.R fit    # modal (svycoxph) + three-step point estimates
##   Rscript code/analysis_model4a.R rep    # JKn standard errors re-running steps 1-3 per replicate
##   Rscript code/analysis_model4a.R boot   # the same with 500 Rao-Wu bootstrap replicates
## Outputs: supporting/RR6_model4a_*.csv / .rds
source("code/three_step_functions.R")
suppressPackageStartupMessages({library(survey); library(poLCA); library(parallel)})
options(survey.lonely.psu = "adjust")
DER <- "data/derived"; OUT <- "supporting"
stage <- commandArgs(TRUE)[1]; if (is.na(stage)) stage <- "fit"
NCORES <- as.integer(Sys.getenv("NCORES", max(1, detectCores() - 1)))

## ---------------------------------------------------------------- data (as in analysis_three_step.R)
F0 <- readRDS(file.path(DER, "analysis_frame.rds"))
f4 <- readRDS(file.path(DER, "lca_fits_cancer.rds"))[[4]]
LEVS <- c("Low", "SleepFatigue", "SomaticDepressive", "High")
post <- f4$posterior[, c(2, 3, 4, 1)]; colnames(post) <- LEVS
PD <- data.frame(SEQN = subset(F0, inAnalysis == 1)$SEQN, post)
A  <- readRDS(file.path(DER, "analysis_frame_primary.rds"))
G4 <- read.csv(file.path(DER, "nhanes_model4_design_frame.csv.gz"))
A <- merge(A, G4[, c("SEQN", "comorb_n_i", "egfr_i", "egfr_m")], by = "SEQN", all.x = TRUE)
A <- merge(A, PD, by = "SEQN", all.x = TRUE)
A$cvd <- as.numeric(A$cvd)
A$lca <- factor(as.character(A$lca),
                levels = c("Low symptom burden", "Insomnia-fatigue", "Hypersomnia-somatic", "High symptom burden"),
                labels = LEVS)
A$W <- as.integer(A$lca)

C3  <- "age+sex+race4+educ+married+pir_i+pir_m+bmi_i+bmi_m+smoke+pa3+htn+dm+cvd+ydx_i+ydx_m+multi_primary+catype"
C4a <- paste(C3, "comorb_n_i+egfr_i+egfr_m", sep = "+")

dat <- subset(A, primary == 1)
dat <- dat[order(dat$SEQN), ]
Dp  <- error_matrix(as.matrix(dat[, LEVS]), dat$W)
cc  <- complete.cases(dat[, c("time", "event", all.vars(as.formula(paste("~", C4a))))])
dat <- dat[cc, ]
sw  <- dat$wt / mean(dat$wt)
A$analysed <- A$SEQN %in% dat$SEQN
des <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt, nest = TRUE, data = A)
P   <- subset(des, primary == 1)
TERMS <- c("Sleep-fatigue vs low", "Somatic-depressive vs low", "High symptom burden vs low",
           "Somatic-depressive vs sleep-fatigue")
theta <- function(e) c(e$coef, SDvsSF = unname(e$coef["clsSomaticDepressive"] - e$coef["clsSleepFatigue"]))
hr <- function(b, se) sprintf("%.2f (%.2f-%.2f)", exp(b), exp(b - 1.96 * se), exp(b + 1.96 * se))

if (stage == "fit") {
  cat("analysed:", nrow(dat), "| deaths:", sum(dat$event), "\n")
  ## modal assignment, design-based (Efron ties, as for Models 3 and 4 in the manuscript)
  sv <- svycoxph(as.formula(paste("Surv(time, event) ~ lca +", C4a)), design = P)
  k <- c("lcaSleepFatigue", "lcaSomaticDepressive", "lcaHigh")
  b <- c(coef(sv)[k], coef(svycontrast(sv, c(lcaSomaticDepressive = 1, lcaSleepFatigue = -1))))
  se <- c(sqrt(diag(vcov(sv)))[k], SE(svycontrast(sv, c(lcaSomaticDepressive = 1, lcaSleepFatigue = -1))))
  modal <- data.frame(model = "Model 4a", assignment = "modal", term = TERMS, logHR = b, se = se,
                      HR = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se),
                      p = 2 * pnorm(-abs(b / se)), n = sv$n, deaths = sv$nevent, terms = length(coef(sv)),
                      row.names = NULL)
  print(transform(modal, HR95 = hr(logHR, se))[, c("term", "HR95", "p")], row.names = FALSE)
  cat("rows used:", sv$n, "| events:", sv$nevent, "| terms:", length(coef(sv)), "\n")
  write.csv(modal, file.path(OUT, "RR6_model4a_modal.csv"), row.names = FALSE)

  ## check: with D = identity the EM reproduces the modal (Breslow, case-weighted) estimate
  e0 <- ml3step(dat, C4a, diag(4), sw, LEVS)
  br <- coxph(as.formula(paste("Surv(time, event) ~ lca +", C4a)), data = dat, weights = sw, ties = "breslow")
  d1 <- max(abs(e0$coef - coef(br)[grep("^lca", names(coef(br)))]))
  cat(sprintf("CHECK (D = I reproduces modal): max |difference in log HR| = %.2e -> %s\n",
              d1, ifelse(d1 < 1e-6, "PASS", "FAIL")))

  ## bias-adjusted fit, with the pre-declared ridge fallback
  e <- ml3step(dat, C4a, Dp, sw, LEVS)
  dec <- 0
  if (e$max_abs_mn > 15 || !e$converged) { e <- ml3step(dat, C4a, Dp, sw, LEVS, decay = 0.01); dec <- 0.01 }
  cat(sprintf("ML three-step: iterations %d, converged %s, max|membership coef| %.1f, decay %g\n",
              e$iter, e$converged, e$max_abs_mn, dec))
  print(round(exp(theta(e)), 3))
  alt <- sapply(1:5, function(r) {                       # local-maximum check, as for Models 3 and 4
    set.seed(100 + r)
    mix <- runif(1, 0.2, 0.8)
    pp <- as.matrix(dat[, LEVS]); pp <- mix * pp + (1 - mix) * matrix(rgamma(length(pp), 1), nrow(pp))
    ee <- ml3step(dat, C4a, Dp, sw, LEVS, init_post = pp, decay = dec)
    c(ee$loglik, ee$coef["clsSomaticDepressive"])
  })
  cat(sprintf("main loglik %.4f | perturbed starts loglik %s | SD log HR %s\n", e$loglik,
              paste(sprintf("%.4f", range(alt[1, ])), collapse = " to "),
              paste(sprintf("%.4f", range(alt[2, ])), collapse = " to ")))
  saveRDS(list(fit = e, decay = dec, check_identity = d1, perturbed = alt), file.path(OUT, "RR6_model4a_fit.rds"))
}

## ---------------------------------------------------------------- replicate standard errors
## Identical to the "rep"/"boot" stages of analysis_three_step.R, for Model 4a only.
if (stage %in% c("rep", "boot")) {
  res <- readRDS(file.path(OUT, "RR6_model4a_fit.rds"))
  ITEMS <- c(paste0("phqi", 1:9), "sleepcat", "slq050b")
  FLCA <- as.formula(paste0("cbind(", paste(ITEMS, collapse = ","), ") ~ 1"))
  LD <- subset(F0, inAnalysis == 1)
  main_probs <- lapply(f4$probs, function(p) p[c(2, 3, 4, 1), , drop = FALSE])
  perms <- as.matrix(expand.grid(1:4, 1:4, 1:4, 1:4))
  perms <- perms[apply(perms, 1, function(r) length(unique(r)) == 4), ]
  match_perm <- function(fit) {
    dist <- apply(perms, 1, function(p) sum(sapply(seq_along(main_probs), function(j)
      sum(abs(fit$probs[[j]][p, ] - main_probs[[j]])))))
    p <- perms[which.min(dist), ]
    list(perm = p, maxdiff = max(sapply(seq_along(main_probs), function(j)
      max(abs(fit$probs[[j]][p, ] - main_probs[[j]])))))
  }
  if (stage == "boot") {
    set.seed(20260919)                                   # same replicate weights as the Models 3/4 bootstrap
    rdes <- as.svrepdesign(des, type = "subbootstrap", replicates = as.integer(Sys.getenv("NBOOT", "500")), mse = TRUE)
  } else rdes <- as.svrepdesign(des, type = "JKn", mse = TRUE)
  SUF <- if (stage == "boot") "_bootstrap" else ""
  rsub <- subset(rdes, analysed)
  RW <- weights(rsub, "analysis")[match(dat$SEQN, rsub$variables$SEQN), , drop = FALSE]
  rall <- subset(rdes, inAnalysis == 1)
  il <- match(LD$SEQN, rall$variables$SEQN)
  FAC <- weights(rall, "analysis")[il, , drop = FALSE] / weights(rall, "sampling")[il]
  stopifnot(!anyNA(RW), !anyNA(FAC), ncol(RW) == ncol(FAC), all(abs(2 * FAC - round(2 * FAC)) < 1e-8))
  YL <- as.matrix(LD[, ITEMS])
  prim <- LD$SEQN %in% A$SEQN[A$primary == 1]
  idx <- match(dat$SEQN, LD$SEQN); stopifnot(!anyNA(idx))
  full <- res$fit
  one <- function(r) {
    f <- FAC[, r]
    fr <- poLCA(FLCA, data = LD[rep(seq_len(nrow(LD)), round(2 * f)), ], nclass = 4, maxiter = 8000, nrep = 1,
                probs.start = f4$probs, verbose = FALSE, calc.se = FALSE)
    mp <- match_perm(fr)
    pr <- poLCA.posterior(fr, YL)[, mp$perm, drop = FALSE]
    Wl <- apply(pr, 1, which.max)
    Dp_r <- error_matrix(pr[prim, , drop = FALSE], Wl[prim], f[prim])
    d_r <- dat; d_r$W <- Wl[idx]
    w <- RW[, r]; w <- w / mean(w[w > 0])
    e <- ml3step(d_r, C4a, Dp_r, w, LEVS, start = list(beta = full$beta, wts = full$wts), init_post = full$post,
                 decay = res$decay, maxit = 3000, tol = 1e-9, ctol = 1e-6, cpat = 20)
    list(theta = theta(e), diag = c(lca_maxdiff = mp$maxdiff, changed_W = sum((d_r$W != dat$W)[w > 0]),
                                    iter = e$iter, stop_rule = match(e$converged_by, c("loglik", "coefficients", "maxit")),
                                    max_abs_mn = e$max_abs_mn))
  }
  t0 <- Sys.time()
  RES <- mclapply(seq_len(ncol(RW)), one, mc.cores = NCORES)
  fail <- which(!sapply(RES, is.list))
  if (length(fail)) { print(RES[[fail[1]]]); stop(length(fail), " replicate(s) failed") }
  TH <- do.call(rbind, lapply(RES, `[[`, "theta")); DG <- do.call(rbind, lapply(RES, `[[`, "diag"))
  saveRDS(list(theta = TH, diag = DG), file.path(OUT, paste0("RR6_model4a_replicates", SUF, ".rds")))
  est <- theta(full)
  V <- svrVar(TH, scale = rsub$scale, rscales = rsub$rscales, mse = rsub$mse, coef = est)
  se <- sqrt(diag(as.matrix(V)))
  tab <- data.frame(model = "Model 4a", assignment = paste0("three-step (", if (stage == "boot") "bootstrap" else "JKn", ")"),
                    term = TERMS, logHR = est, se = se, HR = exp(est), lo = exp(est - 1.96 * se),
                    hi = exp(est + 1.96 * se), p = 2 * pnorm(-abs(est / se)), n = nrow(dat), deaths = sum(dat$event),
                    replicates = nrow(TH), not_converged = sum(DG[, "stop_rule"] == 3),
                    median_iterations = median(DG[, "iter"]), max_membership_coef = max(DG[, "max_abs_mn"]),
                    max_item_prob_difference = max(DG[, "lca_maxdiff"]), row.names = NULL)
  cat(sprintf("%d replicates in %.0f s\n", nrow(TH), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  print(transform(tab, HR95 = hr(logHR, se))[, c("term", "HR95", "p", "not_converged")], row.names = FALSE)
  write.csv(tab, file.path(OUT, paste0("RR6_model4a_three_step", SUF, ".csv")), row.names = FALSE)
}
