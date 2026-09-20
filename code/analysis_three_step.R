## Bias-adjusted (ML three-step) phenotype-mortality estimates, with the validation checks
## and design-based (JKn replicate-weight) standard errors set out in
## review/step2_analysis_plan.md (written before any of these results were seen).
##
## Run from the project root, in stages:
##   Rscript code/analysis_three_step.R fit   # point estimates + checks 1 and local-maximum check
##   Rscript code/analysis_three_step.R sim   # check 3: simulation with a known hazard ratio
##   Rscript code/analysis_three_step.R rep       # JKn standard errors re-running steps 1-3 in every
##                                                 # replicate -> supporting/RR5_three_step_estimates.csv
##   Rscript code/analysis_three_step.R boot      # the same with 500 Rao-Wu bootstrap replicates, a check
##                                                 # for the non-smooth modal step -> *_bootstrap.csv
##   Rscript code/analysis_three_step.R rep3only  # check 2 + the earlier JKn that re-ran step 3 only
##                                                 # (kept for comparison) -> *_step3only.csv
source("code/three_step_functions.R")
suppressPackageStartupMessages({library(survey); library(poLCA); library(parallel)})
options(survey.lonely.psu = "adjust")
set.seed(20260918)
DER <- "data/derived"; OUT <- "supporting"
stage <- commandArgs(TRUE)[1]; if (is.na(stage)) stage <- "fit"
NCORES <- as.integer(Sys.getenv("NCORES", max(1, detectCores() - 1)))

## ---------------------------------------------------------------- data
F0 <- readRDS(file.path(DER, "analysis_frame.rds"))          # row order used to fit the LCA
f4 <- readRDS(file.path(DER, "lca_fits_cancer.rds"))[[4]]
stopifnot(isTRUE(all.equal(apply(f4$posterior, 1, max), subset(F0, inAnalysis == 1)$maxpost,
                           check.attributes = FALSE)))
## poLCA class order is High / Low / Insomnia(=Sleep)-fatigue / Hypersomnia(=Somatic)-somatic
LEVS <- c("Low", "SleepFatigue", "SomaticDepressive", "High")
post <- f4$posterior[, c(2, 3, 4, 1)]; colnames(post) <- LEVS
PD <- data.frame(SEQN = subset(F0, inAnalysis == 1)$SEQN, post)

A  <- readRDS(file.path(DER, "analysis_frame_primary.rds"))
G4 <- read.csv(file.path(DER, "nhanes_model4_design_frame.csv.gz"))
M4KEEP <- c("SEQN", "comorb_n_i", "LBXHGB_i", "LBXHGB_m", "LBXSAL_i", "LBXSAL_m",
            "egfr_i", "egfr_m", "func_lim_i", "n_rx_i", "n_rx_m", "antidep_i")
A <- merge(A, G4[, M4KEEP], by = "SEQN", all.x = TRUE)
A <- merge(A, PD, by = "SEQN", all.x = TRUE)
A$cvd <- as.numeric(A$cvd)
A$bio_m <- pmax(A$LBXSAL_m, A$egfr_m)
A$lca <- factor(as.character(A$lca),
                levels = c("Low symptom burden", "Insomnia-fatigue", "Hypersomnia-somatic", "High symptom burden"),
                labels = LEVS)
A$W <- as.integer(A$lca)

C3 <- "age+sex+race4+educ+married+pir_i+pir_m+bmi_i+bmi_m+smoke+pa3+htn+dm+cvd+ydx_i+ydx_m+multi_primary+catype"
C4 <- paste(C3, "comorb_n_i+LBXHGB_i+LBXHGB_m+LBXSAL_i+egfr_i+bio_m+func_lim_i+n_rx_i+n_rx_m+antidep_i", sep = "+")

dat <- subset(A, primary == 1)
dat <- dat[order(dat$SEQN), ]
Dp  <- error_matrix(as.matrix(dat[, LEVS]), dat$W)                    # step 2, primary cohort (2,569)
## Cox models drop the 5 participants with missing hypertension or diabetes status
## (svycoxph uses na.omit); the three-step models use the same 2,564 participants.
cc  <- complete.cases(dat[, c("time", "event", all.vars(as.formula(paste("~", C4))))])
dat <- dat[cc, ]
sw  <- dat$wt / mean(dat$wt)
A$analysed <- A$SEQN %in% dat$SEQN
des <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt, nest = TRUE, data = A)
P   <- subset(des, primary == 1)
Df  <- error_matrix(as.matrix(A[A$inAnalysis == 1, LEVS]), A$W[A$inAnalysis == 1])  # full cohort (sensitivity)
dimnames(Dp) <- dimnames(Df) <- list(true = LEVS, assigned = LEVS)

hr <- function(b, se) sprintf("%.2f (%.2f-%.2f)", exp(b), exp(b - 1.96 * se), exp(b + 1.96 * se))

if (stage == "fit") {
  cat("primary cohort:", nrow(dat), "| deaths:", sum(dat$event), "\n")
  cat("\nClassification-error matrix D (primary cohort; rows = true class):\n"); print(round(Dp, 3))
  cat("\nFull-cohort D:\n"); print(round(Df, 3))
  saveRDS(list(Dp = Dp, Df = Df), file.path(OUT, "RR5_error_matrices.rds"))

  ## modal assignment: survey-weighted (Efron, as in the manuscript) and Breslow point estimate
  for (m in c("C3", "C4")) {
    sv <- svycoxph(as.formula(paste("Surv(time, event) ~ lca +", get(m))), design = P)
    br <- coxph(as.formula(paste("Surv(time, event) ~ lca +", get(m))), data = dat, weights = sw, ties = "breslow")
    cat(sprintf("\n%s modal  | svycoxph (Efron) SD %.3f  SF %.3f  High %.3f | Breslow SD %.3f\n", m,
                exp(coef(sv)["lcaSomaticDepressive"]), exp(coef(sv)["lcaSleepFatigue"]),
                exp(coef(sv)["lcaHigh"]), exp(coef(br)["lcaSomaticDepressive"])))
  }

  ## check 1: with D = identity the EM must return the modal (Breslow) estimates
  e0 <- ml3step(dat, C3, diag(4), sw, LEVS)
  br3 <- coxph(as.formula(paste("Surv(time, event) ~ lca +", C3)), data = dat, weights = sw, ties = "breslow")
  d1 <- max(abs(e0$coef - coef(br3)[grep("^lca", names(coef(br3)))]))
  cat(sprintf("\nCHECK 1 (D = I reproduces modal): max |difference in log HR| = %.2e -> %s\n",
              d1, ifelse(d1 < 1e-6, "PASS", "FAIL")))

  ## bias-adjusted fits
  res <- list()
  for (m in c("C3", "C4")) {
    t0 <- Sys.time()
    e <- ml3step(dat, get(m), Dp, sw, LEVS, verbose = FALSE)
    dec <- 0
    if (e$max_abs_mn > 15 || !e$converged) {             # pre-declared fallback
      e <- ml3step(dat, get(m), Dp, sw, LEVS, decay = 0.01); dec <- 0.01
    }
    cat(sprintf("\n%s ML three-step: iterations %d, converged %s, max|membership coef| %.1f, decay %g, %.0f s\n",
                m, e$iter, e$converged, e$max_abs_mn, dec, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    print(round(exp(e$coef), 3))
    cat(sprintf("  contrast SD vs SF: %.3f\n", exp(e$coef["clsSomaticDepressive"] - e$coef["clsSleepFatigue"])))
    ## local-maximum check: 5 perturbed starts (start the EM from mixed posteriors)
    alt <- sapply(1:5, function(r) {
      set.seed(100 + r)
      mix <- runif(1, 0.2, 0.8)                            # mix LCA posteriors with random noise
      pp <- as.matrix(dat[, LEVS]); pp <- mix * pp + (1 - mix) * matrix(rgamma(length(pp), 1), nrow(pp))
      ee <- ml3step(dat, get(m), Dp, sw, LEVS, init_post = pp, decay = dec)
      c(ee$loglik, ee$coef["clsSomaticDepressive"])
    })
    cat(sprintf("  main fit loglik %.4f\n", e$loglik))
    cat("  perturbed starts: loglik range", sprintf("%.4f", range(alt[1, ])), "| SD log HR range",
        sprintf("%.4f", range(alt[2, ])), "\n")
    res[[m]] <- list(fit = e, decay = dec)
  }
  ## sensitivity: full-cohort D
  ef <- ml3step(dat, C3, Df, sw, LEVS)
  cat("\nC3 with full-cohort D:", round(exp(ef$coef), 3), "\n")
  res$C3_fullD <- list(fit = ef, decay = 0)
  saveRDS(res, file.path(OUT, "RR5_three_step_fits.rds"))
}

## ---------------------------------------------------------------- check 3: simulation
## Data generated from the fitted model structure with KNOWN hazard ratios; misclassification
## drawn from D. The ML three-step estimator should be approximately unbiased and the modal
## estimator attenuated. Case weights are set to 1 (this checks the estimator, not the design).
if (stage == "sim") {
  res <- readRDS(file.path(OUT, "RR5_three_step_fits.rds"))
  NREP <- 60; TRUE_HR <- c(SleepFatigue = 1.0, SomaticDepressive = 2.0, High = 1.6)
  ## administrative censoring: observed time if censored; for deaths, the median censoring
  ## time of censored participants in the same survey cycle
  cmed <- tapply(dat$time[dat$event == 0], dat$cycle[dat$event == 0], median)
  cens <- ifelse(dat$event == 0, dat$time, cmed[as.character(dat$cycle)])
  one <- function(r, covs, fit, lam, k = 1.2) {
    set.seed(7000 + r)
    n <- nrow(dat)
    pz <- predict(fit$mn, newdata = dat, type = "probs")[, LEVS]
    X <- apply(pz, 1, function(p) sample.int(4, 1, prob = p))
    W <- vapply(X, function(x) sample.int(4, 1, prob = Dp[x, ]), 1L)
    b <- fit$beta; b[paste0("cls", names(TRUE_HR))] <- log(TRUE_HR)
    d1 <- dat; d1$cls <- factor(LEVS[X], levels = LEVS)
    MMx <- model.matrix(as.formula(paste("~ cls +", covs)), d1)[, -1, drop = FALSE]
    eta <- as.vector(MMx[, names(b)] %*% b)
    Tt <- (-log(runif(n)) / (lam * exp(eta)))^(1 / k)
    d2 <- dat; d2$time <- pmin(Tt, cens); d2$event <- as.integer(Tt <= cens)
    d2$W <- W; d2$lca <- factor(LEVS[W], levels = LEVS)
    mo <- coef(coxph(as.formula(paste("Surv(time, event) ~ lca +", covs)), data = d2, ties = "breslow"))
    ml <- ml3step(d2, covs, Dp, rep(1, n), LEVS, maxit = 500, tol = 1e-8)
    c(modal = mo[paste0("lca", names(TRUE_HR))], ml3 = ml$coef[paste0("cls", names(TRUE_HR))],
      deaths = mean(d2$event), conv = ml$converged)
  }
  out <- list()
  for (m in c("C3", "C4")) {
    fit <- res[[m]]$fit
    ## choose the Weibull scale so that about 26% die, as observed
    pz <- predict(fit$mn, newdata = dat, type = "probs")[, LEVS]
    b <- fit$beta; b[paste0("cls", names(TRUE_HR))] <- log(TRUE_HR)
    etas <- sapply(LEVS, function(l) { d1 <- dat; d1$cls <- factor(l, levels = LEVS)
      as.vector(model.matrix(as.formula(paste("~ cls +", get(m))), d1)[, -1][, names(b)] %*% b) })
    f <- function(ll) mean(rowSums(pz * (1 - exp(-exp(ll) * exp(etas) * cens^1.2)))) - mean(dat$event)
    lam <- exp(uniroot(f, c(-40, 10))$root)
    t0 <- Sys.time()
    sims <- do.call(rbind, mclapply(seq_len(NREP), one, covs = get(m), fit = fit, lam = lam,
                                    mc.cores = NCORES))
    cat(sprintf("\n%s simulation: %d replicates, %.0f s, mean deaths %.3f, EM converged %d/%d\n", m, NREP,
                as.numeric(difftime(Sys.time(), t0, units = "secs")), mean(sims[, "deaths"]),
                sum(sims[, "conv"]), NREP))
    tab <- data.frame(phenotype = names(TRUE_HR), true_HR = TRUE_HR,
                      modal_HR = exp(colMeans(sims[, paste0("modal.lca", names(TRUE_HR))])),
                      ml3_HR = exp(colMeans(sims[, paste0("ml3.cls", names(TRUE_HR))])),
                      modal_bias_logHR = colMeans(sims[, paste0("modal.lca", names(TRUE_HR))]) - log(TRUE_HR),
                      ml3_bias_logHR = colMeans(sims[, paste0("ml3.cls", names(TRUE_HR))]) - log(TRUE_HR),
                      ml3_empirical_SE = apply(sims[, paste0("ml3.cls", names(TRUE_HR))], 2, sd),
                      mcse_ml3_bias = apply(sims[, paste0("ml3.cls", names(TRUE_HR))], 2, sd) / sqrt(NREP),
                      row.names = NULL)
    print(format(tab, digits = 3))
    out[[m]] <- cbind(covariates = ifelse(m == "C3", "Model 3", "Model 4"), tab)
  }
  write.csv(do.call(rbind, out), file.path(OUT, "RR5_simulation_check.csv"), row.names = FALSE)
}

## ---------------------------------------------------------------- check 2 + JKn standard errors
if (stage == "rep3only") {
  res <- readRDS(file.path(OUT, "RR5_three_step_fits.rds"))
  rdes <- as.svrepdesign(des, type = "JKn", mse = TRUE)
  rsub <- subset(rdes, analysed)
  RW <- weights(rsub, "analysis"); seqn <- rsub$variables$SEQN
  RW <- RW[match(dat$SEQN, seqn), , drop = FALSE]
  cat("JKn replicates:", ncol(RW), "| analysed:", nrow(RW), "\n")

  ## check 2: modal model - standard errors from the same replicate machinery used below
  ## (weighted coxph refitted per JKn replicate) versus Taylor linearisation (svycoxph)
  Ptl <- subset(des, analysed)
  chk <- list()
  for (m in c("C3", "C4")) {
    fm <- as.formula(paste("Surv(time, event) ~ lca +", get(m)))
    tl <- svycoxph(fm, design = Ptl)
    k <- grep("^lca", names(coef(tl)))
    th <- do.call(rbind, mclapply(seq_len(ncol(RW)), function(r) {
      d <- dat; d$.w <- RW[, r]; d <- d[d$.w > 0, ]        # weights inside data (model.frame scoping)
      coef(coxph(fm, data = d, weights = .w))[names(coef(tl))[k]]
    }, mc.cores = NCORES))
    se_jk <- sqrt(diag(as.matrix(svrVar(th, scale = rsub$scale, rscales = rsub$rscales, mse = rsub$mse,
                                        coef = coef(tl)[k]))))
    se_tl <- sqrt(diag(vcov(tl)))[k]
    cat(sprintf("CHECK 2 %s modal SE (log HR) Taylor: %s | JKn: %s | max relative difference %.1f%%\n", m,
                paste(sprintf("%.3f", se_tl), collapse = "/"), paste(sprintf("%.3f", se_jk), collapse = "/"),
                100 * max(abs(se_jk / se_tl - 1))))
    chk[[m]] <- data.frame(model = m, term = names(se_tl), se_taylor = se_tl, se_jkn = se_jk)
  }
  write.csv(do.call(rbind, chk), file.path(OUT, "RR5_check2_se_comparison.csv"), row.names = FALSE)

  theta <- function(e) c(e$coef, SDvsSF = unname(e$coef["clsSomaticDepressive"] - e$coef["clsSleepFatigue"]))
  rows <- list()
  for (m in c("C3", "C4", "C3_fullD")) {
    full <- res[[m]]$fit; dec <- res[[m]]$decay
    covs <- if (m == "C4") C4 else C3
    Dm <- if (m == "C3_fullD") Df else Dp
    t0 <- Sys.time()
    reps <- mclapply(seq_len(ncol(RW)), function(r) {
      w <- RW[, r]; w <- w / mean(w[w > 0])
      e <- ml3step(dat, covs, Dm, w, LEVS, start = list(beta = full$beta, wts = full$wts),
                   init_post = full$post, decay = dec, maxit = 2000, tol = 1e-8)
      c(theta(e), conv = e$converged, iter = e$iter)
    }, mc.cores = NCORES)
    TH <- do.call(rbind, reps)
    V <- svrVar(TH[, 1:4], scale = rsub$scale, rscales = rsub$rscales, mse = rsub$mse, coef = theta(full))
    se <- sqrt(diag(as.matrix(V)))
    est <- theta(full)
    saveRDS(TH, file.path(OUT, sprintf("RR5_replicates3only_%s.rds", m)))
    cat(sprintf("\n%s: %d replicates in %.0f s; converged (relative change < 1e-8 within 2,000 iterations): %d; median iterations %d\n",
                m, nrow(TH), as.numeric(difftime(Sys.time(), t0, units = "secs")), sum(TH[, "conv"] == 1),
                as.integer(median(TH[, "iter"]))))
    tab <- data.frame(model = m, term = c("Sleep-fatigue vs low", "Somatic-depressive vs low",
                                          "High symptom burden vs low", "Somatic-depressive vs sleep-fatigue"),
                      logHR = est, se = se, HR = exp(est), lo = exp(est - 1.96 * se), hi = exp(est + 1.96 * se),
                      p = 2 * pnorm(-abs(est / se)), row.names = NULL)
    print(transform(tab, HR95 = hr(logHR, se))[, c("term", "HR95", "p")], row.names = FALSE)
    rows[[m]] <- tab
  }
  write.csv(do.call(rbind, rows), file.path(OUT, "RR5_three_step_estimates_step3only.csv"), row.names = FALSE)
}

## ---------------------------------------------------------------- JKn re-running steps 1-3
## In every replicate the measurement model is re-fitted on the retained survivors (unweighted, as
## in the main analysis, starting from the published solution), classes are matched to the main
## solution by their item-response probabilities, modal assignment and both classification-error
## matrices are recomputed, and the step-3 EM is re-run with the replicate weights. The standard
## errors therefore include uncertainty in the measurement model and in the error matrix, which the
## step-3-only jackknife (stage rep3only) treated as known.
if (stage %in% c("rep", "boot")) {
  res <- readRDS(file.path(OUT, "RR5_three_step_fits.rds"))
  ITEMS <- c(paste0("phqi", 1:9), "sleepcat", "slq050b")
  FLCA <- as.formula(paste0("cbind(", paste(ITEMS, collapse = ","), ") ~ 1"))
  LD <- subset(F0, inAnalysis == 1)                                   # LCA sample, original row order
  main_probs <- lapply(f4$probs, function(p) p[c(2, 3, 4, 1), , drop = FALSE])   # in LEVS order
  perms <- as.matrix(expand.grid(1:4, 1:4, 1:4, 1:4))
  perms <- perms[apply(perms, 1, function(r) length(unique(r)) == 4), ]
  match_perm <- function(fit) {
    dist <- apply(perms, 1, function(p) sum(sapply(seq_along(main_probs), function(j)
      sum(abs(fit$probs[[j]][p, ] - main_probs[[j]])))))
    p <- perms[which.min(dist), ]
    list(perm = p, maxdiff = max(sapply(seq_along(main_probs), function(j)
      max(abs(fit$probs[[j]][p, ] - main_probs[[j]])))))
  }
  f0 <- poLCA(FLCA, data = LD, nclass = 4, maxiter = 8000, nrep = 1, probs.start = f4$probs,
              verbose = FALSE, calc.se = FALSE)
  cat(sprintf("Refit of the measurement model on the full sample from the published solution: loglik %.3f (published %.3f); max |posterior difference| %.1e\n",
              f0$llik, f4$llik, max(abs(f0$posterior - f4$posterior))))

  ## The measurement model is fitted without survey weights, so in each replicate it gets the
  ## replicate factors themselves as frequency weights: 0 for the deleted (or unsampled) PSU,
  ## n_h/(n_h - 1) for the rest of its stratum, 1 elsewhere (JKn); bootstrap factors likewise.
  ## poLCA takes no weights, so rows are repeated round(2 * factor) times (factors are multiples of
  ## 1/2, and doubling every row leaves the estimates unchanged). The error matrix uses the same factors.
  if (stage == "boot") {
    set.seed(20260919)
    rdes <- as.svrepdesign(des, type = "subbootstrap", replicates = as.integer(Sys.getenv("NBOOT", "500")), mse = TRUE)
  } else rdes <- as.svrepdesign(des, type = "JKn", mse = TRUE)
  SUF0 <- if (stage == "boot") "_bootstrap" else ""
  rsub <- subset(rdes, analysed)
  RW <- weights(rsub, "analysis")[match(dat$SEQN, rsub$variables$SEQN), , drop = FALSE]
  rall <- subset(rdes, inAnalysis == 1)
  il <- match(LD$SEQN, rall$variables$SEQN)
  FAC <- weights(rall, "analysis")[il, , drop = FALSE] / weights(rall, "sampling")[il]
  stopifnot(!anyNA(RW), !anyNA(FAC), ncol(RW) == ncol(FAC), all(abs(2 * FAC - round(2 * FAC)) < 1e-8))
  YL <- as.matrix(LD[, ITEMS])
  f2 <- poLCA(FLCA, data = LD[rep(seq_len(nrow(LD)), 2), ], nclass = 4, maxiter = 8000, nrep = 1,
              probs.start = f4$probs, verbose = FALSE, calc.se = FALSE)
  cat(sprintf("Every row doubled: loglik/2 %.3f; max |posterior difference| %.1e (posteriors recomputed for the original rows)\n",
              f2$llik / 2, max(abs(poLCA.posterior(f2, YL) - f4$posterior))))
  prim <- LD$SEQN %in% A$SEQN[A$primary == 1]
  theta <- function(e) c(e$coef, SDvsSF = unname(e$coef["clsSomaticDepressive"] - e$coef["clsSleepFatigue"]))
  MODELS <- c("C3", "C4", "C3_fullD")

  idx <- match(dat$SEQN, LD$SEQN); stopifnot(!anyNA(idx))
  one <- function(r) {
    f <- FAC[, r]
    fr <- poLCA(FLCA, data = LD[rep(seq_len(nrow(LD)), round(2 * f)), ], nclass = 4, maxiter = 8000, nrep = 1,
                probs.start = f4$probs, verbose = FALSE, calc.se = FALSE)
    mp <- match_perm(fr)
    post <- poLCA.posterior(fr, YL)[, mp$perm, drop = FALSE]    # every LCA-sample member, original rows
    Wl <- apply(post, 1, which.max)
    Dp_r <- error_matrix(post[prim, , drop = FALSE], Wl[prim], f[prim])
    Df_r <- error_matrix(post, Wl, f)
    d_r <- dat; d_r$W <- Wl[idx]
    w <- RW[, r]; w <- w / mean(w[w > 0])
    th <- c(); dg <- c(lca_maxdiff = mp$maxdiff, lca_loglik = fr$llik / 2,
                       changed_W = sum((d_r$W != dat$W)[w > 0]), D_SD_own = Dp_r[3, 3])
    for (m in MODELS) {
      full <- res[[m]]$fit
      e <- ml3step(d_r, if (m == "C4") C4 else C3, if (m == "C3_fullD") Df_r else Dp_r, w, LEVS,
                   start = list(beta = full$beta, wts = full$wts), init_post = full$post,
                   decay = res[[m]]$decay, maxit = 3000, tol = 1e-9, ctol = 1e-6, cpat = 20)
      tt <- theta(e); th <- c(th, setNames(tt, paste(m, names(tt), sep = ":")))
      dg <- c(dg, setNames(c(e$iter, match(e$converged_by, c("loglik", "coefficients", "maxit")), e$max_abs_mn),
                           paste(m, c("iter", "stop_rule", "max_abs_mn"), sep = ":")))
    }
    list(theta = th, diag = dg)
  }
  NTEST <- as.integer(Sys.getenv("NREP_TEST", "0"))      # > 0: quick test on the first NTEST replicates
  SUF <- paste0(SUF0, if (NTEST > 0) "_TEST" else "")
  t0 <- Sys.time()
  RES <- mclapply(if (NTEST > 0) seq_len(NTEST) else seq_len(ncol(RW)), one, mc.cores = NCORES)
  fail <- which(!sapply(RES, is.list))
  if (length(fail)) { print(RES[[fail[1]]]); stop(length(fail), " replicate(s) failed") }
  TH <- do.call(rbind, lapply(RES, `[[`, "theta")); DG <- do.call(rbind, lapply(RES, `[[`, "diag"))
  saveRDS(list(theta = TH, diag = DG), file.path(OUT, paste0("RR5_jackknife_full", SUF, ".rds")))
  if (NTEST > 0) { print(round(TH, 4)); print(DG); quit(save = "no") }
  cat(sprintf("%d replicates in %.0f s\n", nrow(TH), as.numeric(difftime(Sys.time(), t0, units = "secs"))))

  rows <- list(); dgs <- list()
  for (m in MODELS) {
    full <- res[[m]]$fit; est <- theta(full)
    cols <- paste(m, names(est), sep = ":")
    V <- svrVar(TH[, cols], scale = rsub$scale, rscales = rsub$rscales, mse = rsub$mse, coef = est)
    se <- sqrt(diag(as.matrix(V)))
    tab <- data.frame(model = m, term = c("Sleep-fatigue vs low", "Somatic-depressive vs low",
                                          "High symptom burden vs low", "Somatic-depressive vs sleep-fatigue"),
                      logHR = est, se = se, HR = exp(est), lo = exp(est - 1.96 * se), hi = exp(est + 1.96 * se),
                      p = 2 * pnorm(-abs(est / se)), n = nrow(dat), deaths = sum(dat$event), row.names = NULL)
    cat(sprintf("\n%s\n", m)); print(transform(tab, HR95 = hr(logHR, se))[, c("term", "HR95", "p")], row.names = FALSE)
    rule <- DG[, paste0(m, ":stop_rule")]
    dgs[[m]] <- data.frame(model = m, replicates = nrow(DG), stopped_loglik = sum(rule == 1),
                           stopped_coefficients = sum(rule == 2), not_converged = sum(rule == 3),
                           median_iterations = median(DG[, paste0(m, ":iter")]),
                           max_membership_coef = max(DG[, paste0(m, ":max_abs_mn")]))
    rows[[m]] <- tab
  }
  write.csv(do.call(rbind, rows), file.path(OUT, paste0("RR5_three_step_estimates", SUF, ".csv")), row.names = FALSE)
  lca <- data.frame(model = "measurement model", replicates = nrow(DG),
                    max_item_prob_difference = max(DG[, "lca_maxdiff"]),
                    median_changed_assignments = median(DG[, "changed_W"]),
                    max_changed_assignments = max(DG[, "changed_W"]),
                    somatic_depressive_correctly_assigned_range = sprintf("%.3f to %.3f", min(DG[, "D_SD_own"]), max(DG[, "D_SD_own"])))
  write.csv(do.call(rbind, dgs), file.path(OUT, paste0("RR5_jackknife_diagnostics", SUF, ".csv")), row.names = FALSE)
  write.csv(lca, file.path(OUT, paste0("RR5_jackknife_lca_diagnostics", SUF, ".csv")), row.names = FALSE)
  print(do.call(rbind, dgs)); print(lca)
}
