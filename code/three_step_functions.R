## Functions for the bias-adjusted (ML three-step) latent class - survival analysis,
## bivariate residuals and inclusive pseudo-class draws.
## Plan written before any result was seen: review/step2_analysis_plan.md
##
## ML three-step (Vermunt 2010; Bakk, Tekle and Vermunt 2013) with a Cox step-3 model:
##   step 1  measurement model (symptoms only) -> posteriors p_it
##   step 2  modal class W and classification-error matrix D[t, s] = P(W = s | X = t)
##   step 3  pseudo-likelihood over latent X for (W, T, delta):
##           P(X = t | Z) multinomial logistic, h(t | X, Z) = h0(t) exp(beta_X + gamma'Z),
##           P(W | X) fixed at D; fitted by EM with a Breslow baseline hazard.
suppressPackageStartupMessages({library(survival); library(nnet)})

## classification-error matrix from posteriors (rows = true class t, cols = assigned class s)
error_matrix <- function(post, W = apply(post, 1, which.max), w = rep(1, nrow(post))) {
  K <- ncol(post)                                    # w: frequency weights (replicate factors)
  D <- sapply(seq_len(K), function(s) colSums(w[W == s] * post[W == s, , drop = FALSE]))  # K x K: [t, s]
  D / rowSums(D)
}

## Breslow cumulative baseline hazard at each subject's own time (weighted, uncentred lp)
breslow_H0 <- function(time, event, eta, w) {
  o <- order(time)
  t <- time[o]; d <- event[o] * w[o]; r <- w[o] * exp(eta[o])
  ut <- unique(t[d > 0])
  dN <- tapply(d[d > 0], t[d > 0], sum)[as.character(ut)]
  risk <- sapply(ut, function(u) sum(r[t >= u]))
  jump <- dN / risk
  list(times = ut, jump = as.numeric(jump), H = stepfun(ut, c(0, cumsum(jump))))
}

## Build the expanded (subject x class) data once
expand_classes <- function(dat, K, levs) {
  ex <- dat[rep(seq_len(nrow(dat)), each = K), , drop = FALSE]
  ex$cls <- factor(rep(levs, times = nrow(dat)), levels = levs)
  ex$.id <- rep(seq_len(nrow(dat)), each = K)
  ex
}

## ML three-step EM.
##   dat   : one row per subject with time, event, W (1..K in the order of levs), covariates
##   covs  : right-hand side string of covariates (same in the Cox and membership models)
##   D     : K x K classification-error matrix [true, assigned]
##   sw    : case weights (normalised survey weights)
##   start : optional list(beta, wts) for warm starts
##   init_post : optional n x K starting posterior weights (for local-maximum checks)
##   ctol, cpat: optional second stopping rule - stop when no class log hazard ratio changes by
##               more than ctol for cpat consecutive iterations (the pseudo-log-likelihood can keep
##               creeping when a membership coefficient drifts while the hazard ratios are stable)
ml3step <- function(dat, covs, D, sw, levs, start = NULL, maxit = 1000, tol = 1e-9,
                    decay = 0, verbose = FALSE, init_post = NULL, ctol = NULL, cpat = 20) {
  K <- length(levs); n <- nrow(dat)
  stopifnot(!anyNA(dat[, c("time", "event", "W", all.vars(as.formula(paste("~", covs))))]))
  ex <- expand_classes(dat, K, levs)
  fcox <- as.formula(paste("Surv(time, event) ~ cls +", covs))
  fmn  <- as.formula(paste("cls ~", covs))
  MM <- model.matrix(as.formula(paste("~ cls +", covs)), ex)[, -1, drop = FALSE]
  logD <- log(pmax(D[, dat$W, drop = FALSE], 1e-300))          # K x n : log P(W_i | X = t)

  ## initial posterior weights: modal assignment unless another starting point is given
  ## (W itself is data and is never changed)
  w_it <- matrix(0, n, K); w_it[cbind(seq_len(n), dat$W)] <- 1
  if (!is.null(init_post)) w_it <- init_post / rowSums(init_post)
  mn_wts <- if (!is.null(start)) start$wts else NULL
  ll_old <- -Inf; beta <- NULL; cprev <- NULL; stable <- 0; how <- "maxit"
  for (it in seq_len(maxit)) {
    a <- as.vector(t(w_it)) * rep(sw, each = K)                  # expanded weights
    keep <- a > 1e-10
    init0 <- if (!is.null(beta)) beta else if (!is.null(start)) start$beta else NULL
    cf <- if (is.null(init0)) coxph(fcox, data = ex[keep, ], weights = a[keep], ties = "breslow")
          else coxph(fcox, data = ex[keep, ], weights = a[keep], ties = "breslow", init = init0)
    beta <- coef(cf); beta[is.na(beta)] <- 0
    mn <- if (is.null(mn_wts)) multinom(fmn, data = ex[keep, ], weights = a[keep], trace = FALSE,
                                        maxit = 2000, decay = decay)
          else multinom(fmn, data = ex[keep, ], weights = a[keep], trace = FALSE, maxit = 2000,
                        decay = decay, Wts = mn_wts)
    mn_wts <- mn$wts
    pi_it <- predict(mn, newdata = dat, type = "probs")
    if (K == 2) pi_it <- cbind(1 - pi_it, pi_it)
    pi_it <- pi_it[, levs, drop = FALSE]

    eta <- matrix(MM[, names(beta), drop = FALSE] %*% beta, nrow = n, ncol = K, byrow = TRUE)
    ## baseline hazard from the current posterior-weighted expanded data
    bh <- breslow_H0(ex$time[keep], ex$event[keep], as.vector(MM[keep, names(beta), drop = FALSE] %*% beta), a[keep])
    H0 <- bh$H(dat$time)
    dH <- numeric(n); ev <- dat$event == 1
    dH[ev] <- bh$jump[match(dat$time[ev], bh$times)]

    lp <- log(pmax(pi_it, 1e-300)) + t(logD) + dat$event * eta - H0 * exp(eta)
    m <- apply(lp, 1, max)
    lse <- m + log(rowSums(exp(lp - m)))
    pos <- sw > 0                     # zero replicate weights drop out (their death times may have no jump)
    ll <- sum((sw * (lse + ifelse(ev, log(dH), 0)))[pos])
    w_it <- exp(lp - lse)
    if (verbose && it %% 10 == 1) cat(sprintf("  iter %d  pseudo-loglik %.6f\n", it, ll))
    cc <- beta[grep("^cls", names(beta))]
    if (!is.null(ctol) && !is.null(cprev)) stable <- if (max(abs(cc - cprev)) < ctol) stable + 1 else 0
    cprev <- cc
    if (is.finite(ll_old) && abs(ll - ll_old) < tol * abs(ll_old)) { how <- "loglik"; break }
    if (!is.null(ctol) && stable >= cpat) { how <- "coefficients"; break }
    ll_old <- ll
  }
  list(beta = beta, coef = beta[grep("^cls", names(beta))], loglik = ll, iter = it,
       converged = how != "maxit", converged_by = how, post = w_it, wts = mn_wts,
       max_abs_mn = max(abs(coef(mn))), mn = mn)
}

## Bivariate residuals for all item pairs of a poLCA fit (Pearson X^2 / df)
bvr_table <- function(fit) {
  y <- fit$y; J <- ncol(y); N <- nrow(y); P <- fit$P; pr <- fit$probs
  out <- list()
  for (j in 1:(J - 1)) for (k in (j + 1):J) {
    Kj <- ncol(pr[[j]]); Kk <- ncol(pr[[k]])
    O <- table(factor(y[, j], levels = 1:Kj), factor(y[, k], levels = 1:Kk))
    E <- N * Reduce(`+`, lapply(seq_along(P), function(t) P[t] * outer(pr[[j]][t, ], pr[[k]][t, ])))
    X2 <- sum((O - E)^2 / E); df <- (Kj - 1) * (Kk - 1)
    out[[length(out) + 1]] <- data.frame(item1 = names(pr)[j], item2 = names(pr)[k],
                                         X2 = X2, df = df, BVR = X2 / df)
  }
  do.call(rbind, out)
}

## Rubin's rules for M sets of coefficients (B: M x p) and variances (V: M x p)
rubin <- function(B, V) {
  M <- nrow(B); qb <- colMeans(B); ub <- colMeans(V); bv <- apply(B, 2, var)
  tot <- ub + (1 + 1 / M) * bv
  data.frame(term = colnames(B), logHR = qb, se = sqrt(tot), fmi = (1 + 1 / M) * bv / tot)
}
