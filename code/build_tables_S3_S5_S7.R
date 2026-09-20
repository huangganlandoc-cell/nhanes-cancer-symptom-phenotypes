## =============================================================================
## build_tables_S3_S5_S7.R  --  assemble Supplementary Tables S3, S5 and S7 from saved outputs
##
## Run from the project root:
##   /opt/homebrew/bin/Rscript code/build_tables_S3_S5_S7.R
##       reads supporting/ (and data/derived/), writes supplementary/TableS3_proportional_hazards.csv,
##       supplementary/TableS5_model4_detail.csv and supplementary/TableS7_irreducibility_analyses.csv
##   Options:  --in=DIR   folder holding the analysis outputs (a file missing there is read from
##                        supporting/; every file actually used is printed)
##             --out=DIR  folder the three tables are written to (default supplementary)
##             --legacy   reproduce the conventions of the round-11 tables exactly (P values as typed
##                        into the old helper files, the `events` column, E-values as stored, the
##                        two-row sample-size block). Used only to check this script against the
##                        archived inputs in review/archive_round11/supporting_2557_frame/.
##
## Layout, block names, row labels and display conventions follow the existing tables:
##   HR (95% CI) "%.2f (%.2f-%.2f)"; P "%.4f", or "%.1e" below 1e-4; other numbers as the original
##   Python builders printed them (steps 172, 265, 285 and 289 of the original analysis log);
##   phenotype and term labels passed through the replacements of code/fix_display_labels.py.
##
## Default (non-legacy) mode also applies these corrections:
##   S3  block A: Schoenfeld P values recomputed from the stored chi-square statistic and df
##                (the stored P values were rounded to 4 decimals);
##       block B: "deaths in period" = deaths actually analysed in that period (2,564 survivors with
##                complete Model 3 covariates, from data/derived/analysis_frame_primary.rds);
##   S5  block B: HR, CI and P from the unrounded model outputs (RR_model4.csv, all-cause;
##                HM_model4_cause_specific.csv, cancer and cardiovascular), BH and Bonferroni
##                recomputed from those P values (the old helper file held 3-digit P values);
##       block F: E-value for the confidence limit = 1 whenever the interval includes 1;
##   S7  block C: only the specifications of the existing layout (A1 to C2) are shown;
##       block D: HR, CI and P from the unrounded cause-specific fits (RR_dimensional_affective.csv),
##                BH and Bonferroni recomputed from them;
##       block F: first two rows = bias-adjusted three-step Model 4 estimate with its jackknife and its
##                bootstrap SE (RR5_three_step_estimates.csv and RR5_three_step_estimates_bootstrap.csv,
##                model C4, "Somatic-depressive vs low"),
##                multiplier = (se / (logHR / (qnorm(0.975) + qnorm(0.8))))^2; then the pseudo-class
##                and modal rows of RR_power_confirmatory.csv, relabelled. All three rows use the same
##                "current sample" (2,569 survivors, 676 deaths) as the existing rows.
## =============================================================================
suppressPackageStartupMessages(library(survey))
if (!dir.exists("supporting") || !dir.exists("data/derived")) stop("run from the project root")
invisible(suppressWarnings(Sys.setlocale("LC_CTYPE", "en_US.UTF-8")))   # the "\u00d7" label must stay UTF-8
args <- commandArgs(TRUE)
opt <- function(name, default) {
  a <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(a)) sub(paste0("^--", name, "="), "", a[1]) else default }
IN <- opt("in", "supporting"); OUTD <- opt("out", "supplementary"); LEGACY <- "--legacy" %in% args
dir.create(OUTD, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("mode: %s | inputs: %s | output folder: %s\n", if (LEGACY) "legacy" else "corrected", IN, OUTD))

## ---- input and formatting helpers ------------------------------------------------------------
used <- character(0)
src <- function(f) {
  p <- file.path(IN, f); if (!file.exists(p)) p <- file.path("supporting", f)
  if (!file.exists(p)) stop("input not found: ", f)
  used[f] <<- p; p }
rd <- function(f, ...) read.csv(src(f), stringsAsFactors = FALSE, check.names = FALSE, ...)
fh <- function(HR, lo, hi) sprintf("%.2f (%.2f-%.2f)", HR, lo, hi)
pf <- function(v) ifelse(v >= 1e-4, sprintf("%.4f", v), sprintf("%.1e", v))
py_round <- function(x, d) round(x * 10^d) / 10^d            # numpy round (half to even on x*10^d)
py_str <- function(x) vapply(x, function(v) {                # Python str() of a float
  if (is.na(v)) return("nan"); if (v == 0) return("0.0")
  for (d in 1:17) { s <- sprintf(paste0("%.", d - 1, "e"), v); if (as.numeric(s) == v) break }
  ex <- as.integer(sub(".*e", "", s))
  if (ex >= -4 && ex < 16) sprintf(paste0("%.", max(d - 1 - ex, 1), "f"), v)
  else { m <- sub("e.*", "", s); m <- sub("0+$", "", sub("\\.$", "", m)); sprintf("%se%s%02d", m, ifelse(ex < 0, "-", "+"), abs(ex)) }
}, "")
fmt_int <- function(n) formatC(n, format = "d", big.mark = ",")
REPL <- list(c("lcaInsomnia-fatigue", "Sleep-fatigue"), c("lcaHypersomnia-somatic", "Somatic-depressive"),
             c("lcaHigh symptom burden", "High symptom burden"),
             c("Insomnia-fatigue", "Sleep-fatigue"), c("insomnia-fatigue", "sleep-fatigue"),
             c("Hypersomnia-somatic", "Somatic-depressive"), c("hypersomnia-somatic", "somatic-depressive"),
             c("lca Schoenfeld residuals", "Phenotype Schoenfeld residuals"),
             c("lca x follow-up period interaction", "Phenotype \u00d7 follow-up period interaction"),
             c("Insomnia complaint", "Previously reported sleep trouble"))       # = code/fix_display_labels.py
relabel <- function(s) { for (r in REPL) s <- gsub(r[1], r[2], s, fixed = TRUE); s }
csv_field <- function(s) { s <- as.character(s); s[is.na(s)] <- ""
  q <- grepl('[,"\n\r]', s); s[q] <- paste0('"', gsub('"', '""', s[q]), '"'); s }
write_table <- function(df, name) {                          # pandas to_csv(index=False) style, LF, UTF-8
  df[] <- lapply(df, function(x) relabel(as.character(x)))
  lines <- c(paste(csv_field(names(df)), collapse = ","), do.call(paste, c(lapply(df, csv_field), sep = ",")))
  con <- file(file.path(OUTD, name), "wb"); writeLines(enc2utf8(lines), con, sep = "\n", useBytes = TRUE); close(con)
  cat(sprintf("  wrote %s (%d rows)\n", file.path(OUTD, name), nrow(df))) }
blk <- function(Block, Row, Estimate, Detail, P) data.frame(Block = Block, Row = Row, Estimate = Estimate,
                                                          Detail = Detail, P = P, check.names = FALSE)

## ---- analysed sample (data/derived) ------------------------------------------------------------
AF <- readRDS("data/derived/analysis_frame_primary.rds")
C3v <- c("age","sex","race4","educ","married","pir_i","pir_m","bmi_i","bmi_m","smoke","pa3","htn","dm",
         "cvd","ydx_i","ydx_m","multi_primary","catype")
an <- AF$primary == 1 & complete.cases(AF[, c("time", "event", C3v)])

## ============================================ Table S3 ==========================================
t1 <- rd("ph_tests.csv"); t2 <- rd("ph_period_specific.csv"); t3 <- rd("ph_stratified_baseline.csv")
pA <- t1$p
if (!LEGACY) { k <- !is.na(t1$statistic); pA[k] <- pchisq(t1$statistic[k], t1$df[k], lower.tail = FALSE) }
deaths_note <- if (LEGACY) t2$events else {
  d_by <- c(`0-5 yr` = sum(AF$event[an] == 1 & AF$time[an] <= 5), `>5 yr` = sum(AF$event[an] == 1 & AF$time[an] > 5))
  stopifnot(sum(an) == 2564, sum(d_by) == 675)
  unname(d_by[t2$period]) }
cm <- rd("P_cox_main.csv"); cm <- cm[cm$model == "Model 3" & cm$outcome == "event", ]   # Model 3, primary cohort
cm <- cm[match(relabel(t3$term), relabel(sub("^lca", "", cm$term))), ]; stopifnot(!anyNA(cm$HR))
S3 <- rbind(
  data.frame(Block = "A. Assumption tests", Row = t1$test,
             Estimate = ifelse(is.na(t1$statistic), "joint Wald", sprintf("chi2=%.2f", t1$statistic)),
             df = t1$df, P = pf(pA), Note = t1$note, check.names = FALSE),
  data.frame(Block = "B. Period-specific HR (follow-up split at 5 yr)", Row = paste(t2$period, "|", t2$term),
             Estimate = fh(t2$HR, t2$lo, t2$hi), df = "", P = pf(t2$p),
             Note = sprintf("%d deaths in period", as.integer(deaths_note)), check.names = FALSE),
  data.frame(Block = "C. Baseline hazard stratified on age group x years since diagnosis (16 strata)",
             Row = t3$term, Estimate = fh(t3$HR, t3$lo, t3$hi), df = "", P = pf(t3$p),
             Note = paste0("vs Model 3: ", fh(cm$HR, cm$lo, cm$hi)), check.names = FALSE))
write_table(S3, "TableS3_proportional_hazards.csv")

## ============================================ Table S5 ==========================================
m4 <- rd("RR_model4.csv"); Ev <- rd("RR_model4_evalues.csv"); att <- rd("RR_attenuation_comparison.csv")
pc4 <- rd("RR_model4_pseudoclass.csv"); cn <- rd("RR_contrast_pseudoclass.csv")
if (LEGACY) mult <- rd("RR_model4_multiplicity.csv") else {
  ad <- m4[m4$model == "Model 4 + antidepressant use", ]; cs <- rd("HM_model4_cause_specific.csv")
  mult <- rbind(data.frame(outcome = "All-cause", ad[, c("term","HR","lo","hi","p")]),
                data.frame(outcome = "Cancer",     cs[cs$outcome == "ev_ca",  c("term","HR","lo","hi","p")]),
                data.frame(outcome = "CVD/stroke", cs[cs$outcome == "ev_cvd", c("term","HR","lo","hi","p")]))
  mult$p_bh <- p.adjust(mult$p, "BH"); mult$p_bonf <- p.adjust(mult$p, "bonferroni") }
Eci <- Ev$E_value_CI
if (!LEGACY) Eci[Ev$HR >= 1 & Ev$CI_lower <= 1] <- 1       # interval includes 1 -> E-value for the limit is 1
S5 <- rbind(
  blk("A. Covariate specifications (modal class assignment)", paste(m4$model, "|", m4$term), fh(m4$HR, m4$lo, m4$hi),
      sprintf("phenotype joint P = %.4f; %d parameters", m4$joint_p, as.integer(m4$npar)), pf(m4$p)),
  blk("B. Model 4 primary testing family, multiplicity-corrected", paste(mult$outcome, "|", mult$term),
      fh(mult$HR, mult$lo, mult$hi), sprintf("BH = %.4f; Bonferroni = %.4f", mult$p_bh, mult$p_bonf), pf(mult$p)),
  blk("C. Assignment-uncertainty correction under Model 4 (100 draws)", pc4$term, fh(pc4$HR, pc4$lo, pc4$hi),
      paste("fraction of missing information =", py_str(pc4$fmi)), pf(pc4$p)),
  blk("D. Direct contrast (reference = insomnia-fatigue) under assignment correction", paste(cn$cov, "|", cn$term),
      fh(cn$HR, cn$lo, cn$hi), paste("fraction of missing information =", py_str(cn$fmi)), pf(cn$p)),
  blk("E. Attenuation by comorbidity and frailty adjustment", paste(att$comparison, "|", att$model),
      fh(att$HR, att$lo, att$hi),
      ifelse(is.na(att$attenuation_pct), "", paste0(py_str(att$attenuation_pct), "% of excess risk removed")), ""),
  blk("F. E-values", Ev$estimate, sprintf("HR %.2f (lower CI %.2f)", Ev$HR, Ev$CI_lower),
      sprintf("E-value %s (CI limit %s)", py_str(Ev$E_value_point), py_str(Eci)), ""))
write_table(S5, "TableS5_model4_detail.csv")

## ============================================ Table S7 ==========================================
icA <- rd("lca_fit_all_adults.csv"); pcA <- rd("RR_fullsample_lca_pseudoclass.csv")
dimr <- rd("RR_dimensional.csv"); pw <- rd("RR_power_confirmatory.csv")
efa <- read.csv(src("RR_efa_loadings.csv"), colClasses = "character", check.names = FALSE)
dimr <- dimr[!grepl("^D ", dimr$model), ]                    # keep the existing layout (specifications A1 to C2)
tl <- c("I(tot/5)" = "PHQ-9 total per 5 pts", "I(som/3)" = "somatic dimension per 3 pts",
        "I(aff/3)" = "affective dimension per 3 pts", "som_resid" = "somatic predominance at fixed total")
if (LEGACY) dimf <- rd("RR_dimensional_multiplicity.csv") else {
  da <- rd("RR_dimensional_affective.csv"); da <- da[da$block == "cause-specific", ]
  dimf <- data.frame(outcome = c(event = "All-cause", ev_ca = "Cancer", ev_cvd = "CVD/stroke")[da$model],
                     term = c("I(aff/3)" = "Affective dimension (per 3 pts)", "I(som/3)" = "Somatic dimension (per 3 pts)")[da$term],
                     HR = da$HR, lo = da$lo, hi = da$hi, p = da$p)
  dimf$p_bh <- p.adjust(dimf$p, "BH"); dimf$p_bonf <- p.adjust(dimf$p, "bonferroni") }
## block E: survey-weighted dimension means by phenotype in the 2,569 primary-cohort members (data/derived)
DIM <- read.csv("data/derived/nhanes_dimension_scores.csv")
AD <- merge(AF, DIM[, c("SEQN", "som", "aff", "tot")], by = "SEQN", all.x = TRUE)
Pp <- subset(svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt, nest = TRUE, data = AD), primary == 1)
mn <- lapply(c("som", "aff", "tot"), function(v) svyby(as.formula(paste0("~", v)), ~lca, Pp, svymean)[, 2])
names(mn) <- c("som", "aff", "tot")
## block F
pw_hr <- sub("^.*HR ([0-9.]+).*$", "\\1", pw$target)
base_n <- pw$participants_needed[pw$sample_multiplier == 1]; base_d <- pw$deaths_needed[pw$sample_multiplier == 1]
Fold <- data.frame(Row = pw$target, Estimate = sprintf("log HR %s, SE %s, z %s", py_str(pw$logHR), py_str(pw$SE), py_str(pw$z)),
                   Detail = sprintf("%sx current sample: %s survivors, %s deaths", py_str(pw$sample_multiplier),
                                    fmt_int(pw$participants_needed), fmt_int(pw$deaths_needed)))
if (!LEGACY) {
  ts <- rd("RR5_three_step_estimates.csv"); ts <- ts[ts$model == "C4" & ts$term == "Somatic-depressive vs low", ]
  stopifnot(nrow(ts) == 1)
  ratio <- (ts$se / (ts$logHR / (qnorm(0.975) + qnorm(0.8))))^2
  Fold$Row <- c(sprintf("Indicator-only pseudo-class draws, Model 4 (HR %s)", pw_hr[1]),
                sprintf("Modal assignment, Model 4 (HR %s)", pw_hr[2]))
  tb <- rd("RR5_three_step_estimates_bootstrap.csv"); tb <- tb[tb$model == "C4" & tb$term == "Somatic-depressive vs low", ]
  stopifnot(nrow(tb) == 1)
  f_row <- function(x, lab) {
    r <- (x$se / (x$logHR / (qnorm(0.975) + qnorm(0.8))))^2
    data.frame(Row = sprintf("Bias-adjusted three-step estimate, Model 4 (HR %.2f), %s", x$HR, lab),
               Estimate = sprintf("log HR %s, SE %s, z %s", py_str(py_round(x$logHR, 4)),
                                  py_str(py_round(x$se, 4)), py_str(py_round(x$logHR / x$se, 2))),
               Detail = sprintf("%sx current sample: %s survivors, %s deaths", py_str(py_round(r, 1)),
                                fmt_int(as.integer(trunc(base_n * r))), fmt_int(as.integer(trunc(base_d * r))))) }
  Fold <- rbind(f_row(ts, "jackknife SE"), f_row(tb, "bootstrap SE"), Fold)
  cat(sprintf("  three-step Model 4: log HR %.4f, SE %.4f -> multiplier %.2f\n", ts$logHR, ts$se, ratio)) }
S7 <- rbind(
  blk("A. Measurement model fitted in all adults (n = 34,022)", paste("k =", icA$k), paste("BIC", as.integer(round(icA$BIC))),
      sprintf("entropy %s; mean max posterior %s; smallest class %s%%", py_str(icA$entropy), py_str(icA$mean_maxpost),
              py_str(icA$min_class_pct)), ""),
  blk("B. Survivor subset under the all-adult k=4 fit", pcA$term, fh(pcA$HR, pcA$lo, pcA$hi),
      paste("fraction of missing information", py_str(pcA$fmi)), pf(pcA$p)),
  blk("C. Dimensional models (no class assignment)", paste(dimr$model, "|", tl[dimr$term]), fh(dimr$HR, dimr$lo, dimr$hi),
      "", pf(dimr$p)),
  blk("D. Dimensional family, multiplicity-corrected (Model 4)", paste(dimf$outcome, "|", dimf$term),
      fh(dimf$HR, dimf$lo, dimf$hi), sprintf("BH %.4f; Bonferroni %.4f", dimf$p_bh, dimf$p_bonf), pf(dimf$p)),
  blk("E. Dimension means by phenotype (weighted)", levels(AF$lca),
      sprintf("somatic %.2f / affective %.2f", mn$som, mn$aff), sprintf("total %.2f", mn$tot), ""),
  blk("F. Sample size for a confirmatory study (80% power)", Fold$Row, Fold$Estimate, Fold$Detail, ""),
  blk("G. Exploratory two-factor analysis of the nine PHQ-9 items (all adults, n = 34,022)", efa[[1]],
      sprintf("F1 (somatic) %s / F2 (affective) %s", efa$F1, efa$F2),
      sprintf("dominant %s; conventional subscale %s; %s; |loading difference| %.3f", efa$dominant, efa$literature,
              ifelse(efa$matches == "True", "matches", "differs"), abs(as.numeric(efa$gap))), ""))
write_table(S7, "TableS7_irreducibility_analyses.csv")
cat("inputs used:\n"); print(data.frame(file = names(used), path = unname(used)), row.names = FALSE)
