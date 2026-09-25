## Measurement model without SLQ050 (post hoc, 2026-09-24).
## SLQ050 records whether trouble sleeping was ever raised with a health professional, so it partly
## measures clinical contact, and it is the indicator that most separates the two intermediate
## phenotypes. The ten remaining indicators are refitted with the same settings as the main model,
## classes are matched to the main solution by cross-tabulation, and the phenotype-mortality models
## are refitted with modal assignment, as for the other measurement sensitivity analyses
## (code/reconstructed/measurement_sensitivity_RR4.R).
## Run from the project root: Rscript code/analysis_drop_slq050.R
## Outputs: supporting/RR6_drop_slq050_{fit,profiles,crosstab,cox}.csv
suppressPackageStartupMessages({library(poLCA); library(survey); library(survival)})
options(survey.lonely.psu = "adjust")
DER <- "data/derived"; OUT <- "supporting"

H  <- readRDS(file.path(DER, "analysis_frame_primary.rds"))
G4 <- read.csv(file.path(DER, "nhanes_model4_design_frame.csv.gz"))
keep <- c("SEQN", "comorb_n_i", "LBXHGB_i", "LBXHGB_m", "LBXSAL_i", "LBXSAL_m",
          "egfr_i", "egfr_m", "func_lim_i", "n_rx_i", "n_rx_m", "antidep_i")
H <- merge(H, G4[, keep], by = "SEQN", all.x = TRUE)
H$cvd <- as.numeric(H$cvd); H$bio_m <- pmax(H$LBXSAL_m, H$egfr_m)
H$lca <- factor(as.character(H$lca),
                levels = c("Low symptom burden", "Insomnia-fatigue", "Hypersomnia-somatic", "High symptom burden"),
                labels = c("Low symptom burden", "Sleep-fatigue", "Somatic-depressive", "High symptom burden"))
C3  <- "age+sex+race4+educ+married+pir_i+pir_m+bmi_i+bmi_m+smoke+pa3+htn+dm+cvd+ydx_i+ydx_m+multi_primary+catype"
C4a <- paste(C3, "comorb_n_i+egfr_i+egfr_m", sep = "+")
C4  <- paste(C3, "comorb_n_i+LBXHGB_i+LBXHGB_m+LBXSAL_i+egfr_i+bio_m+func_lim_i+n_rx_i+n_rx_m+antidep_i", sep = "+")

D <- subset(H, inAnalysis == 1)
cat("measurement sample:", nrow(D), "\n")
ent <- function(p) 1 - sum(-p * log(pmax(p, 1e-12))) / (nrow(p) * log(ncol(p)))
F10 <- cbind(phqi1, phqi2, phqi3, phqi4, phqi5, phqi6, phqi7, phqi8, phqi9, sleepcat) ~ 1

## same settings as the main measurement model: 30 random starts, up to 8,000 iterations
set.seed(20260924)
fits <- lapply(3:5, function(k) poLCA(F10, D, nclass = k, maxiter = 8000, nrep = 30, verbose = FALSE, calc.se = FALSE))
fitstat <- do.call(rbind, lapply(fits, function(f) {
  k <- length(f$P); n <- f$N; np <- f$npar
  data.frame(classes = k, loglik = f$llik, BIC = f$bic, cAIC = -2 * f$llik + np * (log(n) + 1),
             entropy = ent(f$posterior), smallest_class_pct = 100 * min(f$P))
}))
print(fitstat, digits = 5)
f <- fits[[2]]                                                   # four classes, as in the main solution
D$new <- apply(f$posterior, 1, which.max)
ct <- table(original = D$lca, refit = D$new)
V <- sqrt(suppressWarnings(chisq.test(ct))$statistic / (sum(ct) * (min(dim(ct)) - 1)))
## label each refitted class by the original phenotype it shares most members with
lab <- rownames(ct)[apply(ct, 2, which.max)]
if (anyDuplicated(lab)) warning("refitted classes do not map one-to-one onto the original phenotypes")
D$new <- factor(lab[D$new], levels = levels(D$lca))
ct2 <- table(original = D$lca, refit = D$new)
agree <- sum(diag(ct2)) / sum(ct2)
cat(sprintf("four classes: entropy %.3f | Cramer's V %.3f | agreement with the main assignment %.1f%%\n",
            ent(f$posterior), V, 100 * agree))
print(ct2)

## item-response profiles (probability of the symptom present; sleep categories)
pr <- f$probs; cls <- lab
prof <- rbind(sapply(1:9, function(j) pr[[j]][, 2]) |> t(), short = pr[[10]][, 1], long = pr[[10]][, 3])
colnames(prof) <- cls; rownames(prof)[1:9] <- paste0("phq", 1:9)
prof <- round(100 * prof[, levels(D$lca)[levels(D$lca) %in% cls]], 1)
print(prof)

## Cox models, modal assignment, primary cohort, reference = the class matched to low symptom burden
Hh <- merge(H, D[, c("SEQN", "new")], by = "SEQN", all.x = TRUE)
des <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt, nest = TRUE, data = Hh)
P <- subset(des, primary == 1 & !is.na(new))
rows <- list()
for (m in c("C3", "C4a", "C4")) {
  sv <- svycoxph(as.formula(paste("Surv(time, event) ~ new +", get(m))), design = P)
  k <- grep("^new", names(coef(sv)))
  b <- coef(sv)[k]; se <- sqrt(diag(vcov(sv)))[k]; nm <- sub("^new", "", names(b))
  ctr <- svycontrast(sv, setNames(c(1, -1), c("newSomatic-depressive", "newSleep-fatigue")))
  b <- c(b, coef(ctr)); se <- c(se, SE(ctr))
  nm <- c(paste(nm, "vs low symptom burden"), "Somatic-depressive vs sleep-fatigue")
  rows[[m]] <- data.frame(model = c(C3 = "Model 3", C4a = "Model 4a", C4 = "Model 4")[m], term = nm,
                          HR = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se),
                          p = 2 * pnorm(-abs(b / se)), n = sv$n, deaths = sv$nevent, row.names = NULL)
}
cox <- do.call(rbind, rows)
print(transform(cox, HR95 = sprintf("%.2f (%.2f-%.2f)", HR, lo, hi))[, c("model", "term", "HR95", "p")], row.names = FALSE)

write.csv(cbind(fitstat, cramer_v_4class = c(NA, V, NA), agreement_4class = c(NA, agree, NA)),
          file.path(OUT, "RR6_drop_slq050_fit.csv"), row.names = FALSE)
write.csv(prof, file.path(OUT, "RR6_drop_slq050_profiles.csv"))
write.csv(as.data.frame.matrix(ct2), file.path(OUT, "RR6_drop_slq050_crosstab.csv"))
write.csv(cox, file.path(OUT, "RR6_drop_slq050_cox.csv"), row.names = FALSE)
cat("analysis_drop_slq050.R done\n")
