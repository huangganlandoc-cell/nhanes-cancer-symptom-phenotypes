## Sleep duration under the Model 3 covariate set without phenotype, the model used for
## Table 3 and Figure S3 (review/step2_analysis_plan.md, section 3).
## Run from the project root:  Rscript code/analysis_sleep_model3.R
suppressPackageStartupMessages({library(survey); library(survival); library(Hmisc)})
options(survey.lonely.psu = "adjust")
A <- readRDS("data/derived/analysis_frame_primary.rds"); A$cvd <- as.numeric(A$cvd)
C3 <- "age+sex+race4+educ+married+pir_i+pir_m+bmi_i+bmi_m+smoke+pa3+htn+dm+cvd+ydx_i+ydx_m+multi_primary+catype"
d0 <- subset(A, primary == 1)
dA <- d0[complete.cases(d0[, all.vars(as.formula(paste("~", C3)))]), ]   # 2,564 analysed in Model 3
kn <- quantile(d0$sleep_h, c(.05, .35, .65, .95), na.rm = TRUE)
B <- rcspline.eval(A$sleep_h, knots = kn, inclx = TRUE); colnames(B) <- paste0("sp", 1:ncol(B))
A <- cbind(A, as.data.frame(B))
P <- subset(svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt, nest = TRUE, data = A), primary == 1)
ms <- svycoxph(as.formula(paste("Surv(time, event) ~ sp1 + sp2 + sp3 +", C3)), design = P)
p_all <- regTermTest(ms, ~ sp1 + sp2 + sp3)$p; p_nl <- regTermTest(ms, ~ sp2 + sp3)$p
grid <- seq(quantile(d0$sleep_h, .01, na.rm = TRUE), quantile(d0$sleep_h, .99, na.rm = TRUE), by = 0.01)
G <- rcspline.eval(grid, knots = kn, inclx = TRUE)
lp <- as.vector(G %*% coef(ms)[c("sp1", "sp2", "sp3")])
nadir <- grid[which.min(lp)]
mc <- svycoxph(as.formula(paste("Surv(time, event) ~ relevel(factor(sleepcat), ref = '2') +", C3)), design = P)
k <- grep("sleepcat", names(coef(mc))); ci <- exp(confint(mc)[k, ])
out <- data.frame(quantity = c("knots (h)", "overall P", "nonlinearity P", "nadir (h)",
                               "short sleep HR (95% CI)", "long sleep HR (95% CI)", "n short / long"),
                  value = c(paste(round(kn, 2), collapse = ", "), signif(p_all, 3), signif(p_nl, 3), nadir,
                            sprintf("%.2f (%.2f-%.2f)", exp(coef(mc)[k[1]]), ci[1, 1], ci[1, 2]),
                            sprintf("%.2f (%.2f-%.2f)", exp(coef(mc)[k[2]]), ci[2, 1], ci[2, 2]),
                            paste(sum(dA$sleepcat == 1), sum(dA$sleepcat == 3), sep = " / ")))   # analysed sample
print(out, row.names = FALSE)
write.csv(out, "supporting/RR5_sleep_model3.csv", row.names = FALSE)
