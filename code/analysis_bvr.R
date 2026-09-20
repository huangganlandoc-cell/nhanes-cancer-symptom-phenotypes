## Local independence of the 4-class measurement model: bivariate residuals (BVR)
## for all 55 indicator pairs = Pearson X^2 of the observed versus model-implied two-way
## table divided by its degrees of freedom (review/step2_analysis_plan.md, section 2).
## Run from the project root:  Rscript code/analysis_bvr.R
source("code/three_step_functions.R")
suppressPackageStartupMessages(library(poLCA))
f4 <- readRDS("data/derived/lca_fits_cancer.rds")[[4]]
LAB <- c(phqi1 = "Anhedonia", phqi2 = "Depressed mood", phqi3 = "Sleep disturbance (PHQ-9)",
         phqi4 = "Fatigue", phqi5 = "Appetite change", phqi6 = "Worthlessness",
         phqi7 = "Concentration", phqi8 = "Psychomotor change", phqi9 = "Suicidal ideation",
         sleepcat = "Sleep duration (3 categories)", slq050b = "Previously reported sleep trouble")
b <- bvr_table(f4)
b$item1 <- LAB[b$item1]; b$item2 <- LAB[b$item2]
b <- b[order(-b$BVR), ]
b$flag <- ifelse(b$BVR > 3.84, "BVR > 3.84", "")
write.csv(transform(b, X2 = round(X2, 2), BVR = round(BVR, 2)), "supporting/RR5_bvr.csv", row.names = FALSE)
cat("pairs:", nrow(b), "| BVR > 3.84:", sum(b$BVR > 3.84), "| BVR > 10:", sum(b$BVR > 10), "\n")
print(head(transform(b, BVR = round(BVR, 1))[, c("item1", "item2", "df", "BVR")], 15), row.names = FALSE)
