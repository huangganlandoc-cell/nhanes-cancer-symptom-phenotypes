## Survivors excluded for incomplete symptom data compared with those included (post hoc, 2026-09-25).
## Of 3,779 survivors eligible for mortality linkage with follow-up and an examination weight, 511 lacked
## one or more PHQ-9 items, sleep duration or SLQ050. Survey-weighted descriptive comparison; no tests.
## Run from the project root (after code/build_cohort.py): Rscript code/analysis_excluded_comparison.R
## Output: supporting/RR7_excluded_comparison.csv
suppressPackageStartupMessages({library(survey)})
options(survey.lonely.psu = "adjust")
DER <- "data/derived"; if (!dir.exists(DER)) stop("run from the project root")
OCSV <- Sys.getenv("OUT", "supporting"); dir.create(OCSV, recursive = TRUE, showWarnings = FALSE)

X <- read.csv(file.path(DER, "eligible_survivors.csv.gz"))
stopifnot(nrow(X) == 3779, sum(X$included) == 3268)
X$group <- factor(ifelse(X$included == 1, "Included", "Excluded"), levels = c("Included", "Excluded"))
X$female <- as.numeric(X$RIAGENDR == 2)
X$race <- factor(c("Mexican American", "Other Hispanic", "NH White", "NH Black", "Other/Multi")[X$RIDRETH1])
X$educ_lt_hs <- ifelse(X$DMDEDUC2 %in% 1:2, 1, ifelse(X$DMDEDUC2 %in% 3:5, 0, NA))
X$htn <- ifelse(X$BPQ020 == 1, 1, ifelse(X$BPQ020 == 2, 0, NA))
X$dm <- ifelse(X$DIQ010 == 1, 1, ifelse(X$DIQ010 %in% 2:3, 0, NA))
X$proxy <- ifelse(X$MIAPROXY == 1, 1, ifelse(X$MIAPROXY == 2, 0, NA))
X$skin_only <- as.numeric(X$only_nms %in% c(TRUE, "True"))
X$yrs_dx <- pmax(X$RIDAGEYR - X$age_dx, 0)
D <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt, nest = TRUE, data = X)

f1 <- function(x, d = 1) formatC(x, format = "f", digits = d, big.mark = ",")
row <- function(label, fn) data.frame(Characteristic = label,
  Included = fn(subset(D, group == "Included")), Excluded = fn(subset(D, group == "Excluded")))
wm  <- function(v) function(d) { m <- svymean(as.formula(paste0("~", v)), d, na.rm = TRUE); sprintf("%s (%s)", f1(coef(m)), f1(SE(m))) }
wp  <- function(v) function(d) f1(100 * coef(svymean(as.formula(paste0("~", v)), d, na.rm = TRUE)))
tab <- rbind(
  row("Participants, n (unweighted)", function(d) f1(nrow(d$variables), 0)),
  row("Age, years, mean (SE)", wm("RIDAGEYR")),
  row("Female, %", wp("female")),
  row("Non-Hispanic White, %", function(d) f1(100 * coef(svymean(~I(race == "NH White"), d))[2])),
  row("Less than high school education, %", wp("educ_lt_hs")),
  row("Income-to-poverty ratio, mean (SE)", wm("INDFMPIR")),
  row("Body mass index, kg/m2, mean (SE)", wm("BMXBMI")),
  row("Hypertension, %", wp("htn")),
  row("Diabetes, %", wp("dm")),
  row("Cardiovascular disease, %", wp("cvd")),
  row("Years since cancer diagnosis, mean (SE)", wm("yrs_dx")),
  row("Skin cancer only, %", wp("skin_only")),
  row("Examination interview by proxy, %", wp("proxy")),
  row("Deaths, n (unweighted)", function(d) f1(sum(d$variables$event), 0)),
  row("Deaths per 1,000 person-years (unweighted)", function(d) f1(1000 * sum(d$variables$event) / sum(d$variables$time))),
  row("Died during follow-up, weighted %", wp("event")))
ex <- subset(X, included == 0)
miss <- data.frame(Characteristic = c("Excluded: all nine PHQ-9 items missing", "Excluded: one to eight PHQ-9 items missing",
                                      "Excluded: PHQ-9 complete, sleep duration or SLQ050 missing"),
                   Included = "",
                   Excluded = f1(c(sum(ex$n_phq_missing == 9), sum(ex$n_phq_missing %in% 1:8),
                                   sum(ex$n_phq_missing == 0 & ex$sleep_missing == 1)), 0))
tab <- rbind(tab, miss)
print(tab, row.names = FALSE)
write.csv(tab, file.path(OCSV, "RR7_excluded_comparison.csv"), row.names = FALSE)
