## =============================================================================
## prior_classification.R  --  our phenotypes versus the Lan et al. (2024) six-cell
##                             classification (PHQ-9 category x sleep complaint)
##
## Regenerates (into review/reproduction/output/):
##   RR_vs_prior_classification.csv  Lan-style classification alone; phenotype adjusted for it
##   RR_crosstab_vs_prior.csv        cross-tabulation phenotype x Lan cell (primary cohort)
##
## Reconstructed from the original analysis log:
##   steps 205/206: setup objects G, C2, C3, mk(), tid() from the
##                  sourced analysis_reviewer_response.R
##   steps 208-216: added helper columns to G (late5, I2-I4, cyc,
##                  .k, catype2) that do not enter these models;
##                  not repeated
##   step 233: a = primary cohort, PHQ-9 category, sleep complaint
##   step 234: Lan-style classification; models mlan and mboth
##   step 235: writes both CSVs
## Run   : /opt/homebrew/bin/Rscript code/reconstructed/prior_classification.R  (from the project root)
## =============================================================================
suppressPackageStartupMessages({library(survey); library(survival); library(rms); library(poLCA)})
if (!file.exists("data/derived/analysis_frame_primary.rds")) stop("run from the project root")
OUT <- "review/reproduction/output"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

## ---- setup lines of analysis_reviewer_response.R (sourced in step 206) ---------------
options(survey.lonely.psu="adjust")
G <- readRDS("data/derived/analysis_frame_primary.rds")
G$cvd <- as.numeric(G$cvd)
levels(G$lca)[levels(G$lca)=="Hypersomnia-somatic"] <- "Somatic-depressive"
C2 <- "age+sex+race4+educ+married+pir_i+pir_m"
C3 <- paste(C2,"bmi_i+bmi_m+smoke+pa3+htn+dm+cvd+ydx_i+ydx_m+multi_primary+catype",sep="+")
mk  <- function(d) svydesign(ids=~SDMVPSU, strata=~SDMVSTRA, weights=~wt, nest=TRUE, data=d)
tid <- function(m,pat){s<-summary(m)$coef;k<-grep(pat,rownames(s));ci<-confint(m)[k,,drop=FALSE]
  data.frame(term=gsub(pat,"",rownames(s)[k]),HR=exp(s[k,"coef"]),lo=exp(ci[,1]),hi=exp(ci[,2]),p=s[k,ncol(s)])}

## ---- step 233 --------------------------------------------------------------------------
a <- subset(G, primary==1)
levels(a$lca)[levels(a$lca)=="Hypersomnia-somatic"] <- "Somatic-depressive"
## Lan 2024 classification: PHQ-9 category x self-perceived sleep disturbance (yes/no)
a$phqcat <- cut(a$phq9_score, c(-1,4,9,27), labels=c("0-4","5-9",">=10"))
a$slpdist <- factor(ifelse(a$insomnia=="Yes" | a$insomnia==1 | a$insomnia=="1","Yes","No"))
cat("insomnia levels:", paste(levels(factor(a$insomnia)),collapse=" / "), "\n")
a$slpdist <- factor(ifelse(as.character(a$insomnia)==levels(factor(a$insomnia))[1],"A","B"))
print(table(a$insomnia, a$slpdist))

## ---- step 234 --------------------------------------------------------------------------
a$slpdist <- factor(as.character(a$insomnia), levels=c("No","Yes"))
a$lan <- interaction(a$phqcat, a$slpdist, sep=" / ", lex.order=TRUE)
a$lan <- relevel(a$lan, ref="0-4 / No")
cat("=== Lan 2024 six cells x our four phenotypes (n=2,569) ===\n")
print(table(a$lca, a$lan, dnn=c("phenotype","PHQ-9 / sleep complaint")))
cat("\ndeaths per cell:\n"); print(tapply(a$event, a$lan, sum))
G$phqcat <- cut(G$phq9_score, c(-1,4,9,27), labels=c("0-4","5-9",">=10"))
G$slpdist <- factor(as.character(G$insomnia), levels=c("No","Yes"))
G$lan <- relevel(interaction(G$phqcat, G$slpdist, sep=" / ", lex.order=TRUE), ref="0-4 / No")
Pl <- subset(mk(G), primary==1)
mlan <- svycoxph(as.formula(paste0("Surv(time,event)~lan+",C3)), design=Pl)
print(transform(tid(mlan,"^lan"), HR95=sprintf("%.2f (%.2f-%.2f)",HR,lo,hi))[,c("term","HR95","p")],
      row.names=FALSE, digits=3)
mboth <- svycoxph(as.formula(paste0("Surv(time,event)~lca+lan+",C3)), design=Pl)
print(transform(tid(mboth,"^lca"), HR95=sprintf("%.2f (%.2f-%.2f)",HR,lo,hi))[,c("term","HR95","p")],
      row.names=FALSE, digits=3)
cat("phenotype joint Wald p =", signif(regTermTest(mboth,~lca)$p,4), "\n")
cat("Lan classification joint Wald p =", signif(regTermTest(mboth,~lan)$p,4),
    " (alone =", signif(regTermTest(mlan,~lan)$p,4), ")\n")

## ---- step 235 --------------------------------------------------------------------------
ct <- table(a$lca, a$lan)
ov <- as.data.frame.matrix(ct); ov$total <- rowSums(ct)
ov$pct_with_complaint <- round(rowSums(ct[,grep("Yes",colnames(ct))])/rowSums(ct)*100,1)
print(ov)
res <- rbind(
 cbind(spec="Lan-style joint classification alone", tid(mlan,"^lan")),
 cbind(spec="phenotype, adjusted for Lan-style classification", tid(mboth,"^lca")))
res$joint_p <- c(rep(signif(regTermTest(mlan,~lan)$p,4), sum(grepl("^Lan",res$spec))),
                 rep(signif(regTermTest(mboth,~lca)$p,4), sum(grepl("^phenotype",res$spec))))
write.csv(res, file.path(OUT, "RR_vs_prior_classification.csv"), row.names=FALSE)
write.csv(ov,  file.path(OUT, "RR_crosstab_vs_prior.csv"))
cat("prior_classification.R done\n")
