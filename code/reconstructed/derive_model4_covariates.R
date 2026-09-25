## =============================================================================
## derive_model4_covariates.R  --  how the Model 4 covariates were built
##
## PART A RUNS ONLY WHEN these 35 NHANES files are in data/raw_nhanes (the original run downloaded
## them but never saved them; code/download_nhanes.py now fetches them):
##   KIQ_U  (kidney conditions, KIQ022 "weak/failing kidneys"):
##          KIQ_U_D.xpt KIQ_U_E.xpt KIQ_U_F.xpt KIQ_U_G.xpt KIQ_U_H.xpt KIQ_U_I.xpt KIQ_U_J.xpt
##   CBC    (complete blood count, LBXHGB haemoglobin):
##          CBC_D.xpt CBC_E.xpt CBC_F.xpt CBC_G.xpt CBC_H.xpt CBC_I.xpt CBC_J.xpt
##   BIOPRO (standard biochemistry, LBXSAL albumin, LBXSCR creatinine):
##          BIOPRO_D.xpt BIOPRO_E.xpt BIOPRO_F.xpt BIOPRO_G.xpt BIOPRO_H.xpt BIOPRO_I.xpt BIOPRO_J.xpt
##   PFQ    (physical functioning: PFQ049/054/057/059/090 and PFQ061A-T):
##          PFQ_D.xpt PFQ_E.xpt PFQ_F.xpt PFQ_G.xpt PFQ_H.xpt PFQ_I.xpt PFQ_J.xpt
##   RXQ_RX (prescription medications in the past 30 days: RXDUSE, RXDDRUG, RXDCOUNT):
##          RXQ_RX_D.xpt RXQ_RX_E.xpt RXQ_RX_F.xpt RXQ_RX_G.xpt RXQ_RX_H.xpt RXQ_RX_I.xpt RXQ_RX_J.xpt
##   (cycles D..J = 2005-2006 .. 2017-2018; source https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/
##    {2005,2007,...,2017}/DataFiles/<file>). download_nhanes.py saves them with the project naming
##    convention (e.g. KIQ_U_D_2005-2006.csv.gz, RXQ_RX_J_2017-2018.csv.gz) and Part A then runs
##    automatically. Nothing is downloaded by this script.
## The analysis scripts (model4.R, dimensional.R, all_adult_lca.R, measurement_sensitivity_RR4.R)
## start from the saved files data/derived/nhanes_model4_covariates.csv and
## data/derived/nhanes_model4_design_frame.csv.gz.
##
## What this script DOES run:
##   Part A-check : the 6 conditions taken from MCQ (which IS in data/raw_nhanes) are re-derived
##                  with Part A's code and compared with the saved nhanes_model4_covariates.csv
##   Part B       : step 251 (imputation + merge into the 70,190-row design frame), run from the
##                  saved covariates file and compared with the saved nhanes_model4_design_frame
## Outputs: review/reproduction/output/input_checks/check_model4_mcq_conditions.csv
##          review/reproduction/output/input_checks/check_model4_design_frame.csv
##          data/derived/nhanes_model4_covariates.csv       only if absent and Part A ran
##          data/derived/nhanes_model4_design_frame.csv.gz  only if absent
##          (both in pandas' to_csv layout, the design-frame columns copied as text from
##          nhanes_cancer_design_frame.csv.gz; an existing file is never overwritten, and the
##          comparisons below then read the file as before)
##
## Reconstructed from the original analysis log:
##   step 244: downloaded KIQ_U, CBC, BIOPRO, PFQ, RXQ_RX (into EXTRA)
##   steps 245/246: inspection of column availability (not repeated)
##   step 247: antidepressant list (30 generic names) and medication counts
##   step 248: conditions, laboratory values, PFQ061 limitations, eGFR (CKD-EPI 2021)
##             -> nhanes_model4_covariates.csv (first version)
##   step 249: PFQ061 missingness by age (not repeated)
##   step 250: all-age functional-limitation items -> func_lim, func_n; file re-written
##   step 251: imputation, merge into the design frame -> nhanes_model4_design_frame.csv
## Faithfulness notes (reproduced, not fixed):
##   * medians for imputing haemoglobin, albumin, eGFR, the condition count and the medication count
##     are taken over ALL 67,203 MCQ respondents, not over the survivor cohort;
##   * a missing condition is imputed as 0 (absent) and a missing functional-limitation or
##     antidepressant value as 0;
##   * albumin and creatinine come from the same file, so LBXSAL_m and egfr_m are identical (the
##     analysis scripts replace them by one indicator, bio_m).
## Run   : /opt/homebrew/bin/Rscript code/reconstructed/derive_model4_covariates.R (from the project root)
## =============================================================================
if (!dir.exists("data/raw_nhanes")) stop("run from the project root")
OUT <- "review/reproduction/output"
dir.create(file.path(OUT, "input_checks"), recursive = TRUE, showWarnings = FALSE)
DER <- "data/derived"; RAWD <- "data/raw_nhanes"
CY <- c(`2005`="_D", `2007`="_E", `2009`="_F", `2011`="_G", `2013`="_H", `2015`="_I", `2017`="_J")
rawfile <- function(mod, y) file.path(RAWD, sprintf("%s%s_%s-%d.csv.gz", mod, CY[[y]], y, as.integer(y) + 1))
left_join <- function(x, y, by = "SEQN") {              # pandas merge(how="left"): keeps x's row order
  i <- match(x[[by]], y[[by]]); out <- cbind(x, y[i, setdiff(names(y), by), drop = FALSE]); rownames(out) <- NULL; out }
yn <- function(v) ifelse(v == 1, 1, ifelse(v == 2, 0, NA))   # 1 = yes, 2 = no, anything else missing
py_num <- function(x) {                                        # Python repr() of a float, as pandas writes it
  fmt <- function(v, k) { e <- as.integer(sub(".*e", "", sprintf("%.*e", k - 1L, v)))
    ifelse(e >= -4 & e < 16, sprintf("%.*f", pmax(k - 1L - e, 1L), v), sprintf("%.*e", k - 1L, v)) }
  s <- ifelse(is.na(x), "", ifelse(x > 0, "inf", "-inf")); i <- which(is.finite(x))
  for (k in 1:17) { if (!length(i)) break
    f <- fmt(x[i], k); ok <- k == 17L | as.numeric(f) == x[i]; s[i[ok]] <- f[ok]; i <- i[!ok] }
  s }
py_lines <- function(df) {                                     # lines of pandas DataFrame.to_csv(index=False)
  cols <- lapply(df, function(col) {
    if (is.logical(col)) ifelse(is.na(col), "", ifelse(col, "True", "False"))
    else if (is.integer(col)) ifelse(is.na(col), "", as.character(col))
    else if (is.numeric(col)) py_num(col)
    else { s <- as.character(col); s[is.na(s)] <- ""; q <- grepl('[,"\n]', s)
           s[q] <- paste0('"', gsub('"', '""', s[q]), '"'); s } })
  c(paste(names(df), collapse = ","), do.call(paste, c(cols, sep = ",")))
}

## =============================== PART A (steps 247, 248, 250) ================================
MODS <- c("KIQ_U", "CBC", "BIOPRO", "PFQ", "RXQ_RX")
need <- unlist(lapply(MODS, function(m) sapply(names(CY), function(y) rawfile(m, y))))
have_all <- all(file.exists(need))
derive_part_A <- function() {
  ## ---- step 247: antidepressants and prescription counts ----
  AD <- c('FLUOXETINE','SERTRALINE','PAROXETINE','CITALOPRAM','ESCITALOPRAM','FLUVOXAMINE',
          'VENLAFAXINE','DESVENLAFAXINE','DULOXETINE','LEVOMILNACIPRAN','MILNACIPRAN',
          'BUPROPION','MIRTAZAPINE','TRAZODONE','NEFAZODONE','VILAZODONE','VORTIOXETINE',
          'AMITRIPTYLINE','NORTRIPTYLINE','IMIPRAMINE','DESIPRAMINE','DOXEPIN','CLOMIPRAMINE',
          'PROTRIPTYLINE','TRIMIPRAMINE','AMOXAPINE','MAPROTILINE',
          'PHENELZINE','TRANYLCYPROMINE','ISOCARBOXAZID')
  RX <- do.call(rbind, lapply(names(CY), function(y) {
    d <- read.csv(rawfile("RXQ_RX", y), stringsAsFactors = FALSE)
    drug <- toupper(trimws(ifelse(is.na(d$RXDDRUG), "", as.character(d$RXDDRUG))))
    data.frame(SEQN = d$SEQN, RXDUSE = d$RXDUSE, RXDCOUNT = d$RXDCOUNT, drug = drug) }))
  RX$is_ad <- Reduce(`|`, lapply(AD, function(a) grepl(a, RX$drug, fixed = TRUE)))   # substring match
  grp_max <- function(v, g) { r <- tapply(v, g, function(z) if (all(is.na(z))) NA else max(z, na.rm = TRUE)); r }
  grp_min <- function(v, g) { r <- tapply(v, g, function(z) if (all(is.na(z))) NA else min(z, na.rm = TRUE)); r }
  adp   <- table(RX$SEQN[RX$is_ad])                                      # n_ad (count of AD records)
  medc  <- grp_max(RX$RXDCOUNT, RX$SEQN)                                 # n_rx
  anyrx <- grp_min(RX$RXDUSE, RX$SEQN)                                   # rxuse
  ## ---- step 248: conditions, labs, PFQ061 limitations, eGFR ----
  blocks <- lapply(names(CY), function(y) {
    mcq <- read.csv(rawfile("MCQ", y)); kiq <- read.csv(rawfile("KIQ_U", y)); cbc <- read.csv(rawfile("CBC", y))
    bio <- read.csv(rawfile("BIOPRO", y)); pfq <- read.csv(rawfile("PFQ", y))
    X <- data.frame(SEQN = mcq$SEQN)
    for (p in list(c("MCQ160B","hf"), c("MCQ160F","stroke"), c("MCQ160G","emphysema"),
                   c("MCQ160K","bronchitis"), c("MCQ160L","liver"), c("MCQ160A","arthritis")))
      X[[p[2]]] <- if (p[1] %in% names(mcq)) yn(mcq[[p[1]]]) else NA
    X <- left_join(X, kiq[, c("SEQN","KIQ022")]); X$ckd <- yn(X$KIQ022); X$KIQ022 <- NULL
    X <- left_join(X, cbc[, c("SEQN","LBXHGB")])
    X <- left_join(X, bio[, c("SEQN","LBXSAL","LBXSCR")])
    p61 <- grep("^PFQ061[A-T]$", names(pfq), value = TRUE)
    lim <- as.matrix(pfq[, p61]); lim[!(lim <= 4) | is.na(lim)] <- NA
    allna <- rowSums(!is.na(lim)) == 0
    pf <- data.frame(SEQN = pfq$SEQN, pfq_n = rowSums(lim >= 2, na.rm = TRUE),
                     pfq_any = as.numeric(rowSums(lim >= 2, na.rm = TRUE) > 0))
    pf$pfq_n[allna] <- NA; pf$pfq_any[allna] <- NA
    X <- left_join(X, pf); X$cycle <- sprintf("%s-%d", y, as.integer(y) + 1); X })
  M4 <- do.call(rbind, blocks)
  M4$n_ad <- as.numeric(adp[as.character(M4$SEQN)])
  M4$n_rx <- as.numeric(medc[as.character(M4$SEQN)]); M4$rxuse <- as.numeric(anyrx[as.character(M4$SEQN)])
  M4$antidep <- ifelse(!is.na(M4$n_ad), 1, ifelse(!is.na(M4$rxuse), 0, NA))
  M4$n_rx[!is.na(M4$n_rx) & M4$n_rx > 90] <- NA
  M4$n_rx[!is.na(M4$rxuse) & M4$rxuse == 2 & is.na(M4$n_rx)] <- 0
  dem <- do.call(rbind, lapply(names(CY), function(y) read.csv(rawfile("DEMO", y))[, c("SEQN","RIAGENDR","RIDAGEYR")]))
  M4 <- left_join(M4, dem)
  f <- M4$RIAGENDR == 2
  k <- ifelse(f, 0.7, 0.9); al <- ifelse(f, -0.241, -0.302)
  scr <- ifelse(M4$LBXSCR >= 0.1 & M4$LBXSCR <= 20, M4$LBXSCR, NA)
  M4$egfr <- 142 * pmin(scr/k, 1)^al * pmax(scr/k, 1)^-1.200 * 0.9938^M4$RIDAGEYR * ifelse(f, 1.012, 1.0)  # CKD-EPI 2021
  M4 <- M4[, setdiff(names(M4), c("RIAGENDR","RIDAGEYR","n_ad","rxuse","LBXSCR"))]
  ## ---- step 250: functional limitation from items asked at all ages ----
  ALT <- do.call(rbind, lapply(names(CY), function(y) {
    pfq <- read.csv(rawfile("PFQ", y))
    cs <- intersect(c("PFQ049","PFQ054","PFQ057","PFQ059","PFQ090"), names(pfq))
    d <- sapply(cs, function(c) yn(pfq[[c]])); d <- matrix(d, ncol = length(cs))
    nn <- rowSums(!is.na(d))
    data.frame(SEQN = pfq$SEQN,
               func_lim = ifelse(nn == 0, NA, suppressWarnings(apply(d, 1, max, na.rm = TRUE))),
               func_n   = ifelse(nn == 0, NA, rowSums(d, na.rm = TRUE))) }))
  left_join(M4, ALT)
}
if (have_all) {
  cat("Part A: all 35 module files found -> deriving nhanes_model4_covariates from raw files\n")
  M4 <- derive_part_A()
  dir.create(file.path(OUT, "intermediate"), showWarnings = FALSE)
  write.csv(M4, file.path(OUT, "intermediate", "nhanes_model4_covariates_rederived.csv"), row.names = FALSE)
  if (!file.exists(file.path(DER, "nhanes_model4_covariates.csv"))) {
    writeLines(py_lines(M4), file.path(DER, "nhanes_model4_covariates.csv"))
    cat("data/derived/nhanes_model4_covariates.csv did not exist: wrote it from Part A\n")
  }
} else {
  cat("Part A NOT RUN: requires modules KIQ_U, CBC, BIOPRO, PFQ and RXQ_RX for all 7 cycles;",
      sum(!file.exists(need)), "of", length(need), "files are missing from data/raw_nhanes.\n",
      "Part B below starts from the saved data/derived/nhanes_model4_covariates.csv instead.\n")
  M4 <- read.csv(file.path(DER, "nhanes_model4_covariates.csv"))
}

## ============== PART A-check: the MCQ-based conditions (MCQ is in data/raw_nhanes) ==============
saved <- read.csv(file.path(DER, "nhanes_model4_covariates.csv"))
mcqc <- do.call(rbind, lapply(names(CY), function(y) {
  mcq <- read.csv(rawfile("MCQ", y)); X <- data.frame(SEQN = mcq$SEQN)
  for (p in list(c("MCQ160B","hf"), c("MCQ160F","stroke"), c("MCQ160G","emphysema"),
                 c("MCQ160K","bronchitis"), c("MCQ160L","liver"), c("MCQ160A","arthritis")))
    X[[p[2]]] <- if (p[1] %in% names(mcq)) yn(mcq[[p[1]]]) else NA
  X$cycle <- sprintf("%s-%d", y, as.integer(y) + 1); X }))
chkA <- data.frame(column = c("SEQN (same rows, same order)", "hf","stroke","emphysema","bronchitis","liver","arthritis","cycle"),
  n_rows = nrow(saved),
  cells_differing = c(sum(mcqc$SEQN != saved$SEQN),
    sapply(c("hf","stroke","emphysema","bronchitis","liver","arthritis"), function(v)
      sum(xor(is.na(mcqc[[v]]), is.na(saved[[v]]))) + sum(mcqc[[v]] != saved[[v]], na.rm = TRUE)),
    sum(mcqc$cycle != saved$cycle)))
cat("\n=== Part A-check: conditions re-derived from raw MCQ vs saved nhanes_model4_covariates.csv ===\n")
print(chkA, row.names = FALSE)
write.csv(chkA, file.path(OUT, "input_checks", "check_model4_mcq_conditions.csv"), row.names = FALSE)

## =============================== PART B (step 251) ==========================================
cond <- c('hf','stroke','emphysema','bronchitis','liver','arthritis','ckd')
nn <- rowSums(!is.na(M4[, cond]))
M4$comorb_n <- ifelse(nn == 0, NA, rowSums(M4[, cond], na.rm = TRUE))      # sum(min_count=1)
for (c in c("LBXHGB","LBXSAL","egfr")) {
  M4[[paste0(c,"_i")]] <- ifelse(is.na(M4[[c]]), median(M4[[c]], na.rm = TRUE), M4[[c]])
  M4[[paste0(c,"_m")]] <- as.integer(is.na(M4[[c]])) }
for (c in cond) M4[[paste0(c,"_i")]] <- ifelse(is.na(M4[[c]]), 0, M4[[c]])
M4$comorb_n_i <- ifelse(is.na(M4$comorb_n), median(M4$comorb_n, na.rm = TRUE), M4$comorb_n)
M4$n_rx_i <- ifelse(is.na(M4$n_rx), median(M4$n_rx, na.rm = TRUE), M4$n_rx); M4$n_rx_m <- as.integer(is.na(M4$n_rx))
M4$func_lim_i <- ifelse(is.na(M4$func_lim), 0, M4$func_lim); M4$antidep_i <- ifelse(is.na(M4$antidep), 0, M4$antidep)
keep <- c('SEQN','comorb_n_i','LBXHGB_i','LBXHGB_m','LBXSAL_i','LBXSAL_m','egfr_i','egfr_m',
          'func_lim_i','func_n','n_rx_i','n_rx_m','antidep_i', paste0(cond, "_i"))
GD <- left_join(read.csv(file.path(DER, "nhanes_cancer_design_frame.csv.gz"), check.names = FALSE), M4[, keep])
for (c in keep[-1]) GD[[c]][is.na(GD[[c]])] <- 0
if (!file.exists(file.path(DER, "nhanes_model4_design_frame.csv.gz"))) {
  txt <- readLines(file.path(DER, "nhanes_cancer_design_frame.csv.gz")); stopifnot(length(txt) == nrow(GD) + 1)
  gz <- gzfile(file.path(DER, "nhanes_model4_design_frame.csv.gz"), "w")
  writeLines(paste(txt, py_lines(GD[, keep[-1]]), sep = ","), gz); close(gz)
  cat("data/derived/nhanes_model4_design_frame.csv.gz did not exist: wrote the rebuilt frame;",
      "the comparison below then checks the written file against the in-memory frame\n")
}
SV <- read.csv(file.path(DER, "nhanes_model4_design_frame.csv.gz"), check.names = FALSE)
chkB <- data.frame(column = names(SV), in_rederived = names(SV) %in% names(GD),
  differing_cells = sapply(names(SV), function(v) {
    if (!v %in% names(GD)) return(NA); a <- SV[[v]]; b <- GD[[v]]
    if (is.numeric(a) && is.numeric(b)) sum(xor(is.na(a), is.na(b))) + sum(abs(a - b) > 1e-9 * pmax(1, abs(a)), na.rm = TRUE)
    else sum(xor(is.na(a), is.na(b))) + sum(as.character(a) != as.character(b), na.rm = TRUE) }),
  max_abs_diff = sapply(names(SV), function(v) {
    a <- SV[[v]]; b <- GD[[v]]; if (is.numeric(a) && is.numeric(b)) max(c(0, abs(a - b)), na.rm = TRUE) else NA }))
cat("\n=== Part B: design frame rebuilt from the saved covariates vs saved nhanes_model4_design_frame ===\n")
cat("rows:", nrow(GD), "vs", nrow(SV), "| columns:", ncol(GD), "vs", ncol(SV),
    "| columns with any differing cell:", sum(chkB$differing_cells > 0, na.rm = TRUE), "\n")
print(chkB[chkB$column %in% keep, ], row.names = FALSE)
write.csv(chkB, file.path(OUT, "input_checks", "check_model4_design_frame.csv"), row.names = FALSE)
cat("derive_model4_covariates.R done\n")
