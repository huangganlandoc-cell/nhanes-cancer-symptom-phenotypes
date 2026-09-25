"""Build the 70,190-row survey design frame for the cancer-survivor cohort.

Pools NHANES 2005-2018 (7 cycles), flags self-reported cancer survivors (MCQ220), links the
NCHS public-use mortality file, derives the PHQ-9 items, sleep indicators, covariates, cancer
type and years since diagnosis, and marks the analytic subpopulation (inAnalysis) and the
skin-cancer-only survivors excluded from the primary cohort (only_nms). Prints the participant
flow of Additional file 1: Fig. S1.

Inputs : data/raw_nhanes/{DEMO,MCQ,DPQ,SLQ,BMX,SMQ,BPQ,DIQ,PAQ}_{D..J}_{cycle}.csv.gz
         data/mortality/nchs_linked_mortality_2019.csv.gz
Output : <OUT>/nhanes_cancer_design_frame.csv.gz   (OUT defaults to data/derived)
Run from the project root:  python3 code/build_cohort.py [OUT]     or  OUT=dir python3 code/build_cohort.py

Notes
* The raw files are read with exact (round-trip) float parsing, so values equal those of the
  original XPORT files. The frame is then written, re-read with pandas' default parser and
  written again with the skin-cancer flags, as in the original two-step build; this fixes the
  last digit of about 8,800 WTMEC2YR, wt and time cells in the saved file.
* The age at diagnosis for 2005-2016 is read from the MCQ240 column of each participant's first
  reported cancer. Codes 37-39 (thyroid, uterus, other) are MCQ240BB/CC/DD; the original
  interactive build used MCQ240AB/AC/AD, which do not exist, leaving yrs_dx missing for 270
  primary-cohort survivors with a recorded age. Corrected on 2026-09-25 and every analysis rerun.
"""
import io
import os
import sys

import numpy as np
import pandas as pd

RAWDIR, MORT = "data/raw_nhanes", "data/mortality/nchs_linked_mortality_2019.csv.gz"
OUT = sys.argv[1] if len(sys.argv) > 1 else os.environ.get("OUT", "data/derived")
if not os.path.isdir(RAWDIR):
    sys.exit("run from the project root")
os.makedirs(OUT, exist_ok=True)

CY = [(2005, "_D"), (2007, "_E"), (2009, "_F"), (2011, "_G"), (2013, "_H"), (2015, "_I"), (2017, "_J")]
MODULES = ["DEMO", "MCQ", "DPQ", "SLQ", "BMX", "SMQ", "BPQ", "DIQ", "PAQ"]
RAW = {(y, s): {f: pd.read_csv(f"{RAWDIR}/{f}{s}_{y}-{y + 1}.csv.gz", float_precision="round_trip")
              for f in MODULES} for y, s in CY}
mo2 = pd.read_csv(MORT)
for c in ["ELIGSTAT", "MORTSTAT", "UCOD_LEADING", "PERMTH_INT", "PERMTH_EXM"]:
    mo2[c] = pd.to_numeric(mo2[c], errors="coerce")

DPQ9 = [f"DPQ0{i}0" for i in range(1, 10)]
TYPE = {14: "Breast", 30: "Prostate", 16: "Colorectal", 31: "Colorectal", 25: "Melanoma/skin",
        32: "Melanoma/skin", 33: "Melanoma/skin", 15: "Gynecologic", 28: "Gynecologic", 38: "Gynecologic",
        11: "Hematologic", 21: "Hematologic", 24: "Hematologic", 23: "Lung"}
LETTERS = {10: "A", 11: "B", 12: "C", 13: "D", 14: "E", 15: "F", 16: "G", 17: "H", 18: "I", 19: "J",
           20: "K", 21: "L", 22: "M", 23: "N", 24: "O", 25: "P", 26: "Q", 27: "R", 28: "S", 29: "T",
           30: "U", 31: "V", 32: "W", 33: "X", 34: "Y", 35: "Z", 36: "AA", 37: "BB", 38: "CC", 39: "DD"}


def build2(y, s, d):
    dem = d["DEMO"]; mcq = d["MCQ"]
    base = dem[["SEQN", "RIAGENDR", "RIDAGEYR", "RIDRETH1", "DMDEDUC2", "DMDMARTL", "INDFMPIR",
                "SDMVPSU", "SDMVSTRA", "WTMEC2YR"]].copy()
    acols = [c for c in ["MCQ230A", "MCQ230B", "MCQ230C", "MCQ230D"] if c in mcq.columns]
    m = mcq[["SEQN", "MCQ220"] + acols].copy()
    m["surv"] = m.MCQ220 == 1
    m["n_ca"] = m[acols].apply(lambda r: r.between(10, 39).sum(), axis=1)

    def firsttype(r):
        v = [x for x in r if pd.notna(x) and 10 <= x <= 39]
        if not v:
            return np.nan
        pri = [t for t in v if t in TYPE]
        return TYPE[pri[0]] if pri else "Other"
    m["ca_type"] = m[acols].apply(firsttype, axis=1)
    if "MCQ240E" in mcq.columns:
        ad = np.full(len(mcq), np.nan)
        first = mcq[acols].apply(lambda r: next((x for x in r if pd.notna(x) and 10 <= x <= 39), np.nan), axis=1)
        for code, L in LETTERS.items():
            col = f"MCQ240{L}"
            if col in mcq.columns:
                hit = (first == code).values
                ad[hit] = mcq[col].values[hit]
        m["age_dx"] = ad
    else:
        pairs = [(a, a.replace("MCQ230", "MCD240")) for a in acols if a.replace("MCQ230", "MCD240") in mcq.columns]
        ad = np.full(len(mcq), np.nan)
        for a, dc in pairs:
            hit = mcq[a].between(10, 39).values & np.isnan(ad)
            ad[hit] = mcq[dc].values[hit]
        m["age_dx"] = ad
    m.loc[m.age_dx >= 777, "age_dx"] = np.nan
    cc = [c for c in ["MCQ160B", "MCQ160C", "MCQ160E", "MCQ160F"] if c in mcq.columns]
    m["cvd"] = mcq[cc].eq(1).any(axis=1).astype(float).values
    X = base.merge(m[["SEQN", "surv", "n_ca", "ca_type", "age_dx", "cvd"]], on="SEQN")
    dpq = d["DPQ"][["SEQN"] + DPQ9].copy()
    for c in DPQ9:
        dpq[c] = dpq[c].where(dpq[c] <= 3)
    X = X.merge(dpq, on="SEQN", how="left")
    sl = d["SLQ"].copy(); h = "SLD010H" if "SLD010H" in sl.columns else "SLD012"
    sl["sleep_h"] = sl[h].where(sl[h] <= 24); sl["slq050"] = sl.SLQ050.where(sl.SLQ050.isin([1, 2]))
    X = X.merge(sl[["SEQN", "sleep_h", "slq050"]], on="SEQN", how="left")
    X = X.merge(d["BMX"][["SEQN", "BMXBMI"]], on="SEQN", how="left")
    X = X.merge(d["SMQ"][["SEQN", "SMQ020"] + (["SMQ040"] if "SMQ040" in d["SMQ"].columns else [])], on="SEQN", how="left")
    X = X.merge(d["BPQ"][["SEQN", "BPQ020"]], on="SEQN", how="left")
    X = X.merge(d["DIQ"][["SEQN", "DIQ010"]], on="SEQN", how="left")
    pa = d["PAQ"].copy()
    if "PAQ650" in pa.columns:
        def mm(q, dc, met):
            if q not in pa.columns or dc not in pa.columns:
                return 0
            dd = pa[dc].where(pa[dc] <= 1000)
            return np.where(pa[q] == 1, np.nan_to_num(dd.values) * met, 0)
        pa["met"] = (mm("PAQ605", "PAD615", 8) + mm("PAQ620", "PAD630", 4) + mm("PAQ635", "PAD645", 4)
                     + mm("PAQ650", "PAD660", 8) + mm("PAQ665", "PAD675", 4))
        X = X.merge(pa[["SEQN", "met"]], on="SEQN", how="left")
    else:
        X["met"] = np.nan
    X["cycle"] = f"{y}-{y + 1}"
    return X


FULL = pd.concat([RAW[(y, s)]["DEMO"][["SEQN", "SDMVPSU", "SDMVSTRA", "WTMEC2YR", "RIDAGEYR", "RIAGENDR"]]
                  .assign(cycle=f"{y}-{y + 1}") for y, s in CY], ignore_index=True)
D2 = pd.concat([build2(y, s, RAW[(y, s)]) for y, s in CY], ignore_index=True)
S = D2[D2.surv].copy()

acols_all = ["MCQ230A", "MCQ230B", "MCQ230C", "MCQ230D"]
codes = {"Breast": [14], "Prostate": [30], "Colorectal": [16, 31], "Melanoma/skin": [25, 32, 33],
         "Gynecologic": [15, 28, 38], "Hematologic": [11, 21, 24], "Lung": [23]}
MC = pd.concat([RAW[(y, s)]["MCQ"][["SEQN"] + [c for c in acols_all if c in RAW[(y, s)]["MCQ"].columns]]
                for y, s in CY], ignore_index=True)
for k, v in codes.items():
    S = S.merge(MC.assign(**{f"has_{k}": MC[[c for c in acols_all if c in MC.columns]].isin(v).any(axis=1)})
                [["SEQN", f"has_{k}"]], on="SEQN", how="left")

E2 = S.merge(mo2[["SEQN", "ELIGSTAT", "MORTSTAT", "UCOD_LEADING", "PERMTH_EXM", "PERMTH_INT"]], on="SEQN", how="left")
E2 = E2[E2.ELIGSTAT == 1]
E2["time"] = E2.PERMTH_EXM.fillna(E2.PERMTH_INT) / 12
E2 = E2[(E2.time > 0) & E2.time.notna() & (E2.WTMEC2YR > 0)]
B = E2.dropna(subset=DPQ9 + ["sleep_h", "slq050"]).copy()
B["wt"] = B.WTMEC2YR / 7
B["event"] = (B.MORTSTAT == 1).astype(int)
B["ev_ca"] = ((B.MORTSTAT == 1) & (B.UCOD_LEADING == 2)).astype(int)
B["ev_cvd"] = ((B.MORTSTAT == 1) & (B.UCOD_LEADING.isin([1, 5]))).astype(int)
B["phq9_score"] = B[DPQ9].sum(axis=1)
for i, c in enumerate(DPQ9, 1):
    B[f"phqi{i}"] = (B[c] >= 2).astype(int) + 1
B["sleepcat"] = np.select([B.sleep_h < 6, B.sleep_h >= 9], [1, 3], default=2)
B["slq050b"] = np.where(B.slq050 == 1, 1, 2)
B["dep_cat"] = pd.cut(B.phq9_score, [-1, 4, 9, 27], labels=["None/minimal (0-4)", "Mild (5-9)", "Moderate+ (>=10)"])
B["race"] = B.RIDRETH1.map({1: "Mexican American", 2: "Other Hispanic", 3: "NH White", 4: "NH Black", 5: "Other/Multi"})
B["educ"] = pd.cut(B.DMDEDUC2.where(B.DMDEDUC2 <= 5), [0, 2, 3, 5], labels=["<High school", "High school", ">High school"])
B["married"] = np.select([B.DMDMARTL.isin([1, 6]), B.DMDMARTL.isin([2, 3, 4, 5])], ["Married/partnered", "Not partnered"], default=None)
B["smoke"] = np.select([B.SMQ020 == 2, B.SMQ040.isin([1, 2]), (B.SMQ020 == 1) & (B.SMQ040 == 3)], ["Never", "Current", "Former"], default=None)
B["htn"] = np.where(B.BPQ020 == 1, 1, np.where(B.BPQ020 == 2, 0, np.nan))
B["dm"] = np.where(B.DIQ010 == 1, 1, np.where(B.DIQ010.isin([2, 3]), 0, np.nan))
B["pa_active"] = np.where(B.met.isna(), np.nan, (B.met >= 600).astype(float))
B["yrs_dx"] = (B.RIDAGEYR - B.age_dx).clip(lower=0)
B["female"] = (B.RIAGENDR == 2).astype(int)
B["multi_primary"] = (B.n_ca > 1).astype(int)

keep2 = (["SEQN", "time", "event", "ev_ca", "ev_cvd", "phq9_score", "dep_cat", "sleep_h", "sleepcat", "slq050b",
          "race", "educ", "married", "smoke", "htn", "dm", "cvd", "pa_active", "yrs_dx", "multi_primary", "female",
          "ca_type", "INDFMPIR", "BMXBMI", "met", "RIDAGEYR"]
         + [f"phqi{i}" for i in range(1, 10)] + [f"has_{k}" for k in codes])
G = FULL[["SEQN", "SDMVPSU", "SDMVSTRA", "WTMEC2YR", "cycle"]].copy(); G["wt"] = G.WTMEC2YR / 7
G = G.merge(B[keep2].rename(columns={"RIDAGEYR": "age"}), on="SEQN", how="left")
G["inAnalysis"] = G.SEQN.isin(B.SEQN).astype(int)
assert G.inAnalysis.sum() == len(B)

cols = [c for c in acols_all if c in MC.columns]
mel = MC[cols].eq(25).any(axis=1)
nms = MC[cols].isin([32, 33]).any(axis=1)
other = MC[cols].apply(lambda r: any(pd.notna(x) and 10 <= x <= 39 and x not in (25, 32, 33) for x in r), axis=1)
flag = MC[["SEQN"]].assign(only_nms=(nms & ~mel & ~other).values, any_mel=mel.values)
buf = io.StringIO(); G.to_csv(buf, index=False); buf.seek(0)
gd = pd.read_csv(buf).merge(flag, on="SEQN", how="left")
out = os.path.join(OUT, "nhanes_cancer_design_frame.csv.gz")
gd.to_csv(out, index=False, compression={"method": "gzip", "mtime": 0})

sub = gd[gd.inAnalysis == 1]
skin = sub.only_nms.astype(bool)
c32 = sub.SEQN.isin(MC.SEQN[MC[cols].eq(32).any(axis=1)])
c33 = sub.SEQN.isin(MC.SEQN[MC[cols].eq(33).any(axis=1)])
print("Participant flow (Additional file 1: Fig. S1)")
print(f"  NHANES 2005-2018 participants, pooled             {len(FULL):6,d}")
print(f"  self-reported cancer (MCQ220 = 1)                 {len(S):6,d}")
print(f"  linkage-eligible, follow-up > 0, MEC weight > 0   {len(E2):6,d}")
print(f"  complete PHQ-9 items, sleep duration and SLQ050   {len(B):6,d}   (inAnalysis = 1)")
print(f"  excluded, skin cancer only                        {int(skin.sum()):6,d}")
print(f"      non-melanoma skin cancer only (code 32)       {int((skin & c32 & ~c33).sum()):6,d}")
print(f"      skin cancer of unknown type only (code 33)    {int((skin & c33 & ~c32).sum()):6,d}")
print(f"      both codes 32 and 33                          {int((skin & c32 & c33).sum()):6,d}")
print(f"  primary cohort                                    {int((~skin).sum()):6,d}")
print(f"wrote {out}  {gd.shape}")

# Survivors eligible before the symptom-completeness requirement (3,779), with their missingness and a few
# baseline variables, for the comparison of included and excluded survivors
# (code/analysis_excluded_comparison.R). The design frame above does not depend on this block.
import glob
prox = pd.concat([pd.read_csv(f, usecols=["SEQN", "MIAPROXY"]) for f in sorted(glob.glob(os.path.join(RAWDIR, "DEMO_*.csv.gz")))])
X = E2.merge(prox, on="SEQN", how="left")
X["included"] = X.SEQN.isin(B.SEQN).astype(int)
X["n_phq_missing"] = X[DPQ9].isna().sum(axis=1)
X["sleep_missing"] = (X.sleep_h.isna() | X.slq050.isna()).astype(int)
X["event"] = (X.MORTSTAT == 1).astype(int)
X["wt"] = X.WTMEC2YR / 7
X = X.merge(flag, on="SEQN", how="left")
keep3 = ["SEQN", "SDMVPSU", "SDMVSTRA", "wt", "cycle", "time", "event", "included", "n_phq_missing",
         "sleep_missing", "MIAPROXY", "RIDAGEYR", "RIAGENDR", "RIDRETH1", "DMDEDUC2", "INDFMPIR", "BMXBMI",
         "BPQ020", "DIQ010", "cvd", "age_dx", "ca_type", "only_nms"]
out2 = os.path.join(OUT, "eligible_survivors.csv.gz")
X[keep3].to_csv(out2, index=False, compression={"method": "gzip", "mtime": 0})
print(f"wrote {out2}  {X[keep3].shape}; excluded {int((X.included == 0).sum())}")

