"""Three derived inputs built from the raw NHANES files, as the original analysis built them.

  data/derived/nhanes_all_adults_lca_input.csv.gz  adults aged 20+ with complete PHQ-9 items, sleep duration
                                                   and SLQ050 (the measurement model refitted in all adults:
                                                   code/reconstructed/all_adult_lca.R, dimensional.R)
  data/derived/nhanes_dimension_scores.csv         somatic and affective dimension scores of the 3,268
                                                   survivors and the somatic residual at a fixed total score
                                                   (dimensional.R, the harmonised scripts, Table S7)
  data/derived/nhanes_dpq_raw.csv                  raw 0-3 PHQ-9 item scores of the 3,268 survivors
                                                   (code/reconstructed/measurement_sensitivity_RR4.R)

The original analysis wrote these files with pandas in its interactive session; the code below repeats
those steps line for line (steps 276, 278 and 407 of the analysis log). The raw files are read with
exact (round-trip) float parsing, so the values equal those of the XPORT files the original read, and
the somatic residual comes from the same statsmodels OLS fit. An existing file is never overwritten: it
is compared with the re-derived table instead (text of every line, then values), and the result printed.

Run from the project root after code/build_cohort.py:  python3 code/build_derived_inputs.py
Requires pandas, numpy and statsmodels.
"""
import gzip
import io
import os

import numpy as np
import pandas as pd
import statsmodels.api as sm

RAWDIR, DER = "data/raw_nhanes", "data/derived"
CY = [(2005, "_D"), (2007, "_E"), (2009, "_F"), (2011, "_G"), (2013, "_H"), (2015, "_I"), (2017, "_J")]
RAW = {(y, s): {f: pd.read_csv(f"{RAWDIR}/{f}{s}_{y}-{y + 1}.csv.gz", float_precision="round_trip")
                for f in ("DEMO", "DPQ", "SLQ")} for y, s in CY}
GD = pd.read_csv(f"{DER}/nhanes_cancer_design_frame.csv.gz")
DPQ9 = [f'DPQ0{i}0' for i in range(1, 10)]
SOM = ['DPQ030', 'DPQ040', 'DPQ050', 'DPQ080']                  # sleep, fatigue, appetite, psychomotor
AFF = ['DPQ010', 'DPQ020', 'DPQ060', 'DPQ070', 'DPQ090']        # anhedonia, mood, worthlessness, concentration, ideation

# ---- step 276: all adults with complete symptom data
rows = []
for y, s in CY:
    dem = RAW[(y, s)]['DEMO']; dpq = RAW[(y, s)]['DPQ']; slq = RAW[(y, s)]['SLQ']
    d = dem[['SEQN', 'RIDAGEYR', 'RIAGENDR', 'SDMVPSU', 'SDMVSTRA', 'WTMEC2YR']].copy()
    q = dpq[['SEQN'] + DPQ9].copy()
    for c in DPQ9: q[c] = q[c].where(q[c] <= 3)
    h = 'SLD010H' if 'SLD010H' in slq.columns else 'SLD012'
    sl = slq[['SEQN', h, 'SLQ050']].rename(columns={h: 'sleep_h'})
    sl['sleep_h'] = sl.sleep_h.where(sl.sleep_h <= 24); sl['slq050'] = sl.SLQ050.where(sl.SLQ050 <= 2)
    m = d.merge(q, on='SEQN', how='inner').merge(sl[['SEQN', 'sleep_h', 'slq050']], on='SEQN', how='inner')
    m['cycle'] = f"{y}-{y + 1}"; rows.append(m)
ALLAD = pd.concat(rows, ignore_index=True)
ALLAD = ALLAD[(ALLAD.RIDAGEYR >= 20) & ALLAD[DPQ9].notna().all(axis=1) & ALLAD.sleep_h.notna() & ALLAD.slq050.notna()].copy()
ALLAD['phq9_score'] = ALLAD[DPQ9].sum(axis=1)
for i, c in enumerate(DPQ9, 1): ALLAD[f'phqi{i}'] = (ALLAD[c] >= 2).astype(int) + 1
ALLAD['sleepcat'] = np.select([ALLAD.sleep_h < 6, ALLAD.sleep_h >= 9], [1, 3], default=2)
ALLAD['slq050b'] = np.where(ALLAD.slq050 == 1, 1, 2)
ALLAD['wt'] = ALLAD.WTMEC2YR / 7
surv = set(GD.query('inAnalysis==1').SEQN)
ALLAD['is_survivor'] = ALLAD.SEQN.isin(surv).astype(int)

# ---- step 278, first half: dimension scores and the somatic residual at a fixed total score
sub = GD[GD.inAnalysis == 1][['SEQN']]
dq = []
for y, s in CY:
    d = RAW[(y, s)]['DPQ'][['SEQN'] + DPQ9].copy()
    for c in DPQ9: d[c] = d[c].where(d[c] <= 3)
    dq.append(d)
DQ = pd.concat(dq, ignore_index=True)
DQ['som'] = DQ[SOM].sum(axis=1, min_count=4); DQ['aff'] = DQ[AFF].sum(axis=1, min_count=5)
DQ['tot'] = DQ[DPQ9].sum(axis=1, min_count=9)
DIM = sub.merge(DQ[['SEQN', 'som', 'aff', 'tot']], on='SEQN', how='left').dropna()
mod = sm.OLS(DIM.som, sm.add_constant(np.column_stack([DIM.tot, DIM.tot ** 2]))).fit()
DIM['som_resid'] = mod.resid

# ---- step 407: raw item scores of the survivors
DQ2 = DQ[['SEQN'] + DPQ9].copy()
for i, c in enumerate(DPQ9, 1): DQ2[f'dpq{i}_raw'] = DQ2[c]
RAWSC = sub.merge(DQ2[['SEQN'] + [f'dpq{i}_raw' for i in range(1, 10)]], on='SEQN', how='left')


def put(df, path):
    """Write df as pandas does, unless the file exists; then compare it with df instead."""
    text = df.to_csv(index=False)
    if not os.path.exists(path):
        with (gzip.open(path, "wt", newline="") if path.endswith(".gz") else open(path, "w", newline="")) as fh:
            fh.write(text)
        print(f"wrote {path}  {df.shape}")
        return
    with (gzip.open(path, "rt", newline="") if path.endswith(".gz") else open(path, newline="")) as fh:
        saved = fh.read()
    if saved == text:
        print(f"{path} exists and is identical, line for line, to the re-derived table {df.shape}")
        return
    old = pd.read_csv(io.StringIO(saved), float_precision="round_trip")
    new = pd.read_csv(io.StringIO(text), float_precision="round_trip")
    same = list(old.columns) == list(new.columns) and old.shape == new.shape
    diff = {c: int((~np.isclose(old[c], new[c], rtol=1e-9, atol=0, equal_nan=True)).sum())
            if pd.api.types.is_numeric_dtype(old[c]) else int((old[c].astype(str) != new[c].astype(str)).sum())
            for c in old.columns} if same else "shape or columns differ"
    nlines = sum(a != b for a, b in zip(saved.splitlines(), text.splitlines()))
    print(f"{path} exists and was kept; it differs from the re-derived table in {nlines} lines; "
          f"cells beyond 1e-9 relative by column: {diff}")


put(ALLAD, f"{DER}/nhanes_all_adults_lca_input.csv.gz")
put(DIM, f"{DER}/nhanes_dimension_scores.csv")
put(RAWSC, f"{DER}/nhanes_dpq_raw.csv")
print(f"all adults {len(ALLAD):,} (survivors {int(ALLAD.is_survivor.sum()):,}) | dimension scores {len(DIM):,} "
      f"(R2 of som on tot {mod.rsquared:.3f}) | raw item scores {len(RAWSC):,}")
