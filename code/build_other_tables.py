"""Build the eight display tables that are formatted directly from saved outputs (no model fitting).

  tables/Table1_weighted_baseline.csv                   P_table1_weighted_baseline.csv
  tables/Table3_single_symptom_exposures.csv            P_multiplicity_secondary.csv
  supplementary/TableS1_lca_fit_indices.csv             lca_fit_indices_cancer.csv (copied)
  supplementary/TableS2_subgroups_full.csv              P_cox_subgroups.csv (copied)
  supplementary/TableS4_reviewer_response_analyses.csv  RR_events_per_phenotype.csv, RR_contrast_vs_insomnia_fatigue.csv,
                                                        RR_landmark_full.csv, RR_interaction_time_since_dx.csv,
                                                        RR_cycle_sensitivity.csv; deaths by time since
                                                        diagnosis from P_cox_subgroups.csv
  supplementary/TableS6_vs_prior_classification.csv     RR_crosstab_vs_prior.csv, RR_vs_prior_classification.csv
  supplementary/TableS8_all_adult_lca_profiles.csv      RR_fullsample_lca_profiles.csv
  supplementary/TableS9_incremental_value.csv           RR_incremental_value.csv, RR_incremental_joint.csv
All sources are in supporting/. The tables keep the labels of the R outputs; code/fix_display_labels.py
is then run on the output directory, as for every other table (with the default OUTDIR it also
re-checks the other tables in the project; it is idempotent).
Numbers are read with pandas' default parser, as when the tables were first made: the full-precision
joint P values in Table S9 are printed as pandas parses them.

Run from the project root:  python3 code/build_other_tables.py [OUTDIR]
  OUTDIR (argument, or environment variable OUTDIR; default ".") receives tables/ and supplementary/.
"""
import csv
import os, shutil, subprocess, sys
import pandas as pd

S = "supporting/"
OUT = sys.argv[1] if len(sys.argv) > 1 else os.environ.get("OUTDIR", ".")
for d in ("tables", "supplementary"):
    os.makedirs(os.path.join(OUT, d), exist_ok=True)
rd = lambda f, **kw: pd.read_csv(S + f, **kw)
hr = lambda r: f"{r.HR:.2f} ({r.lo:.2f}-{r.hi:.2f})"
p4 = lambda v: f"{v:.4f}" if v >= 1e-4 else f"{v:.1e}"
p3 = lambda v: f"{v:.3f}" if v >= 0.001 else f"{v:.1e}"
def save(df, path, index=False):
    df.to_csv(os.path.join(OUT, path), index=index, lineterminator="\n")
    print("wrote", path, flush=True)

# ---------------------------------------------------------------- Table 1
V = {"age": "Age, years", "phq9_score": "PHQ-9 total score", "sleep_h": "Sleep duration, hours",
     "bmi_i": "Body mass index, kg/m2", "ydx_i": "Years since cancer diagnosis", "sex": "Sex",
     "race4": "Race and ethnicity", "educ": "Education", "married": "Marital status",
     "smoke": "Smoking status", "pa3": "Physical activity", "catype": "First primary cancer site",
     "htn": "Hypertension", "dm": "Diabetes", "cvd": "Cardiovascular disease",
     "multi_primary": "More than one primary cancer"}
t1 = rd("P_table1_weighted_baseline.csv")
t1["variable"] = t1.variable.map(lambda v: V.get(v, v) if isinstance(v, str) else "")
t1["p"] = t1.p.map(lambda v: "" if pd.isna(v) else p3(v))
save(t1.rename(columns={"variable": "Characteristic", "level": "Level", "p": "P"}),
     "tables/Table1_weighted_baseline.csv")

# ---------------------------------------------------------------- Table 3
LBL = {"dep_catMild (5-9)": ("Mild depressive symptoms (PHQ-9 5 to 9)", "No depressive symptoms (PHQ-9 < 5)"),
       "dep_catModerate+ (>=10)": ("Moderate or severe depressive symptoms (PHQ-9 ≥ 10)",
                                   "No depressive symptoms (PHQ-9 < 5)"),
       "sleep3<6 h": ("Short sleep (< 6 h)", "6 to < 9 h"), "sleep3>=9 h": ("Long sleep (≥ 9 h)", "6 to < 9 h"),
       "insomniaYes": ("Previously reported sleep trouble", "No"),
       "I(phq9_score/5)": ("PHQ-9 total score, per 5 points", "per 5-point increment")}
t3 = rd("P_multiplicity_secondary.csv")
save(pd.DataFrame({"Exposure (entered as a single variable)": t3.term.map(lambda t: LBL[t][0]),
                   "Reference": t3.term.map(lambda t: LBL[t][1]),
                   "Hazard ratio (95% CI)": [f"{r.HR:.2f} ({r.lo:.2f}–{r.hi:.2f})" for r in t3.itertuples()],
                   "P": t3.p.map(p3), "BH-adjusted P": t3.p_bh.map(lambda v: f"{v:.3f}")}),
     "tables/Table3_single_symptom_exposures.csv")

# ---------------------------------------------------------------- Tables S1, S2
for src, dst in [("lca_fit_indices_cancer.csv", "supplementary/TableS1_lca_fit_indices.csv"),
                 ("P_cox_subgroups.csv", "supplementary/TableS2_subgroups_full.csv")]:
    shutil.copyfile(S + src, os.path.join(OUT, dst)); print("wrote", dst, flush=True)

# ---------------------------------------------------------------- Table S4
ev, ct = rd("RR_events_per_phenotype.csv"), rd("RR_contrast_vs_insomnia_fatigue.csv")
lnd, itx, cyc = rd("RR_landmark_full.csv"), rd("RR_interaction_time_since_dx.csv"), rd("RR_cycle_sensitivity.csv")
sg = rd("P_cox_subgroups.csv").drop_duplicates("subgroup").set_index("subgroup").events
dx = f"{sg['<5 yr since dx']} deaths <5 yr, {sg['>=5 yr since dx']} deaths >=5 yr"
s4 = pd.concat([
    pd.DataFrame({"Block": "A. Events and person-years by phenotype", "Row": ev.phenotype,
                  "Estimate": [f"{int(r.deaths)} deaths / {int(r.person_years)} py" for r in ev.itertuples()],
                  "Detail": [f"crude {r.crude_rate_per_1000py}/1000 py; weighted {r.weighted_pct_dead}% dead (SE {r.se})"
                             for r in ev.itertuples()], "P": ""}),
    pd.DataFrame({"Block": "B. Direct contrast, insomnia-fatigue as reference", "Row": ct.term,
                  "Estimate": ct.apply(hr, axis=1), "Detail": "fully adjusted (Model 3)",
                  "P": ct.p.map(lambda v: f"{v:.4f}")}),
    pd.DataFrame({"Block": "C. Landmark exclusion series",
                  "Row": lnd.exclude_yr.astype(str) + " yr | " + lnd.term.str.replace("lca", "", regex=False),
                  "Estimate": lnd.apply(hr, axis=1),
                  "Detail": [f"n={int(r.n)}, {int(r.events)} deaths, {int(r.person_years)} py" for r in lnd.itertuples()],
                  "P": lnd.p.map(p4)}),
    pd.DataFrame({"Block": "D. Phenotype x time since diagnosis (>=5 yr vs <5 yr)", "Row": itx.term,
                  "Estimate": [f"{r.ratio:.2f} ({r.lo:.2f}-{r.hi:.2f})" for r in itx.itertuples()],
                  "Detail": f"joint Wald P = {itx.joint_p.iloc[0]:.3f}; {dx}",
                  "P": itx.p.map(lambda v: f"{v:.4f}")}),
    pd.DataFrame({"Block": "E. Survey cycle", "Row": cyc.spec + " | " + cyc.term.str.replace("lca", "", regex=False),
                  "Estimate": cyc.apply(hr, axis=1),
                  "Detail": "cycle not identifiable as a covariate (nested in design strata)", "P": cyc.p.map(p4)}),
], ignore_index=True)
save(s4, "supplementary/TableS4_reviewer_response_analyses.csv")

# ---------------------------------------------------------------- Table S6
ov, vp = rd("RR_crosstab_vs_prior.csv", index_col=0), rd("RR_vs_prior_classification.csv")
cells = [c for c in ov.columns if "/" in c]
blocks = [pd.DataFrame({"Block": "A. Cross-tabulation against the prior six-cell classification", "Row": ov.index,
                        "Estimate": [" | ".join(f"{c}: {int(ov.loc[i, c])}" for c in cells) for i in ov.index],
                        "Detail": [f"total {int(ov.loc[i, 'total'])}; {ov.loc[i, 'pct_with_complaint']}% report a sleep complaint"
                                   for i in ov.index], "P": ""})]
jp = {}
for spec, lab in [("Lan-style joint classification alone", "B. Prior classification fitted alone"),
                  ("phenotype, adjusted for Lan-style classification", "C. Phenotype adjusted for the prior classification")]:
    g = vp[vp.spec == spec]
    jp[spec] = g.joint_p.iloc[0]
    blocks.append(pd.DataFrame({"Block": lab, "Row": g.term.str.replace("lan", "", regex=False).str.replace("lca", "", regex=False),
                                "Estimate": g.apply(hr, axis=1), "Detail": f"design-based joint Wald P = {g.joint_p.iloc[0]}",
                                "P": g.p.map(lambda v: f"{v:.4f}")}))
LAN_WITH_PHENOTYPE = float(next(r["joint_p"] for r in csv.DictReader(open("supporting/RR_prior_joint_tests.csv"))
                                if r["test"] == "Lan-style classification with phenotype"))   # code/reconstructed/prior_classification.R
blocks.append(pd.DataFrame({"Block": "D. Joint tests with both blocks in one model",
                            "Row": ["phenotype block", "prior classification block"], "Estimate": ["", ""],
                            "Detail": ["absorbs the prior classification", "loses significance when phenotype is present"],
                            "P": [f"{jp['phenotype, adjusted for Lan-style classification']:.4f}",
                                  f"{LAN_WITH_PHENOTYPE:.4f} ({jp['Lan-style joint classification alone']:.4f} when fitted alone)"]}))
save(pd.concat(blocks, ignore_index=True), "supplementary/TableS6_vs_prior_classification.csv")

# ---------------------------------------------------------------- Table S8
save(rd("RR_fullsample_lca_profiles.csv", index_col=0), "supplementary/TableS8_all_adult_lca_profiles.csv", index=True)

# ---------------------------------------------------------------- Table S9
inc = rd("RR_incremental_value.csv")
inc["term"] = inc.term.str.replace("lca", "", regex=False)
s9 = inc.merge(rd("RR_incremental_joint.csv"), on="model", how="left")
s9["HR (95% CI)"] = s9.apply(hr, axis=1)
s9["P"] = s9.p.map(p4)
save(s9[["model", "term", "HR (95% CI)", "P", "joint_p_phenotype", "joint_p_11_indicators"]],
     "supplementary/TableS9_incremental_value.csv")

# ---------------------------------------------------------------- display labels
subprocess.run([sys.executable, os.path.join(os.path.dirname(os.path.abspath(__file__)), "fix_display_labels.py")],
               cwd=OUT, check=True)
