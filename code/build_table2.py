"""Build the main Table 2 and Supplementary Table S13 from the saved outputs (no model fitting).

Table 2 (key estimates): all-cause mortality by phenotype under the two covariate sets
(Model 3, specified; Model 4, post hoc) and the two ways of assigning class membership
(modal assignment; bias-adjusted three-step correction for misclassification), plus the
direct somatic-depressive versus sleep-fatigue contrast.
Table S13 (model progression): Models 1-4 with modal assignment, all-cause and
cause-specific, with Benjamini-Hochberg P values across the nine Model 3 tests of the
primary family (the former Table 2).

Sources (supporting/):
  P_cox_main.csv                         Models 1-3, three outcomes (analysis_primary_exclNMS.R)
  RR_model4.csv                          Model 4 = Model 3 + comorbidity/frailty + antidepressant use
  RR_contrast_vs_insomnia_fatigue.csv    Model 3 contrast, modal (analysis_reviewer_response.R)
  RR_model4_contrast.csv                 Model 4 contrast, modal
  RR5_three_step_estimates.csv           bias-adjusted three-step, Models 3 and 4 (analysis_three_step.R)
Run from the project root:  python3 code/build_table2.py
"""
import csv

S = "supporting/"
def rows(f): return list(csv.DictReader(open(S + f, encoding="utf-8")))
LAB = {"lcaInsomnia-fatigue": "Sleep-fatigue", "Insomnia-fatigue": "Sleep-fatigue",
       "lcaHypersomnia-somatic": "Somatic-depressive", "Somatic-depressive": "Somatic-depressive",
       "lcaHigh symptom burden": "High symptom burden", "High symptom burden": "High symptom burden"}
OUT = {"event": "All-cause", "ev_ca": "Cancer", "ev_cvd": "Cardiovascular or stroke"}
PHEN = ["Sleep-fatigue", "Somatic-depressive", "High symptom burden"]
hr = lambda r: f"{float(r['HR']):.2f} ({float(r['lo']):.2f}–{float(r['hi']):.2f})"
pfmt = lambda p: f"{p:.1e}" if p < 0.001 else f"{p:.3f}"


def bh(ps):
    m = len(ps); order = sorted(range(m), key=lambda i: ps[i], reverse=True)
    adj, prev = [0.0] * m, 1.0
    for rank, i in zip(range(m, 0, -1), order):
        prev = min(prev, ps[i] * m / rank); adj[i] = prev
    return adj


main = {(r["model"], OUT[r["outcome"]], LAB[r["term"]]): r for r in rows("P_cox_main.csv")}
m4 = {LAB[r["term"]]: r for r in rows("RR_model4.csv") if r["model"] == "Model 4 + antidepressant use"}
c3 = {r["term"]: r for r in rows("RR_contrast_vs_insomnia_fatigue.csv")}["Somatic-depressive"]
c4 = {r["term"]: r for r in rows("RR_model4_contrast.csv")}["Somatic-depressive"]
T = {(r["model"], r["term"]): r for r in rows("RR5_three_step_estimates.csv")}
TT = {"Sleep-fatigue": "Sleep-fatigue vs low", "Somatic-depressive": "Somatic-depressive vs low",
      "High symptom burden": "High symptom burden vs low"}

# ---------------------------------------------------------------- Table 2
t2 = [["Phenotype (reference: low symptom burden)",
       "Model 3, modal: HR (95% CI)", "P", "Model 3, corrected: HR (95% CI)", "P",
       "Model 4, modal: HR (95% CI)", "P", "Model 4, corrected: HR (95% CI)", "P"]]
for p in PHEN:
    a, b = main[("Model 3", "All-cause", p)], T[("C3", TT[p])]
    c, d = m4[p], T[("C4", TT[p])]
    t2.append([p, hr(a), pfmt(float(a["p"])), hr(b), pfmt(float(b["p"])),
               hr(c), pfmt(float(c["p"])), hr(d), pfmt(float(d["p"]))])
b, d = T[("C3", "Somatic-depressive vs sleep-fatigue")], T[("C4", "Somatic-depressive vs sleep-fatigue")]
t2.append(["Somatic-depressive vs sleep-fatigue", hr(c3), pfmt(float(c3["p"])), hr(b), pfmt(float(b["p"])),
           hr(c4), pfmt(float(c4["p"])), hr(d), pfmt(float(d["p"]))])
with open("tables/Table2_main_cox.csv", "w", newline="", encoding="utf-8") as fh:
    csv.writer(fh).writerows(t2)

# ---------------------------------------------------------------- Table S13 (former Table 2)
keys = [(o, p) for o in OUT.values() for p in PHEN]
bh3 = dict(zip(keys, bh([float(main[("Model 3", o, p)]["p"]) for o, p in keys])))
s13 = [["Outcome", "Phenotype", "Model 1 HR (95% CI)", "Model 2 HR (95% CI)", "Model 3 HR (95% CI)",
        "Model 3 P", "Model 3 BH-adjusted P", "Model 4 HR (95% CI)", "Model 4 P"]]
for o in OUT.values():
    s13.append([o, "Low symptom burden (reference)", "1.00", "1.00", "1.00", "", "", "1.00" if o == "All-cause" else "", ""])
    for p in PHEN:
        r1, r2, r3 = (main[(f"Model {k}", o, p)] for k in (1, 2, 3))
        row = ["", p, hr(r1), hr(r2), hr(r3), pfmt(float(r3["p"])), pfmt(bh3[(o, p)])]
        row += [hr(m4[p]), f"{float(m4[p]['p']):.3f}"] if o == "All-cause" else ["", ""]
        s13.append(row)
with open("supplementary/TableS13_model_progression.csv", "w", newline="", encoding="utf-8") as fh:
    csv.writer(fh).writerows(s13)

# values quoted in the text must still be what the table shows
chk = {r[1]: r for r in s13[1:]}
assert main[("Model 3", "All-cause", "Somatic-depressive")]["HR"].startswith("1.80")
assert pfmt(bh3[("Cardiovascular or stroke", "Somatic-depressive")]) == "0.053"
print("wrote tables/Table2_main_cox.csv and supplementary/TableS13_model_progression.csv")
for r in t2:
    print("  " + " | ".join(r))
