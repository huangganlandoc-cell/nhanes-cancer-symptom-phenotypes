"""Build tables/Table4_robustness.csv and supplementary/TableS10_robustness_grid.csv
directly from the saved analysis outputs in supporting/ (no model fitting, no copied numbers).

All Model 4 rows use the harmonised 2,564-person frame (code/reconstructed/harmonised/).
Sources, by row:
  modal, Model 3 vs low          P_cox_main.csv (Model 3, all-cause)
  modal, Model 4 vs low          RR_model4.csv ("Model 4 + antidepressant use")
  modal, Model 3 vs sleep-fatigue RR_contrast_vs_insomnia_fatigue.csv
  modal, Model 4 vs sleep-fatigue RR_model4_contrast.csv
  bias-adjusted three-step       RR5_three_step_estimates.csv   (code/analysis_three_step.R)
  inclusive draws                RR5_inclusive_draws_score.csv  (code/analysis_inclusive_draws_score.R)
  indicator-only draws           sens_pseudoclass_M100.csv (Model 3), RR_model4_pseudoclass.csv,
                                 RR_contrast_pseudoclass.csv
  PHQ-9 sleep item removed       RR4_localdep_cox.csv, class L4 (the somatic-depressive profile)
  lower item threshold           RR4_threshold_cox.csv (no class associated; P range)
  classes fitted in all adults   RR_fullsample_lca_pseudoclass.csv, class C3 (somatic-depressive profile)
Run from the project root:  python3 code/build_table4.py
"""
import csv

S = "supporting/"
def rows(f): return list(csv.DictReader(open(S + f, encoding="utf-8")))
def pick(f, **kw):
    hit = [r for r in rows(f) if all(r[k] == v for k, v in kw.items())]
    assert len(hit) == 1, (f, kw, len(hit)); return hit[0]
fp = lambda p: f"{p:.1e}" if p < 0.001 else f"{p:.4f}"
def cell(r, dash="–"):
    hr, lo, hi, p = (float(r[k]) for k in ("HR", "lo", "hi", "p"))
    return f"{hr:.2f} ({lo:.2f}{dash}{hi:.2f})", fp(p), "yes" if (lo > 1 or hi < 1) else "no"

MODAL = "Modal assignment (uncorrected)"
THREE = "Bias-adjusted three-step (classes defined by symptoms only)"
INCL = "Inclusive draws (membership model includes outcome, cumulative hazard and a covariate prognostic score)"
IND = "Indicator-only draws (posteriors from the symptoms alone; attenuate towards the null)"
LOW, SF = "vs low symptom burden", "vs sleep-fatigue"
src = {
    ("modal", "Model 3", LOW): pick("P_cox_main.csv", model="Model 3", outcome="event", term="lcaHypersomnia-somatic"),
    ("modal", "Model 4", LOW): pick("RR_model4.csv", model="Model 4 + antidepressant use", term="Somatic-depressive"),
    ("modal", "Model 3", SF): pick("RR_contrast_vs_insomnia_fatigue.csv", term="Somatic-depressive"),
    ("modal", "Model 4", SF): pick("RR_model4_contrast.csv", term="Somatic-depressive"),
    ("three", "Model 3", LOW): pick("RR5_three_step_estimates.csv", model="C3", term="Somatic-depressive vs low"),
    ("three", "Model 3", SF): pick("RR5_three_step_estimates.csv", model="C3", term="Somatic-depressive vs sleep-fatigue"),
    ("three", "Model 4", LOW): pick("RR5_three_step_estimates.csv", model="C4", term="Somatic-depressive vs low"),
    ("three", "Model 4", SF): pick("RR5_three_step_estimates.csv", model="C4", term="Somatic-depressive vs sleep-fatigue"),
    ("incl", "Model 3", LOW): pick("RR5_inclusive_draws_score.csv", model="C3", term="Somatic-depressive vs low"),
    ("incl", "Model 4", LOW): pick("RR5_inclusive_draws_score.csv", model="C4", term="Somatic-depressive vs low"),
    ("ind", "Model 3", LOW): pick("sens_pseudoclass_M100.csv", cohort="primary", term="Somatic-depressive"),
    ("ind", "Model 4", LOW): pick("RR_model4_pseudoclass.csv", term="Somatic-depressive"),
    ("ind", "Model 3", SF): pick("RR_contrast_pseudoclass.csv", cov="Model 3", term="Somatic-depressive"),
    ("ind", "Model 4", SF): pick("RR_contrast_pseudoclass.csv", cov="Model 4 + antidepressant", term="Somatic-depressive"),
}
sleep = pick("RR4_localdep_cox.csv", term="L4")
allad = pick("RR_fullsample_lca_pseudoclass.csv", term="C3")
thr_p = [float(r["p"]) for r in rows("RR4_threshold_cox.csv")]

H4 = ["Domain", "Specification", "Covariate set", "Comparison", "Hazard ratio (95% CI)", "P", "Interval excludes 1"]
t4 = []
def add(domain, spec, cov, comp, r):
    h, p, x = cell(r); t4.append([domain, spec, cov, comp, h, p, x])
for cov, comp in [("Model 3", LOW), ("Model 4", LOW), ("Model 3", SF), ("Model 4", SF)]:
    add("Class assignment", MODAL, cov, comp, src[("modal", cov, comp)])
for cov in ("Model 3", "Model 4"):
    for comp in (LOW, SF):
        add("Class assignment", THREE, cov, comp, src[("three", cov, comp)])
for cov in ("Model 3", "Model 4"):
    add("Class assignment", INCL, cov, LOW, src[("incl", cov, LOW)])
for cov, comp in [("Model 3", LOW), ("Model 4", LOW), ("Model 4", SF)]:
    add("Class assignment", IND, cov, comp, src[("ind", cov, comp)])
add("Measurement decisions", "PHQ-9 sleep item removed (overlaps the two sleep indicators)", "Model 4", "vs largest class", sleep)
t4.append(["Measurement decisions", 'Items dichotomised at "several days" rather than "more than half the days"', "Model 4",
           "vs largest class", "no class associated", f"{min(thr_p):.3f} to {max(thr_p):.3f}", "no"])
add("Measurement decisions", "Measurement model refitted in all 34,022 adults", "Model 4", "vs largest class", allad)
with open("tables/Table4_robustness.csv", "w", newline="", encoding="utf-8") as fh:
    csv.writer(fh).writerows([H4] + t4)
print(f"Table 4: {len(t4)} rows, {sum(r[6] == 'no' for r in t4)} with an interval including one")

# ---------------------------------------------------------------- Table S10 grid
H10 = ["Covariate set", "Class assignment", "Comparison", "HR (95% CI)", "P", "Interval excludes 1"]
lab = {"modal": "modal class assignment", "ind": "100 pseudo-class draws (symptom-based posteriors)",
       "three": "bias-adjusted three-step estimator"}
covlab = {"Model 3": "Model 3 (demographic, behavioural, clinical)",
          "Model 4": "Model 4 (+ comorbidity, frailty, medication, antidepressant)"}
t10 = []
for comp in (LOW, SF):
    for cov in ("Model 3", "Model 4"):
        for k in ("modal", "ind", "three"):
            h, p, x = cell(src[(k, cov, comp)], dash="-")
            t10.append([covlab[cov], lab[k], comp, h, p, x])
with open("supplementary/TableS10_robustness_grid.csv", "w", newline="", encoding="utf-8") as fh:
    csv.writer(fh).writerows([H10] + t10)
print(f"Table S10: {len(t10)} rows")
for r in t4:
    print(f"  {r[1][:34]:34s} {r[2]} {r[3][:22]:22s} {r[4]:20s} {r[5]}")
