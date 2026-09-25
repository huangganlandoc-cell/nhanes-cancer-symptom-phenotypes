"""Supplementary Table S14: Model 4a and the measurement model without SLQ050 (post hoc, 2026-09-24).

Reads the outputs of code/analysis_model4a.R and code/analysis_drop_slq050.R and writes
supplementary/TableS14_model4a_drop_slq050.csv. Run from the project root after both scripts.
"""
import csv

S = "supporting/"
fp = lambda p: f"{p:.1e}" if p < 0.001 else f"{p:.4f}"
hr = lambda r: f"{float(r['HR']):.2f} ({float(r['lo']):.2f}-{float(r['hi']):.2f})"
rows = lambda f: list(csv.DictReader(open(S + f, encoding="utf-8")))

out = []
def add(block, cov, assign, comp, r):
    out.append({"Analysis": block, "Covariate set": cov, "Class assignment": assign, "Comparison": comp,
                "HR (95% CI)": hr(r), "P": fp(float(r["p"])),
                "Interval excludes 1": "yes" if float(r["lo"]) > 1 or float(r["hi"]) < 1 else "no"})

M4A = "Model 4a (Model 3 + condition count and eGFR)"
comp = {"Sleep-fatigue vs low": "sleep-fatigue vs low symptom burden",
        "Somatic-depressive vs low": "somatic-depressive vs low symptom burden",
        "High symptom burden vs low": "high symptom burden vs low symptom burden",
        "Somatic-depressive vs sleep-fatigue": "somatic-depressive vs sleep-fatigue"}
for f, assign in [("RR6_model4a_modal.csv", "modal class assignment"),
                  ("RR6_model4a_three_step.csv", "bias-adjusted three-step estimator (JKn)"),
                  ("RR6_model4a_three_step_bootstrap.csv", "bias-adjusted three-step estimator (bootstrap)")]:
    for r in rows(f):
        add("A. Model 4a, main measurement model", M4A, assign, comp[r["term"]], r)

covname = {"Model 3": "Model 3 (demographic, behavioural, clinical)", "Model 4a": M4A,
           "Model 4": "Model 4 (+ comorbidity, frailty, medication, antidepressant)"}
for r in rows("RR6_drop_slq050_cox.csv"):
    add("B. Measurement model refitted without SLQ050 (ten indicators, four classes)", covname[r["model"]],
        "modal class assignment", r["term"].replace("Sleep-fatigue", "sleep-fatigue")
        .replace("Somatic-depressive", "somatic-depressive").replace("High symptom burden", "high symptom burden"), r)

fit = rows("RR6_drop_slq050_fit.csv")
f4 = next(r for r in fit if r["classes"] == "4")
note = (f"B: ten-indicator solutions with 3, 4 and 5 classes had BIC "
        + ", ".join(f"{float(r['BIC']):.0f}" for r in fit)
        + f"; four-class entropy {float(f4['entropy']):.3f}, Cramer's V with the main assignment "
        f"{float(f4['cramer_v_4class']):.3f}, same class for {100 * float(f4['agreement_4class']):.1f}% of survivors.")
with open("supplementary/TableS14_model4a_drop_slq050.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(out[0]))
    w.writeheader(); w.writerows(out)
print(f"wrote supplementary/TableS14_model4a_drop_slq050.csv ({len(out)} rows)")
print(note)
