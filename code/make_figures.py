"""Figures 1 and 2 and Additional file 1: Figs. S1-S3, drawn from the saved outputs in supporting/.

Figure 1: item-endorsement profiles and weighted prevalence of the four phenotypes
          (supporting/lca_profiles_cancer.csv, supporting/P_lca_prevalence.csv).
Figure 2: survey-weighted Kaplan-Meier curves (supporting/P_km_curves.csv).
Fig. S1 : participant flow. Counts as printed by code/build_cohort.py; deaths and person-years of
          the primary cohort from supporting/RR_landmark_full.csv (no follow-up excluded).
Fig. S2 : phenotype hazard ratios within subgroups (supporting/P_cox_subgroups.csv); black outline
          where the BH-adjusted P < 0.05.
Fig. S3 : restricted cubic splines (supporting/P_rcs_curves.csv), drawn between the 1st and 99th
          percentiles of the primary cohort (data/derived/nhanes_cancer_design_frame.csv.gz);
          the curves are also written out as its source data, FigureS3_source_data.csv.
Figure 3 is drawn by code/make_figure3.py.

Outputs: figures/Figure1_symptom_phenotypes.png, figures/Figure2_survival_curves.png,
         supplementary/FigureS1_participant_flow.png, supplementary/FigureS2_subgroups.png,
         supplementary/FigureS3_dose_response.png, supplementary/FigureS3_source_data.csv
         (all six into OUT when OUT is given)
Run from the project root:  python3 code/make_figures.py [OUT]     or  OUT=dir python3 code/make_figures.py
"""
import os
import shutil
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
import pandas as pd

S = "supporting/"
if not os.path.isdir(S):
    sys.exit("run from the project root")
OUT = sys.argv[1] if len(sys.argv) > 1 else os.environ.get("OUT")
FIG, SUP = (OUT, OUT) if OUT else ("figures", "supplementary")
for d in {FIG, SUP}:
    os.makedirs(d, exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8,
    "legend.fontsize": 7, "xtick.labelsize": 6, "ytick.labelsize": 6, "axes.linewidth": 0.6,
    "xtick.direction": "out", "ytick.direction": "out", "xtick.major.size": 3, "ytick.major.size": 3,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "axes.spines.top": False,
    "axes.spines.right": False, "axes.spines.left": True, "axes.spines.bottom": True,
    "axes.grid": False, "legend.frameon": False, "figure.dpi": 200, "savefig.dpi": 300,
    "savefig.bbox": "tight", "axes.titleweight": "normal", "axes.titlelocation": "left",
    "axes.labelweight": "normal", "lines.linewidth": 1.2, "patch.linewidth": 0.6,
    "pdf.fonttype": 42, "ps.fonttype": 42})
META_GREY = "#888888"
LB = ["Low symptom burden", "Sleep-fatigue", "Somatic-depressive", "High symptom burden"]
CL = dict(zip(LB, ["#9aa3ad", "#2c7fb8", "#c0392b", "#e08214"]))
PH3 = ["Sleep-fatigue", "Somatic-depressive", "High symptom burden"]
RN = {"Insomnia-fatigue": "Sleep-fatigue", "Hypersomnia-somatic": "Somatic-depressive"}

pf = pd.read_csv(S + "lca_profiles_cancer.csv"); pv = pd.read_csv(S + "P_lca_prevalence.csv")
km = pd.read_csv(S + "P_km_curves.csv"); sg = pd.read_csv(S + "P_cox_subgroups.csv")
lm = pd.read_csv(S + "RR_landmark_full.csv"); rc = pd.read_csv(S + "P_rcs_curves.csv")
pf["class"] = pf["class"].replace(RN); pv["class"] = pv["class"].replace(RN)
km["lca"] = km["lca"].replace(RN)
sg["ph"] = sg.term.str.replace("lca", "", regex=False).replace(RN)

# ------------------------------------------------------------------ Figure 1
pf["item"] = pf["item"].replace({"Insomnia complaint": "Ever reported sleep trouble"})
order = ["Anhedonia", "Depressed mood", "Worthlessness", "Concentration", "Psychomotor change",
         "Suicidal ideation", "Appetite change", "Fatigue", "Sleep disturbance",
         "Short sleep (<6 h)", "Long sleep (>=9 h)", "Ever reported sleep trouble"]
PV = {r["class"]: (r.weighted_pct, r.lo, r.hi) for _, r in pv.iterrows()}

fig1, axes = plt.subplots(1, 2, figsize=(7.2, 4.0), gridspec_kw={"width_ratios": [2.4, 1], "wspace": 0.42})
ax = axes[0]; y = np.arange(len(order))[::-1]
for L in LB:
    v = pf[pf["class"] == L].set_index("item").reindex(order)["prob"].values
    ax.plot(v, y, "-o", color=CL[L], ms=3.6, lw=1.7 if L == "Somatic-depressive" else 1.3,
            label=L, zorder=4 if L == "Somatic-depressive" else 2,
            alpha=1.0 if L in ("Somatic-depressive", "Sleep-fatigue") else .75)
ax.axhline(y[8] - 0.5, color=META_GREY, lw=.7, ls=":")
ax.set_yticks(y); ax.set_yticklabels(order); ax.set_xlim(-3, 103); ax.margins(y=0.05)
ax.set_xlabel("Item endorsement probability (%)")
ax.set_title("Four symptom phenotypes among US cancer survivors", loc="left")
ax.legend(frameon=False, ncol=2, handlelength=1.3, columnspacing=1.2,
          loc="upper center", bbox_to_anchor=(0.5, -0.155))
ax2 = axes[1]
for i, L in enumerate(LB[::-1]):
    v, lo, hi = PV[L]
    ax2.plot([lo, hi], [i, i], color=CL[L], lw=1.4); ax2.plot(v, i, "o", color=CL[L], ms=5)
    ax2.text(v + 4, i, f"{v:.1f}%", ha="left", va="center", color=CL[L])
ax2.set_yticks(range(4)); ax2.set_yticklabels([]); ax2.set_xlim(0, 95); ax2.set_ylim(-0.5, 3.5)
ax2.set_xlabel("Weighted prevalence (%)")
ax2.set_title("Population prevalence", loc="left")
fig1.savefig(f"{FIG}/Figure1_symptom_phenotypes.png", dpi=300, bbox_inches="tight")

# ------------------------------------------------------------------ Figure 2
fig2, ax = plt.subplots(figsize=(3.7, 3.2))
for L in LB:
    k = km[km.lca == L].sort_values("time")
    ax.step(k.time, k.surv * 100, where="post", color=CL[L],
            lw=2.0 if L == "Somatic-depressive" else 1.4,
            label=L, zorder=4 if L == "Somatic-depressive" else 2)
ax.set_xlabel("Years since examination"); ax.set_ylabel("Survival (%)")
ax.set_xlim(0, 15); ax.set_ylim(40, 101); ax.margins(0.02)
ax.set_title("The somatic-depressive phenotype has\nthe lowest survival", loc="left")
ax.legend(frameon=False, loc="lower left", handlelength=1.4, labelspacing=.35)
ax.text(14.6, 99.5, "survey-weighted", ha="right", va="top", color=META_GREY, fontsize=6)
fig2.savefig(f"{FIG}/Figure2_survival_curves.png", dpi=300, bbox_inches="tight")

# ------------------------------------------------------------------ Fig. S1
# pooled participants, self-reported cancer, linkage-eligible with weight > 0, complete items, primary cohort
FLOW = [70190, 3782, 3779, 3268, 2569]
base = lm[lm.exclude_yr == 0].iloc[0]
assert int(base.n) == FLOW[-1]
figA, ax = plt.subplots(figsize=(6.6, 5.6)); ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")
main = [(8.9, "NHANES 2005–2018 participants", f"n = {FLOW[0]:,}"),
        (7.0, "Self-reported prior cancer diagnosis\n(MCQ220 = 1)", f"n = {FLOW[1]:,}"),
        (5.1, "Eligible for mortality linkage, with\nfollow-up time and MEC weight > 0", f"n = {FLOW[2]:,}"),
        (3.2, "Complete PHQ-9 (9 items), sleep duration\nand previously reported sleep trouble", f"n = {FLOW[3]:,}")]
exc = [(7.95, f"Excluded: no cancer history\nor missing (n = {FLOW[0] - FLOW[1]:,})"),
       (6.05, f"Excluded: ineligible for linkage\nor weight = 0 (n = {FLOW[1] - FLOW[2]:,})"),
       (4.15, f"Excluded: incomplete PHQ-9\nor sleep items (n = {FLOW[2] - FLOW[3]:,})"),
       (2.25, f"Excluded: reported only\nnon-melanoma / unspecified\nskin cancer (n = {FLOW[3] - FLOW[4]:,})")]
for y, t, n in main:
    ax.text(3.0, y, f"{t}\n{n}", ha="center", va="center", fontsize=7.2,
            bbox=dict(boxstyle="round,pad=0.42", fc="white", ec="#4a4a4a", lw=1.0))
for y, t in exc:
    ax.text(7.3, y, t, ha="left", va="center", fontsize=6.6, color="#3a3a3a",
            bbox=dict(boxstyle="round,pad=0.34", fc="#f6f6f6", ec="#b8b8b8", lw=0.7))
for y0, y1 in [(8.9, 7.0), (7.0, 5.1), (5.1, 3.2), (3.2, 1.15)]:
    ax.annotate("", xy=(3.0, y1 + 0.62), xytext=(3.0, y0 - 0.62),
                arrowprops=dict(arrowstyle="-|>", color="#4a4a4a", lw=1.0))
for y in [7.95, 6.05, 4.15, 2.25]:
    ax.annotate("", xy=(7.15, y), xytext=(3.0, y),
                arrowprops=dict(arrowstyle="-|>", color="#9a9a9a", lw=0.8))
ax.text(3.0, 1.15, f"PRIMARY ANALYTIC COHORT\nn = {FLOW[4]:,}  ·  {int(base.events)} deaths\n"
        f"{int(base.person_years):,} person-years",
        ha="center", va="center", fontsize=7.4, fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.46", fc="#eaf2f8", ec="#2c7fb8", lw=1.5))
ax.set_title("Participant flow", loc="left", fontsize=8)
figA.savefig(f"{SUP}/FigureS1_participant_flow.png", dpi=300, bbox_inches="tight")

# ------------------------------------------------------------------ Fig. S2
SGO = ["Age < 65 yr", "Age >= 65 yr", "Female", "Male", "<5 yr since dx", ">=5 yr since dx",
       "Breast cancer", "Prostate cancer", "Colorectal cancer", "Gynecologic cancer", "Melanoma"]
sgb = sg[sg.subgroup.isin(SGO)].copy()
bh = set(zip(sgb[sgb.p_bh < 0.05].subgroup, sgb[sgb.p_bh < 0.05].ph))
off = {PH3[0]: 0.26, PH3[1]: 0.0, PH3[2]: -0.26}
XS2 = (0.35, 11)
figS, ax = plt.subplots(figsize=(5.6, 6.0))
for i, s_ in enumerate(SGO[::-1]):
    g = sgb[sgb.subgroup == s_]
    for L in PH3:
        r = g[g.ph == L]
        if not len(r): continue
        r = r.iloc[0]
        if not np.isfinite(r.hi): continue
        y0 = i + off[L]; sigb = (s_, L) in bh
        ax.plot([max(r.lo, XS2[0]), min(r.hi, XS2[1])], [y0, y0], color=CL[L], lw=1.2, solid_capstyle="butt")
        if r.hi > XS2[1]: ax.plot(XS2[1], y0, ">", color=CL[L], ms=3.4, mew=0, clip_on=False)   # interval runs off the axis
        if r.lo < XS2[0]: ax.plot(XS2[0], y0, "<", color=CL[L], ms=3.4, mew=0, clip_on=False)
        ax.plot(r.HR, y0, "o", color=CL[L], ms=4.4, mec="k" if sigb else CL[L], mew=.7 if sigb else 0)
ax.axvline(1, color=META_GREY, lw=.9, ls="--")
ax.set_xscale("log"); ax.set_xticks([0.5, 1, 2, 4, 8])
ax.xaxis.set_major_formatter(mtick.FixedFormatter(["0.5", "1", "2", "4", "8"]))
ax.xaxis.set_minor_locator(mtick.NullLocator()); ax.set_xlim(*XS2)
ax.set_yticks(list(range(len(SGO))))
ax.set_yticklabels([f"{s}  (n={int(sgb[sgb.subgroup == s].n.iloc[0]):,}, {int(sgb[sgb.subgroup == s].events.iloc[0])} deaths)"
                    for s in SGO[::-1]])
ax.set_ylim(-0.75, len(SGO) - 0.25)
ax.set_xlabel("Hazard ratio for all-cause death (95% CI)")
ax.set_title("Phenotype and mortality within subgroups", loc="left")
for L in PH3: ax.plot([], [], "o-", color=CL[L], ms=4.4, lw=1.2, label=L)
ax.legend(frameon=False, ncol=3, handlelength=1.3, columnspacing=1.1,
          loc="upper center", bbox_to_anchor=(0.5, -0.085),
          title=f"vs low symptom burden (black outline: survives BH correction, {len(sg)} tests; arrow: interval runs off the axis)", title_fontsize=6)
figS.savefig(f"{SUP}/FigureS2_subgroups.png", dpi=300, bbox_inches="tight")

# ------------------------------------------------------------------ Fig. S3
B = pd.read_csv("data/derived/nhanes_cancer_design_frame.csv.gz")
sub = B[B.inAnalysis == 1]
sub = sub[~sub["only_nms"].astype(str).isin(["True", "TRUE", "1", "1.0"])]
assert len(sub) == FLOW[-1]
lims = {v: (sub[v].quantile(.01), sub[v].quantile(.99)) for v in ["phq9_score", "sleep_h"]}
g = rc[rc["var"] == "sleep_h"]
nad = g.loc[g.HR.idxmin()]
figB, axes = plt.subplots(1, 2, figsize=(6.6, 2.9), sharey=True, gridspec_kw={"wspace": 0.07})
for ax, (v, xl, ttl) in zip(axes, [("phq9_score", "PHQ-9 depression score", "Depressive symptoms: graded increase"),
                                   ("sleep_h", "Self-reported sleep duration (h)", "Sleep duration: U-shaped")]):
    lo_, hi_ = lims[v]
    d = rc[(rc["var"] == v) & (rc.x >= lo_) & (rc.x <= hi_)].sort_values("x")
    ax.fill_between(d.x, d.lo, d.hi, color="#2c7fb8", alpha=.15, lw=0)
    ax.plot(d.x, d.HR, color="#2c7fb8", lw=1.9)
    ax.axhline(1, color=META_GREY, lw=.8, ls="--")
    ax.set_xlabel(xl); ax.set_title(ttl, loc="left"); ax.margins(x=0.02)
axes[0].set_ylim(0.45, 3.5); axes[0].set_ylabel("Hazard ratio for all-cause death")
axes[0].text(0.03, 0.955, "reference = cohort median", transform=axes[0].transAxes,
             fontsize=6, color=META_GREY, va="top")
axes[1].annotate(f"nadir {nad.x:.2f} h", xy=(round(nad.x, 2), round(nad.HR, 2)), xytext=(8.1, 0.62), fontsize=6, color=META_GREY,
                 arrowprops=dict(arrowstyle="-", color=META_GREY, lw=.7))
figB.savefig(f"{SUP}/FigureS3_dose_response.png", dpi=300, bbox_inches="tight")
shutil.copyfile(S + "P_rcs_curves.csv", f"{SUP}/FigureS3_source_data.csv")

print(f"wrote Figure1, Figure2 to {FIG}/ and FigureS1-S3 (+ FigureS3_source_data.csv) to {SUP}/")
