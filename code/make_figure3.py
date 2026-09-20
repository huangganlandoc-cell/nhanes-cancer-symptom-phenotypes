"""Figure 3: robustness of the somatic-depressive association (two stacked panels, 170 mm wide).

Left : stepped landmark exclusion, all three phenotypes (supporting/RR_landmark_full.csv).
Right: somatic-depressive vs low symptom burden under alternative class-assignment and
       measurement specifications, read from tables/Table4_robustness.csv so that the
       figure and Table 4 cannot disagree.
Panel titles describe content only. Colours match Figures 1 and S2
(validated with the dataviz palette checker: CVD and normal-vision separation pass).

Run from the project root:  python3 code/make_figure3.py
"""
import csv, math, re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, NullFormatter, FixedFormatter

COL = {"Sleep-fatigue": "#2c7fb8", "Somatic-depressive": "#c0392b", "High symptom burden": "#e08214"}
INK, MUTED, GRID = "#222222", "#6b6b6b", "#e6e6e6"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": "#9a9a9a", "axes.linewidth": 0.8,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.labelcolor": INK})

# ------------------------------------------------------------------ data
lm = list(csv.DictReader(open("supporting/RR_landmark_full.csv", encoding="utf-8")))
for r in lm:
    r["term"] = r["term"].replace("Insomnia-fatigue", "Sleep-fatigue")
years = sorted({int(r["exclude_yr"]) for r in lm})
deaths = {int(r["exclude_yr"]): int(r["events"]) for r in lm}

t4 = list(csv.DictReader(open("tables/Table4_robustness.csv", encoding="utf-8")))
H = list(t4[0].keys())


def get(spec_start, cov, comp="vs low symptom burden"):
    for r in t4:
        if r[H[1]].startswith(spec_start) and r[H[2]] == cov and r[H[3]] in (comp, "vs largest class"):
            m = re.match(r"([\d.]+) \(([\d.]+)–([\d.]+)\)", r[H[4]])
            return (float(m.group(1)), float(m.group(2)), float(m.group(3))) if m else None
    raise KeyError((spec_start, cov))


RIGHT = [("Model 3", None),
         ("Modal assignment", get("Modal assignment", "Model 3")),
         ("Bias-adjusted three-step", get("Bias-adjusted three-step", "Model 3")),
         ("Inclusive draws", get("Inclusive draws", "Model 3")),
         ("Model 4", None),
         ("Modal assignment", get("Modal assignment", "Model 4")),
         ("Bias-adjusted three-step", get("Bias-adjusted three-step", "Model 4")),
         ("Inclusive draws", get("Inclusive draws", "Model 4")),
         ("Indicator-only draws", get("Indicator-only draws", "Model 4")),
         ("Measurement choices", None),
         ("PHQ-9 sleep item removed", get("PHQ-9 sleep item removed", "Model 4")),
         ("Lower item threshold", "none"),
         ("Classes fitted in all adults", get("Measurement model refitted", "Model 4"))]

# ------------------------------------------------------------------ figure
# 170 mm (6.7 in) wide, the BMC full-page width; axes placed by hand so that the direct labels of
# panel a and the estimate annotations of panel b have room without shrinking the plotting areas
fig = plt.figure(figsize=(6.7, 8.4))
a1 = fig.add_axes([0.13, 0.63, 0.59, 0.30])
a2 = fig.add_axes([0.27, 0.06, 0.55, 0.43])

# left: landmark series
for ph in ["Sleep-fatigue", "High symptom burden", "Somatic-depressive"]:
    rs = sorted((r for r in lm if r["term"] == ph), key=lambda r: int(r["exclude_yr"]))
    x = [int(r["exclude_yr"]) for r in rs]
    y = [float(r["HR"]) for r in rs]
    lo = [float(r["lo"]) for r in rs]; hi = [float(r["hi"]) for r in rs]
    a1.fill_between(x, lo, hi, color=COL[ph], alpha=0.10, lw=0)
    a1.plot(x, y, color=COL[ph], lw=2, solid_capstyle="round", zorder=3)
    a1.plot(x, y, "o", ms=6, color=COL[ph], mec="white", mew=1.5, zorder=4, label=ph)
    a1.annotate(ph, (x[-1], y[-1]), xytext=(8, 0), textcoords="offset points",
                va="center", fontsize=8.5, color=INK)
a1.axhline(1, color="#8c8c8c", lw=0.9, ls=(0, (4, 3)), zorder=1)
a1.set_yscale("log")
a1.set_yticks([0.6, 1, 2, 4]); a1.set_yticklabels(["0.6", "1", "2", "4"])
a1.set_ylim(0.55, 5.8)
a1.set_xticks(years); a1.set_xlim(-0.4, 7.9)
a1.set_xlabel("Follow-up excluded after examination (years)")
a1.set_ylabel("Hazard ratio for all-cause death (log scale)")
a1.set_title("a  Stepped landmark exclusion", loc="left", fontsize=10, color=INK, pad=22)
for yr in years:
    a1.annotate(f"{deaths[yr]}", (yr, 5.8), xytext=(0, 3), textcoords="offset points",
                ha="center", va="bottom", fontsize=7.5, color=MUTED, annotation_clip=False)
a1.annotate("deaths remaining", (years[-1] + 0.62, 5.8), xytext=(0, 3), textcoords="offset points",
            ha="left", va="bottom", fontsize=7.5, color=MUTED, annotation_clip=False)
a1.grid(axis="y", color=GRID, lw=0.8); a1.set_axisbelow(True)
for s in ("top", "right"):
    a1.spines[s].set_visible(False)
# no legend: the three series are labelled directly at their last point

# right: specification forest plot (single series: somatic-depressive vs low symptom burden)
n = len(RIGHT)
ys = list(range(n))[::-1]
labels = []
for yv, (lab, val) in zip(ys, RIGHT):
    if val is None:                                   # group header
        labels.append((yv, lab, True)); continue
    labels.append((yv, lab, False))
    if val == "none":
        a2.annotate("no class associated", (1.03, yv), va="center", fontsize=8, color=MUTED, style="italic")
        continue
    hr, lo, hi = val
    a2.plot([lo, hi], [yv, yv], color=COL["Somatic-depressive"], lw=2, solid_capstyle="round")
    a2.plot(hr, yv, "o", ms=6.5, color=COL["Somatic-depressive"], mec="white", mew=1.5, zorder=3)
    a2.annotate(f"{hr:.2f} ({lo:.2f}–{hi:.2f})", (1.02, yv), xycoords=("axes fraction", "data"),
                va="center", ha="left", fontsize=7.8, color=INK, annotation_clip=False)
a2.axvline(1, color="#8c8c8c", lw=0.9, ls=(0, (4, 3)))
a2.set_xscale("log")
vals = [v for _, v in RIGHT if isinstance(v, tuple)]           # axis range from the data, so no interval is clipped
XL = (min(0.55, 0.9 * min(v[1] for v in vals)), max(4.6, 1.1 * max(v[2] for v in vals)))
TICKS = [t for t in (0.25, 0.5, 1, 2, 4, 8) if XL[0] <= t <= XL[1]]
a2.xaxis.set_major_locator(FixedLocator(TICKS)); a2.xaxis.set_major_formatter(FixedFormatter([f"{t:g}" for t in TICKS]))
a2.xaxis.set_minor_formatter(NullFormatter())
a2.set_xlim(*XL)
a2.set_yticks([y for y, _, _ in labels])
a2.set_yticklabels([("" if head else "   ") + lab for _, lab, head in labels])
for tick, (_, _, head) in zip(a2.get_yticklabels(), labels):
    tick.set_color(INK if head else MUTED)
    tick.set_fontweight("bold" if head else "normal")
a2.tick_params(axis="y", length=0)
a2.set_ylim(-0.7, n - 0.3)
a2.set_xlabel("Hazard ratio vs low symptom burden or largest class (95% CI, log scale)")
a2.set_title("b  Somatic-depressive phenotype, alternative specifications", loc="left", fontsize=10, color=INK, pad=22)
a2.grid(axis="x", color=GRID, lw=0.8); a2.set_axisbelow(True)
for s in ("top", "right", "left"):
    a2.spines[s].set_visible(False)

fig.savefig("figures/Figure3_sensitivity.png", dpi=300, facecolor="white")
print("wrote figures/Figure3_sensitivity.png")
