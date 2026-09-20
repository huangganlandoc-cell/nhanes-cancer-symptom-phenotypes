"""Normalise labels in the display tables (labels only; no number is changed).

  * phenotype names: the R outputs still carry the pre-renaming labels
    (Insomnia-fatigue -> Sleep-fatigue, Hypersomnia-somatic -> Somatic-depressive)
    and raw model term names (lcaXxx, "lca ..."); SLQ050 is "previously reported
    sleep trouble" throughout, as defined in Methods
  * Table 1: the unlabelled first-primary category is the explicit "Unknown"
    category described in Methods (27 participants); 0/1 levels read No/Yes

Idempotent: safe to re-run after any table is regenerated.
Run from the project root:  python3 code/fix_display_labels.py
"""
import csv, glob, re

REPL = [("lcaInsomnia-fatigue", "Sleep-fatigue"), ("lcaHypersomnia-somatic", "Somatic-depressive"),
        ("lcaHigh symptom burden", "High symptom burden"),
        ("Insomnia-fatigue", "Sleep-fatigue"), ("insomnia-fatigue", "sleep-fatigue"),
        ("Hypersomnia-somatic", "Somatic-depressive"), ("hypersomnia-somatic", "somatic-depressive"),
        ("lca Schoenfeld residuals", "Phenotype Schoenfeld residuals"),
        ("lca x follow-up period interaction", "Phenotype × follow-up period interaction"),
        ("Insomnia complaint", "Previously reported sleep trouble")]   # SLQ050, named as in Methods

changed = {}
for path in sorted(glob.glob("tables/*.csv") + glob.glob("supplementary/Table*.csv")):
    txt = open(path, encoding="utf-8").read()
    new = txt
    for a, b in REPL:
        new = new.replace(a, b)
    if new != txt:
        open(path, "w", encoding="utf-8").write(new)
        changed[path] = sum(txt.count(a) for a, _ in REPL)

# ---- Table 1 levels
p = "tables/Table1_weighted_baseline.csv"
rows = list(csv.reader(open(p, encoding="utf-8")))
out, block, unknown = [], None, None
for r in rows:
    if r[0]:
        block = r[0]
    if block == "First primary cancer site" and r[0] == "" and r[1] == "":
        unknown = ["", "Unknown"] + r[2:]          # hold back; append at the end of the block
        continue
    if block in ("Hypertension", "Diabetes", "Cardiovascular disease",
                 "More than one primary cancer") and r[1] in ("0", "1"):
        r[1] = {"0": "No", "1": "Yes"}[r[1]]
    out.append(r)
if unknown:
    start = next(i for i, r in enumerate(out) if r[0] == "First primary cancer site")
    end = start + 1
    while end < len(out) and out[end][0] == "":     # block continues on rows with no characteristic
        end += 1
    out.insert(end, unknown)
    changed[p] = changed.get(p, 0) + 1
with open(p, "w", newline="", encoding="utf-8") as fh:
    csv.writer(fh).writerows(out)

left = [f for f in glob.glob("tables/*.csv") + glob.glob("supplementary/Table*.csv")
        if re.search(r"Insomnia|Hypersomnia|insomnia-fatigue|hypersomnia", open(f, encoding="utf-8").read())]
print("files changed:", changed)
print("old phenotype names remaining in:", left or "none")
