"""Assemble Manuscript.md + figures + tables into a single .docx.

Figures and tables are placed immediately after the paragraph that first cites
them; the markdown "Display items" manifest is dropped (its legends are reused
as captions). Supplementary items follow the references.

Run from the project root:  python3 code/build_docx.py
"""
import re, csv, datetime
from pathlib import Path
from docx import Document
from docx.shared import Pt, Inches, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

FIG_FILES = {
    "Figure 1": "figures/Figure1_symptom_phenotypes.png",
    "Figure 2": "figures/Figure2_survival_curves.png",
    "Figure 3": "figures/Figure3_sensitivity.png",
    "Figure S1": "supplementary/FigureS1_participant_flow.png",
    "Figure S2": "supplementary/FigureS2_subgroups.png",
    "Figure S3": "supplementary/FigureS3_dose_response.png",
}
TBL_FILES = {
    "Table 1": "tables/Table1_weighted_baseline.csv",
    "Table 2": "tables/Table2_main_cox.csv",
    "Table 3": "tables/Table3_single_symptom_exposures.csv",
    "Table 4": "tables/Table4_robustness.csv",
}
FIG_W = {"Figure 1": 6.4, "Figure 2": 3.7, "Figure 3": 6.5,
         "Figure S1": 5.9, "Figure S2": 5.4, "Figure S3": 6.2}


# ---------------------------------------------------------------- legends
def parse_legends(md):
    """Pull `| Figure 1 | `path` | Legend |` rows out of the Display items table."""
    out = {}
    for m in re.finditer(r'^\|\s*((?:Figure|Table)\s+S?\d+)\s*\|([^|]*)\|(.+?)\|\s*$', md, re.M):
        out[m.group(1).strip()] = m.group(3).strip()
    return out


# ---------------------------------------------------------------- inline runs
def add_runs(par, text, base=10):
    """Render **bold**, *italic*, `code`, and superscript digits after ×10."""
    text = text.replace('\u2212', '-')
    for tok in re.split(r'(\*\*[^*]+\*\*|\*[^*]+\*|`[^`]+`)', text):
        if not tok:
            continue
        if tok.startswith('**') and tok.endswith('**'):
            r = par.add_run(tok[2:-2]); r.bold = True
        elif tok.startswith('*') and tok.endswith('*'):
            r = par.add_run(tok[1:-1]); r.italic = True
        elif tok.startswith('`') and tok.endswith('`'):
            r = par.add_run(tok[1:-1]); r.font.name = 'Consolas'; r.font.size = Pt(base - 1)
        else:
            # split out unicode superscripts so they render as real superscript
            for piece in re.split(r'([\u2070\u00b9\u00b2\u00b3\u2074-\u2079\u207b]+)', tok):
                if not piece:
                    continue
                if re.fullmatch(r'[\u2070\u00b9\u00b2\u00b3\u2074-\u2079\u207b]+', piece):
                    tr = str.maketrans('\u2070\u00b9\u00b2\u00b3\u2074\u2075\u2076\u2077\u2078\u2079\u207b',
                                       '0123456789-')
                    r = par.add_run(piece.translate(tr)); r.font.superscript = True
                else:
                    par.add_run(piece)
    for r in par.runs:
        if r.font.size is None:
            r.font.size = Pt(base)


def para(doc, text, base=10, space_after=6, align=None, style=None):
    p = doc.add_paragraph(style=style)
    add_runs(p, text, base)
    p.paragraph_format.space_after = Pt(space_after)
    p.paragraph_format.line_spacing = 1.5 if base >= 10 else 1.15
    if align:
        p.alignment = align
    return p


def caption(doc, label, legend):
    p = doc.add_paragraph()
    r = p.add_run(label + '. '); r.bold = True; r.font.size = Pt(9)
    add_runs(p, legend, base=9)
    for r in p.runs:
        r.font.size = Pt(9)
    p.paragraph_format.space_before = Pt(4)
    p.paragraph_format.space_after = Pt(14)
    p.paragraph_format.line_spacing = 1.15
    return p


def add_numbering(doc):
    """Arabic page number in the footer and continuous line numbers, both asked for by Nature Portfolio."""
    for sec in doc.sections:
        sec.footer.is_linked_to_previous = False
        p = sec.footer.paragraphs[0]
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        r = p.add_run(); r.font.size = Pt(9)
        for tag, attr, txt in (('w:fldChar', ('w:fldCharType', 'begin'), None),
                               ('w:instrText', ('xml:space', 'preserve'), ' PAGE '),
                               ('w:fldChar', ('w:fldCharType', 'end'), None)):
            el = OxmlElement(tag); el.set(qn(attr[0]), attr[1])
            if txt:
                el.text = txt
            r._r.append(el)
        ln = OxmlElement('w:lnNumType')
        for k, v in (('w:countBy', '1'), ('w:start', '1'), ('w:restart', 'continuous')):
            ln.set(qn(k), v)
        sec._sectPr.insert_element_before(ln, 'w:pgNumType', 'w:cols', 'w:formProt', 'w:vAlign',
                                          'w:noEndnote', 'w:titlePg', 'w:textDirection', 'w:bidi',
                                          'w:rtlGutter', 'w:docGrid', 'w:printerSettings', 'w:sectPrChange')


AUTHOR = 'Yongyong Bao'   # corresponding author; written into every generated file's properties


def set_properties(doc, out):
    """Set the author, title and dates in the document properties; python-docx otherwise leaves
    author 'python-docx', a 'generated by python-docx' comment and a 2013 creation date."""
    cp = doc.core_properties
    cp.author = cp.last_modified_by = AUTHOR
    cp.comments = cp.keywords = cp.subject = cp.category = ''
    cp.title = Path(out).stem
    cp.created = cp.modified = datetime.datetime.now().replace(microsecond=0)
    cp.revision = 1


# ---------------------------------------------------------------- inserters
def insert_figure(doc, key, legends, label=None):
    doc.add_picture(FIG_FILES[key], width=Inches(FIG_W.get(key, 6.2)))
    doc.paragraphs[-1].alignment = WD_ALIGN_PARAGRAPH.CENTER
    doc.paragraphs[-1].paragraph_format.space_before = Pt(10)
    doc.paragraphs[-1].paragraph_format.space_after = Pt(2)
    caption(doc, label or key, legends.get(key, ''))


# ---------------------------------------------------------------- display formatting of table cells
# The CSVs are analysis outputs and stay as written; only what is printed changes: P values to three decimals
# (below 0.001 as a x 10^-n), confidence limits joined by an en dash, code-style headers and terms spelled out.
_SUP = str.maketrans("0123456789-", "\u2070\u00b9\u00b2\u00b3\u2074\u2075\u2076\u2077\u2078\u2079\u207b")
_ENOTE = re.compile(r"(?<![\w.])(\d+(?:\.\d+)?)[eE]([-+]?\d+)(?![\w.])")
HEADERS = {"k": "Classes", "logLik": "Log-likelihood", "entropy": "Entropy", "min_class_pct": "Smallest class (%)",
           "subgroup": "Subgroup", "events": "Deaths", "p": "P", "p_bh": "BH-adjusted P",
           "joint_p_phenotype": "Joint P, phenotype", "joint_p_11_indicators": "Joint P, 11 indicators",
           "Pearson X2": "Pearson \u03c7\u00b2", "BVR": "Bivariate residual"}
TABLE_HEADERS = {"Table S2": {"model": "Covariate set", "term": "Phenotype"},
                 "Table S8": {"": "Indicator"},
                 "Table S9": {"model": "Model", "term": "Term"}}
S9_MODELS = {"1 phenotype only": "A. Phenotype only", "2 PHQ-9 continuous only": "B. PHQ-9 score only (linear)",
             "3 phenotype + PHQ-9 continuous": "C. Phenotype + PHQ-9 score (linear)",
             "4 phenotype + PHQ-9 spline": "D. Phenotype + PHQ-9 score (spline)",
             "6 all 11 indicators + phenotype": "E. Phenotype + all 11 indicators"}
CELL_TEXT = [("I(phq9_score/5)", "PHQ-9 score (per 5 points)"), ("JKn", "jackknife"), ("chi2=", "\u03c7\u00b2 = "),
             ("kg/m2", "kg/m\u00b2"), ("NH White", "Non-Hispanic White"), ("NH Black", "Non-Hispanic Black"),
             ("report a sleep complaint", "with previously reported sleep trouble"), (" since dx", " since diagnosis"),
             ("subscale AFF", "subscale affective"), ("subscale SOM", "subscale somatic"),
             ("CVD/stroke", "Cardiovascular or cerebrovascular"), ("Cardiovascular or stroke", "Cardiovascular or cerebrovascular"),
             ("time-varying phenotype effect", "time-varying phenotype association"),
             ("exposure satisfies PH", "phenotype meets the proportional hazards assumption"), (">=", "\u2265"), ("<=", "\u2264")]


def _sci(x, digits=2):
    mant, exp = f"{x:.{digits - 1}e}".split("e")
    return f"{mant} \u00d7 10{str(int(exp)).translate(_SUP)}"


def _sig(s):
    """Significant digits written in the source string (0.0006 -> 1, 6.3e-04 -> 2), at most 2."""
    m = s.lower().split("e")[0].replace("-", "").replace(".", "").lstrip("0")
    return max(1, min(2, len(m)))


def _is_p(header):
    return header in ("P", "p", "p_bh") or header.endswith(" P") or header.startswith(("Joint P", "joint_p"))


def _p1(s):
    x = float(s)
    return _sci(x, _sig(s)) if 0 < x < 0.001 else (f"{x:.3f}" if 0 <= x <= 1 else s)


_PNUM = r"(\d\.\d{4,}|\d(?:\.\d+)?[eE]-\d+)"


def _p(cell):
    try:
        return _p1(cell)
    except ValueError:                           # e.g. "0.5393 (0.0348 when fitted alone)"
        return re.sub(r"(?<![\d.])" + _PNUM + r"(?![\d])", lambda m: _p1(m.group(1)), cell)


S5_LABELS = [("Model 4: + comorbidity and frailty block", "Model 4 without antidepressant use"),
             ("Model 4b: individual conditions + antidepressant", "Model 4 with the seven conditions entered individually"),
             ("Model 4b individual conditions", "Model 4 with the seven conditions entered individually"),
             ("Model 4 + antidepressant use", "Model 4"), ("Model 4 + antidepressant (new primary)", "Model 4"),
             ("Model 4 + antidepressant", "Model 4"), ("Model 4 + antidep", "Model 4"),
             ("Model 3 (previous primary)", "Model 3"),
             ("Model 4 + pseudo-class (most conservative)", "Model 4, symptom-based pseudo-class draws"),
             ("Assignment-uncertainty correction under Model 4 (100 draws)",
              "Symptom-based pseudo-class draws under Model 4 (100 draws)"),
             ("under assignment correction", "with symptom-based pseudo-class draws")]


def _text(cell, key=None):
    if key == "Table S5":                        # labels as the manuscript defines Model 4 (antidepressant use included)
        for a, b in S5_LABELS:
            cell = cell.replace(a, b)
        cell = re.sub(r"(-?\d+(?:\.\d+)?)% of excess risk removed",
                      lambda m: (f"excess hazard (HR \u2212 1) {abs(float(m.group(1))):g}% "
                                 + ("smaller" if float(m.group(1)) >= 0 else "larger")), cell)
    for a, b in CELL_TEXT:
        cell = cell.replace(a, b)
    cell = re.sub(r"\b((?:joint |design-based joint Wald |joint Wald |Wald )?P|BH|Bonferroni)( =)? " + _PNUM,
                  lambda m: f"{m.group(1)}{m.group(2) or ''} {_p1(m.group(3))}", cell)
    cell = re.sub(r"(?<=[A-Za-z]) x (?=[A-Za-z])", " \u00d7 ", cell)       # age group x years -> age group × years
    cell = re.sub(r"(?<=\d)x(?= )", "\u00d7", cell)                        # 66.6x current sample -> 66.6× ...
    cell = re.sub(r"\b([kn])=(\d)", r"\1 = \2", cell)
    cell = _ENOTE.sub(lambda m: _sci(float(m.group(0))) if float(m.group(0)) < 0.001 else m.group(0), cell)
    cell = re.sub(r"(?<=\d) py\b", " person-years", cell)
    cell = re.sub(r"(?<=\d)-(?=\d)", "\u2013", cell)            # 0.78-1.39 -> 0.78\u20131.39, 0-4 -> 0\u20134
    cell = re.sub(r"(?<![\w.])nan(?![\w.])", "\u2013", cell)
    return "\u2013" if cell.strip() == "NA" else cell


def display_rows(key, rows):
    head = [TABLE_HEADERS.get(key, {}).get(h, HEADERS.get(h, h)) for h in rows[0]]
    body = [list(r) for r in rows[1:]]
    if key == "Table S2" and {"HR", "lo", "hi"} <= set(rows[0]):   # one "HR (95% CI)" column, as in the other tables
        i, j, k = (rows[0].index(c) for c in ("HR", "lo", "hi"))
        for r in body:
            r[i] = f"{float(r[i]):.2f} ({float(r[j]):.2f}\u2013{float(r[k]):.2f})" if r[i] not in ("", "NA") else r[i]
        keep = [c for c in range(len(head)) if c not in (j, k)]
        head = ["HR (95% CI)" if c == i else head[c] for c in keep]
        body = [[r[c] for c in keep] for r in body]
    if key == "Table S9":
        mi = rows[0].index("model")
        for r in body:
            r[mi] = S9_MODELS.get(r[mi], r[mi])
    out = [head]
    for r in body:
        out.append([_p(c.strip()) if _is_p(head[n]) else _text(c, key) for n, c in enumerate(r)])
    return out


def insert_table(doc, key, legends, label=None):
    with open(TBL_FILES[key], newline='', encoding='utf-8') as fh:
        rows = list(csv.reader(fh))
    if not rows:
        return
    rows = display_rows(key, rows)
    caption(doc, label or key, legends.get(key, ''))
    ncol = len(rows[0])
    size = 7.5 if ncol <= 7 else 7.0                 # wide tables (e.g. Table 2, nine columns) a little smaller
    t = doc.add_table(rows=len(rows), cols=ncol)
    t.style = 'Table Grid'
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    mar = OxmlElement('w:tblCellMar')                # narrow cell margins (0.03 in) so numbers need not wrap
    for side in ('left', 'right'):
        el = OxmlElement(f'w:{side}'); el.set(qn('w:w'), '43'); el.set(qn('w:type'), 'dxa'); mar.append(el)
    t._tbl.tblPr.append(mar)
    for j, cell in enumerate(rows[0]):
        c = t.cell(0, j); c.text = ''
        p = c.paragraphs[0]; r = p.add_run(cell); r.bold = True; r.font.size = Pt(size)
        p.paragraph_format.space_after = Pt(0)
    for i, row in enumerate(rows[1:], start=1):
        for j, cell in enumerate(row):
            c = t.cell(i, j); c.text = ''
            p = c.paragraphs[0]
            add_runs(p, '' if cell in ('nan', 'None') else cell, base=size)   # x 10^-n as real superscript
            p.paragraph_format.space_after = Pt(0)
    # column widths: numeric columns (HR, P) and very short ones get the width of their longest cell so they
    # never wrap; text columns share what is left in proportion to their length, never below their longest word
    body = rows[1:] or rows
    longest = [max(len(r[j]) for r in body if j < len(r)) for j in range(ncol)]
    word = [max(len(w) for r in rows if j < len(r) for w in (r[j].split() or [''])) for j in range(ncol)]
    inch = lambda n: n * 0.55 * size / 72 + 0.10     # ~0.55 em per character, plus margins
    def numeric(j):
        cells = [r[j] for r in body if j < len(r) and r[j].strip()]
        return bool(cells) and sum(any(ch.isdigit() for ch in c) for c in cells) >= 0.8 * len(cells)
    total, short = 6.5, [j for j in range(ncol) if longest[j] <= 12 or numeric(j)]
    width = {j: inch(max(longest[j], word[j])) for j in short}
    rest = [j for j in range(ncol) if j not in short]
    free = total - sum(width.values())          # may be <= 0 in a very wide table; every column still gets a width
    for j in rest:
        share = free * longest[j] / sum(longest[k] for k in rest) if free > 0 else 0
        width[j] = max(inch(word[j]), share)
    scale = total / sum(width.values())
    t.autofit = False
    for j, col in enumerate(t.columns):              # grid widths (read by LibreOffice) and cell widths (Word)
        col.width = Inches(width[j] * scale)
    for row in t.rows:
        for j, c in enumerate(row.cells):
            c.width = Inches(width[j] * scale)
    doc.add_paragraph().paragraph_format.space_after = Pt(10)


# ---------------------------------------------------------------- main build
def render_tail(doc, text):
    """Render a block that belongs after the references (e.g. Springer's Statements and Declarations)."""
    for block in [b for b in text.split('\n') if b.strip() != '---']:
        ln = block.rstrip()
        if not ln:
            continue
        if ln.startswith('## '):
            h = doc.add_heading(ln[3:].strip(), level=1)
            for r in h.runs:
                r.font.size = Pt(13); r.font.color.rgb = RGBColor(0, 0, 0)
        elif ln.startswith('### '):
            h = doc.add_heading(ln[4:].strip(), level=2)
            for r in h.runs:
                r.font.size = Pt(11); r.font.color.rgb = RGBColor(0, 0, 0)
        else:
            para(doc, ln.strip())


def build(md_path='manuscript/Manuscript.md', out='Manuscript_with_figures.docx', tail_md=None,
          numbering=False):
    md = open(md_path, encoding='utf-8').read()
    legends = parse_legends(md)

    doc = Document()
    st = doc.styles['Normal']
    st.font.name = 'Calibri'; st.font.size = Pt(10)
    st.element.rPr.rFonts.set(qn('w:eastAsia'), 'Calibri')
    for s in doc.sections:
        s.top_margin = s.bottom_margin = Inches(0.9)
        s.left_margin = s.right_margin = Inches(0.9)

    body, _, refs = md.partition('## References')
    body = body.split('## Display items')[0]

    placed = set()
    lines = body.split('\n')
    i = 0
    while i < len(lines):
        ln = lines[i].rstrip()

        if not ln or ln == '---':
            i += 1
            continue

        if ln.startswith('# ') and not ln.startswith('## '):
            h = doc.add_heading(ln[2:].strip(), level=0)
            for r in h.runs:
                r.font.size = Pt(16); r.font.color.rgb = RGBColor(0, 0, 0)
            i += 1
            continue

        if ln.startswith('## '):
            h = doc.add_heading(ln[3:].strip(), level=1)
            for r in h.runs:
                r.font.size = Pt(13); r.font.color.rgb = RGBColor(0, 0, 0)
            i += 1
            continue

        if ln.startswith('### '):
            h = doc.add_heading(ln[4:].strip(), level=2)
            for r in h.runs:
                r.font.size = Pt(11); r.font.color.rgb = RGBColor(0, 0, 0)
            i += 1
            continue

        if ln.startswith('> '):
            para(doc, ln[2:], base=9, space_after=8)
            i += 1
            continue

        if ln.lstrip().startswith('|'):          # stray markdown table -> skip block
            while i < len(lines) and lines[i].lstrip().startswith('|'):
                i += 1
            continue

        # ordinary paragraph: gather until blank line
        buf = []
        while i < len(lines) and lines[i].strip() and not lines[i].startswith(('#', '|', '> ', '---')):
            buf.append(lines[i].strip())
            i += 1
        text = ' '.join(buf)
        if not text:
            continue
        para(doc, text)

        # place any display item this paragraph is the first to cite
        for key in list(FIG_FILES) + list(TBL_FILES):
            if key in placed or key.startswith(('Figure S', 'Table S')):
                continue
            if re.search(r'\b' + key.replace(' ', r'\s+') + r'\b', text):
                if key in FIG_FILES:
                    insert_figure(doc, key, legends)
                else:
                    insert_table(doc, key, legends)
                placed.add(key)

    # anything never cited in prose still belongs in the document
    for key in list(FIG_FILES) + list(TBL_FILES):
        if key not in placed and not key.startswith(('Figure S', 'Table S')):
            h = doc.add_heading(key, level=2)
            for r in h.runs:
                r.font.size = Pt(11); r.font.color.rgb = RGBColor(0, 0, 0)
            (insert_figure if key in FIG_FILES else insert_table)(doc, key, legends)
            placed.add(key)

    # ---- references
    h = doc.add_heading('References', level=1)
    for r in h.runs:
        r.font.size = Pt(13); r.font.color.rgb = RGBColor(0, 0, 0)
    for ln in refs.split('\n'):
        s_ = ln.strip()
        if re.match(r'^\d+\.\s', s_):
            p = para(doc, s_, base=9, space_after=3)
            p.paragraph_format.line_spacing = 1.15
        elif s_.startswith('> '):
            para(doc, s_[2:], base=9, space_after=8)

    if tail_md:
        doc.add_page_break()
        render_tail(doc, tail_md)
    if numbering:
        add_numbering(doc)
    set_properties(doc, out)
    doc.save(out)
    return out, sorted(placed), legends


def build_additional_file(md_path='manuscript/Manuscript.md', out='Additional_file_1.docx',
                          title='Additional file 1', label_prefix='', intro=None):
    """One supplementary file: figures S1 to S3 then tables S1 to S15, each with its manuscript legend.
    BMC calls it Additional file 1; Nature Portfolio wants Supplementary Information and requires the
    word "Supplementary" on every item label, which is what label_prefix supplies."""
    import glob
    md = open(md_path, encoding='utf-8').read()
    legends = parse_legends(md)
    doc = Document()
    st = doc.styles['Normal']
    st.font.name = 'Calibri'; st.font.size = Pt(10)
    st.element.rPr.rFonts.set(qn('w:eastAsia'), 'Calibri')
    for sec in doc.sections:
        sec.top_margin = sec.bottom_margin = Inches(0.9)
        sec.left_margin = sec.right_margin = Inches(0.9)
    h = doc.add_heading(title, level=1)
    for r in h.runs:
        r.font.size = Pt(14); r.font.color.rgb = RGBColor(0, 0, 0)
    para(doc, intro or 'Supplementary figures S1 to S3 and supplementary tables S1 to S15 for '
         '"Symptom phenotypes and mortality in US cancer survivors: a cohort study".',
         base=10, space_after=12)
    for key in ['Figure S1', 'Figure S2', 'Figure S3']:
        insert_figure(doc, key, legends, label=label_prefix + key)
    keys = sorted(glob.glob('supplementary/TableS*.csv'),
                  key=lambda f: int(re.search(r'TableS(\d+)', f).group(1)))
    for f in keys:
        key = 'Table S' + re.search(r'TableS(\d+)', f).group(1)
        TBL_FILES[key] = f
        doc.add_page_break()
        insert_table(doc, key, legends, label=label_prefix + key)
    set_properties(doc, out)
    doc.save(out)
    return out, len(keys)


if __name__ == '__main__':
    out, placed, legends = build()
    print('wrote', out)
    print('inline items:', placed)
    print('legends found:', len(legends))
    add, n = build_additional_file()
    print('wrote', add, 'with 3 figures and', n, 'tables')
