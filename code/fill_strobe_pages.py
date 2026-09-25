"""Fill the Page column of the STROBE checklist from a built manuscript file.

The checklist ships with that column blank because page numbers only exist once the manuscript has
been laid out. This paginates the built .docx with LibreOffice, finds the page each heading, figure
caption and table caption falls on, and resolves every "Where reported" entry against that index.

    python3 code/fill_strobe_pages.py submission/Scientific_Reports/02_Manuscript.docx out.md
"""
import re, subprocess, sys, tempfile
from pathlib import Path

SI = "Supplementary Information"


def page_index(docx_path):
    """{heading or caption label -> first page it appears on}, 1-based, from the PDF rendering."""
    from pypdf import PdfReader
    with tempfile.TemporaryDirectory() as tmp:
        subprocess.run(["soffice", "--headless", "--convert-to", "pdf", "--outdir", tmp, str(docx_path)],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        pdf = Path(tmp) / (Path(docx_path).stem + ".pdf")
        pages = [p.extract_text() for p in PdfReader(str(pdf)).pages]
    index = {}
    for n, text in enumerate(pages, start=1):
        for line in (l.strip() for l in text.split("\n")):
            cap = re.match(r"^((?:Figure|Table) \d+)\.", line)
            key = cap.group(1) if cap else line
            if key and len(key) < 80:
                index.setdefault(key.lower(), n)
    return index


def resolve(where, index):
    """'Methods: Statistical analysis; Discussion, final paragraph' -> '5, 9'."""
    if where.strip() in {"-", ""}:
        return "-"
    pages, supplementary = [], False
    for token in where.split(";"):
        t = re.sub(r"^(Methods|Results|Discussion|Declarations):\s*", "", token.strip())
        t = t.split(",")[0].strip()
        t = re.sub(r"\s*\(.*\)$", "", t)
        if re.search(r"(?:Table|Figure|Fig\.)s?\s+S\d", t):
            supplementary = True
            continue
        if t.lower() == "title":
            pages.append(1)
            continue
        hit = index.get(t.lower())
        if hit is None:                       # headings are cited by their opening words
            cands = [v for k, v in index.items() if k.startswith(t.lower())]
            hit = min(cands) if cands else None
        if hit is None and token.strip().split(":")[0] in {"Methods", "Results", "Discussion"}:
            hit = index.get(token.strip().split(":")[0].lower())
        if hit is not None:
            pages.append(hit)
    out = ", ".join(str(p) for p in sorted(set(pages)))
    return f"{out}; {SI}" if supplementary and out else (SI if supplementary else out)


def fill(md_text, docx_path):
    index = page_index(docx_path)
    out, unresolved = [], []
    for line in md_text.split("\n"):
        cells = [c.strip() for c in line.strip().strip("|").split("|")] if line.startswith("|") else None
        if cells and len(cells) == 5 and cells[0] not in {"Item", "---"}:
            cells[4] = resolve(cells[3], index)
            if not cells[4]:
                unresolved.append(cells[0])
            line = "| " + " | ".join(cells) + " |"
        out.append(line)
    assert not unresolved, f"no page found for STROBE items {unresolved}"
    text = "\n".join(out)
    return text.replace(" and should be filled in once it is paginated for submission", "")


if __name__ == "__main__":
    src = Path("supporting/STROBE_checklist_EN.md").read_text(encoding="utf-8")
    text = fill(src, sys.argv[1])
    Path(sys.argv[2]).write_text(text, encoding="utf-8")
    print(f"wrote {sys.argv[2]}")
    for line in text.split("\n"):
        if line.startswith("|"):
            c = [x.strip() for x in line.strip().strip("|").split("|")]
            if len(c) == 5:
                print(f"  {c[0]:>4}  {c[3][:58]:<58}  {c[4]}")
