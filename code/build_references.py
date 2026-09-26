"""Re-emit the manuscript's reference list in a chosen journal style.

The bibliographic fields come from `supporting/reference_metadata.json`, which was built from CrossRef
and PubMed on 2026-09-19 (journal names are the NLM/Index Medicus abbreviations). The order and the
numbers come from the manuscript itself, so this script never renumbers anything: it matches each
listed reference by its DOI and rewrites only the formatting. References without a DOI (the two data
sources) are passed through unchanged apart from style-specific punctuation.

Styles
  bmc       Vancouver as BMC uses it: up to 6 authors then "et al", Journal. Year;vol(issue):pages. DOI
  nature    Nature/Scientific Reports: up to 5 authors then "et al.", Journal vol, pages (year).
  springer  Springer basic (e.g. Supportive Care in Cancer): up to 3 authors then "et al", (year) after
            the authors, Journal vol:pages, DOI as a URL

Run from the project root:
    python3 code/build_references.py bmc            # print the list
    python3 code/build_references.py nature --apply manuscript/Manuscript.md
"""
import argparse, difflib, json, re, sys
from pathlib import Path

META = "supporting/reference_metadata.json"


def load_manuscript(path):
    md = Path(path).read_text(encoding="utf-8")
    head, sep, refs = md.partition("## References")
    if not sep:
        sys.exit(f"no References section in {path}")
    entries = []
    for ln in refs.splitlines():
        m = re.match(r"^(\d+)\.\s+(.*)$", ln.strip())
        if m:
            entries.append((int(m.group(1)), m.group(2)))
    return head + sep, entries


def fmt_authors(authors, style):
    cap = {"bmc": 6, "nature": 5, "springer": 3}[style]
    if not authors:
        return ""
    shown = authors[:cap]
    if style == "nature":
        name = lambda a: f"{a['family']}, {' '.join(c + '.' for c in a['initials'])}".strip()
        if len(authors) >= 6:                  # Nature style: six or more authors -> first author et al.
            return name(authors[0]) + " et al."
        names = [name(a) for a in shown]
        if len(names) == 1:
            return names[0]
        return ", ".join(names[:-1]) + " & " + names[-1]
    names = [f"{a['family']} {a['initials']}".strip() for a in shown]
    tail = ", et al" if len(authors) > cap else ""
    return ", ".join(names) + tail


def pages(pg):
    """Full page ranges with an en dash, one form throughout: 1797-810 -> 1797\u20131810, e458-e466 -> e458\u2013e466."""
    m = re.fullmatch(r"([A-Za-z]*)(\d+)[-\u2013]([A-Za-z]*)(\d+)(.*)", pg or "")
    if not m:
        return pg
    p1, a, p2, b, rest = m.groups()
    if not p2 and len(b) < len(a):
        b = a[:len(a) - len(b)] + b
    return f"{p1}{a}\u2013{p2}{b}{rest}"


def fmt(rec, style):
    au = fmt_authors(rec["authors"], style)
    vol, iss, pg, yr = rec["volume"], rec["issue"], rec["pages"], rec["year"]
    pg = pages(pg)
    if style == "bmc":
        loc = yr + (f";{vol}" + (f"({iss})" if iss else "") + (f":{pg}" if pg else "") if vol else (f":{pg}" if pg else ""))
        doi = f" https://doi.org/{rec['doi']}" if rec["doi"] else ""
        return f"{au}. {rec['title']}. {rec['journal']}. {loc}.{doi}"
    if style == "nature":
        bits = [b for b in (f"**{vol}**" if vol else "", pg) if b]
        loc = ", ".join(bits)
        doi = f" https://doi.org/{rec['doi']}" if rec["doi"] else ""
        return f"{au} {rec['title']}. *{rec['journal']}* {loc} ({yr}).{doi}".replace("  ", " ")
    loc = (f"{vol}:{pg}" if vol and pg else vol or pg)
    doi = f". https://doi.org/{rec['doi']}" if rec["doi"] else ""
    return f"{au} ({yr}) {rec['title']}. {rec['journal']}{' ' + loc if loc else ''}{doi}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("style", choices=["bmc", "nature", "springer"])
    ap.add_argument("--apply", metavar="MANUSCRIPT", help="rewrite the reference list of this file in place")
    ap.add_argument("--source", default="manuscript/Manuscript.md", help="manuscript to read the list from")
    a = ap.parse_args()
    records = json.loads(Path(META).read_text(encoding="utf-8"))
    meta = {r["doi"]: r for r in records if r["doi"]}
    norm = lambda t: re.sub(r"[^a-z0-9 ]", " ", t.lower()).split()
    by_title = {" ".join(norm(r["title"])): r for r in records}
    head, entries = load_manuscript(a.apply or a.source)
    out, untouched = [], 0
    for n, text in entries:
        m = re.search(r"https?://doi\.org/(\S+?)\.?$", text)
        rec = meta.get(m.group(1)) if m else None
        if rec is None:                          # records without a DOI are matched on their title
            best = difflib.get_close_matches(" ".join(norm(text)), by_title, n=1, cutoff=0.45)
            rec = by_title[best[0]] if best else None
        if rec:
            out.append(f"{n}. {fmt(rec, a.style)}")
        else:                                    # data sources and anything not in the metadata file
            untouched += 1
            out.append(f"{n}. {text}")
    body = "\n" + "\n".join(out) + "\n"
    if a.apply:
        Path(a.apply).write_text(head + body, encoding="utf-8")
        print(f"rewrote {len(out)} references in {a.apply} ({a.style} style; {untouched} passed through)")
    else:
        print(body.strip())


if __name__ == "__main__":
    main()
