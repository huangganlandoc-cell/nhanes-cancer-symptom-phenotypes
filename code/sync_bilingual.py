"""Keep manuscript/Manuscript_中英对照.md in step with the source manuscript.

The bilingual file pairs every paragraph of manuscript/Manuscript.md (English above) with its Chinese
translation (below). After the source changes, this script re-derives the English paragraphs from the
source, keeps the existing Chinese wherever the English is unchanged, and lists the paragraphs that
need a new translation. Nothing is translated here.

    python3 code/sync_bilingual.py plan  OUT.json       # paragraphs needing translation (new or changed)
    python3 code/sync_bilingual.py build IN.json        # IN.json: {key: chinese}; writes the .md and .docx
"""
import json, re, subprocess, sys, datetime
from pathlib import Path

SRC = Path("manuscript/Manuscript.md")
BIL = Path("manuscript/Manuscript_中英对照.md")
PREFIX = {"Abstract": "A", "Introduction": "I", "Methods": "M", "Results": "R", "Discussion": "D",
          "Conclusion": "C", "Supplementary information": "S", "List of abbreviations": "X", "Declarations": "X"}


def old_blocks():
    """(english, chinese) pairs and English->bilingual heading map from the current bilingual file."""
    t = BIL.read_text(encoding="utf-8")
    head = t[:t.index("## Title, authors and affiliations")]
    pairs, heads = {}, {}
    for m in re.finditer(r"^(#{2,3}) (.+?)（(.+?)）\s*$", t, re.M):
        heads[m.group(2).strip()] = f"{m.group(2).strip()}（{m.group(3).strip()}）"
    for m in re.finditer(r"^\*\*([A-Z]-\d+)\*\*[^\n]*\n\n(.+?)\n\n(.+?)(?=\n\n(?:\*\*[A-Z]-\d+\*\*|---|#)|\n*\Z)", t, re.M | re.S):
        pairs[m.group(2).strip()] = m.group(3).strip()
    return head, pairs, heads


def new_units():
    """Ordered units from the source: ('h2'|'h3', text) headings and ('p', id, english, header) paragraphs."""
    src = SRC.read_text(encoding="utf-8")
    body, _, _refs = src.partition("\n## References")
    units, sec, n = [], "Title", {}
    def nid(p):
        n[p] = n.get(p, 0) + 1
        return f"{p}-{n[p]}"
    lines = body.split("\n")
    in_disp = False
    for ln in lines:
        s = ln.strip()
        if not s or s == "---":
            continue
        if s.startswith("# ") and not s.startswith("## "):
            units.append(("p", nid("T"), s[2:].strip(), "")); continue
        if s.startswith("## "):
            sec = s[3:].strip(); in_disp = sec == "Display items"
            units.append(("h2", sec)); continue
        if s.startswith("### "):
            units.append(("h3", s[4:].strip())); continue
        if in_disp:
            m = re.match(r"^\|\s*((?:Figure|Table)\s+S?\d+)\s*\|\s*`([^`]+)`\s*\|(.+)\|\s*$", s)
            if m:
                units.append(("p", nid("L"), m.group(3).strip(), f"{m.group(1)} — `{m.group(2)}`"))
            continue
        units.append(("p", nid(PREFIX.get(sec, "T")), s, ""))
    return units


def plan(out):
    _, pairs, _ = old_blocks()
    todo = {}
    for u in new_units():
        if u[0] == "p" and u[2] not in pairs:
            # the closest old paragraph, to reuse its wording where the English changed only in part
            import difflib
            best = difflib.get_close_matches(u[2], list(pairs), n=1, cutoff=0.5)
            todo[u[1]] = {"english": u[2], "old_english": best[0] if best else "",
                          "old_chinese": pairs[best[0]] if best else ""}
    Path(out).write_text(json.dumps(todo, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"{len(todo)} paragraph(s) need translation -> {out}")
    for k, v in todo.items():
        print(f"  {k:6s} {'changed' if v['old_english'] else 'new':8s} {v['english'][:90]}")


def build(tr_path):
    head, pairs, heads = old_blocks()
    tr = json.loads(Path(tr_path).read_text(encoding="utf-8")) if tr_path else {}
    today = datetime.date.today().isoformat()
    head = re.sub(r"本稿于 \d{4}-\d{2}-\d{2} 由", f"本稿于 {today} 由", head, count=1)
    out, missing = [head.rstrip("\n"), ""], []
    for u in new_units():
        if u[0] in ("h2", "h3"):
            if u[0] == "h2" and out[-1] != "---" and not out[-2:] == ["", ""]:
                out += ["---", ""]
            out += [f"{'##' if u[0] == 'h2' else '###'} {heads.get(u[1], u[1])}", ""]
            continue
        _, pid, en, hdr = u
        zh = pairs.get(en) or tr.get(pid)
        if not zh:
            missing.append(pid); zh = "（待译）"
        out += [f"**{pid}**" + (f" {hdr}" if hdr else ""), "", en, "", zh, ""]
    out += ["---", "", f"## {heads.get('References', 'References（参考文献）')}", "",
            "参考文献未作任何改动，完整列表见源文件 `manuscript/Manuscript.md` 的 References 部分（编号 1 至 48），本对照稿不再重复。", ""]
    assert not missing, f"no Chinese for {missing}"
    BIL.write_text("\n".join(out), encoding="utf-8")
    subprocess.run(["pandoc", "-f", "markdown+east_asian_line_breaks", "-t", "docx", str(BIL), "-o",
                    str(BIL.with_suffix(".docx"))], check=True)
    sys.path.insert(0, "code")
    import build_docx
    from docx import Document
    d = Document(str(BIL.with_suffix(".docx"))); build_docx.set_properties(d, BIL.with_suffix(".docx"))
    d.save(str(BIL.with_suffix(".docx")))
    print(f"wrote {BIL} and {BIL.with_suffix('.docx').name}")


if __name__ == "__main__":
    {"plan": lambda: plan(sys.argv[2]), "build": lambda: build(sys.argv[2] if len(sys.argv) > 2 else None)}[sys.argv[1]]()
