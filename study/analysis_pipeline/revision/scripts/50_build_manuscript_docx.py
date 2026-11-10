#!/usr/bin/env python3
"""Build the manuscript .docx from markdown and fix table layout.

Usage: python 50_build_manuscript_docx.py MANUSCRIPT.md REFERENCE.docx OUT.docx
"""
import os
import re
import shutil
import subprocess
import sys
import zipfile
import docx
from docx.shared import Inches, Pt
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

# Column widths (inches) keyed by column count; must sum to <= 6.3 in.
WIDTHS = {
    6: [0.98, 1.05, 0.78, 2.05, 0.88, 0.56],
    4: [1.35, 1.85, 1.60, 1.50],
}


def border(el, edge, sz=6, val="single"):
    tag = f"w:{edge}"
    e = OxmlElement(tag)
    e.set(qn("w:val"), val)
    if val != "nil":
        e.set(qn("w:sz"), str(sz))
        e.set(qn("w:space"), "0")
        e.set(qn("w:color"), "000000")
    el.append(e)


def style_table(t):
    n = len(t.columns)
    w = WIDTHS.get(n)
    if w is None:
        return
    t.autofit = False
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    tblPr = t._tbl.tblPr

    total = OxmlElement("w:tblW")
    total.set(qn("w:type"), "dxa")
    total.set(qn("w:w"), str(int(sum(w) * 1440)))
    tblPr.append(total)

    layout = OxmlElement("w:tblLayout")
    layout.set(qn("w:type"), "fixed")
    tblPr.append(layout)

    borders = OxmlElement("w:tblBorders")
    border(borders, "top", 8)
    border(borders, "bottom", 8)
    border(borders, "left", 0, "nil")
    border(borders, "right", 0, "nil")
    border(borders, "insideH", 0, "nil")
    border(borders, "insideV", 0, "nil")
    tblPr.append(borders)

    for row in t.rows:
        trPr = row._tr.get_or_add_trPr()
        trPr.append(OxmlElement("w:cantSplit"))
        for i, cell in enumerate(row.cells):
            if i < len(w):
                cell.width = Inches(w[i])
            for p in cell.paragraphs:
                p.paragraph_format.space_before = Pt(2)
                p.paragraph_format.space_after = Pt(2)
                for r in p.runs:
                    r.font.size = Pt(9)

    # rule under the header row, bold header text, repeating header on page break
    for cell in t.rows[0].cells:
        tcPr = cell._tc.get_or_add_tcPr()
        b = OxmlElement("w:tcBorders")
        border(b, "bottom", 6)
        tcPr.append(b)
        for p in cell.paragraphs:
            for r in p.runs:
                r.bold = True
    trPr = t.rows[0]._tr.get_or_add_trPr()
    th = OxmlElement("w:tblHeader")
    th.set(qn("w:val"), "true")
    trPr.append(th)


def prune_media(path):
    """Drop media left over from the reference doc that the document never uses."""
    tmp = path + ".tmp"
    with zipfile.ZipFile(path) as zin:
        doc_xml = zin.read("word/document.xml").decode("utf-8")
        rels = zin.read("word/_rels/document.xml.rels").decode("utf-8")
        used_ids = set(re.findall(r'r:(?:embed|id)="(rId\d+)"', doc_xml))
        entries = re.findall(r"<Relationship [^>]*/>", rels)
        keep, drop_targets = [], set()
        for e in entries:
            rid = re.search(r'Id="(rId\d+)"', e).group(1)
            target = re.search(r'Target="([^"]+)"', e).group(1)
            if rid in used_ids or not target.startswith("media/"):
                keep.append(e)
            else:
                drop_targets.add("word/" + target)
        new_rels = re.sub(r"(<Relationships[^>]*>).*(</Relationships>)",
                          lambda m: m.group(1) + "".join(keep) + m.group(2),
                          rels, flags=re.S)
        with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as zout:
            for item in zin.infolist():
                if item.filename in drop_targets:
                    continue
                data = zin.read(item.filename)
                if item.filename == "word/_rels/document.xml.rels":
                    data = new_rels.encode("utf-8")
                zout.writestr(item, data)
    shutil.move(tmp, path)


def main(md, ref, out):
    base = os.path.dirname(os.path.abspath(md))
    resource = os.pathsep.join([base, os.path.dirname(base)])
    subprocess.run(
        ["pandoc", md, "-o", out, f"--reference-doc={ref}",
         f"--resource-path={resource}"],
        check=True,
    )
    doc = docx.Document(out)
    for t in doc.tables:
        style_table(t)
    for p in doc.paragraphs:
        txt = p.text.strip()
        if txt.startswith("Table 1.") or txt.startswith("Table 2."):
            p.paragraph_format.page_break_before = True
            p.paragraph_format.keep_with_next = True
            p.paragraph_format.space_after = Pt(6)
    doc.save(out)
    prune_media(out)
    print(f"built {out}")


if __name__ == "__main__":
    main(*sys.argv[1:4])
