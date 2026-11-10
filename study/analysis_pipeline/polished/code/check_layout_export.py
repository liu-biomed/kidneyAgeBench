"""Verify content preservation and physical sizes against the preceding figure set."""
from collections import Counter
from pathlib import Path
from xml.etree import ElementTree as ET
import csv
import hashlib
import math
import re
import subprocess

from PIL import Image
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT.parent / "genome_medicine_spaced"

def signatures(folder):
    return {p.relative_to(folder): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in folder.rglob("*") if p.is_file() and not p.name.startswith(".")}

def styled_words(svg):
    return Counter((word, node.get("style", "")) for node in svg.findall(".//{*}text")
                   for word in "".join(node.itertext()).split())

def styles_of_art(svg):
    return Counter((node.tag.split("}")[-1], node.get("style", ""))
                   for node in svg.iter()
                   if node.tag.split("}")[-1] in {"path","circle","polygon","line","polyline","rect"})

def run():
    for name in ("source_data","text","supplement"):
        assert signatures(BASE/name) == signatures(ROOT/name), name
    result = []
    for i in range(1,6):
        prefix = Path(f"figures/Figure{i}")
        new = ET.parse(ROOT/prefix.with_suffix(".svg"))
        old = ET.parse(BASE/prefix.with_suffix(".svg"))
        assert styled_words(new) == styled_words(old), f"word, font or color difference in {i}"
        assert styles_of_art(new) == styles_of_art(old), f"mark or style difference in {i}"
        assert not new.findall(".//{*}image"), f"embedded raster in {i}"
        pdf = ROOT/prefix.with_suffix(".pdf")
        page = PdfReader(pdf).pages[0]
        sizes = [float(re.search(r"font-size: ([.0-9]+)px", node.get("style"))[1])
                 for node in new.findall(".//{*}text")]
        assert min(sizes)>=7 and max(sizes)==11
        xml = ET.fromstring(subprocess.check_output(["pdftotext","-bbox",pdf,"-"]))
        pg = xml.find(".//{*}page")
        width, height = (float(pg.get(k)) for k in ("width","height"))
        for word in pg.findall("{*}word"):
            assert 0<=float(word.get("xMin"))<=float(word.get("xMax"))<=width
            assert 0<=float(word.get("yMin"))<=float(word.get("yMax"))<=height
        font_info = subprocess.check_output(["pdffonts",pdf],text=True).splitlines()[2:]
        assert len(font_info)>0 and all(ln.split()[-5] == "yes" for ln in font_info)
        pixel_size = Image.open(ROOT/prefix.with_suffix(".png")).size
        expected = tuple(math.ceil(float(x)*300/72) for x in (page.mediabox.width,page.mediabox.height))
        assert pixel_size == expected, (i,pixel_size,expected)
        result.append((i,len(sizes),min(sizes),max(sizes),*pixel_size,
                       "identical","identical","inside page","embedded"))
    with (ROOT/"audit/layout_export_checks.tsv").open("w",newline="") as out:
        writer = csv.writer(out,delimiter=chr(9),lineterminator=chr(10))
        writer.writerow(("figure","svg_text_nodes","smallest_pt","largest_pt","png_width_px",
                         "png_height_px","text_and_styles_vs_prior","art_counts_and_styles_vs_prior",
                         "pdf_word_bounds","pdf_fonts"))
        writer.writerows(result)
    return result

if __name__=="__main__":
    for row in run():
        print(*row)
