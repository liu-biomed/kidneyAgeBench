"""Assemble the five approved vector figures at native size, with separate Markdown legends."""
from pathlib import Path
import csv
import hashlib
import re
import subprocess
import tempfile
from pypdf import PdfReader
from PIL import Image

root=Path(__file__).resolve().parents[1]
files=[root/f"figures/Figure{i}.pdf" for i in range(1,6)]
dest=root/"preview/Genome_Medicine_figures_only.pdf"
legend=(root/"text/Figure_legends.md").read_text()
man=(root/"text/manuscript_redesign.md").read_text()
blocks=re.findall(r"^\*\*Figure (\d)\. (.+?)\*\* (.+)$",legend,re.M)
assert [n for n,_,_ in blocks]==list("12345")
for n,title,body in blocks:
    assert f"**Figure {n}. {title}** {body}" in man
    for link in re.findall(r"\]\(([^)]+)\)",body):
        assert (root/"text"/link).is_file()

with tempfile.TemporaryDirectory(dir=dest.parent) as temp:
    work=Path(temp)
    pdf=work/"proof.pdf"
    subprocess.run(["pdfunite",*map(str,files),str(pdf)],check=True)
    doc=PdfReader(pdf)
    assert len(doc.pages)==5
    data=[]
    thumb=Image.new("RGB",(1400,1170),"#E8EBEA")
    for i,(page,original) in enumerate(zip(doc.pages,files),1):
        source=PdfReader(original).pages[0]
        assert tuple(page.mediabox)==tuple(source.mediabox)
        w=float(page.mediabox.width)*25.4/72
        h=float(page.mediabox.height)*25.4/72
        assert abs(w-170)<.01 and h<=225
        assert page.extract_text()==source.extract_text()
        title=blocks[i-1][1]
        assert title not in page.extract_text()
        prefix=work/f"page-{i}"
        subprocess.run(["pdftoppm","-singlefile","-f",str(i),"-l",str(i),"-png","-r","300",str(pdf),str(prefix)],check=True)
        with Image.open(prefix.with_suffix(".png")) as out, Image.open(original.with_suffix(".png")) as src:
            pixels=out.convert("RGB")
            assert pixels.size==src.size and pixels.tobytes()==src.convert("RGB").tobytes()
            pixels.thumbnail((435,544),Image.Resampling.LANCZOS)
            x=(18,480,942)[i-1] if i<=3 else (249,711)[i-4]
            thumb.paste(pixels,(x,15 if i<=3 else 587))
        data.append((i,round(w,2),round(h,2),"yes","identical"))
    pdf.replace(dest)
    thumb.save(root/"preview/Five_figures_only_overview.png")
with (root/"audit/figures_only_fit.tsv").open("w",newline="") as f:
    wr=csv.writer(f,delimiter="\t",lineterminator="\n")
    wr.writerow(("figure","width_mm","figure_height_mm","within_225_mm","vs_approved_png_300dpi"))
    wr.writerows(data)
print(*data,sep="\n")
