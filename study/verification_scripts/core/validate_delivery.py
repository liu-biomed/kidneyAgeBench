from pathlib import Path
from lxml import etree as E
from pypdf import PdfReader
from PIL import Image
import math,hashlib,re,subprocess,tempfile,csv,zipfile
from collections import Counter
root=Path(__file__).resolve().parents[1];f=root/'Figures_final'
assert (root/'RAAG_provenance_audit.csv').read_bytes()==(root/'raag/RAAG_provenance_audit.csv').read_bytes()
doc=PdfReader(f/'Figures_1-5_preview.pdf');assert len(doc.pages)==5
rows=[]
for i,page in enumerate(doc.pages,1):
 path=f/f'Figure{i}'; orig=PdfReader(path.with_suffix('.pdf')).pages[0]
 w,h=map(float,(page.mediabox.width,page.mediabox.height))
 assert abs(w*25.4/72-170)<.02 and page.extract_text()==orig.extract_text()
 svg=E.parse(path.with_suffix('.svg'));assert len(svg.xpath('//*[local-name()="image"]'))==0
 before=E.parse(Path("/Users/zheranliu/Kidney_Manuscript_Packages/Genome_Medicine_终稿资料包_20261106")/'03_Publication_figures'/f'Figure{i}.svg')
 def point_geometry(tree):
  return Counter((node.tag,tuple(sorted((k,v) for k,v in node.attrib.items() if k!='style'))) for node in tree.xpath('//*[local-name()="circle" or local-name()="polygon"]'))
 assert point_geometry(before)==point_geometry(svg)

 png=Image.open(path.with_suffix('.png')); assert png.size==(math.ceil(w*300/72),math.ceil(h*300/72))
 with tempfile.TemporaryDirectory() as tmp:
  base=Path(tmp)/'page'
  subprocess.run(['pdftoppm','-singlefile','-r','300','-f',str(i),'-l',str(i),'-png',f/'Figures_1-5_preview.pdf',base],check=True)
  assert Image.open(base.with_suffix('.png')).convert('RGB').tobytes()==png.convert('RGB').tobytes()
 fonts=svg.xpath('//*[local-name()="text"]/@style')
 sizes=[float(re.search(r'font-size: ([.0-9]+)px',x)[1]) for x in fonts]
 assert min(sizes)>=7 and max(sizes)==11
 rows.append((i,round(w*25.4/72,2),round(h*25.4/72,2),*png.size,min(sizes),'identical','vector'))
with (root/'audit/final_output_check.csv').open('w',newline='') as file:
 wr=csv.writer(file);wr.writerow('figure width_mm height_mm png_width png_height min_pt preview_vs_png_300dpi svg'.split());wr.writerows(rows)
for p in [root/'Manuscript_audit.md',root/'Supplementary/Supplementary_index.md',root/'Figure_legends.md']:
 for x in re.findall(r'\]\(([^)]+)\)',p.read_text()):
  if not x.startswith('https://'):assert (p.parent/x).is_file(),(p,x)
for i in range(1,5):
 with zipfile.ZipFile(root/f'Supplementary/TableS{i}.xlsx') as z:assert z.testzip() is None
for i in range(1,11):assert len(PdfReader(root/f'Supplementary/FigureS{i}.pdf').pages)==1
for row in rows:print(row)

