# Preserve quantitative SVG elements. Move only text and text-side endpoints of callouts.
from pathlib import Path
from xml.etree import ElementTree as ET
from collections import Counter
import re,subprocess
src=Path("/Users/zheranliu/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030/benchmark/figures/v7/Figure8_tissue_source.svg")
s=src.read_text()
a="x='227.28' y='223.46'";assert s.count(a)==1
s=s.replace(a,"x='305.31' y='223.46' text-anchor='end'")
s=s.replace(" textLength='78.03px' lengthAdjust='spacingAndGlyphs'>RAAG 2019 CKD-concordant",">CKD-prog. RAAGs")
for label,y,old_end in [('GenAge',299,'286.52'),('CellAge',289,'286.48'),('SenMayo',279,'286.43')]:
 pat=r"<text x='340.09' y='289.04'([^>]+)>"+label+r"</text>"
 s,n=re.subn(pat,lambda m:f"<text x='340.09' y='{y}'"+m[1]+">"+label+"</text>",s)
 assert n==1
 s=s.replace(f"x1='338.29' y1='{old_end}'",f"x1='338.29' y1='{y-2.2}'")
points=lambda x:Counter(ET.tostring(node) for node in x.iter() if node.tag.endswith(('circle','polygon')))
assert points(ET.fromstring(s))==points(ET.parse(src).getroot())
p=Path('Supplementary/FigureS9.svg');p.write_text(s)
subprocess.run(['rsvg-convert','-f','pdf','-o',p.with_suffix('.pdf'),p],check=True)
subprocess.run(['pdftoppm','-singlefile','-png','-r','160',p.with_suffix('.pdf'),'audit/pages/FigureS9'],check=True)

