from pathlib import Path
from xml.etree import ElementTree as ET
from zipfile import ZipFile
import hashlib, csv, re
P=Path("/Users/zheranliu/Kidney_Manuscript_Packages/Genome_Medicine_终稿资料包_20261106")
D=Path("/Users/zheranliu/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030")
def digest(file):
 h=hashlib.sha256()
 with file.open('rb') as stream:
  for chunk in iter(lambda:stream.read(2**20),b''):h.update(chunk)
 return h.hexdigest()
def save(path, cols, rows):
 with open(path,"w",newline="") as f:
  w=csv.writer(f);w.writerow(cols);w.writerows(rows)
def pars(path):
 if path.suffix=='.docx':
  with ZipFile(path) as z:
   doc=ET.fromstring(z.read('word/document.xml'))
   return [''.join(t.text or '' for t in n.findall('.//{*}t')) for n in doc.findall('.//{*}p')]
 return [l.strip('# *') for l in path.read_text().splitlines() if l.strip()]
files=[(P/'01_Manuscript/Manuscript.docx','Working title/content-matched DOCX; exact identity with user-named suffix file unproven'),(P/'01_Manuscript/Manuscript.md','paired working source'),(D/'genome_medicine_polished/text/manuscript_rewritten.md','latest matching edited text'),(D/'genome_medicine_polished/text/manuscript_redesign.md','prior main text'),(D/'genome_medicine_polished/text/Figure_legends.md','earlier captions'),(P/'01_Manuscript/Figure_legends.md','matching captions'),(Path.home()/'Downloads/manuscript(1).docx','name match only; excluded'),(Path.home()/'Library/Mobile Documents/com~apple~CloudDocs/sandbox/杂事/杂事归档/20260607/高琨师妹论文/修订版/Manuscript(1).docx','name match only; excluded')]
versions=[]
for path,role in files:
 s=pars(path);title=next((v for v in s if v.startswith('title:')),s[0])
 def sample(header):
  if 'excluded' in role:return ''
  pos=next((i for i,l in enumerate(s) if l==header),None)
  return s[pos+1][:210] if pos is not None and pos+1<len(s) else ''
 fig=next((v for v in s if v.startswith('Figure 1.')),'')
 versions.append([str(path),role,title[:180],sample('Abstract'),sample('Discussion'),fig[:170],digest(path)])
save("version/version_table.csv",['path','role','title_or_first_text','abstract_opening','discussion_opening','figure1_legend','sha256'],versions)
# Compare physical contents to the full preceding release manifest.
known=list(csv.DictReader((P/'05_Other/Inventory/package_files.tsv').open(),delimiter='\t'))
locked=[]
for row in known:
 file=P/row['relative_path'];h=digest(file);assert h==row['sha256'],file
 locked.append((str(file),row['size_bytes'],h))
for f in list((D/'scripts').glob('*.R'))+list((D/'benchmark/results').glob('*.csv'))+list((Path.home()/'kidney_aging_data/harmonized').glob('*.rds')):
 locked.append((str(f),f.stat().st_size,digest(f)))
save("version/input_lock.csv",['path','bytes','sha256'],locked)
for i in range(1,6):
 for ext in ('pdf','svg','png'):
  assert digest(P/f'03_Publication_figures/Figure{i}.{ext}')==digest(D/f'genome_medicine_polished/figures/Figure{i}.{ext}')

