"""Read local GEO deposits, downloaded public-record headers and supplementary table inventories."""
import csv,gzip,hashlib,re,zipfile
from pathlib import Path
from collections import Counter
from xml.etree import ElementTree as E
from os import environ
d=Path(environ.get('KIDNEY_DATA',Path.home()/'kidney_aging_data'))
sources=Path('sources'); fields=[]; allrows=[]
def rows(p):
 opener=gzip.open if p.name.endswith('.gz') else open
 with opener(p,'rt') as f:
  dt={}; traits=[]
  for v in csv.reader(f,delimiter='\t'):
   if not v:continue
   if v[0]=='!series_matrix_table_begin':break
   if v[0]=='!Sample_characteristics_ch1':traits.append(v[1:])
   elif v[0] in ('!Sample_title','!Sample_geo_accession','!Sample_source_name_ch1'):dt[v[0]]=v[1:]
  rec=[]
  for i,gsm in enumerate(dt['!Sample_geo_accession']):
   attrs=dict(v.split(': ',1) for col in traits for v in [col[i]] if ': ' in v)
   rec.append(dict(gsm=gsm,title=dt['!Sample_title'][i],**attrs))
  return rec
def save(fn,rec):
 with open(fn,'w') as out:
  w=csv.DictWriter(out,fieldnames=list(dict.fromkeys(k for r in rec for k in r)));w.writeheader();w.writerows(rec)
for p in sorted([*d.glob('geo/GSE104954*'),*d.glob('geo/GSE147451*'),*d.glob('geo/GSE175759*'),*sources.glob('*_headers.txt')]):
 rr=rows(p)
 keys=list(dict.fromkeys(k for r in rr for k in r if k not in ('gsm','title')))
 fields.append(dict(source=p.name,records=len(rr),characteristic_keys='|'.join(keys),n_linkable_patient_age=sum('age' in r and r['age'].strip() not in ('NA','') for r in rr),n_linkable_patient_sex=sum('sex' in r and r['sex'].strip() not in ('NA','') for r in rr)))
 if '175759' in p.name:
  save('GSE175759_all90_records.csv',rr);save('GSE175759_three_excluded.csv',[r for r in rr if r['technical outlier']!='No'])
 if '147451' in p.name:
  link=list(csv.DictReader(open('transplant_linkage_audit.csv'))); lk={k:r for r in link for k in r['gsm'].split('+')}
  for r in rr:
   k=r['gsm']; group=lk.get(k);age=r.get('donor age','')
   r['audit_status']=('adult_outcome_available' if group and group['first_outcome'] else 'adult_outcome_missing') if group else ('donor_under18' if age not in ('','NA') and float(age)<18 else 'other_exclusion')
   r['inferred_group_gsms']='' if not group else group['gsm'];r['retained_label']='' if not group else group['first_outcome']
  save('GSE147451_192_raw_eligibility.csv',rr)
  print('GSE147451',Counter(r['audit_status'] for r in rr),Counter((r['audit_status'],r.get('gfr group','missing')) for r in rr))
save('GEO_sample_metadata_fields.csv',fields)
b=E.parse(sources/'GSE175759_biosamples.xml');bio=b.findall('.//BioSample');ax=Counter(a.get('attribute_name') for r in bio for a in r.findall('.//Attribute'))
print('GSE175759 BioSamples',len(bio),ax)
supp=[]
for p in sorted(sources.rglob('*.xlsx')):
 with zipfile.ZipFile(p) as z:
  names=[x.get('name') for x in E.fromstring(z.read('xl/workbook.xml')).iter() if x.tag.endswith('}sheet')]
  txt=' '.join(E.fromstring(z.read('xl/sharedStrings.xml')).itertext()) if 'xl/sharedStrings.xml' in z.namelist() else ''
  supp.append(dict(file=p.relative_to(sources),sheets=' | '.join(names),GSM_strings=len(re.findall(r'GSM\d{5,}',txt)),ERCB_tubule_title_strings=len(re.findall(r'H\d+-Tub-',txt))))
save('publication_workbook_inventory.csv',supp)
inputs=[*d.glob('geo/GSE104954*'),*d.glob('geo/GSE147451*'),*d.glob('geo/GSE175759*'),d/'harmonized/clinical_bundle.rds',d/'harmonized/harmonized_bundle_corrected.rds']
revision=Path(environ.get('PROJECT_REV',Path.home()/'Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030'))
inputs.extend(revision/'scripts'/n for n in ['01b_clinical.R','01c_corrected_set.R','04_locked_definitions.R','21_benchmark.R','51_genomemed_upgrades.R']);inputs.extend([revision/'kidneyAgeBench/data/kab_library.rda',revision/'benchmark/results/REVISION_clinical_CKD_and_transplant.csv'])
save('clinical_input_sha256.csv',[dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),bytes=p.stat().st_size) for p in inputs])
