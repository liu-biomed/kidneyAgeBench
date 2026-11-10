from pathlib import Path
import csv, shutil
from openpyxl import load_workbook
from openpyxl.styles import Font, PatternFill
from openpyxl.utils import get_column_letter
P=Path("/Users/zheranliu/Kidney_Manuscript_Packages/Genome_Medicine_终稿资料包_20261106")/'02_Supplementary'
out=Path('Supplementary')
def rows(path):return list(csv.DictReader(open(path)))
def asvalue(x):
 if x in ('','NA'):return None
 if x in ('TRUE','FALSE'):return x=='TRUE'
 try:return int(x)
 except ValueError:
  try:return float(x)
  except ValueError:return x
def extra(w,title,source,caption):
 s=w.create_sheet(title);records=rows(source);cols=list(records[0]);s.append(cols)
 for r in records:s.append([asvalue(r[k]) for k in cols])
 for k,col in enumerate(cols,1):
  cell=s.cell(1,k);cell.font=Font(bold=True,color='FFFFFF');cell.fill=PatternFill('solid',fgColor='315A67')
  s.column_dimensions[get_column_letter(k)].width=min(52,max(13,len(col)+2))
 s.freeze_panes='B2';s.auto_filter.ref=s.dimensions
 w['README'].append([title,caption])
def readmefix(w,tab,msg):
 for row in w['README']:
  if row[0].value==tab:row[1].value=msg;return
 w['README'].append([tab,msg])
w=load_workbook(P/'TableS1.xlsx')
audit=rows('raag/RAAG_provenance_audit.csv'); by={r['current_gene']:r for r in audit}
s=w['S1b_gene_lists']
for i,head in enumerate(('A09_age_CKD_relationship_Fig3','A09_source_CKD_FC','A09_age_CKD_concordant','A09_direction_basis'),7):s.cell(1,i,head)
n=0
for row in list(s.values)[1:]:
 if row[0]=='A09_Perco2019_RAAG_CKD':
  n+=1
assert n==31
for row in range(2,s.max_row+1):
 if s.cell(row,1).value!='A09_Perco2019_RAAG_CKD':continue
 r=by[s.cell(row,5).value];assert int(r['current_signed_score_direction'])==s.cell(row,6).value
 for i,k in enumerate(('figure3_age_ckd_relationship','source_ckd_fc','is_source_concordant23','source_age_direction_basis'),7):s.cell(row,i,asvalue(r[k]))
readmefix(w,'S1b_gene_lists','Original and mapped members. dir is the source-defined scoring polarity (+1 or −1; 0 = unsigned). For A09 it is the CKD progression direction, not the age direction. Additional A09 columns trace the age-direction category to Figure 3 and the progression fold change to Table 3 of the source article.')
extra(w,'RAAG_member_audit','raag/RAAG_provenance_audit.csv','A09, all 31 members, scoring directions, source age-direction categories and measured coverage by cohort.')
w.save(out/'TableS1.xlsx')
w=load_workbook(P/'TableS2.xlsx')
for row in w['GenePools_Comp_and_CKD']:
 for c in row:
  if c.value=='RAAG, CKD-concordant':c.value='CKD progression-associated RAAGs'
  if c.value=='RAAG, non-CKD':c.value='Other RAAGs'
extra(w,'VIF_null_calibration','calibration/false_positive_by_set.csv','Post hoc competitive-null calibration; 1,000 experiments per scenario, 2,000 matched draws. These are simulated rejection frequencies.')
w.save(out/'TableS2.xlsx')
w=load_workbook(P/'TableS3.xlsx');s=w['KPMP_Deconv_Concordance']
methods=['base','marker5','MuSiC','Bisque5','Bisque_class','NNLS']
v={(r['sig_id'],r['method']):r for r in rows('deconvolution/pooled_effects.csv') if r['scope']=='common_five'}
for i in range(2,s.max_row+1):
 id=s.cell(i,1).value
 for j,method in enumerate(methods,2):s.cell(i,j).value=float(v[(id,method)]['rho'])
 s.cell(i,8).value=int(v[(id,'base')]['donors_actually_used'])
s.cell(1,8,'common_profiles')
readmefix(w,'KPMP_Deconv_Concordance','Corrected same-profile comparison; actual linked n and model rank in weights. See before/after sheet for the old results.')
extra(w,'Deconv_before_after','deconvolution/before_after.csv','Historical vs corrected results; estimated fraction inputs averaged for inferred bilateral units.')
extra(w,'Deconv_by_scope','deconvolution/pooled_effects.csv','Matched five-cohort and maximum-available method estimates; historical fits retained with scope labels.')
extra(w,'Clinical_audit','clinical/all_clinical_effects.csv','Primary and post hoc fits, CKD-only, overlap trimming, arm threshold and paired-donor exclusions. q = covariate count; q_BH15 = adjusted P across 15 scores in the named family.')
extra(w,'Gene_proxy_overlap','clinical/overlap_all48_by_cohort.csv','All 48 sets, overlap with each proxy, by cohort.')
extra(w,'Reference_overlap','calibration/reference_overlap_pooled.csv','Reference-tissue post hoc check removing shared genes from signatures, with fixed proxies.')
w.save(out/'TableS3.xlsx')
shutil.copyfile(P/'TableS4.xlsx',out/'TableS4.xlsx')
shutil.copyfile('raag/RAAG_provenance_audit.csv','RAAG_provenance_audit.csv')

