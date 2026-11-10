# Long table of main-figure and sensitivity results; observational units and interval types stay separate.
from pathlib import Path
import csv
from collections import Counter
def read(p): return list(csv.DictReader(open(p)))
result=[]
ids={r['sig_id']:r['label'] for r in read('source_data/signed_key.csv')}
def put(fig,src,key,condition,metric,n,unit,value,lo='',hi='',interval='',p='',bh='',status='frozen plot input'):
 result.append([fig,src,key,ids.get(key,key),condition,metric,n,unit,value,lo,hi,interval,p,bh,status])
def src(name):return 'source_data/'+name+'.csv'
for r in read(src('F1_reference_n')):put('1A',src('F1_reference_n'),r['Var1'],'reference','count',r['Freq'],'profiles',r['Freq'])
for r in read(src('F1_folds')):put('1B',src('F1_folds'),r['signature'],r['heldout']+'; '+r['sig_rule'],'Spearman rho',r['n'],'held-out profiles',r['rho'])
for r in read(src('F1_pooled_KAS')):put('1C',src('F1_pooled_KAS'),r['signature'],r['analysis'],'pooled rho',r['k'],'cohorts',r['rho'],r['lo'],r['hi'],'95% effect CI',r['p'])
for r in read(src('F1_GTEx_estimate')):put('1D',src('F1_GTEx_estimate'),'KAS-94','GTEx','rho',r['n'],'profiles',r['spearman_rho'],r['boot_lo'],r['boot_hi'],'95% bootstrap CI',r['p'])
perm=read(src('F1_sex_perm'));obs=read(src('F1_concordance_observed'))[0]
put('1E',src('F1_sex_perm'),'KAS-94','within-sex permutation','concordant genes',obs['n_genes'],'measured genes',obs['concordant'],p=(1+sum(int(x['permutation_concordant'])>=int(obs['concordant']) for x in perm))/(len(perm)+1))
for r in read(src('F1_background')):
 for m in ['naive_ratio','matched_ratio']:
  matched=m.startswith('matched')
  put('1F',src('F1_background'),r['pool'],r['config'],m,r['n_genes'],'genes',r[m],r['lo'] if matched else '',r['hi'] if matched else '','matched-null range' if matched else '',r['emp_p'] if matched else '')
co=read(src('F2_signed_cohorts'))
for r in read(src('F2_48_sets')):put('2A',src('F2_48_sets'),r['sig_id'],'primary benchmark','enrichment',r['n_measured'],'measured genes',r['enrichment'],p=r['emp_p'],bh=r['fdr_all'])
for r in read(src('F2_F3_signed_intervals')):
 n=sum(int(float(c['n'])) for c in co if c['sig_id']==r['sig_id'] and c['use']=='Observed' and c['dataset']!='GTEx')
 put('2B;3C' if r['model']=='base' else '3C',src('F2_F3_signed_intervals'),r['sig_id'],r['model'],'pooled score-age rho',n,'eligible reference profiles across contributing cohorts',r['rho'],r['lo'],r['hi'],'95% meta-analytic CI',r['pooled_p'],status='reproduced from expression')
for r in co:
 if r['use']=='Observed':put('2C',src('F2_signed_cohorts'),r['sig_id'],r['dataset']+'; rank covariate df '+r['q'],'rho',r['n'],'profiles',r['rho'],r['lo'],r['hi'],'Fisher interval',r['p'])
for r in read(src('F2_VIF')):
 for m in ['fdr','fdr_vif']:put('2D',src('F2_VIF'),r['sig_id'],m,'BH q',48,'tested sets',r[m],bh=r[m])
for r in read(src('F2_power')):put('2E',src('F2_power'),'set size '+r['k'],'true ratio '+r['true_ER'],'power',1000,'simulated sets',r['power'])
primary_n={r['sig_id']:r['n_measured'] for r in read(src('F2_48_sets'))}
alt_file=Path("/Users/zheranliu/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030/benchmark/results/REVISION_gene_level_four_configs.csv")
alt_n={r['sig_id']:r['n'] for r in read(alt_file) if r['config']=='noCPTAC_evalbins'}
for r in read(src('F2_without_CPTAC')):
 for m in ['enrichment_primary','enrichment_noCPTAC']:put('2F',src('F2_without_CPTAC'),r['sig_id'],m,'ratio',primary_n[r['sig_id']] if m=='enrichment_primary' else alt_n.get(r['sig_id'],''),'measured genes',r[m])
for r in read(src('F3_age_proxies')):put('3A',src('F3_age_proxies'),r['proxy'],r['ds']+'; rank covariate df '+r['q'],'rho',r['n'],'profiles',r['rho'],r['lo'],r['hi'],'Fisher interval',r['p'])
for r in read(src('F3E_descriptive_rho')):put('3E',src('F3E_descriptive_rho'),r['sig_id'],r['cohort']+'; '+r['proxy'],'raw Spearman rho',r['n'],'inferred donors',r['rho'])
for r in read(src('F3_pools')):put('3D',src('F3_pools'),r['pool'],r['analysis'],'enrichment',r['n'],'genes',r['ER'],r['lo'],r['hi'],'matched-null range',r['emp_p'])
for r in read('clinical/all_clinical_effects.csv'):
 fig={'poor_function_group':'5D-F','eGFR':'4C,F','CKD_reference_difference':'4D,E','score_vs_PT_injury':'4B'}.get(r['endpoint'],'4; sensitivity')
 put(fig,'clinical/all_clinical_effects.csv',r['sig_id'],'; '.join(r[k] for k in ('cohort','population','model','variant','min_arm')),r['endpoint'],r['n'],'inferred donor' if r['cohort']=='GSE147451' else 'biopsy',r['estimate'],r['lo'],r['hi'],'pointwise t CI' if r['endpoint']=='CKD_reference_difference' else 'pointwise Fisher CI',r['p'],r['q_BH15'],'recomputed; population and variant identify primary vs sensitivity')
for r in read('deconvolution/pooled_effects.csv'):put('S10;3', 'deconvolution/pooled_effects.csv',r['sig_id'],r['scope']+'; '+r['method'],'pooled score-age rho',r['donors_actually_used'],'profiles across cohorts',r['rho'],r['lo'],r['hi'],'meta-analytic CI',r['p'],status='recomputed')
refrows=read('calibration/reference_overlap_by_cohort.csv')
for r in read('calibration/reference_overlap_pooled.csv'):
 key=r['sig_id'];included={x['dataset'] for x in co if x['sig_id']==key and x['use']=='Observed' and x['dataset']!='GTEx'};n=sum(int(x['n']) for x in refrows if x['id']==key and x['ds'] in included and x['model']==r['model'] and x['rule']==r['rule'])
 put('3C; sensitivity','calibration/reference_overlap_pooled.csv',key,r['model']+'; '+r['rule'],'pooled score-age rho',n,'profiles; cohort-specific n in reference_overlap_by_cohort.csv',r['rho'],r['lo'],r['hi'],'95% effect CI',r['p'],status='recomputed')
for r in read('calibration/false_positive_by_set.csv'):put('2D; calibration','calibration/false_positive_by_set.csv',r['sig_id'],r['scenario']+'; '+r['test'],'null rejection rate',r['replicates'],'simulated experiments',r['rate'],status='post hoc calibration')

# Counts and marginal patient examples are keyed to the exact rows shown, not reconstructed from graphics.
import statistics
donors=read(src('F3_donor_samples'))
for k in ['immune','PT_content']:
 put('3B',src('F3_donor_samples'),k,'GSE147451','marginal score-age Spearman rho',len(donors),'inferred donors',statistics.correlation([float(r['age']) for r in donors],[float(r[k]) for r in donors],method='ranked'))
for cohort,file in [('ERCB','F4_ercb_biopsies'),('GSE175759','F4_gse175759_biopsies')]:
 profiles=read(src(file))
 for k,v in Counter(r['diagnosis'] for r in profiles).items():
  put('4A' if cohort=='ERCB' else '4C',src(file),k,cohort,'profiles',v,'biopsies',v)
allunits=read(src('F5_donor_linkage'));graft=read(src('F5_donors'))
for unit,label,n in [('kidney profile','adult profiles',sum(int(r['kidneys']) for r in allunits)),('inferred donor','inferred donors',len(allunits)),('inferred donor','classified',len(graft)),('inferred donor','unclassified',len(allunits)-len(graft))]:
 put('5A',src('F5_donor_linkage'),label,'GSE147451','count',n,unit,n)
for k in ('good','poor'):
 group=[float(r['age']) for r in graft if r['outcome']==k]
 put('5B',src('F5_donors'),k,'GSE147451','median age; years',len(group),'inferred donors',statistics.median(group))
for r in read('raag/RAAG_verified_results.csv'):
 put(r['panel'],'raag/RAAG_verified_results.csv',r['sig_id'],r['model']+'; '+r['cohort'],r['endpoint']+'; '+r['unit'],r['n'],r['n_unit'],r['est'],r['lo'],r['hi'],r['interval'],r['p'],r['fdr'],'31-member provenance / original numerical source check')
with open('verified_results.csv','w',newline='') as f:
 w=csv.writer(f);w.writerow('figure source key display_name model_or_population measure n unit value lower upper interval_definition nominal_p bh_q evidence'.split());w.writerows(result)

