from pathlib import Path
import csv
R=Path("/Users/zheranliu/Kidney_Manuscript_Packages/Genome_Medicine_终稿资料包_20261106")
source=R/'04_Source_data/c_Datasets/03_Figure_source_tables'
for nm in ('F2_48_sets','sets_key','signed_key','F3_pools'):
 with (source/f'{nm}.csv').open() as f:
  rd=csv.DictReader(f); rows=list(rd);keys=rd.fieldnames
 for r in rows:
  if r.get('sig_id')=='A09_Perco2019_RAAG_CKD':
   field='label' if nm=='signed_key' else 'short_label'
   r[field]='CKD-prog. RAAGs'
   if 'full_set_label' in r:r['full_set_label']='CKD progression-associated RAAGs (31; signed by CKD progression)'
  if 'pool' in r:r['pool']={'RAAG, CKD-concordant':'CKD progression-associated RAAGs','RAAG, non-CKD':'Other RAAGs'}.get(r['pool'],r['pool'])
 with Path('source_data',f'{nm}.csv').open('w',newline='') as f:
  wr=csv.DictWriter(f,keys);wr.writeheader();wr.writerows(rows)

