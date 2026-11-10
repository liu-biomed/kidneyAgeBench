from pathlib import Path
import shutil,re,subprocess,csv
rev=Path("/Users/zheranliu/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030")
pack=Path("/Users/zheranliu/Kidney_Manuscript_Packages/Genome_Medicine_终稿资料包_20261106")
P=Path('Supplementary');P.mkdir(exist_ok=True)
for i in [1,2,3,4,6,7,9]:
 shutil.copyfile(pack/f'02_Supplementary/FigureS{i}.pdf',P/f'FigureS{i}.pdf')
for i,name in [(5,'Figure3_scorecard'),(8,'Figure7_signed_replication')]:
 original=rev/f'benchmark/figures/v7/{name}.svg'
 text=original.read_text()
 assert text.count('RAAG 2019 CKD-concordant')==1
 text=re.sub(r"( textLength='[^']+' lengthAdjust='spacingAndGlyphs')>RAAG 2019 CKD-concordant</text>",'>CKD-prog. RAAGs</text>',text)
 assert 'CKD-concordant' not in text
 path=P/f'FigureS{i}.svg';path.write_text(text)
 subprocess.run(['rsvg-convert','-f','pdf','-o',path.with_suffix('.pdf'),path],check=True)
# Preserve panel A and B's drawing functions and source order from the originating supplement.
layout=rev/'redesign_v1/code/layout.R'
shutil.copyfile(layout,'code/layout_supplement.R')
original=(rev/'redesign_v1/code/supplement.R').read_text()
s=original[original.index('s10 <- function()'):]
s=s.replace('sdfile("S10_pools")','read.csv("source_data/S10_pools.csv")')
s=s.replace('sdfile("S10_deconvolution")','read.csv("source_data/S10_deconvolution.csv")')
s=s.replace('"RAAG, CKD subset"','"CKD-prog. RAAGs"').replace('Five cohorts contributed 414 donor records to the marker analysis.','The five-cohort comparison uses the same profiles for each method.')
s=s.replace('Deconvolution matched 403 records. Eleven bilateral identifiers did not match; pooling retained full cohort sizes.','414 profiles contributed for 12 signatures; 261 for tubular MHC II. Weights use actual n and model rank.')
s=s.replace('"supplement/FigureS10"','"Supplementary/FigureS10"')
Path('code/figureS10.R').write_text('source("code/layout_supplement.R")\n'+s)
pools=(rev/'redesign_v1/source_data/S10_pools.csv').read_text().replace('RAAG, CKD-concordant','CKD progression-associated RAAGs').replace('RAAG, non-CKD','Other RAAGs')
Path('source_data/S10_pools.csv').write_text(pools)
rows=[r for r in csv.DictReader(open('deconvolution/pooled_effects.csv')) if r['scope']=='common_five']
methods={'G6_marker5':'marker5','G6_KPMP_MuSiC':'MuSiC','G6_KPMP_Bisque_5sub':'Bisque5'}
with open('source_data/S10_deconvolution.csv','w',newline='') as f:
 wr=csv.writer(f);wr.writerow(['sig_id',*methods])
 for id in sorted({r['sig_id'] for r in rows}):
  wr.writerow([id,*[next(r['rho'] for r in rows if r['sig_id']==id and r['method']==m) for m in methods.values()]])
subprocess.run(['Rscript','code/figureS10.R'],check=True)
for i in (5,8,10):
 subprocess.run(['pdftoppm','-singlefile','-png','-r','160',P/f'FigureS{i}.pdf',f'audit/pages/FigureS{i}'],check=True)


subprocess.run(["python3","code/supplement_s9.py"],check=True)
