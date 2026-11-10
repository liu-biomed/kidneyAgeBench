"""Cross-check key numbers quoted in Manuscript.md and Table 2 against Additional file 2.

Each check: (description, value read from Additional file 2, value quoted in text/Table 2, tolerance)."""
from pathlib import Path
import numpy as np
import pandas as pd

import argparse
parser = argparse.ArgumentParser(description='Check 77 published numerical assertions against Supplementary Tables S1-S6.')
parser.add_argument('workbook', type=Path)
F = parser.parse_args().workbook.resolve()
S1 = pd.read_excel(F, 'Table S1', header=3)
S2 = pd.read_excel(F, 'Table S2', header=3)
S3 = pd.read_excel(F, 'Table S3', header=3)
S4 = pd.read_excel(F, 'Table S4', header=3)
S5 = pd.read_excel(F, 'Table S5', header=3)
S6 = pd.read_excel(F, 'Table S6', header=3)

s3 = S3[S3.Unit == 'Gene set'].set_index('Archived set ID')
pool = S3[S3.Unit == 'Pooled gene source'].set_index('Signature name')
P = S4[S4['Dataset or pooled analysis'] == 'Pooled (G6, TCGA, CPTAC)']
pr = lambda sid, adj: P[(P['Archived set ID'] == sid) & (P.Adjustment.str.startswith(adj))].iloc[0]
gtex_r = lambda sid, adj='Sex': S4[(S4['Archived set ID'] == sid) & (S4['Dataset or pooled analysis'] == 'GTEx') & (S4.Adjustment.str.startswith(adj))].iloc[0]

# Deconvolution subset
DEC = S4[S4['Analysis block'].str.startswith('5.')]
dec_piv = DEC.pivot(index='Archived set ID', columns='Adjustment', values='ρ')

def c(cohort, sid, ep, pop, model):
    m_cond = (S5.Model == model) if model == 'Adjusted for donor age' else S5.Model.str.startswith(model)
    x = S5[(S5.Cohort.str.startswith(cohort)) & (S5['Archived set ID'] == sid) &
           (S5.Endpoint.str.startswith(ep)) & (S5.Population.str.startswith(pop)) & m_cond]
    assert len(x) == 1, (cohort, sid, ep, pop, model, len(x))
    return x.iloc[0]

A = s3[s3.Class == 'Kidney age-derived']
B = s3[s3.Class.str.startswith('General ageing')]
C = s3[s3.Class == 'Disease-derived']

# Attenuation checks
base_all = P[P.Adjustment.str.startswith('Sex')].set_index('Archived set ID')['ρ']
comp_all = P[P.Adjustment.str.startswith('Adjusted')].set_index('Archived set ID')['ρ']
att_all = 100 * (base_all - comp_all) / base_all

att_marker_g6 = 100 * (dec_piv['Unadjusted'] - dec_piv['Five marker proxies']) / dec_piv['Unadjusted']
att_music = 100 * (dec_piv['Unadjusted'] - dec_piv['KPMP deconvolution, MuSiC']) / dec_piv['Unadjusted']
att_nnls = 100 * (dec_piv['Unadjusted'] - dec_piv['KPMP deconvolution, non-negative least squares']) / dec_piv['Unadjusted']
att_b5 = 100 * (dec_piv['Unadjusted'] - dec_piv['KPMP deconvolution, Bisque (five subclasses)']) / dec_piv['Unadjusted']
att_bc = 100 * (dec_piv['Unadjusted'] - dec_piv['KPMP deconvolution, Bisque (broad classes)']) / dec_piv['Unadjusted']

checks = [
    # Screening & library
    ('Table S2 unique screened records', len(S2), 268, 0),
    ('Table S2 included disease sets', (S2['Included in library'] == 'Yes').sum(), 29, 0),
    # Gene-level benchmark
    ('kidney sets FDR<0.05', (A.FDR < 0.05).sum(), 9, 0),
    ('general ageing sets FDR<0.05', (B.FDR < 0.05).sum(), 1, 0),
    ('disease sets FDR<0.05', (C.FDR < 0.05).sum(), 5, 0),
    ('second draw kidney FDR<0.05', (A['Second draw (sensitivity): FDR'] < 0.05).sum(), 8, 0),
    ('second draw general FDR<0.05', (B['Second draw (sensitivity): FDR'] < 0.05).sum(), 2, 0),
    ('second draw disease FDR<0.05', (C['Second draw (sensitivity): FDR'] < 0.05).sum(), 5, 0),
    ('kidney median ER', A['Enrichment ratio'].median(), 1.23, .005),
    ('disease median ER', C['Enrichment ratio'].median(), 1.05, .005),
    ('disease ER<1', (C['Enrichment ratio'] < 1).sum(), 14, 0),
    ('disease median genes', C['Measured genes'].median(), 5, 0),
    ('VIF Rodwell', s3.loc['A01_Rodwell2004', 'VIF'], 68, .5),
    ('KAS-94 VIF FDR', s3.loc['A11_KAS94_original', 'VIF-adjusted FDR'], .028, .0005),
    ('VIF-adjusted significant sets', (s3['VIF-adjusted FDR'] < .05).sum(), 1, 0),
    ('Rowland rep ER', s3.loc['A06_Rowland2019_rep19', 'Enrichment ratio'], 1.23, .005),
    ('Rowland rep FDR', s3.loc['A06_Rowland2019_rep19', 'FDR'], .15, .005),
    ('Rowland disc comp ER', s3.loc['A05_Rowland2019_37', 'Composition-adjusted enrichment ratio'], 1.81, .005),
    ('Rowland rep comp ER', s3.loc['A06_Rowland2019_rep19', 'Composition-adjusted enrichment ratio'], 1.97, .005),
    ('comp-adjusted FDR<0.05 sets', (s3['Composition-adjusted FDR'] < .05).sum(), 5, 0),
    ('comp-adjusted disease FDR<0.05 sets', (C['Composition-adjusted FDR'] < .05).sum(), 0, 0),
    ('no-signal matched mean rejection', s3['Null rejection rate, no signal, matched test'].mean(), .093, .0005),
    ('no-signal VIF mean rejection', s3['Null rejection rate, no signal, VIF test'].mean(), .031, .0006),
    # Gene pools
    ('disease pool naive', pool.loc['Disease-derived genes', 'Pool: ratio versus whole genome'], 1.50, .005),
    ('disease pool matched', pool.loc['Disease-derived genes', 'Enrichment ratio'], 1.16, .005),
    ('disease pool comp', pool.loc['Disease-derived genes', 'Composition-adjusted enrichment ratio'], 1.05, .005),
    ('RAAG pool', pool.loc['CKD progression-associated RAAGs', 'Enrichment ratio'], 1.40, .005),
    ('RAAG pool comp', pool.loc['CKD progression-associated RAAGs', 'Composition-adjusted enrichment ratio'], 1.62, .005),
    ('RAAG pool null lower', pool.loc['CKD progression-associated RAAGs', 'Pool: matched-null range, lower'], 1.17, .005),
    ('RAAG pool null upper', pool.loc['CKD progression-associated RAAGs', 'Pool: matched-null range, upper'], 1.73, .005),
    # Score-level & composition attenuation
    ('min signed rho', base_all.min(), .13, .005),
    ('max signed rho', base_all.max(), .30, .005),
    ('median marker attenuation (%)', att_all.median(), 57, .6),
    ('G6 deconv marker attenuation (%)', att_marker_g6.median(), 57, .6),
    ('G6 deconv MuSiC attenuation (%)', att_music.median(), 48, .6),
    ('G6 deconv NNLS attenuation (%)', att_nnls.median(), 16, .6),
    ('G6 deconv Bisque5 attenuation (%)', att_b5.median(), 8, .6),
    ('G6 deconv Bisque class attenuation (%)', att_bc.median(), 9, .6),
    ('Rowland rep rho', pr('A06_Rowland2019_rep19', 'Sex')['ρ'], .28, .005),
    ('Rowland rep comp', pr('A06_Rowland2019_rep19', 'Adjusted')['ρ'], .22, .005),
    ('Rowland disc GTEx rho', gtex_r('A05_Rowland2019_37')['ρ'], .23, .005),
    ('Rowland rep GTEx rho', gtex_r('A06_Rowland2019_rep19')['ρ'], .21, .005),
    ('KAS-94 GTEx rho', gtex_r('A11_KAS94_original')['ρ'], .04, .005),
    ('Rodwell rho', pr('A01_Rodwell2004', 'Sex')['ρ'], .23, .005),
    ('Rodwell comp', pr('A01_Rodwell2004', 'Adjusted')['ρ'], .14, .005),
    ('SenMayo rho', pr('B05_SenMayo', 'Sex')['ρ'], .20, .005),
    ('SenMayo comp', pr('B05_SenMayo', 'Adjusted')['ρ'], .07, .005),
    ('RAAG rho comp', pr('A09_Perco2019_RAAG_CKD', 'Adjusted')['ρ'], .16, .005),
    # Y chromosome check in Table S4
    ('CPTAC Y-chr rho', S4[(S4.Code == 'CHR_Y') & (S4['Dataset or pooled analysis'] == 'CPTAC')]['ρ'].iloc[0], -.75, .005),
    ('TCGA Y-chr rho', S4[(S4.Code == 'CHR_Y') & (S4['Dataset or pooled analysis'] == 'TCGA')]['ρ'].iloc[0], .03, .005),
    ('GTEx Y-chr rho', S4[(S4.Code == 'CHR_Y') & (S4['Dataset or pooled analysis'] == 'GTEx')]['ρ'].iloc[0], -.01, .005),
    # CKD analyses
    ('ERCB Rodwell base', c('ERCB', 'A01_Rodwell2004', 'CKD', 'All', 'Baseline')['Estimate (all members)'], 2.86, .005),
    ('ERCB Rodwell adj', c('ERCB', 'A01_Rodwell2004', 'CKD', 'All', 'Adjusted')['Estimate (all members)'], .72, .005),
    ('ERCB Rowland rep adj', c('ERCB', 'A06_Rowland2019_rep19', 'CKD', 'All', 'Adjusted')['Estimate (all members)'], .95, .005),
    ('ERCB SenMayo adj', c('ERCB', 'B05_SenMayo', 'CKD', 'All', 'Adjusted')['Estimate (all members)'], -.78, .005),
    ('GSE175759 Rodwell base', c('GSE175759', 'A01_Rodwell2004', 'CKD', 'All', 'Baseline')['Estimate (all members)'], .82, .005),
    ('GSE175759 Rodwell adj', c('GSE175759', 'A01_Rodwell2004', 'CKD', 'All', 'Adjusted')['Estimate (all members)'], -.19, .005),
    ('GSE175759 de Magalhaes adj', c('GSE175759', 'B01_deMagalhaes2009', 'CKD', 'All', 'Adjusted')['Estimate (all members)'], .30, .005),
    ('GSE175759 RAAG adj', c('GSE175759', 'A09_Perco2019_RAAG_CKD', 'CKD', 'All', 'Adjusted')['Estimate (all members)'], .01, .005),
    ('ERCB RAAG injury', c('ERCB', 'A09_Perco2019_RAAG_CKD', 'Correlation with proximal tubular injury', 'All', 'Unadjusted')['Estimate (all members)'], .89, .005),
    ('ERCB SenMayo injury', c('ERCB', 'B05_SenMayo', 'Correlation with proximal tubular injury', 'All', 'Unadjusted')['Estimate (all members)'], .54, .005),
    ('GSE175759 SenMayo injury', c('GSE175759', 'B05_SenMayo', 'Correlation with proximal tubular injury', 'All', 'Unadjusted')['Estimate (all members)'], .85, .005),
    ('ERCB Rodwell injury trimmed', c('ERCB', 'A01_Rodwell2004', 'Correlation with proximal tubular injury', 'All', 'Unadjusted')['Estimate (shared genes removed)'], .845, .0005),
    ('GSE175759 SenMayo eGFR raw P', c('GSE175759', 'B05_SenMayo', 'Correlation with eGFR', 'All', 'Unadjusted')['Nominal P'], .056, .001),
    ('GSE175759 SenMayo eGFR adj', c('GSE175759', 'B05_SenMayo', 'Correlation with eGFR', 'All', 'Adjusted')['Estimate (all members)'], .02, .005),
    ('GSE175759 RAAG eGFR raw', c('GSE175759', 'A09_Perco2019_RAAG_CKD', 'Correlation with eGFR', 'All', 'Unadjusted')['Estimate (all members)'], -.36, .005),
    ('GSE175759 RAAG eGFR adj', c('GSE175759', 'A09_Perco2019_RAAG_CKD', 'Correlation with eGFR', 'All', 'Adjusted')['Estimate (all members)'], -.51, .005),
    ('GSE175759 RAAG eGFR CKD-only adj', c('GSE175759', 'A09_Perco2019_RAAG_CKD', 'Correlation with eGFR', 'CKD', 'Adjusted')['Estimate (all members)'], -.50, .005),
    # Transplant analyses
    ('TX Rodwell age-adj rho', c('GSE147451', 'A01_Rodwell2004', 'Correlation with lower', '151', 'Adjusted for donor age')['Estimate (all members)'], .25, .005),
    ('TX Rodwell age-adj q', c('GSE147451', 'A01_Rodwell2004', 'Correlation with lower', '151', 'Adjusted for donor age')['BH q (15 scores)'], .034, .0005),
    ('TX Rodwell noblood age-adj q', c('GSE147451', 'A02_Rodwell2004_noblood', 'Correlation with lower', '151', 'Adjusted for donor age')['BH q (15 scores)'], .040, .0005),
    ('TX MHCII age-adj q', c('GSE147451', 'A10_Sinning2023_MHCII', 'Correlation with lower', '151', 'Adjusted for donor age')['BH q (15 scores)'], .040, .0005),
    ('TX Rodwell full', c('GSE147451', 'A01_Rodwell2004', 'Correlation with lower', '151', 'Adjusted for donor age and')['Estimate (all members)'], .20, .005),
    ('TX Rodwell full q', c('GSE147451', 'A01_Rodwell2004', 'Correlation with lower', '151', 'Adjusted for donor age and')['BH q (15 scores)'], .174, .0005),
    ('TX Rowland rep full', c('GSE147451', 'A06_Rowland2019_rep19', 'Correlation with lower', '151', 'Adjusted for donor age and')['Estimate (all members)'], .06, .005),
    ('TX Kainz raw', c('GSE147451', 'A07_Kainz2007_donor', 'Correlation with lower', '151', 'Unadjusted')['Estimate (all members)'], .02, .005),
    # GTEx multi-tissue
    ('GTEx CellAge directional FDR<0.05 tissues', ((S6['Gene set'].str.startswith('CellAge inducers')) & (S6.Adjustment.str.startswith('Sex, death')) & (S6['Signed FDR'] < .05)).sum(), 29, 0),
]

bad = 0
for d, v, t, tol in checks:
    ok = abs(float(v) - t) <= tol + 1e-9
    bad += not ok
    if not ok:
        print(f"BAD {d}: workbook {float(v):.4f} vs text {t}")
print(f"Checked {len(checks)} quantitative assertions against Additional file 2; mismatches: {bad}")
assert bad == 0, f"{bad} mismatches found!"
