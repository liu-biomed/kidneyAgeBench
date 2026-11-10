# Locked addendum 5: benchmark of published kidney-ageing gene sets

Version 1.0, 2026-10-31. Locked after the signature library was assembled (`benchmark/signatures/signature_library_long.csv`, 34 gene sets) and **before** any signature in this library was scored against age in any dataset.

Earlier work already scored some of these sets against age in some datasets:
- Rodwell, de Magalhães and Rowland, in the discovery cohorts (`R1_published_signatures.csv`);
- SenMayo, in all datasets;
- KAS-94 and KAS-82, in all datasets.

These prior results are disclosed. Everything else is new.

## Datasets and leakage rules

Every dataset is analysed **separately** with the same code.

| Dataset | n | Covariates |
|---|---|---|
| G6: six eligible GEO discovery cohorts (corrected set) | 441 donors | sex where available; pooled by random effects |
| TCGA kidney Solid Tissue Normal | 129 | sex, project |
| CPTAC ccRCC NAT RNA | 75 | sex |
| GTEx Kidney Cortex | 85 | sex; reported separately, as a post-mortem tissue |

**Leakage-free evaluation set per signature.** A dataset is excluded from a signature's evaluation if it was used to derive that signature:
- KAS-94 and KAS-82: exclude G6.
- Rowland 2019 (A05, A06): exclude TCGA, which was its discovery cohort.
- C-class models built on TCGA tumours: TCGA normals are a different tissue type, so they are retained, but this is flagged.
- No other overlaps are known; Rodwell and Melk data are not in any evaluation dataset.

**Primary pooled estimate.** Random-effects pooling over the leakage-free datasets among G6 (as one pooled estimate), TCGA and CPTAC. GTEx is reported separately.

## Metrics

**M-signed (signed sets only: up and/or down genes).**
- Score = mean within-dataset z of up genes minus mean z of down genes. One-directional sets use the mean only.
- Statistic: partial Spearman ρ with age, as in `pspear`.
- Two versions: base, and composition-adjusted (locked marker proxies from `04_locked_definitions.R`).
- A signature is "replicates" if the primary pooled base ρ > 0 with P < 0.05.

**M-gene (all 34 sets, signed or not).** Tests whether member genes are more age-associated than chance.
- Gene-level age statistic in each leakage-free dataset: Fisher z of the partial correlation with age.
- Pooled across the leakage-free datasets by inverse-variance weighting.
- Statistic: mean |z| of measured members.
- Null: 2,000 random gene sets of equal size, matched to members on expression mean and variance deciles (computed in CPTAC NAT RNA).
- Outputs: empirical P and enrichment ratio (observed / null mean). FDR (BH) across the 34 sets.

**M-direction (signed sets).** Proportion of measured members whose pooled validation z has the claimed sign, with a one-sided binomial test.

## Pre-specified expectations (to be reported as confirmed or not)

- **E1.** Among C-class sets (16 disease-model "ageing/senescence-related" gene sets), the median M-gene enrichment ratio is < 1.2, and fewer than half have FDR < 0.05.
- **E2.** Human-kidney age-derived A-class sets with signed directions (A01–A06, A12) replicate (M-signed) more often than C-class sets are age-informative (M-gene FDR < 0.05).
- **E3.** For signed sets, the composition-adjusted ρ is attenuated relative to base, with median attenuation > 30%.

## Controls

- Negative control: 1,000 random 50-gene sets (M-signed with random directions) give a null ρ distribution.
- Housekeeping control: Eisenberg housekeeping genes are not used, because they are not available offline. The matched-null procedure serves this role.
