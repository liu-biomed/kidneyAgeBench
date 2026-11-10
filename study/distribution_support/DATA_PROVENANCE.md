# Data provenance

The package contains summary statistics and annotations, not individual-level matrices. Source access and processing are described in Additional file 1 and the complete study code archive (Additional file 3).

| Object | Content and source |
|---|---|
| kab_gene_z | Gene-level age statistics for six GEO donor cohorts, TCGA KIRC/KIRP/KICH normal kidney, CPTAC ccRCC adjacent tissue and GTEx v8 cortex |
| kab_gene_z_comp | Corresponding age statistics adjusted for five marker-based composition proxies |
| kab_bins | CPTAC expression mean-by-SD matching bins |
| kab_bins_coprimary | Matching bins derived from G6 and TCGA |
| kab_library | 48 published gene sets; PMID, class, source label, gene and direction fields |
| kab_power_grid | Precomputed simulation grid used for interpolation |

GEO series: GSE36895, GSE60807, GSE147451, GSE66272, GSE93069 and GSE167093. Reference cohorts contribute 441 inferred GEO donors, 129 TCGA profiles and 75 CPTAC profiles. GTEx provides 85 separately analysed cortex samples. The GSE147451 donor identities include inferred bilateral pairs. Cohort overlap is not fully resolved for every historical signature; users must review the source publications and the manuscript's exclusions before interpreting a test as external validation.

The MIT code license does not replace the conditions attached to original public datasets or third-party gene-set databases. The historical HPA/GTEx overlap question and ERCB/RAAG independence limitation remain scientific provenance issues, not solved by software packaging.
