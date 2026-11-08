# kidneyAgeBench

Leakage-free benchmark of whether a gene set labelled "kidney ageing" or "senescence-related"
actually tracks chronological age in independent normal human kidney transcriptomes.

```r
# install.packages("kidneyAgeBench_0.2.0.tar.gz", repos = NULL, type = "source")
library(kidneyAgeBench)
kab_benchmark(c("CCR2", "VCAM1", "CSF1R", "ITGAM"))               # unsigned gene set
kab_benchmark(c("CXCL13", "CCL19", "EGF"), direction = c(1, 1, -1)) # signed signature
kab_benchmark(genes, exclude = "TCGA")   # drop datasets used to derive the signature
res <- kab_benchmark_library()           # re-run the full 48-set benchmark
```

Output: `enrichment` = mean |age z| of member genes divided by that of expression-matched random
genes (1 = no better than random); `emp_p` = empirical P; `pct_dir_correct` = share of signed
genes whose age association has the claimed direction.

Datasets bundled as gene-level statistics only (no individual-level data): six GEO donor cohorts
(441 donors), TCGA solid tissue normal (129), CPTAC ccRCC normal adjacent RNA (75), GTEx v8 kidney
cortex (85; post-mortem, excluded from the default pool).
