# kidneyAgeBench

Reproducible benchmarking of kidney ageing gene sets against chronological age in human kidney reference transcriptomes.

kidneyAgeBench provides reusable R functions, reference summary statistics, examples and tests for benchmarking kidney ageing gene sets. The package is developed and maintained by Menghan Liu.

## Install

The study package is **kidneyAgeBench 0.2.0**, written in R (R >= 4.1). The package uses base/recommended R functionality and has no additional declared runtime package dependencies. A local check was run with R 4.4.2 on macOS arm64.

Install the fixed manuscript version directly from GitHub:

```r
install.packages("remotes") # only needed if remotes is not already installed
remotes::install_github(
  "liu-biomed/kidneyAgeBench",
  subdir = "kidneyAgeBench",
  ref = "manuscript-20261108"
)
```

The R package lives in the `kidneyAgeBench/` subdirectory. Alternatively, download the source installation archive from [Releases](https://github.com/liu-biomed/kidneyAgeBench/releases) and install it with `install.packages("kidneyAgeBench_0.2.0.tar.gz", repos = NULL, type = "source")`.

To build and test a cloned repository, run from the repository root:

```sh
R CMD build kidneyAgeBench
R CMD INSTALL kidneyAgeBench_0.2.0.tar.gz
Rscript examples/quick_start.R
Rscript kidneyAgeBench/tests/test-reproduce.R
Rscript tests/check_supported_configurations.R
```

The source under `kidneyAgeBench/` is identical to the study package. The README, examples, tests and workflow at repository root are distribution support files; they do not change the study algorithms, gene sets or stored results. The embedded package README and DESCRIPTION are historical metadata; the scope and limitations below clarify their broad wording.

## What the package does

- Tests a user-supplied HGNC gene list against bundled gene-level age statistics.
- Constructs expression-matched random gene sets and reports enrichment and empirical P values.
- Reports estimates before and after adjustment for five marker-based composition proxies.
- Allows exclusion of reference groups used to derive a signature; known KAS and Rowland exclusions are applied by the library helper.
- Provides directional concordance and interpolated power estimates from a precomputed simulation grid.

```r
library(kidneyAgeBench)
kab_benchmark(c("CCR2", "VCAM1", "CSF1R", "ITGAM"), n_null = 2000, seed = 1)
kab_benchmark(c("CXCL13", "CCL19", "EGF"), direction = c(1, 1, -1))
kab_benchmark_library(n_null = 2000, config = "locked")
```

The library contains 48 gene sets. The default configuration pools six GEO donor cohorts, TCGA and CPTAC. GTEx is supplied separately and is excluded from this default pool. Only gene-level statistics and gene annotations are bundled; the repository does not contain individual-level expression matrices or clinical records.

## Interpretation and known limits

This is a benchmarking resource, not a clinical biological-age predictor. It produces no validated patient-level age prediction, clinical cut-off or treatment recommendation. Group exclusions reduce known discovery overlap; they do not establish that every historical cohort is independent.

`config = "coprimary"` is the historical API name for the sensitivity analysis without CPTAC. It should not be interpreted as a newly designated co-primary analysis in the current manuscript. It still includes TCGA tumour-adjacent normal tissue. The historical DESCRIPTION phrase “without tumour-adjacent tissue” is therefore too broad.

The pooling rule requires measurements in at least two eligible datasets. Consequently, after excluding G6 for KAS and CPTAC through `config = "coprimary"`, only TCGA remains: KAS-94 and KAS-82 are not evaluated. The library helper returns 48 sets for `locked` and 46 for `coprimary`, and omits sets that cannot be evaluated. Inspect row counts and IDs. Do not treat a missing row as a null association. Direct single-group use of `kab_pooled_z()` can return no finite statistics under this rule.

Power and minimum-detectable-effect values are interpolated from a stored grid (set sizes 2–100); they are not fresh simulations for a new dataset. Values outside that size range are clamped by the underlying interpolation. The simple binomial directional P value assumes independent genes. Expression matching alone does not correct within-set co-expression.

**Co-expression/VIF calibration, deconvolution, individual-sample signature scoring and clinical analyses are in the study analysis scripts, not functions implemented by this package.** The complete study scripts, figure inputs and data-access instructions are supplied as Additional file 3. The GitHub package is a reusable component of that study archive.

## Reproduction and versioning

The bundled reproduction test compares four representative enrichment ratios with fixed reference results using the original 5% Monte Carlo tolerance. This is a numerical regression test, not a full rerun from raw expression matrices. The root configuration test checks library coverage and documents the exclusion behavior.

The fixed manuscript snapshot uses release tag `manuscript-20261108` and contains package version `0.2.0`. The tag identifies the repository snapshot; it does not denote a change to the study algorithms. A version-specific archival DOI will be added after Zenodo archiving is complete.

## License and citation

The original package declares **MIT + file LICENSE**. Its copyright record is preserved under `kidneyAgeBench/LICENSE`; the full MIT text is provided as `LICENSE`. Derived gene-level data and gene annotations retain the attribution and source conditions of their original datasets and publications; the code license does not supersede third-party rights. See `DATA_PROVENANCE.md`.

Use `CITATION.cff` to cite the software. Contact: Menghan Liu, lmhdoc@163.com, as recorded in the original package DESCRIPTION.
