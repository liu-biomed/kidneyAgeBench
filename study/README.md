# Study source code and reproducibility materials

Supporting code for **Ageing or injury? What transcriptomic ageing scores measure in kidney biopsies: a benchmark of 48 signatures**.

File names: this archive is Supplementary File S3. In the text below and in the dated scripts and validation logs, "Additional file 1", "Additional file 2" and "Additional file 3" denote Supplementary File S1 (methods and figures), Supplementary File S2 (tables, `Supplementary_File_S2_Tables.xlsx`) and Supplementary File S3 (this archive). Historical script headers that name an earlier target journal are retained unchanged for provenance.

## Published software

The kidneyAgeBench package (0.2.0) is available at https://github.com/liu-biomed/kidneyAgeBench. The fixed manuscript snapshot `manuscript-20261108` (commit `fd35b57f823ce6badbfebecab422eae57ebf2306`) is archived at https://doi.org/10.5281/zenodo.24046928. This version DOI refers to the reusable R package repository. The wider study analysis scripts and figure inputs are supplied in this Additional file 3 archive.

GitHub Actions with R 4.4.2 on Linux returned `R CMD check: Status: OK`, and installation, examples and supported-configuration tests passed. The 32 files in the Zenodo archive were checked against the fixed repository tag. See `distribution_support/publication_record.json` for exact publication and validation links.

## Reproduction routes

| Task | Entry point | Required inputs |
|---|---|---|
| Install the reusable R package | `R CMD INSTALL kidneyAgeBench_0.2.0.tar.gz` | R >= 4.1 |
| Run four representative benchmark checks | `Rscript kidneyAgeBench/tests/test-reproduce.R` | Installed package |
| Check primary and no-CPTAC coverage | `Rscript distribution_support/check_supported_configurations.R` | Installed package |
| Recreate Figures 1–5 | `python3 analysis_pipeline/rerender_figures.py NEW_OUTPUT_DIRECTORY` | Bundled figure inputs; R, svglite, Arial, rsvg-convert, pdftoppm |
| Check 77 reported values | `python3 verification_scripts/check_reported_numbers.py /path/to/Supplementary_File_S2_Tables.xlsx` | Additional file 2; pandas, numpy, openpyxl |
| Recompute from raw expression data | Original processing scripts and verification scripts | Retrieved raw data; original dataset layout; path configuration; relevant R packages |

Use a previously nonexistent output directory for figure recreation. The figure entry point resolves its own bundled paths. The input workbook is supplied separately as Additional file 2 and is not silently embedded in this archive.

## Included material

- `analysis_plan/`: recorded plan, original checksums, search records and three locked gene-library snapshots. `PACKAGED_FILE_MAP.json` maps historical names to bundled copies.
- `analysis_pipeline/`: original data processing and analyses, plus the portable main-figure entry point.
- `verification_scripts/`: clinical sensitivity, deconvolution, RAAG provenance and co-expression/null-calibration analyses.
- `source_data/`: fixed plotting inputs for the main figures and Figure S10. Some files contain de-identified sample-level values derived from public datasets. These are not new clinical enrollment records.
- `kidneyAgeBench/` and `kidneyAgeBench_0.2.0.tar.gz`: the unchanged study package and source archive.
- `distribution_support/`: current installation, scope and example documentation from the published GitHub repository.
- `validation/`: dated historical verification records. New checks are recorded separately under `validation/20261108/`.

## Software scope

The package tests gene lists using bundled gene-level summary statistics. It implements expression matching, composition-adjusted gene-level comparisons, directional tests and power-grid interpolation. Co-expression/VIF calibration, sample-level scoring, deconvolution and clinical analyses are implemented in the accompanying scripts, not in exported package functions. Refer to `distribution_support/README.md` for the exact 48-set/46-set coverage and the historical `coprimary` option.

## Raw-data reproduction

Original matrices are not bundled. Retrieve GEO, GDC/TCGA, CPTAC/PDC, GTEx v8 and KPMP data using Additional file 1 and the download scripts. The reference GEO series are GSE36895, GSE60807, GSE147451, GSE66272, GSE93069 and GSE167093. CKD series are GSE104954 and GSE175759; transplantation reuses GSE147451. Review discovery exclusions, inferred donor pairing and historical cohort overlap before analysis.

Historical analysis scripts retain their original project and input paths to preserve provenance. They are not advertised as a one-command raw-data pipeline. `PATH_CONFIGURATION.tsv` lists script files containing local path references. Set each script's project, input and output roots to a locally reconstructed data layout; the clinical verification accepts `PROJECT_REV` and `KIDNEY_DATA` environment overrides. Scripts referring to legacy results require those prerequisite outputs to be generated first. `ENVIRONMENT.md` records versions and dependencies without inventing a fully resolved environment lockfile.

## Licenses

The unchanged kidneyAgeBench package declares MIT + file LICENSE. Third-party datasets and gene-set resources retain their original source terms. Historical analysis scripts do not all carry a project-wide license; no new blanket license has been imposed on them. Authors should choose a license covering their own analysis scripts before publishing a separate analysis repository. This does not change the existing MIT license of kidneyAgeBench.

## Verification boundary

This distribution includes local package installation/checks, a four-set numerical regression test, library-coverage checks, 77 stored-output comparisons and figure recreation checks. These do not constitute a complete rerun from raw expression matrices and do not resolve source cohort identities or missing outcomes.
