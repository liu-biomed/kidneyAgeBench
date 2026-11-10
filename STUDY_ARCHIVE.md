# Complete study source code and reproducibility materials

Materials for **Ageing or injury? What transcriptomic ageing scores measure in kidney biopsies: a benchmark of 48 signatures**.

The complete study material is available in [`study/`](study/) and as a single ZIP under [release study-20261110](https://github.com/liu-biomed/kidneyAgeBench/releases/tag/study-20261110). The same ZIP is archived on [Zenodo](https://doi.org/10.5281/zenodo.24114198). This release expands the distributed study material; it does not change the analysis algorithms or reported results.

## Contents and use

The 202 original files include analysis and verification scripts, the recorded analysis plan, locked gene-library snapshots, fixed figure inputs, validation records, and kidneyAgeBench 0.2.0. They are preserved byte-for-byte from the complete study archive. The existing reusable R package remains in `kidneyAgeBench/` for installation using the repository's original instructions.

Run commands from the `study/` directory. See [`study/README.md`](study/README.md) for reproduction entry points, [`study/ENVIRONMENT.md`](study/ENVIRONMENT.md) for dependencies, and [`study/PATH_CONFIGURATION.tsv`](study/PATH_CONFIGURATION.tsv) for historical paths that require local configuration. Raw expression matrices are not bundled and must be retrieved from the source repositories. The tables workbook, Supplementary File S2, accompanies the manuscript separately and is required for the stored-number verification script.

Some figure inputs contain de-identified sample-level values derived from public datasets. These are not new clinical enrollment records. The R package itself contains gene-level summary statistics. Package checks and stored-output comparisons do not constitute a complete rerun from raw expression matrices.

## Version and historical labels

The full-study version is **study-20261110**, archived at https://doi.org/10.5281/zenodo.24114198. Its release ZIP is named `kidneyAgeBench-study-20261110.zip`; Zenodo retains the original filename `Supplementary_File_S3_Source_Code.zip`. Both contain the identical archive bytes (SHA256 `4069abd6bd4b010f8d0b4755633ba5ac346aaafeed99edb7dffa158835c72811`).

Internal references to “Additional file 3” or “Supplementary File S3” are historical labels for this study archive, retained with the original files for provenance. The current manuscript cites the public repository archive instead of listing S3 as a journal attachment. “Additional file 1” and “Additional file 2” refer to the manuscript's Supplementary Files S1 and S2.

The earlier `manuscript-20261108` release and DOI [10.5281/zenodo.24046928](https://doi.org/10.5281/zenodo.24046928) remain the fixed **R-package repository snapshot**. That earlier archive does not contain the full study scripts and plotting inputs. Its original tag and files are unchanged.

## Licenses and rights

Existing component-specific licensing is preserved. The kidneyAgeBench R package is distributed under MIT + file LICENSE. The repository-root MIT license continues to cover the originally licensed software and distribution support; it does not newly license the contents of `study/` as a whole. Third-party datasets, derived data and gene-set resources retain their original source terms and attribution requirements. No additional blanket reuse license is granted for historical scripts or other materials that lack an explicit license. See the included component license files and the Licenses section of `study/README.md`.

## Citation

Liu M, Luo Y, Liu B, Li G (2026). *Complete study source code and reproducibility materials for Ageing or injury? What transcriptomic ageing scores measure in kidney biopsies: a benchmark of 48 signatures*. Version study-20261110. Zenodo. https://doi.org/10.5281/zenodo.24114198
