# Manuscript snapshot: 8 November 2026

This release contains kidneyAgeBench **0.2.0**, the R package used in the accompanying kidney transcriptomic ageing benchmark.

- Preserves the study package's 23 source, data, documentation and test files without changes.
- Provides installation instructions, quick-start examples, configuration checks, citation metadata and a GitHub Actions workflow.
- Documents expression-matched enrichment, composition-adjusted estimates, directional concordance and stored power-grid interpolation.
- Clarifies that co-expression/VIF calibration, deconvolution and clinical analyses are supplied in the separate study analysis archive.
- Documents the historical `coprimary` option, the two-dataset minimum, and the resulting 48-set/46-set coverage.

Local validation with R 4.4.2 on macOS arm64 returned 0 errors, 0 warnings and one installed-size NOTE (5.2 MB). Four representative reproduction tests passed. Independent Linux validation with R 4.4.2 completed successfully: `R CMD check` returned `Status: OK`; package installation, quick-start examples and the 48-set/46-set configuration tests passed. [GitHub Actions validation](https://github.com/liu-biomed/kidneyAgeBench/actions/runs/42437757263).

The release tag identifies a fixed manuscript snapshot. Package version 0.2.0 and the underlying study results are unchanged.

Archived manuscript snapshot: [doi:10.5281/zenodo.24046928](https://doi.org/10.5281/zenodo.24046928). The archive contains package version 0.2.0 and the repository files at commit `fd35b57f823ce6badbfebecab422eae57ebf2306`.

Documentation correction: removed an unpublished article title and author list from the repository description and archived README. Package source, data and version 0.2.0 are unchanged.
