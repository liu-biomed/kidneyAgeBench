# Analysis environments

The recorded study environment used R 4.4.2, metafor 4.8.0, limma 3.62.2, DESeq2 1.46.0, edgeR 4.4.0, msigdbr 24.1.0 and org.Hs.eg.db 3.20.0. These are reported historical versions, not a claim that all scripts were re-run under a newly reconstructed environment.

The reusable kidneyAgeBench package declares R >= 4.1 and no extra runtime package dependencies. The new package checks were run with R 4.4.2 on macOS arm64. Existing full session records are in `verification_scripts/clinical/R_session.txt` and `verification_scripts/deconvolution/sessionInfo.txt`.

The following package names were statically detected in R library/require calls (including historical and exploratory scripts; this is not a minimal or fully versioned dependency lock):

AnnotationDbi, DESeq2, Matrix, data.table, edgeR, fgsea, ggbeeswarm, ggnewscale, ggplot2, ggrepel, grid, jsonlite, kidneyAgeBench, limma, metafor, mgcv, msigdbr, n_null, org.Hs.eg.db, pROC, parallel, patchwork, readxl, xml2

Figure rendering additionally requires svglite, Arial fonts, rsvg-convert and pdftoppm. The numerical workbook checks require Python with pandas, numpy and openpyxl. Other Python dependencies should be read from the script used. No synthetic renv.lock file has been created; a fully pinned raw-data environment would require running the complete analysis under a freshly resolved installation.
