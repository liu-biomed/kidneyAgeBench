# PubMed screening log (2026-10-31), C-class "ageing/senescence-related" kidney gene models

Search: PubMed E-utilities, query stored in `~/kidney_aging_data/lit/pubmed_search.json` (key "C"); 2010-01-01 to 2026-10-31, English.
- Records identified: PubMed 200; plus Europe PMC/OpenAlex records screened earlier (Europe PMC queries in `/tmp/epmc*.py`, archived in benchmark/).
- Automated pre-screen (`screen.py`): removed reviews/editorials, non-kidney, animal-only, ncRNA-only, other cell-death/metabolic labels (autophagy, anoikis, ferroptosis, cuproptosis, pyroptosis, m6A, glycolysis) → 133 records retained for manual screening.
- Manual screening criteria: (i) primary study; (ii) gene set explicitly framed/selected as ageing- or senescence-related; (iii) derived from human kidney tissue transcriptomics (kidney disease, transplant, or RCC); (iv) final gene list (≥2 genes) reported in abstract/main text.
- Included overall: 29 C-class gene sets (16 initial + 9 extension 1 + 4 extension 2).
- Excluded at full-text/abstract stage with reason (examples): 42219284 and 41005317 (derived from blood/sepsis or SLE PBMC data, not kidney tissue); 40243723 (single gene); 37800815, 41811638, 41770477, 38605343, 37503331, 39747640, 38304805, 36189255-type cluster-only studies where final genes only in figures (36189255 included after full-text extraction); 37708239 (lncRNA); 41620921 (telomere-labelled, not ageing-labelled); clinical/proteomic "biological age" models without gene lists.
- Databases: PubMed (primary, systematic) and Europe PMC (supplementary). Web of Science, Scopus, Embase and CNKI were not searched (author decision, 2026-10-31); stated as a limitation in the manuscript.
