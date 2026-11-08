# Distribution checks; original package implementation is unchanged.
library(kidneyAgeBench)
primary <- kab_benchmark_library(n_null = 100, config = "locked")
no_cptac <- kab_benchmark_library(n_null = 100, config = "coprimary")
stopifnot(nrow(primary) == 48L, nrow(no_cptac) == 46L)
stopifnot(setequal(setdiff(primary$sig_id, no_cptac$sig_id),
                  c("A11_KAS94_original", "A12_KAS82_corrected")))
stopifnot(all(is.finite(primary$enrichment)),
          all(primary$emp_p > 0 & primary$emp_p <= 1))
stopifnot(all(!is.finite(kab_pooled_z(groups = "TCGA"))))
cat("PASS: primary 48 sets; no-CPTAC 46 sets; excluded KAS IDs and minimum-dataset rule verified.\n")
