# Reproduces the published benchmark enrichment ratios within Monte Carlo error
library(kidneyAgeBench)
ref <- read.csv(system.file("extdata", "benchmark_reference_results.csv", package = "kidneyAgeBench"))
L <- kab_library
for (id in c("A01_Rodwell2004", "B03_GenAge", "C07_PMID40319112", "A12_KAS82_corrected")) {
  s <- L[L$sig_id == id, ]; ex <- if (grepl("^A1[12]_", id)) "G6" else character()
  r <- kab_benchmark(s$gene, s$dir, exclude = ex, n_null = 1000)
  e0 <- ref$enrichment[ref$sig_id == id]
  stopifnot(abs(r$enrichment - e0) / e0 < 0.05)
  cat(id, round(r$enrichment, 3), "vs", round(e0, 3), "\n")
}
