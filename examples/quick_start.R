library(kidneyAgeBench)
cat("kidneyAgeBench", as.character(packageVersion("kidneyAgeBench")), "\n")
print(kab_benchmark(c("CCR2", "VCAM1", "CSF1R", "ITGAM"), n_null = 2000, seed = 1))
print(kab_benchmark(c("CXCL13", "CCL19", "EGF"), direction = c(1, 1, -1), n_null = 2000, seed = 1))
