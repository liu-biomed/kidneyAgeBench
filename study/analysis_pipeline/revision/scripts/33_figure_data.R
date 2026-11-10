# Descriptive data for the redesigned figures (no new inferential endpoints).
# Rebuilds the benchmark datasets exactly as 21_benchmark.R and extracts: donor ages per dataset, composition-proxy
# associations with age, per-sample Y-chromosome scores in men (post hoc field-effect check), and gene-level provenance.
suppressMessages({library(metafor)}); REV <- path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030")
source(file.path(REV, "scripts", "04_locked_definitions.R")); set.seed(20261031)
src <- readLines(file.path(REV, "scripts", "21_benchmark.R")); eval(parse(text = src[1:grep("^Z <- lapply", src) - 1]))
lib <- rbind(read.csv(file.path(REV, "benchmark/signatures/signature_library_long.csv")),
             read.csv(file.path(REV, "benchmark/signatures/signature_library_extension1.csv")),
             read.csv(file.path(REV, "benchmark/signatures/signature_library_extension2.csv")))
sexof <- function(d) if (is.null(d$base)) rep(NA, length(d$age)) else toupper(substr(as.character(d$base$sex), 1, 1))
ages <- do.call(rbind, lapply(names(D), function(k) data.frame(dataset = k, grp = D[[k]]$grp, age = D[[k]]$age, sex = sexof(D[[k]]))))
comp <- do.call(rbind, lapply(names(D), function(k) { d <- D[[k]]; do.call(rbind, lapply(names(LOCK$proxy), function(p) {
  v <- mz(d$e, LOCK$proxy[[p]]); cbind(dataset = k, grp = d$grp, proxy = p, pspear(d$age, v, d$base)) })) }))
Yg <- c("RPS4Y1","DDX3Y","EIF1AY","KDM5D","UTY","USP9Y","ZFY","NLGN4Y","TXLNGY")
ysc <- do.call(rbind, lapply(c("CPTAC","TCGA","GTEx"), function(k) { d <- D[[k]]; m <- which(sexof(d) == "M")
  data.frame(dataset = k, age = d$age[m], yscore = mz(d$e[, m], intersect(Yg, rownames(d$e)), 3)) }))
yst <- do.call(rbind, lapply(split(ysc, ysc$dataset), function(x) cbind(dataset = x$dataset[1], pspear(x$age, x$yscore, NULL))))
saveRDS(list(ages = ages, comp = comp, ysc = ysc, yst = yst, lib = lib), path.expand("~/kidney_aging_data/figure_data_v2.rds"))
print(yst[, c("dataset","n","rho","p")]); print(aggregate(age ~ dataset, ages, function(a) c(n = length(a), min = min(a), med = median(a), max = max(a))))
print(format(comp[comp$grp != "G6", c("dataset","proxy","rho","p")], digits = 2))
