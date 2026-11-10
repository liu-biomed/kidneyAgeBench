# Supplementary Figure S1. kidneyAgeBench reproduces the gene-level test for a new gene set in one call.
# Run with R_LIBS pointing to an installed kidneyAgeBench. (The earlier panel a moved to Figure 7c.)
source(file.path(path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030"), "scripts", "40_fig_common.R"))
suppressMessages(library(kidneyAgeBench))
## (b) package demo: one call per gene set; null rebuilt with the package's own data and procedure
demo <- list(list(id = "A05_Rowland2019_37", ex = "TCGA", nm = "Rowland 2019 (37 genes)"), list(id = "C07_PMID40319112", ex = character(), nm = "CKD ageing hub genes [40319112]"))
L <- kab_library; pz <- kab_pooled_z(); bins <- kab_bins
D <- do.call(rbind, lapply(demo, function(d) { s <- L[L$sig_id == d$id, ]; r <- kab_benchmark(s$gene, s$dir, exclude = d$ex, n_null = 2000, seed = 1)
  pz <- kab_pooled_z(exclude = d$ex); ok <- intersect(names(pz)[is.finite(pz)], names(bins)); g <- intersect(unique(toupper(s$gene)), ok); pool <- split(ok, bins[ok]); set.seed(1)
  nul <- replicate(2000, mean(abs(pz[vapply(g, function(x) { p <- pool[[bins[[x]]]]; p[sample.int(length(p), 1)] }, "")])))
  stopifnot(abs(mean(nul) - r$null_mean) < 1e-12)
  data.frame(set = d$nm, null = nul, obs = r$mean_abs_z, ratio = r$enrichment, p = r$emp_p, n = r$n_measured) }))
DS <- unique(D[, c("set","obs","ratio","p","n")]); D$set <- factor(D$set, levels = sapply(demo, `[[`, "nm")); DS$set <- factor(DS$set, levels = levels(D$set))
slab <- setNames(sprintf("%s\nn = %d; ratio %.2f; P = %.3f", DS$set, DS$n, DS$ratio, DS$p), as.character(DS$set))
pb <- ggplot(D, aes(null)) + geom_histogram(bins = 45, fill = "grey78") + geom_vline(data = DS, aes(xintercept = obs), colour = c(cls_col[["A"]], cls_col[["C"]]), linewidth = .6) +
  facet_wrap(~ set, ncol = 1, scales = "free", labeller = as_labeller(slab)) + scale_y_continuous(expand = expansion(c(0, .08))) + labs(x = "Mean |age z|: member genes (line) vs\n2,000 expression-matched random sets", y = "Random sets")
## (c) usage and recommended checks
code <- c("library(kidneyAgeBench)", "", "# test a proposed gene set", "kab_benchmark(", "  genes = c(\"SOD2\", \"FOS\", ...),", "  direction = NULL, # or +1/-1", "  exclude = \"TCGA\",", "  n_null = 2000)", "", "# re-run all 48 library sets", "kab_benchmark_library()")
chk <- c("Test in normal kidney not used\nto select the genes", "Compare with expression-\nmatched random genes", "Report composition-adjusted\nestimates", "State the tissue source; do not\nrely on adjacent or post-mortem\ntissue alone")
pc <- ggplot() + coord_cartesian(xlim = c(0, 10), ylim = c(0, 20), expand = FALSE, clip = "off") + theme_void() + theme(plot.tag = element_text(size = 9, face = "bold")) +
  annotate("rect", xmin = .1, xmax = 9.9, ymin = 11.2, ymax = 19.9, fill = "#F3F5F8", colour = "grey70", linewidth = .25) +
  annotate("text", x = .4 + .235 * nchar(sub("^( *).*", "\\1", code)), y = 19.5 - (seq_along(code) - 1) * .75, label = trimws(code, "left"), hjust = 0, vjust = 1, size = TX, family = "Courier New", colour = ifelse(grepl("^#", code), "grey45", "black")) +
  annotate("text", x = .2, y = 10.5, label = "Before calling a gene set a\nkidney ageing readout:", hjust = 0, vjust = 1, size = TXL, fontface = "bold", lineheight = .9) +
  annotate("text", x = .2, y = 8.6 - c(0, 2.1, 4.2, 6.3), label = 1:4, hjust = 0, vjust = 1, size = TX, fontface = "bold") +
  annotate("text", x = .75, y = 8.6 - c(0, 2.1, 4.2, 6.3), label = chk, hjust = 0, vjust = 1, size = TX, lineheight = .9)
f <- pb + pc + plot_layout(widths = c(1, 1)) + plot_annotation(tag_levels = "a")
saveRDS(list(pb = pb, pc = pc), file.path(FD, "panels_figS1.rds"))
save_fig(f, "FigureS1_kidneyAgeBench", 140, 80)
