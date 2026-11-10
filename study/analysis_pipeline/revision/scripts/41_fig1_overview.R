# Figure 1. Study design: a leakage-free benchmark of 48 published gene sets in normal human kidney transcriptomes.
source(file.path(path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030"), "scripts", "40_fig_common.R"))
## (a) schematic (placeholder until the BioRender version replaces it)
bx <- function(x, y, w, h, lab, fill, col = "black", sz = TX, ff = "plain", tc = "black") list(
  annotate("rect", xmin = x - w / 2, xmax = x + w / 2, ymin = y - h / 2, ymax = y + h / 2, fill = fill, colour = col, linewidth = .25),
  annotate("text", x = x, y = y, label = lab, size = sz, lineheight = .9, fontface = ff, colour = tc))
hd <- function(x, lab) annotate("text", x = x, y = 10.55, label = lab, size = fs(7), fontface = "bold")
ar <- function(x0, x1, y = 5.2) annotate("segment", x = x0, xend = x1, y = y, yend = y, arrow = arrow(length = unit(1.6, "mm"), type = "closed"), linewidth = .45, colour = "grey35")
lt <- function(c) grDevices::adjustcolor(c, alpha.f = .16)
pa <- ggplot() + coord_cartesian(xlim = c(0, 20), ylim = c(0, 11), expand = FALSE) + theme_void() + theme(plot.tag = element_text(size = 9, face = "bold")) +
  hd(2.6, "Gene-set library (n = 48)") + hd(7.7, "Leakage-free normal kidney data") + hd(12.75, "Evaluation metrics") + hd(17.6, "Outputs") +
  bx(2.6, 8.6, 5.0, 2.2, "A  Kidney age-derived (12)\nRodwell, Melk, Rowland, Kainz,\nRAAG, MHC-II, KAS-94/82", lt(cls_col[["A"]]), cls_col[["A"]]) +
  bx(2.6, 5.6, 5.0, 2.2, "B  General ageing (7)\nde Magalhães, Palmer, Sood,\nGenAge, CellAge, SenMayo, Fridman", lt(cls_col[["B"]]), cls_col[["B"]]) +
  bx(2.6, 2.4, 5.0, 2.7, "C  'Ageing-related' disease\nmodels (29): DKD, IgAN, lupus,\ntransplant, RCC; 268 records\nscreened → 29 full gene lists", lt(cls_col[["C"]]), cls_col[["C"]]) +
  bx(7.7, 8.9, 4.7, 1.6, "6 GEO donor cohorts: 441 donors\n(9 series, 498 samples reviewed)", lt(src_col[["G6"]]), src_col[["G6"]]) +
  bx(7.7, 6.95, 4.7, 1.2, "TCGA solid tissue normal: 129", lt(src_col[["TCGA"]]), src_col[["TCGA"]]) +
  bx(7.7, 5.45, 4.7, 1.2, "CPTAC tumour-adjacent: 75", lt(src_col[["CPTAC"]]), src_col[["CPTAC"]]) +
  bx(7.7, 3.95, 4.7, 1.2, "GTEx post-mortem cortex: 85\n(reported separately)", lt(src_col[["GTEx"]]), src_col[["GTEx"]]) +
  bx(7.7, 1.7, 4.7, 2.1, "Leakage rule: a dataset used to\nderive a set is never used to\nevaluate it (KAS: no GEO;\nRowland: no TCGA)", "white", "grey40", ff = "italic") +
  bx(12.75, 8.75, 4.8, 1.9, "Gene-level enrichment vs\n2,000 expression-matched\nrandom gene sets", "grey96", "grey40") +
  bx(12.75, 6.45, 4.8, 1.9, "Signed score vs age;\nrandom-effects pooling;\n500 random signatures", "grey96", "grey40") +
  bx(12.75, 4.25, 4.8, 1.8, "Composition adjustment\n(immune, stroma, endothelium,\ndistal and proximal tubule)", "grey96", "grey40") +
  bx(12.75, 2.1, 4.8, 1.8, "Sensitivity without CPTAC;\npower by set size;\n47 GTEx tissues", "grey96", "grey40") +
  bx(17.6, 8.3, 4.1, 2.4, "Benchmark scorecard\nand power (Figs 3, 4)", "white", "black", sz = TXL) +
  bx(17.6, 5.3, 4.1, 2.4, "What tracks kidney\nage; composition;\ntissue source (Figs 5, 6)", "white", "black", sz = TXL) +
  bx(17.6, 2.3, 4.1, 2.4, "kidneyAgeBench\nR package (Fig. S1)", "grey25", "grey25", sz = TXL, ff = "bold", tc = "white") +
  ar(5.1, 5.35) + ar(10.05, 10.35) + ar(15.15, 15.5)
## (b) donor ages per dataset; n and % women sit outside the panel so they never overlap points
A <- FDv$ages; dsl <- c(G6_GSE36895 = "GSE36895", G6_GSE60807 = "GSE60807", G6_GSE147451 = "GSE147451", G6_GSE66272 = "GSE66272", G6_GSE93069 = "GSE93069", G6_GSE167093 = "GSE167093", TCGA = "TCGA", CPTAC = "CPTAC", GTEx = "GTEx")
A$ds <- factor(dsl[A$dataset], levels = rev(dsl)); an <- aggregate(age ~ ds + grp, A, length); fw <- aggregate(sex ~ ds, A, function(s) mean(s == "F", na.rm = TRUE)); an$f <- fw$sex[match(an$ds, fw$ds)]
pb <- ggplot(A, aes(age, ds, colour = grp)) + geom_violin(aes(fill = grp), colour = NA, alpha = .22, scale = "width", width = .8) +
  geom_quasirandom(size = .35, alpha = .7, width = .3, orientation = "y") + stat_summary(fun = median, geom = "point", shape = 124, size = 3.2, colour = "black") +
  geom_text(data = an, aes(x = Inf, y = ds, label = ifelse(is.na(f), sprintf("%d (n/a)", age), sprintf("%d (%.0f%%)", age, 100 * f))), hjust = -.15, size = TX, colour = "black") +
  annotate("text", x = Inf, y = length(dsl) + .75, label = "n (% women)", hjust = -.1, size = TX, fontface = "bold") +
  scale_colour_manual(values = src_col, guide = "none") + scale_fill_manual(values = src_col, guide = "none") + scale_x_continuous(breaks = seq(20, 80, 20), limits = c(15, 92)) +
  coord_cartesian(clip = "off") + labs(x = "Age (years)", y = NULL) + theme(plot.margin = margin(5.5, 52, 5.5, 5.5))
## (c) eligibility review (Table 1); retained donors take the GEO-cohort colour
E <- data.frame(s = c("GSE36895","GSE60807","GSE147451","GSE66272","GSE93069","GSE167093","GSE232249","GSE21785","GSE217427"),
  prev = c(22, 37, 177, 27, 36, 153, 12, 12, 22), kept = c(22, 37, 166, 27, 36, 153, 0, 0, 0),
  why = c("", "", "11 bilateral donors averaged", "", "", "", "m6A-IP / supernatant fractions", "Tubule/glomerulus microdissections", "Bilateral samples; 8/16 donors with AKI"))
E$s <- factor(E$s, levels = rev(E$s)); E2 <- rbind(data.frame(s = E$s, part = "Retained (one profile per donor)", n = E$kept), data.frame(s = E$s, part = "Removed", n = E$prev - E$kept))
pc <- ggplot(E2, aes(n, s, fill = part)) + geom_col(width = .7, colour = "grey45", linewidth = .2) + geom_text(data = E, aes(x = prev + 3, y = s, label = why), inherit.aes = FALSE, hjust = 0, size = TX, colour = "grey25") +
  scale_fill_manual(values = c("Retained (one profile per donor)" = src_col[["G6"]], Removed = "grey90"), name = NULL) + scale_x_continuous(limits = c(0, 330), breaks = c(0, 50, 100, 150), expand = c(0, 0)) +
  labs(x = "Samples", y = NULL) + theme(legend.position = "top", legend.justification = "left", legend.margin = margin(0, 0, 0, 0))
f <- pa / (pb + pc + plot_layout(widths = c(1, 1.05))) + plot_layout(heights = c(1.05, 1)) + plot_annotation(tag_levels = "a")
saveRDS(list(pa = pa, pb = pb, pc = pc), file.path(FD, "panels_fig1.rds"))
save_fig(f, "Figure1_overview", 183, 140)
