# Supplementary Figure S2. Outside the kidney (47 GTEx tissues), the CellAge inducer/inhibitor annotation tracks age in
# most tissues, as do expression-derived ageing signatures; in kidney (GTEx cortex and the kidney evaluation data) the
# CellAge direction does not, whereas the expression-derived signatures do. Unsigned tests are insensitive in GTEx.
source(file.path(path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030"), "scripts", "40_fig_common.R"))
R <- read.csv(file.path(RR, "REVISION_gtex_multitissue_sets.csv")); R <- R[R$adjustment == "full", ]; K <- read.csv(file.path(RR, "REVISION_kidney_signed_dbtests.csv"))
lab <- c("CellAge inducers (+) / inhibitors (-)" = "CellAge\ndirection", "de Magalhães 2009" = "de Magalhães\n2009", "Palmer 2021" = "Palmer\n2021", SenMayo = "SenMayo", "Fridman 2008" = "Fridman\n2008", GenAge = "GenAge", CellAge = "CellAge")
K$analysis <- c("Locked primary datasets" = "Kidney, all eligible datasets", "Co-primary datasets" = "Kidney, without CPTAC")[K$analysis]
kcol <- c("Kidney, all eligible datasets" = src_col[["G6"]], "Kidney, without CPTAC" = src_col[["TCGA"]])
## (a) signed test: excess of the claimed direction over matched random genes
sg <- names(lab)[1:5]; A <- R[R$set %in% sg, ]; A$s <- factor(lab[A$set], levels = lab[sg]); A$kid <- A$tissue == "kidney_cortex"
KA <- K[K$set %in% sg, ]; KA$s <- factor(lab[KA$set], levels = lab[sg])
pa <- ggplot(A, aes(s, D_excess_z)) + geom_hline(yintercept = c(0, 1.96), linetype = c(1, 2), linewidth = .3, colour = c("grey40", "grey55")) +
  geom_quasirandom(data = A[!A$kid, ], size = .9, width = .28, colour = "grey60", alpha = .85) +
  geom_point(data = A[A$kid, ], shape = 21, size = 2.3, fill = src_col[["GTEx"]], colour = "black", stroke = .5) +
  geom_point(data = KA, aes(fill = analysis), shape = 23, size = 2.4, stroke = .35, position = position_nudge(x = .32)) +
  scale_fill_manual(values = kcol, name = NULL) + labs(x = NULL, y = "Claimed direction vs matched genes (z)") +
  theme(legend.position = "top", legend.margin = margin(0, 0, 0, 0))
## (b) unsigned test across tissues
us <- c("GenAge", "CellAge", "de Magalhães 2009", "Palmer 2021"); B <- R[R$set %in% us, ]; B$s <- factor(lab[B$set], levels = lab[us]); B$kid <- B$tissue == "kidney_cortex"
KB <- K[K$set %in% us, ]; KB$s <- factor(lab[KB$set], levels = lab[us])
pb <- ggplot(B, aes(s, ER)) + geom_hline(yintercept = 1, linetype = 2, linewidth = .3, colour = "grey40") +
  geom_quasirandom(data = B[!B$kid, ], size = .9, width = .28, colour = "grey60", alpha = .85) +
  geom_point(data = B[B$kid, ], shape = 21, size = 2.3, fill = src_col[["GTEx"]], colour = "black", stroke = .5) +
  geom_point(data = KB, aes(fill = analysis), shape = 23, size = 2.4, stroke = .35, position = position_nudge(x = .32), show.legend = FALSE) +
  scale_fill_manual(values = kcol) + scale_y_continuous(trans = "log2", breaks = c(.8, .9, 1, 1.1, 1.25)) +
  labs(x = NULL, y = "Unsigned enrichment ratio")
## (c) CellAge direction by tissue
C <- R[R$set == "CellAge inducers (+) / inhibitors (-)", ]; C$tissue <- gsub("_", " ", C$tissue); C$tissue <- factor(C$tissue, levels = C$tissue[order(C$D_excess_z)])
C$sig <- ifelse(C$fdr_signed < .05, "FDR < 0.05", "n.s.")
pc <- ggplot(C, aes(D_excess_z, tissue, fill = sig)) + geom_col(width = .75, colour = "grey30", linewidth = .15) + geom_vline(xintercept = 0, linewidth = .3) +
  scale_fill_manual(values = c("FDR < 0.05" = "#2B2B2B", "n.s." = "white"), name = NULL) +
  labs(x = "CellAge direction vs matched genes (z)", y = NULL) + theme(axis.text.y = element_text(size = 6, colour = ifelse(levels(C$tissue) == "kidney cortex", warn_col, "black")), legend.position = "inside", legend.position.inside = c(.98, .02), legend.justification = c(1, 0))
f <- (pa / pb + plot_layout(heights = c(1.15, 1))) | pc
f <- f + plot_layout(widths = c(1.25, 1)) + plot_annotation(tag_levels = "a")
saveRDS(list(pa = pa, pb = pb, pc = pc), file.path(FD, "panels_figS2.rds"))
save_fig(f, "FigureS2_gtex_multitissue", 183, 150)
write.csv(R, file.path(FD, "SourceData_FigureS2.csv"), row.names = FALSE)
