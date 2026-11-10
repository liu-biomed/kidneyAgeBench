# Figure 7. Tissue source shapes apparent ageing signals: a Y-chromosome expression loss with age only in tumour-adjacent
# tissue (post hoc), signed-score associations by source, and gene-level enrichment with and without CPTAC.
source(file.path(path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030"), "scripts", "40_fig_common.R"))
## (a) Y-chromosome gene score in men; shared y axis, points beyond the window drawn as edge triangles
Y <- FDv$ysc; Y$dataset <- factor(Y$dataset, levels = c("TCGA", "GTEx", "CPTAC")); ys <- FDv$yst; ys$dataset <- factor(ys$dataset, levels = levels(Y$dataset))
yl <- c(-2.6, 2.6); Yo <- Y[Y$yscore < yl[1] | Y$yscore > yl[2], ]; Yo$yc <- pmin(pmax(Yo$yscore, yl[1]), yl[2])
lab <- c(TCGA = "TCGA solid tissue normal", GTEx = "GTEx cortex (post-mortem)", CPTAC = "CPTAC tumour-adjacent")
pa <- ggplot(Y, aes(age, yscore, colour = dataset)) + geom_point(size = .8, alpha = .75) + geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = .45, fill = "grey85") +
  geom_point(data = Yo, aes(age, yc), shape = 25, fill = "white", size = 1.4, stroke = .4) +
  geom_text(data = Yo, aes(age, yc, label = sprintf("%.1f", yscore)), nudge_x = 3, hjust = 0, vjust = 0, size = TX, colour = "grey30") +
  facet_wrap(~ dataset, nrow = 1, labeller = as_labeller(lab)) + scale_colour_manual(values = src_col, guide = "none") +
  geom_text(data = ys, aes(x = -Inf, y = Inf, label = sprintf("ρ = %.2f, n = %d", rho, n)), hjust = -.08, vjust = 1.4, size = TX, colour = "black", inherit.aes = FALSE) +
  coord_cartesian(ylim = yl) + labs(x = "Age (years)", y = "Y-chromosome gene score (men)")
## (b) signed-score rho with age by tissue source
SM <- read.csv(file.path(RR, "BENCHMARK_summary.csv")); SM <- SM[SM$sig_id %in% F$sig_id[is.finite(F$pooled_rho)], ]
E <- rbind(data.frame(src = "G6", rho = SM$G6_rho), data.frame(src = "TCGA", rho = SM$TCGA_rho), data.frame(src = "CPTAC", rho = SM$CPTAC_rho), data.frame(src = "GTEx", rho = SM$GTEx_rho)); E <- E[is.finite(E$rho), ]
E$src <- factor(E$src, levels = c("G6", "TCGA", "CPTAC", "GTEx")); em <- aggregate(rho ~ src, E, median)
pb <- ggplot(E, aes(src, rho, colour = src)) + geom_hline(yintercept = 0, linewidth = .3, colour = "grey40") +
  geom_quasirandom(size = 1.1, width = .25, alpha = .85) + geom_crossbar(data = em, aes(ymin = rho, ymax = rho), width = .55, linewidth = .3, colour = "black") +
  geom_text(data = em, aes(y = .62, label = sprintf("%.2f", rho)), colour = "black", size = TX) +
  scale_colour_manual(values = src_col, guide = "none") + scale_x_discrete(labels = c(G6 = "GEO", TCGA = "TCGA", CPTAC = "CPTAC", GTEx = "GTEx")) +
  scale_y_continuous(limits = c(-.15, .66)) + labs(x = NULL, y = "Signed score ρ with age")
## (c) gene-level enrichment in all eligible datasets against the analysis without CPTAC (null re-matched)
X <- F[is.finite(F$enrichment_cop), ]; X <- X[order(X$class, decreasing = TRUE), ]
hl <- X$sig_id %in% c("B03_GenAge", "B04_CellAge", "B05_SenMayo", "B01_deMagalhaes2009", "A01_Rodwell2004", "A09_Perco2019_RAAG_CKD")
rs <- cor(X$enrichment, X$enrichment_cop, method = "spearman"); up <- X[hl & X$enrichment >= 1.05, ]; lo <- X[hl & X$enrichment < 1.05, ]
pc <- ggplot(X, aes(enrichment, enrichment_cop)) + annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, fill = "grey93") +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = .3, colour = "grey40") + geom_hline(yintercept = 1, linewidth = .25, colour = "grey60") + geom_vline(xintercept = 1, linewidth = .25, colour = "grey60") +
  geom_point(aes(fill = class, shape = fdr_cop < .05), size = 1.9, stroke = .3) +
  geom_text_repel(data = up, aes(label = lab), size = TX, segment.size = .2, segment.colour = "grey40", min.segment.length = 0, xlim = log2(c(.15, .8)), ylim = log2(c(1.3, 2.5)), direction = "y", hjust = 0, seed = 3) +
  geom_text_repel(data = lo, aes(label = lab), size = TX, segment.size = .2, segment.colour = "grey40", min.segment.length = 0, xlim = log2(c(1.3, 2.5)), ylim = log2(c(.3, .75)), direction = "y", hjust = 0, seed = 3) +
  scale_fill_manual(values = cls_col, labels = cls_name, name = NULL) + scale_shape_manual(values = c(`TRUE` = 24, `FALSE` = 21), breaks = c("TRUE", "FALSE"), labels = c(`TRUE` = "FDR < 0.05 without CPTAC", `FALSE` = "n.s. without CPTAC"), name = NULL) +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 2), ncol = 1), shape = guide_legend(ncol = 1)) +
  scale_x_continuous(trans = "log2", breaks = c(.25, .5, 1, 2), limits = c(.13, 2.6)) + scale_y_continuous(trans = "log2", breaks = c(.25, .5, 1, 2), limits = c(.13, 2.6)) +
  annotate("text", x = 2.5, y = .14, label = sprintf("Spearman ρ = %.2f, %d sets\nShaded: below random in both", rs, nrow(X)), hjust = 1, vjust = 0, size = TX, lineheight = .9) +
  labs(x = "Enrichment, all eligible datasets (G6 + TCGA + CPTAC)", y = "Enrichment, without CPTAC (G6 + TCGA)") +
  theme(legend.position = "right", legend.box = "vertical", legend.margin = margin(0, 0, 0, 0)) + coord_equal()
f <- pa / (pb + pc + plot_layout(widths = c(.8, 1.6))) + plot_layout(heights = c(.85, 1.1)) + plot_annotation(tag_levels = "a")
saveRDS(list(pa = pa, pb = pb, pc = pc), file.path(FD, "panels_fig7.rds"))
save_fig(f, "Figure7_tissue_source", 183, 150)
write.csv(Y, file.path(FD, "SourceData_Figure7a.csv"), row.names = FALSE); write.csv(E, file.path(FD, "SourceData_Figure7b.csv"), row.names = FALSE)
write.csv(X[, c("sig_id", "class", "enrichment", "fdr_all", "enrichment_cop", "fdr_cop")], file.path(FD, "SourceData_Figure7c.csv"), row.names = FALSE)
cat("Spearman all vs without CPTAC:", round(rs, 2), "n =", nrow(X), "\n"); print(em)
