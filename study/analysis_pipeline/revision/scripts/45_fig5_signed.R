# Figure 6. Claim: signed signatures beat size-matched random signatures but with modest effects, the signal is uneven
# across cohorts, and about half of it cannot be separated from age-related shifts in tissue composition.
source(file.path(path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030"), "scripts", "40_fig_common.R"))
G <- F[is.finite(F$pooled_rho), ]; G <- G[order(-match(G$class, c("A","B")), G$pooled_rho), ]; G$labf <- factor(G$lab, levels = G$lab)
G$cl <- factor(c(A = "Kidney age-derived", B = "General ageing")[G$class], levels = c("Kidney age-derived", "General ageing"))
fac <- facet_grid(cl ~ ., scales = "free_y", space = "free_y")
## (a) observed pooled rho against its own random-signature null
NL <- do.call(rbind, lapply(G$sig_id, function(s) data.frame(sig_id = s, rho = NUL[[s]]))); NL <- merge(NL, G[, c("sig_id","labf","cl")])
pa <- ggplot(G, aes(y = labf)) + geom_vline(xintercept = 0, linewidth = .3, colour = "grey40") +
  geom_violin(data = NL, aes(rho, labf), fill = "grey80", colour = NA, scale = "width", width = .85) +
  geom_errorbar(aes(xmin = pooled_lo, xmax = pooled_hi, colour = class), width = 0, orientation = "y", linewidth = .5) + scale_colour_manual(values = cls_col, guide = "none") +
  geom_point(aes(pooled_rho, fill = emp_p_null < .05), shape = 23, size = 2, stroke = .35) +
  geom_point(aes(rho_noCPTAC, shape = "Without CPTAC (G6 + TCGA)"), size = 1.5, stroke = .45, colour = "grey30") + scale_shape_manual(values = 1, name = NULL) + fac +
  scale_fill_manual(values = c(`TRUE` = "#2B2B2B", `FALSE` = "white"), breaks = c("TRUE", "FALSE"), labels = c(`TRUE` = "P < 0.05 vs random signatures", `FALSE` = "n.s."), name = NULL) +
  labs(x = "Pooled partial Spearman ρ with age", y = NULL) +
  theme(strip.text.y = element_text(angle = -90), legend.position = "bottom", legend.box = "vertical", legend.spacing.y = unit(0, "mm"), legend.margin = margin(0, 0, 0, 0))
## (b) every cohort, base model
dsl <- c(G6_GSE36895 = "GSE36895", G6_GSE60807 = "GSE60807", G6_GSE147451 = "GSE147451", G6_GSE66272 = "GSE66272", G6_GSE93069 = "GSE93069", G6_GSE167093 = "GSE167093", TCGA = "TCGA", CPTAC = "CPTAC", GTEx = "GTEx")
B <- MS[MS$model == "base" & MS$sig_id %in% G$sig_id, ]; B <- merge(B, G[, c("sig_id","labf","cl")])
B$excl <- (grepl("^A1[12]_", B$sig_id) & B$grp == "G6") | (grepl("^A0[56]_", B$sig_id) & B$grp == "TCGA")
full <- expand.grid(sig_id = G$sig_id, dataset = names(dsl), stringsAsFactors = FALSE); B <- merge(full, B, all.x = TRUE); B <- merge(B[, setdiff(names(B), c("labf","cl"))], G[, c("sig_id","labf","cl")])
B$ds <- factor(dsl[B$dataset], levels = dsl); B$v <- ifelse(B$excl %in% TRUE, NA, B$rho)
nds <- aggregate(n ~ dataset, MS[MS$model == "base", ], max); dlab <- setNames(sprintf("%s (%d)", dsl[nds$dataset], nds$n), dsl[nds$dataset])
pb <- ggplot(B, aes(ds, labf, fill = v)) + geom_tile(colour = "white", linewidth = .35) +
  geom_text(aes(label = ifelse(B$excl %in% TRUE, "×", ifelse(is.na(v), "", ifelse(p < .05, "•", "")))), size = ifelse(B$excl %in% TRUE, fs(7), fs(8)), colour = ifelse(B$excl %in% TRUE, "grey45", "black")) + fac +
  annotate("rect", xmin = .5, xmax = 6.5, ymin = -Inf, ymax = Inf, fill = NA, colour = src_col[["G6"]], linewidth = .3) +
  scale_fill_gradientn(colours = c("#B2182B", "#F4A582", "white", "#92C5DE", "#2166AC", "#053061"), values = scales::rescale(c(-.5, -.2, 0, .2, .4, .6)), limits = c(-.5, .6), oob = scales::squish, na.value = "grey93", name = "ρ") +
  scale_x_discrete(labels = dlab) + guides(fill = guide_colourbar(barwidth = unit(2.2, "mm"), barheight = unit(18, "mm"))) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank(), strip.text = element_blank(), axis.text.x = element_text(angle = 50, hjust = 1))
## (c) composition adjustment
P <- rbind(data.frame(labf = G$labf, cl = G$cl, model = "Unadjusted", rho = G$pooled_rho), data.frame(labf = G$labf, cl = G$cl, model = "Adjusted", rho = G$pooled_rho_comp))
P$model <- factor(P$model, levels = c("Unadjusted", "Adjusted")); att <- 1 - G$pooled_rho_comp / G$pooled_rho
hl <- G$labf[G$sig_id %in% c("A11_KAS94_original","A05_Rowland2019_37","B05_SenMayo","A07_Kainz2007_donor","A04_Melk2005_PAM")]
pc <- ggplot(P, aes(model, rho, group = labf)) + geom_hline(yintercept = 0, linewidth = .3, colour = "grey40") +
  geom_line(aes(colour = cl), linewidth = .35, alpha = .7) + geom_point(aes(colour = cl), size = 1.2) +
  geom_text_repel(data = P[P$model != "Unadjusted" & P$labf %in% hl, ], aes(label = c("KAS-94 (this study, original)" = "KAS-94", "Rowland 2019 (37)" = "Rowland 2019", "Melk 2005 (PAM)" = "Melk PAM", "SenMayo" = "SenMayo", "Kainz 2007 donor quality" = "Kainz 2007")[as.character(labf)]), size = TX, nudge_x = .3, direction = "y", hjust = 0, segment.size = .2, min.segment.length = 0, xlim = c(2.2, 4.2), box.padding = .15) +
  scale_colour_manual(values = c("Kidney age-derived" = cls_col[["A"]], "General ageing" = cls_col[["B"]]), name = NULL) +
  scale_x_discrete(expand = expansion(add = c(.3, 1.6))) +
  annotate("text", x = 2.3, y = .335, label = sprintf("Median attenuation %.0f%%\n(IQR %.0f–%.0f%%)", 100 * median(att), 100 * quantile(att, .25), 100 * quantile(att, .75)), size = TX, lineheight = .9) +
  labs(x = NULL, y = "Pooled ρ with age") + theme(legend.position = "bottom", legend.direction = "vertical", legend.margin = margin(0, 0, 0, 0))
## (d) composition proxies themselves change with age
cp <- FDv$comp; cp$ds <- factor(dsl[cp$dataset], levels = dsl)
cp$proxy <- factor(c(immune = "Immune", stroma = "Stroma", endothelial = "Endothelial", distal_tub = "Distal tubule", PT_content = "Proximal tubule")[cp$proxy], levels = c("Immune", "Stroma", "Endothelial", "Distal tubule", "Proximal tubule"))
pd <- ggplot(cp, aes(ds, proxy, fill = rho)) + geom_tile(colour = "white", linewidth = .35) + geom_text(aes(label = sub("^(-?)0\\.", "\\1.", sprintf("%.2f", rho)), fontface = ifelse(p < .05, "bold", "plain")), size = TX) +
  scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC", limits = c(-.45, .45), oob = scales::squish, name = "ρ") + guides(fill = guide_colourbar(barwidth = unit(2.2, "mm"), barheight = unit(14, "mm"))) +
  labs(x = NULL, y = NULL) + theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_text(angle = 40, hjust = 1))
# panel (e) of the earlier version (tissue source) moved to Figure 7b
f <- (pa + pb + pc + plot_layout(widths = c(1.25, 1.05, .95))) / (pd + plot_spacer() + plot_layout(widths = c(1.6, .25))) + plot_layout(heights = c(1, .42)) + plot_annotation(tag_levels = "a")
saveRDS(list(pa = pa, pb = pb, pc = pc, pd = pd), file.path(FD, "panels_fig6.rds"))
save_fig(f, "Figure6_signed_replication", 183, 150)
write.csv(B[, c("sig_id","dataset","rho","p","excl")], file.path(FD, "SourceData_Figure6b.csv"), row.names = FALSE); write.csv(cp, file.path(FD, "SourceData_Figure6d.csv"), row.names = FALSE)
