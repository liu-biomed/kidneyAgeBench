# Figure 5. Claim: 'ageing-related' disease-model genes look age-associated only against an unmatched genome background,
# because they are highly expressed and variable; against expression-matched genes the advantage largely disappears,
# GenAge/CellAge genes are not age-associated, and the most re-used 'ageing' genes do not track age. Descriptive (post hoc).
source(file.path(path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030"), "scripts", "40_fig_common.R"))
load(file.path(REV, "kidneyAgeBench/data/kab_bins.rda")); L <- FDv$lib; pz <- BI$PZ$none
ok <- intersect(names(pz)[is.finite(pz)], names(kab_bins)); pz <- pz[ok]; bins <- kab_bins[ok]
mud <- as.integer(sub(" .*", "", bins)); sdd <- as.integer(sub(".* ", "", bins)); names(mud) <- names(sdd) <- ok
ga <- unique(L$gene[L$sig_id == "B03_GenAge"]); ca <- unique(L$gene[L$sig_id == "B04_CellAge"])
Aext <- unique(L$gene[L$class == "A" & !grepl("^A1[12]_|^A0[56]_", L$sig_id)])   # class A sets with no derivation overlap with G6/TCGA/CPTAC
Cg <- unique(L$gene[L$class == "C"])
pools <- list("All measured genes" = ok, "Class A genes\n(external sets)" = Aext, "GenAge" = ga, "CellAge" = ca, "Class C genes" = Cg)
pools <- lapply(pools, intersect, ok); pc <- c("#7F7F7F", cls_col[["A"]], "#A9D3CE", "#2E6E69", cls_col[["C"]]); names(pc) <- names(pools)
## (a) |z| distributions against the genome background
D1 <- do.call(rbind, lapply(names(pools), function(k) data.frame(pool = k, az = abs(pz[pools[[k]]]))))
nl <- setNames(paste0(names(pools), " (n = ", sapply(pools, length), ")"), names(pools)); nl["Class A genes\n(external sets)"] <- sprintf("Class A genes, external\nsets (n = %d)", length(pools[[2]]))
D1$pool <- factor(D1$pool, levels = rev(names(pools)))
md <- aggregate(az ~ pool, D1, median); nn <- aggregate(az ~ pool, D1, length)
pa <- ggplot(D1, aes(az, pool, fill = pool)) + geom_violin(scale = "width", width = .8, linewidth = .2, colour = NA, alpha = .85) +
  geom_boxplot(width = .14, outlier.shape = NA, linewidth = .25, fill = "white") + geom_vline(xintercept = median(abs(pz)), linetype = 2, linewidth = .3, colour = "grey30") +
  scale_fill_manual(values = pc, guide = "none") + scale_y_discrete(labels = nl) + coord_cartesian(xlim = c(0, 7)) + labs(x = "|Pooled age z| per gene (G6 + TCGA + CPTAC)", y = NULL)
## (b) why: |z| depends on expression level and variability; class C genes sit in high-variance bins
G <- data.frame(g = ok, mu = mud, sd = sdd, az = abs(pz)); bz <- aggregate(az ~ mu + sd, G, mean); bn <- aggregate(az ~ mu + sd, G, length); bz <- bz[bn$az >= 20, ]
cc <- as.data.frame(table(mu = mud[pools[["Class C genes"]]], sd = sdd[pools[["Class C genes"]]])); cc$mu <- as.integer(as.character(cc$mu)); cc$sd <- as.integer(as.character(cc$sd)); cc <- cc[cc$Freq > 0, ]
pb <- ggplot(bz, aes(mu, sd)) + geom_tile(aes(fill = az), colour = "white", linewidth = .2) +
  geom_point(data = cc, aes(size = Freq), shape = 21, fill = cls_col[["C"]], colour = "black", stroke = .2, alpha = .9) +
  scale_fill_gradient(low = "white", high = "#4A4A4A", name = "Mean |z|\n(all genes)") + scale_size_area(max_size = 3.2, breaks = c(1, 3, 6), name = "Class C\ngenes") +
  scale_x_continuous(breaks = c(1, 5, 10), expand = c(0, 0)) + scale_y_continuous(breaks = c(1, 5, 10), expand = c(0, 0)) + coord_equal() +
  labs(x = "Expression mean decile", y = "Expression SD decile") +
  theme(legend.key.height = unit(3, "mm"), legend.key.width = unit(2.4, "mm"), axis.line = element_blank())
## (d) provenance of class C sets vs their benchmark result
Cs <- F[F$class == "C", ]; Cs$from <- sapply(Cs$sig_id, function(s) { g <- unique(L$gene[L$sig_id == s]); mean(g %in% union(ga, ca)) })
rs <- cor.test(Cs$from, Cs$enrichment, method = "spearman", exact = FALSE)
pd <- ggplot(Cs, aes(100 * from, enrichment)) + geom_hline(yintercept = 1, linetype = 2, linewidth = .3, colour = "grey40") +
  geom_point(aes(fill = status), shape = 21, size = 1.9, stroke = .35, position = position_jitter(width = 1.5, height = 0, seed = 2)) +
  scale_fill_manual(values = st_col, breaks = names(st_col), name = NULL) + scale_y_continuous(trans = "log2", breaks = c(.25, .5, 1, 2), limits = c(.13, 4)) +
  annotate("text", x = 0, y = .2, hjust = 0, size = TX, label = sprintf("%d/29 sets: ≥ 50%% of\ngenes from GenAge/CellAge\nSpearman ρ = %.2f\nP = %.2f", sum(Cs$from >= .5), rs$estimate, rs$p.value), lineheight = .9) +
  labs(x = "Genes also in GenAge or CellAge (%)", y = "Enrichment vs matched genes") + theme(legend.position = "top", legend.direction = "vertical", legend.margin = margin(0, 0, 0, 0))
## (e) most re-used 'ageing' genes across class C models
Cl <- L[L$class == "C", ]; Cl <- Cl[!duplicated(paste(Cl$sig_id, Cl$gene)), ]; tb <- sort(table(Cl$gene), decreasing = TRUE); tb <- tb[tb >= 2]
E <- data.frame(gene = names(tb), k = as.integer(tb), z = pz[names(tb)], src = ifelse(names(tb) %in% intersect(ga, ca), "GenAge + CellAge", ifelse(names(tb) %in% ga, "GenAge", ifelse(names(tb) %in% ca, "CellAge", "Neither"))))
E$gene <- factor(E$gene, levels = E$gene[order(E$z)])
pe <- ggplot(E, aes(z, gene)) + annotate("rect", xmin = -1.96, xmax = 1.96, ymin = -Inf, ymax = Inf, fill = "grey92") + geom_vline(xintercept = 0, linewidth = .3) +
  geom_segment(aes(x = 0, xend = z, yend = gene), linewidth = .35, colour = "grey55") + geom_point(aes(size = k, fill = src), shape = 21, stroke = .25) +
  scale_size_area(max_size = 3.5, breaks = c(2, 3, 6), name = "Models using\nthe gene") + scale_fill_manual(values = c("GenAge + CellAge" = cls_col[["B"]], GenAge = "#A9D3CE", CellAge = "#2E6E69", Neither = "white"), name = "Database") + guides(fill = guide_legend(override.aes = list(size = 2.2))) +
  labs(x = "Pooled age z (grey band: |z| < 1.96)", y = NULL) + theme(axis.text.y = element_text(face = "italic"))
# panel (c) of the earlier version (pool ratios) moved to Figure 4c with both analyses
f <- (pa + pb + plot_layout(widths = c(1.1, 1))) / (free(pd) + pe + plot_layout(widths = c(1, 1.1))) + plot_layout(heights = c(1, 1.05)) + plot_annotation(tag_levels = "a")
saveRDS(list(pa = pa, pb = pb, pd = pd, pe = pe), file.path(FD, "panels_fig5.rds"))
save_fig(f, "Figure5_why_labels", 183, 150)
write.csv(Cs[, c("sig_id","from","enrichment","status")], file.path(FD, "SourceData_Figure5c.csv"), row.names = FALSE); write.csv(E, file.path(FD, "SourceData_Figure5d.csv"), row.names = FALSE)
print(rs$estimate); print(rs$p.value); print(E)
