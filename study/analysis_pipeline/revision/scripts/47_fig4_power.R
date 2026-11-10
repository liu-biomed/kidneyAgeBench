# Figure 4. What a non-significant result means for small gene sets (post hoc). The gene-level test has little power at
# the sizes of disease-model sets, the number of significant class C sets matches what a common modest enrichment would
# produce, and as a pool their genes are as age-associated as genes from kidney age signatures, unlike GenAge/CellAge.
source(file.path(path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030"), "scripts", "40_fig_common.R"))
P1 <- PW[PW$config == "primary_locked" & PW$true_ER > 1, ]; P1$er <- factor(sprintf("%.1f", P1$true_ER))
er_col <- setNames(c("#BDBDBD", "#969696", "#666666", "#383838", "#000000"), levels(P1$er))
Cn <- F[F$class == "C", ]; An <- F[F$class == "A", ]
## (a) power by set size
pa <- ggplot(P1, aes(k, power, colour = er)) + geom_hline(yintercept = c(.05, .8), linetype = c(3, 2), linewidth = .3, colour = "grey50") +
  geom_line(linewidth = .55) + geom_point(size = .9) +
  geom_rug(data = Cn, aes(x = n_measured), inherit.aes = FALSE, sides = "b", colour = cls_col[["C"]], linewidth = .35, length = unit(2, "mm")) +
  scale_colour_manual(values = er_col, name = "True\nenrichment") + scale_x_log10(breaks = c(2, 5, 10, 20, 50, 100)) +
  scale_y_continuous(labels = scales::percent, breaks = c(0, .05, .25, .5, .8, 1)) + coord_cartesian(clip = "off") +
  labs(x = "Genes measured in the set (log scale); ticks = class C sets", y = "Power (empirical P < 0.05)") +
  theme(legend.position = "right")
## (b) observed enrichment of class C sets against the smallest enrichment detectable with 80% power
mde <- do.call(rbind, lapply(sort(unique(P1$k)), function(k) { s <- P1[P1$k == k, ]; s <- s[order(s$true_ER), ]
  data.frame(k = k, mde = if (max(s$power) < .8) NA else if (s$power[1] >= .8) s$true_ER[1] else approx(s$power, s$true_ER, .8, ties = "ordered")$y) }))
mde <- mde[is.finite(mde$mde), ]
pb <- ggplot(Cn, aes(n_measured, enrichment)) + geom_ribbon(data = mde, aes(x = k, ymin = 1, ymax = mde), inherit.aes = FALSE, fill = "grey92") +
  geom_line(data = mde, aes(k, mde), inherit.aes = FALSE, linewidth = .45, linetype = 2) + geom_hline(yintercept = 1, linewidth = .3, colour = "grey40") +
  geom_point(aes(fill = status), shape = 21, size = 2, stroke = .35, position = position_jitter(width = .03, height = 0, seed = 4)) +
  annotate("text", x = 14, y = 1.75, label = "Detectable with\n80% power", size = TX, hjust = 0, lineheight = .9) +
  annotate("text", x = 30, y = 1.09, label = "Power < 80%", size = TX, hjust = .5, colour = "grey30") +
  scale_fill_manual(values = st_col, breaks = names(st_col), name = NULL) + scale_x_log10(breaks = c(2, 5, 10, 20, 50, 100), limits = c(1.8, 110)) +
  scale_y_continuous(trans = "log2", breaks = c(.25, .5, 1, 1.5, 2, 3), limits = c(.13, 3.3)) +
  labs(x = "Genes measured (log scale)", y = "Enrichment ratio, class C sets") + theme(legend.position = "inside", legend.position.inside = c(.98, .02), legend.justification = c(1, 0), legend.background = element_blank())
## (c) gene pools against expression-matched genes, with and without CPTAC
PL <- read.csv(file.path(RR, "REVISION_gene_pools_four_configs.csv")); PL <- PL[PL$config %in% c("primary_locked", "noCPTAC_evalbins"), ]
PL$cf <- factor(c(primary_locked = "All eligible datasets", noCPTAC_evalbins = "Without CPTAC")[PL$config], levels = c("All eligible datasets", "Without CPTAC"))
PL$pool <- factor(PL$pool, levels = rev(c("Class A, external sets", "GenAge", "CellAge", "Class C"))); pcol <- c("Class A, external sets" = cls_col[["A"]], GenAge = "#A9D3CE", CellAge = "#2E6E69", "Class C" = cls_col[["C"]])
plab <- setNames(sprintf("%s\n(%d genes)", levels(PL$pool), PL$n_genes[PL$config == "primary_locked"][match(levels(PL$pool), PL$pool[PL$config == "primary_locked"])]), levels(PL$pool))
pc <- ggplot(PL, aes(matched_ratio, pool, shape = cf)) + geom_vline(xintercept = 1, linetype = 2, linewidth = .3, colour = "grey40") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, linewidth = .45, position = position_dodge(width = .55), orientation = "y") +
  geom_point(aes(fill = pool), size = 2.1, stroke = .35, position = position_dodge(width = .55)) +
  scale_shape_manual(values = c(21, 24), name = NULL) + scale_fill_manual(values = pcol, guide = "none") + guides(shape = guide_legend(override.aes = list(fill = "grey60"))) +
  scale_y_discrete(labels = plab) + scale_x_continuous(breaks = c(.9, 1, 1.1, 1.2, 1.3)) +
  labs(x = "Mean |age z| relative to matched genes (95% null interval)", y = NULL) + theme(legend.position = "top", legend.margin = margin(0, 0, 0, 0))
## (d) observed number of class C sets with P < 0.05 against the number expected if every set had the same true enrichment
RV4 <- read.csv(file.path(RR, "REVISION_gene_level_four_configs.csv"))
ex <- do.call(rbind, lapply(c("primary_locked", "noCPTAC_evalbins"), function(cf) { d <- RV4[RV4$config == cf & RV4$class == "C", ]
  do.call(rbind, lapply(c(1, 1.1, 1.2, 1.3, 1.5), function(er) data.frame(config = cf, er = er, expected = sum(sapply(d$n, pw_at, er = er, cf = cf)), observed = sum(d$emp_p < .05)))) }))
ex$cf <- factor(c(primary_locked = "All eligible datasets", noCPTAC_evalbins = "Without CPTAC")[ex$config], levels = c("All eligible datasets", "Without CPTAC"))
ob <- unique(ex[, c("cf", "observed")])
pd <- ggplot(ex, aes(er, expected, shape = cf, linetype = cf)) + geom_line(linewidth = .45) + geom_point(size = 2, fill = "white", stroke = .4) +
  geom_hline(data = ob, aes(yintercept = observed, linetype = cf), colour = cls_col[["C"]], linewidth = .55) +
  geom_text(data = ob, aes(x = 1.0, y = observed, label = sprintf("observed: %d", observed)), vjust = -.5, hjust = 0, size = TX, colour = cls_col[["C"]], inherit.aes = FALSE) +
  scale_shape_manual(values = c(21, 24), name = NULL) + scale_linetype_manual(values = c(1, 2), name = NULL) +
  scale_x_continuous(breaks = c(1, 1.1, 1.2, 1.3, 1.5)) + scale_y_continuous(limits = c(0, 17), breaks = seq(0, 15, 5)) +
  labs(x = "Assumed true enrichment, same for all 29 sets", y = "Class C sets with P < 0.05") + theme(legend.position = "top", legend.margin = margin(0, 0, 0, 0))
f <- (pa + pb + plot_layout(widths = c(1, 1))) / (pc + pd + plot_layout(widths = c(1.15, 1))) + plot_layout(heights = c(1, 1)) + plot_annotation(tag_levels = "a")
saveRDS(list(pa = pa, pb = pb, pc = pc, pd = pd), file.path(FD, "panels_fig4.rds"))
save_fig(f, "Figure4_power", 183, 140)
write.csv(P1, file.path(FD, "SourceData_Figure4a.csv"), row.names = FALSE); write.csv(mde, file.path(FD, "SourceData_Figure4b_mde.csv"), row.names = FALSE)
write.csv(PL, file.path(FD, "SourceData_Figure4c.csv"), row.names = FALSE); write.csv(ex, file.path(FD, "SourceData_Figure4d.csv"), row.names = FALSE)
print(mde); print(ex)
