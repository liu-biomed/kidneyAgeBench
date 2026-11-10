# Figure 3: benchmark scorecard of 48 gene sets in all eligible datasets and without CPTAC, with the power of the
# gene-level test at each set's size; most disease-model sets are too small for a non-significant result to be informative.
source(file.path(path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030"), "scripts", "40_fig_common.R"))
X <- F[order(match(F$class, c("C","B","A")), F$enrichment), ]
X$y <- seq_len(nrow(X)) + 1.2 * (match(X$class, c("C","B","A")) - 1)
top <- max(X$y)
## column geometry (x units)
xb <- c(.45, 2.0); en <- c(2.75, 6.55); lr <- c(-3, 1.2); ex <- function(r) en[1] + (log2(r) - lr[1]) / diff(lr) * diff(en)
xpw <- 7.5; xt <- c(pnom = 8.8, pdir = 10.05, rho = 11.5, rhoc = 12.75, gtex = 14.0); xn <- 15.2
X$xn_bar <- xb[1] + log10(X$n_measured) / log10(1000) * diff(xb)
X$rnd <- ifelse(is.na(X$emp_p_null), NA, ifelse(X$emp_p_null < .05, "Beats random", "Not better"))
tl <- function(col, x) data.frame(y = X$y, x = x, v = X[[col]])
P1 <- rbind(tl("pct_nominal", xt["pnom"]), tl("pct_dir_correct", xt["pdir"])); R1 <- rbind(tl("pooled_rho", xt["rho"]), tl("pooled_rho_comp", xt["rhoc"]), tl("GTEx_rho", xt["gtex"]))
band <- aggregate(y ~ class, X, range); band <- data.frame(class = band$class, lo = band$y[, 1] - .5, hi = band$y[, 2] + .5)
hdr <- data.frame(x = c(mean(xb), mean(en), xpw, xt, xn), lab = c("Genes\nmeasured", "Gene-level age enrichment\nvs matched random genes", "Power,\nER 1.3\n(%)", "Nominal\nP < 0.05\n(%)", "Correct\nsign\n(%)", "Pooled\nρ", "Comp.-\nadj. ρ", "GTEx\nρ", "Random\nsigna-\ntures"))
grp <- data.frame(x = c(mean(xt[1:2]), mean(xt[3:5])), lab = c("Member genes", "Signed score vs age"), x0 = c(xt[1] - .5, xt[3] - .5), x1 = c(xt[2] + .5, xt[5] + .5))
ticks <- c(.25, .5, 1, 2)
pa <- ggplot() +
  geom_rect(data = band, aes(xmin = -.05, xmax = .2, ymin = lo, ymax = hi, fill = class), inherit.aes = FALSE) + scale_fill_manual(values = cls_col, guide = "none") +
  geom_text(data = band, aes(x = -.25, y = (lo + hi) / 2, label = c(A = "A", B = "B", C = "C")[class], colour = class), fontface = "bold", size = fs(8)) + scale_colour_manual(values = cls_col, guide = "none") +
  geom_rect(data = X[seq(1, nrow(X), 2), ], aes(xmin = .3, xmax = 15.75, ymin = y - .5, ymax = y + .5), fill = "grey97", inherit.aes = FALSE) +
  # genes bar
  geom_segment(data = X, aes(x = xb[1], xend = xn_bar, y = y, yend = y), linewidth = 1.3, colour = "grey60") +
  geom_text(data = X, aes(x = xn_bar + .07, y = y, label = n_measured), hjust = 0, size = TX, colour = "grey30") +
  # enrichment
  geom_segment(data = data.frame(x = ex(ticks)), aes(x = x, xend = x, y = .4, yend = top + .5), colour = ifelse(ticks == 1, "grey30", "grey88"), linetype = ifelse(ticks == 1, 2, 1), linewidth = .25) +
  geom_segment(data = X, aes(x = ex(enrichment), xend = ex(enrichment_cop), y = y, yend = y), colour = "grey60", linewidth = .3, na.rm = TRUE) +
  geom_point(data = X, aes(x = ex(enrichment_cop), y = y), shape = 124, colour = "grey35", size = 1.6, na.rm = TRUE) +
  new_scale_fill() + geom_point(data = X, aes(x = ex(enrichment), y = y, fill = status), shape = 21, colour = "black", size = 1.9, stroke = .3) +
  scale_fill_manual(values = st_col, breaks = names(st_col), name = "Enrichment", guide = guide_legend(order = 1, ncol = 1, override.aes = list(size = 2))) +
  annotate("text", x = ex(ticks), y = .0, label = ticks, size = TX, vjust = 1) + annotate("text", x = mean(en), y = -1.3, label = "Ratio (log scale); circle = all eligible datasets, tick = without CPTAC", size = TX, vjust = 1) +
  # power of the gene-level test at the set's size (greyscale, so it never reads as an effect)
  new_scale_fill() + geom_tile(data = data.frame(x = xpw, y = X$y, v = X$power13), aes(x, y, fill = v), width = 1.05, height = .9) +
  scale_fill_gradient(low = "white", high = "grey25", limits = c(0, 100), name = "Power (%)", guide = guide_colourbar(order = 2, barheight = unit(1.8, "mm"), barwidth = unit(16, "mm"))) +
  geom_text(data = data.frame(x = xpw, y = X$y, v = X$power13), aes(x, y, label = sprintf("%.0f", v)), colour = ifelse(X$power13 > 60, "white", "black"), size = TX) +
  # tiles
  new_scale_colour() + scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
  new_scale_fill() + geom_tile(data = P1, aes(x, y, fill = v), width = 1.05, height = .9) + geom_text(data = P1[!is.na(P1$v), ], aes(x, y, label = sprintf("%.0f", v), colour = v > 60), size = TX) +
  scale_fill_gradientn(colours = c("white", "#DEEBF7", "#6BAED6", "#08306B"), limits = c(0, 100), na.value = NA, name = "Genes (%)", guide = guide_colourbar(order = 4, barheight = unit(1.8, "mm"), barwidth = unit(16, "mm"))) +
  new_scale_fill() + geom_tile(data = R1, aes(x, y, fill = v), width = 1.05, height = .9) + geom_text(data = R1[!is.na(R1$v), ], aes(x, y, label = sub("^(-?)0\\.", "\\1.", sprintf("%.2f", v)), colour = abs(v) > .22), size = TX) +
  scale_fill_gradientn(colours = c("#B2182B", "#F4A582", "white", "#92C5DE", "#2166AC"), limits = c(-.35, .35), oob = scales::squish, na.value = NA, name = "Spearman ρ", guide = guide_colourbar(order = 4, barheight = unit(1.8, "mm"), barwidth = unit(16, "mm"))) +
  new_scale_colour() + geom_point(data = X[!is.na(X$rnd), ], aes(xn, y, shape = rnd, colour = rnd), size = 1.6, stroke = .6) +
  scale_shape_manual(values = c(16, 4), name = "Signed score vs random", guide = guide_legend(order = 5, ncol = 1)) + scale_colour_manual(values = c("#2B2B2B", warn_col), name = "Signed score vs random", guide = guide_legend(order = 5, ncol = 1)) +
  # headers
  geom_text(data = hdr, aes(x, top + 1.1, label = lab), vjust = 0, size = TXL, lineheight = .85) +
  geom_segment(data = grp, aes(x = x0, xend = x1, y = top + 5.6, yend = top + 5.6), linewidth = .3) + geom_text(data = grp, aes(x, top + 5.9, label = lab), vjust = 0, size = fs(7), fontface = "bold") +
  scale_y_continuous(breaks = X$y, labels = X$lab, expand = expansion(c(.005, .005))) + scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(-.45, 15.8), ylim = c(-2.6, top + 8.4), clip = "off") + labs(x = NULL, y = NULL) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(size = 6, margin = margin(r = 0)),
        legend.position = "bottom", legend.box = "horizontal", legend.title.position = "top", legend.justification = "right", legend.spacing.x = unit(1, "mm"), legend.margin = margin(0, 4, 0, 4))
# class summary in both analyses
lv <- c("A  Kidney age-derived", "B  General ageing", "C  Disease models")
sm <- rbind(data.frame(cl = X$class, an = "All eligible datasets (G6 + TCGA + CPTAC)", st = X$status), data.frame(cl = X$class, an = "Without CPTAC (G6 + TCGA)", st = X$status_cop))
sm$an <- factor(sm$an, levels = unique(sm$an))
sm <- sm[!is.na(sm$st), ]; sm$st <- factor(sm$st, levels = rev(names(st_col))); sm$cl <- factor(setNames(lv, c("A","B","C"))[sm$cl], levels = rev(lv))
lab <- aggregate(st ~ cl + an, sm, function(z) sprintf("%d/%d", sum(z == "FDR < 0.05"), length(z)))
pb <- ggplot(sm, aes(y = cl, fill = st)) + geom_bar(position = "fill", width = .72, colour = "grey30", linewidth = .2) + facet_grid(. ~ an) +
  geom_text(data = lab, aes(x = 1.03, y = cl, label = st), inherit.aes = FALSE, hjust = 0, size = TX) +
  scale_fill_manual(values = st_col, guide = "none") + scale_x_continuous(labels = scales::percent, breaks = c(0, .5, 1), expand = expansion(c(0, .16))) +
  labs(x = "Share of gene sets (label: number with FDR < 0.05)", y = NULL) + theme(panel.spacing.x = unit(8, "mm"))
f <- pa / pb + plot_layout(heights = c(1, .14)) + plot_annotation(tag_levels = "a")
saveRDS(list(pa = pa, pb = pb), file.path(FD, "panels_fig3.rds"))
save_fig(f, "Figure3_scorecard", 183, 230)
write.csv(X[, c("sig_id","class","lab","n_measured","enrichment","fdr_all","enrichment_cop","fdr_cop","power13","pct_nominal","pct_dir_correct","pooled_rho","pooled_rho_comp","GTEx_rho","emp_p_null")], file.path(FD, "SourceData_Figure3.csv"), row.names = FALSE)
