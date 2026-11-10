# Figures for manuscript v7 (Genome Medicine). Rebuilt from saved results and source data only; no cohort is reprocessed.
#  - Figure 1: overview with v7 figure numbers and set-class names
#  - Figure 5: gene-level composition adjustment of pooled gene sources + KPMP single-nucleus deconvolution concordance
#              (previously Figure 8a and 8b)
#  - Figure 9: clinical evaluation in CKD and transplant biopsies (previously Figure 8c-e)
#  - Figures 2-4, 6-8 and Supplementary Figures S1-S2 are copied from v3 under their v7 numbers.
# The deconvolution concordance labels are recomputed (Spearman, 13 signatures); the v6 labels did not match the saved results.

REV <- path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030")
source(file.path(REV, "scripts", "40_fig_common.R"))
FD3 <- FD
FD <- file.path(REV, "benchmark", "figures", "v7"); dir.create(FD, showWarnings = FALSE, recursive = TRUE)
RR <- file.path(REV, "benchmark", "results")

## ---- Figure 5: composition ----
pa_df <- read.csv(file.path(FD3, "SourceData_Figure8a.csv"))
an_new <- c("Normal kidney age (unadjusted)" = "Normal kidney age, unadjusted",
            "Normal kidney age (composition-adjusted)" = "Normal kidney age, composition-adjusted",
            "CKD vs reference biopsies (n = 277)" = "CKD vs reference biopsies (n = 277)")
pa_df$an_lab <- factor(an_new[pa_df$an_lab], levels = unname(an_new))
pool_new <- c("Class A (external age sets)" = "Kidney age-derived sets (A)",
              "RAAG, CKD-concordant" = "RAAG, CKD-concordant",
              "RAAG, non-CKD" = "RAAG, other",
              "Class C ('ageing' disease sets)" = "Disease-derived sets (C)",
              "CellAge" = "CellAge", "GenAge" = "GenAge")
pa_df$pool <- factor(pool_new[pa_df$pool], levels = rev(unname(pool_new)))
p5a <- ggplot(pa_df, aes(ER, pool, colour = an_lab, shape = an_lab)) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = .3, colour = "grey40") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, linewidth = .4, orientation = "y", position = position_dodge(width = .62)) +
  geom_point(size = 1.9, stroke = .45, fill = "white", position = position_dodge(width = .62)) +
  scale_colour_manual(values = setNames(c("#1D63B8", "#008B7D", "#C87028"), levels(pa_df$an_lab)), name = NULL) +
  scale_shape_manual(values = setNames(c(16, 21, 17), levels(pa_df$an_lab)), name = NULL) +
  labs(x = "Enrichment vs expression-matched genes (95% null interval)", y = NULL) +
  theme(legend.position = "top", legend.direction = "vertical", legend.margin = margin(0, 0, 0, 0))

dc <- read.csv(file.path(RR, "REVISION_deconv_method_concordance.csv"))
r_mu <- cor(dc$G6_marker5, dc$G6_KPMP_MuSiC, method = "spearman")
r_bi <- cor(dc$G6_marker5, dc$G6_KPMP_Bisque_5sub, method = "spearman")
lab_mu <- sprintf("MuSiC (Spearman ρ = %.2f)", r_mu)
lab_bi <- sprintf("Bisque, five subclasses (Spearman ρ = %.2f)", r_bi)
pe_df <- rbind(data.frame(sig_id = dc$sig_id, marker5 = dc$G6_marker5, deconv = dc$G6_KPMP_MuSiC, method = lab_mu),
               data.frame(sig_id = dc$sig_id, marker5 = dc$G6_marker5, deconv = dc$G6_KPMP_Bisque_5sub, method = lab_bi))
pe_df$method <- factor(pe_df$method, levels = c(lab_mu, lab_bi))
p5b <- ggplot(pe_df, aes(marker5, deconv, colour = method, shape = method)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = .3, colour = "grey50") +
  geom_hline(yintercept = 0, linetype = 3, linewidth = .25, colour = "grey65") +
  geom_vline(xintercept = 0, linetype = 3, linewidth = .25, colour = "grey65") +
  geom_point(size = 1.9, stroke = .45, fill = "white") +
  scale_colour_manual(values = setNames(c("#1D63B8", "#C87028"), c(lab_mu, lab_bi)), name = NULL) +
  scale_shape_manual(values = setNames(c(16, 21), c(lab_mu, lab_bi)), name = NULL) +
  labs(x = "Age ρ adjusted for five marker scores",
       y = "Age ρ adjusted for KPMP\nsingle-nucleus deconvolution") +
  theme(legend.position = "top", legend.direction = "vertical", legend.margin = margin(0, 0, 0, 0))

fig5 <- p5a + p5b + plot_layout(widths = c(1.3, 1)) + plot_annotation(tag_levels = "a")
save_fig(fig5, "Figure5_composition", 183, 88)
write.csv(pa_df, file.path(FD, "SourceData_Figure5a.csv"), row.names = FALSE)
write.csv(pe_df, file.path(FD, "SourceData_Figure5b.csv"), row.names = FALSE)

## ---- Figure 9: clinical ----
pb_long <- read.csv(file.path(FD3, "SourceData_Figure8b.csv"))
show_lab <- c("Rodwell 2004", "Melk 2005 (PAM)", "Rowland 2019 (rep)", "KAS-82", "de Magalhães 2009", "SenMayo")
pb_agg_inj <- unique(pb_long[, c("diag", "PT_injury")])
pb_long$sig_lab <- factor(pb_long$sig_lab, levels = show_lab)
end_df <- do.call(rbind, lapply(split(pb_long, pb_long$sig_lab), function(d) {
  mx <- max(pb_agg_inj$PT_injury); f <- lm(score_SD ~ PT_injury, data = d)
  data.frame(sig_lab = d$sig_lab[1], cls = d$cls[1], PT_injury = mx, score_SD = as.numeric(predict(f, newdata = data.frame(PT_injury = mx))))
}))
end_df <- end_df[order(end_df$score_SD), ]
end_df$score_end <- end_df$score_SD
for (i in 2:nrow(end_df)) if (end_df$score_SD[i] - end_df$score_SD[i - 1] < 0.36) end_df$score_SD[i] <- end_df$score_SD[i - 1] + 0.36
end_df$y_lab <- end_df$score_SD - c(0, 0, 0.14, 0.14, 0.14, 0.14)
diag_ticks <- data.frame(PT_injury = c(-0.28, 0.33, 0.68, 1.58, 1.89, 2.20, 2.82, 4.01),
                         lab = c("Living\ndonor", "TMD /\nMCD", "MN", "IgAN /\nHTN", "LN", "FSGS", "DKD", "ANCA"))
cls_new <- c("A  Kidney age-derived" = "Kidney age-derived (A)", "B  General ageing / senescence" = "General ageing resource (B)")
pb_long$cls <- factor(cls_new[pb_long$cls], levels = unname(cls_new)); end_df$cls <- factor(cls_new[end_df$cls], levels = unname(cls_new))
p9a <- ggplot(pb_long, aes(PT_injury, score_SD, group = sig_lab, colour = cls)) +
  geom_vline(data = diag_ticks, aes(xintercept = PT_injury), linetype = 3, linewidth = .22, colour = "grey78") +
  geom_hline(yintercept = 0, linetype = 2, linewidth = .25, colour = "grey55") +
  geom_smooth(method = "lm", se = FALSE, linewidth = .4, formula = y ~ x) +
  geom_point(size = 1.35, alpha = .85) +
  geom_text(data = diag_ticks, aes(x = PT_injury, y = 5.25, label = lab), inherit.aes = FALSE, size = TX * 0.85, vjust = 1, lineheight = 0.85, colour = "grey30") +
  geom_segment(data = end_df, aes(x = PT_injury + 0.03, xend = PT_injury + 0.10, y = score_end, yend = y_lab), inherit.aes = FALSE, linewidth = .22, colour = "grey60") +
  geom_text(data = end_df, aes(x = PT_injury + 0.12, y = y_lab, label = sig_lab), inherit.aes = FALSE, colour = "grey15", hjust = 0, size = TX * 0.92, show.legend = FALSE) +
  scale_colour_manual(values = setNames(c("#1D63B8", "#C87028"), unname(cls_new)), name = NULL) +
  scale_x_continuous(limits = c(-0.45, 5.05), breaks = 0:4) +
  scale_y_continuous(limits = c(-0.8, 5.40), breaks = 0:4) +
  labs(x = "Median proximal tubular injury score by diagnosis (SD from living donors; ERCB, n = 190)",
       y = "Median signature score\n(SD from living donors)") +
  theme(legend.position = "top", legend.direction = "horizontal", legend.margin = margin(0, 0, 0, 0))

clin <- read.csv(file.path(RR, "REVISION_clinical_CKD_and_transplant.csv"))
clin$lab <- F$lab[match(clin$sig_id, F$sig_id)]
pc_long <- rbind(data.frame(lab = clin$lab, model = "Unadjusted", beta = clin$ERCB_beta_raw),
                 data.frame(lab = clin$lab, model = "Adjusted for composition and tubular injury", beta = clin$ERCB_beta_adj))
pc_long$lab <- factor(pc_long$lab, levels = clin$lab[order(clin$ERCB_beta_raw)])
pc_long$model <- factor(pc_long$model, levels = c("Unadjusted", "Adjusted for composition and tubular injury"))
p9b <- ggplot(pc_long, aes(beta, lab, shape = model, fill = model)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = .3, colour = "grey40") +
  geom_line(aes(group = lab), linewidth = .35, colour = "grey55") +
  geom_point(size = 1.8, stroke = .45, colour = "#1F5A8A") +
  scale_shape_manual(values = c(16, 21), name = NULL) +
  scale_fill_manual(values = c("#1F5A8A", "white"), name = NULL) +
  labs(x = "Apparent ageing shift, CKD vs living donors (SD; ERCB)", y = NULL) +
  theme(legend.position = "top", legend.direction = "vertical", legend.margin = margin(0, 0, 0, 0))

pd_long <- rbind(data.frame(lab = clin$lab, model = "Unadjusted", rho = clin$TX_poorGFR_rho_raw),
                 data.frame(lab = clin$lab, model = "Adjusted for donor age", rho = clin$TX_poorGFR_rho_adj_age),
                 data.frame(lab = clin$lab, model = "Adjusted for donor age and composition", rho = clin$TX_poorGFR_rho_adj_age_comp))
pd_long$lab <- factor(pd_long$lab, levels = clin$lab[order(clin$TX_poorGFR_rho_raw)])
pd_long$model <- factor(pd_long$model, levels = c("Unadjusted", "Adjusted for donor age", "Adjusted for donor age and composition"))
p9c <- ggplot(pd_long, aes(rho, lab, colour = model, shape = model)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = .3, colour = "grey40") +
  geom_vline(xintercept = 0.301, linetype = 3, linewidth = .35, colour = "#B2182B") +
  geom_line(aes(group = lab), linewidth = .3, colour = "grey70") +
  geom_point(size = 1.7, stroke = .45, fill = "white") +
  annotate("text", x = 0.295, y = 2.2, label = "Donor age\nalone (ρ = 0.30)", hjust = 1, size = TX * 0.9, colour = "grey20") +
  scale_colour_manual(values = c("#1D63B8", "#C87028", "#008B7D"), name = NULL) +
  scale_shape_manual(values = c(16, 17, 21), name = NULL) +
  labs(x = "Correlation with poor graft function (ρ)", y = NULL) +
  theme(legend.position = "top", legend.direction = "vertical", legend.margin = margin(0, 0, 0, 0))

fig9 <- p9a / (p9b + p9c + plot_layout(widths = c(1, 1))) + plot_layout(heights = c(0.9, 1.25)) + plot_annotation(tag_levels = "a")
save_fig(fig9, "Figure9_clinical", 183, 165)
write.csv(pb_long, file.path(FD, "SourceData_Figure9a.csv"), row.names = FALSE)
write.csv(pc_long, file.path(FD, "SourceData_Figure9b.csv"), row.names = FALSE)
write.csv(pd_long, file.path(FD, "SourceData_Figure9c.csv"), row.names = FALSE)

## ---- Figure 1: overview with v7 numbering ----
f1_pan <- readRDS(file.path(FD3, "panels_fig1.rds"))
bx <- function(x, y, w, h, lab, fill, col = "black", sz = TX, ff = "plain", tc = "black") list(
  annotate("rect", xmin = x - w / 2, xmax = x + w / 2, ymin = y - h / 2, ymax = y + h / 2, fill = fill, colour = col, linewidth = .25),
  annotate("text", x = x, y = y, label = lab, size = sz, lineheight = .9, fontface = ff, colour = tc))
hd <- function(x, lab) annotate("text", x = x, y = 10.55, label = lab, size = fs(7), fontface = "bold")
ar <- function(x0, x1, y = 5.2) annotate("segment", x = x0, xend = x1, y = y, yend = y, arrow = arrow(length = unit(1.6, "mm"), type = "closed"), linewidth = .45, colour = "grey35")
lt <- function(c) grDevices::adjustcolor(c, alpha.f = .16)
pa1 <- ggplot() + coord_cartesian(xlim = c(0, 20), ylim = c(0, 11), expand = FALSE) + theme_void() + theme(plot.tag = element_text(size = 9, face = "bold")) +
  hd(2.6, "Gene-set library (n = 48)") + hd(7.7, "Normal & clinical kidney data") + hd(12.75, "Evaluation & controls") + hd(17.6, "Outputs") +
  bx(2.6, 8.6, 5.0, 2.2, "A  Kidney age-derived (12)\nRodwell, Melk, Rowland, Kainz,\nRAAG, MHC-II, KAS-94/82", lt(cls_col[["A"]]), cls_col[["A"]]) +
  bx(2.6, 5.6, 5.0, 2.2, "B  General ageing (7)\nde Magalhães, Palmer, Sood,\nGenAge, CellAge, SenMayo, Fridman", lt(cls_col[["B"]]), cls_col[["B"]]) +
  bx(2.6, 2.4, 5.0, 2.7, "C  Disease-derived 'ageing'\nsets (29): DKD, IgAN, lupus,\ntransplant, RCC; 268 records\nscreened → 29 full gene lists", lt(cls_col[["C"]]), cls_col[["C"]]) +
  bx(7.7, 9.05, 4.7, 1.45, "6 GEO donor cohorts: 441 donors\n(9 series, 498 samples reviewed)", lt(src_col[["G6"]]), src_col[["G6"]]) +
  bx(7.7, 7.25, 4.7, 1.35, "TCGA normal (129); CPTAC\nadjacent (75); GTEx cortex (85)", lt(src_col[["TCGA"]]), src_col[["TCGA"]]) +
  bx(7.7, 5.35, 4.7, 1.65, "Clinical biopsies (n = 428):\nERCB 9 CKD diseases (190),\nGSE175759 eGFR (87),\nGSE147451 graft outcome (151)", lt("#C87028"), "#C87028") +
  bx(7.7, 2.35, 4.7, 2.6, "Leakage rule: a dataset used to\nderive a set is never used to\nevaluate it (KAS: no GEO;\nRowland: no TCGA)", "white", "grey40", ff = "italic") +
  bx(12.75, 8.75, 4.8, 1.9, "Gene-level enrichment vs\n2,000 expression-matched\nrandom sets", "grey96", "grey40") +
  bx(12.75, 6.45, 4.8, 1.9, "Signed score vs age;\nrandom-effects pooling;\n500 random signatures", "grey96", "grey40") +
  bx(12.75, 4.25, 4.8, 1.8, "Gene & score composition\nadjustment + KPMP snRNA-seq\ndeconvolution (MuSiC, Bisque)", "grey96", "grey40") +
  bx(12.75, 2.1, 4.8, 1.8, "Co-primary without CPTAC;\npower by set size;\n47 GTEx tissues", "grey96", "grey40") +
  bx(17.6, 8.3, 4.1, 2.4, "Benchmark, power,\ncomposition & tissue\nsource (Figs 2–8)", "white", "black", sz = TXL) +
  bx(17.6, 5.3, 4.1, 2.4, "Clinical evaluation in\nCKD & transplantation\n(Fig. 9)", "white", "black", sz = TXL) +
  bx(17.6, 2.3, 4.1, 2.4, "kidneyAgeBench v0.2.0\nR package (Table 1,\nFig. S2)", "grey25", "grey25", sz = TXL, ff = "bold", tc = "white") +
  ar(5.1, 5.35) + ar(10.05, 10.35) + ar(15.15, 15.5)
fig1 <- pa1 / (f1_pan$pb + f1_pan$pc + plot_layout(widths = c(1, 1.05))) + plot_layout(heights = c(1.05, 1)) + plot_annotation(tag_levels = "a")
save_fig(fig1, "Figure1_overview", 183, 140)

## ---- unchanged figures under v7 numbers ----
cp <- c(Figure2_signature_audit = "Figure2_signature_audit", Figure3_scorecard = "Figure3_scorecard", Figure4_power = "Figure4_power",
        Figure5_why_labels = "Figure6_why_labels", Figure6_signed_replication = "Figure7_signed_replication",
        Figure7_tissue_source = "Figure8_tissue_source", FigureS2_gtex_multitissue = "FigureS1_gtex_multitissue",
        FigureS1_kidneyAgeBench = "FigureS2_kidneyAgeBench")
for (old in names(cp)) for (ext in c("png", "pdf", "svg", "tiff"))
  file.copy(file.path(FD3, paste0(old, ".", ext)), file.path(FD, paste0(cp[[old]], ".", ext)), overwrite = TRUE)
cat(sprintf("deconvolution concordance: MuSiC %.3f, Bisque five-subclass %.3f (Spearman, %d signatures)\n", r_mu, r_bi, nrow(dc)))
cat("v7 figures saved to", FD, "\n")
