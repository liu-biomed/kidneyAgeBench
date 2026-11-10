# Figure set for the kidney-journal version (manuscript v5): panels built by scripts 41-49 are re-assembled around the
# argument "transcriptomic kidney age resembles chronic tubular injury", and two panels are added (biological
# programmes in independent cohorts; renal age genes split by concordance with CKD).
source(file.path(path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030"), "scripts", "40_fig_common.R"))
# rebuild each panel set in its own environment (panels refer to objects defined in their script), without saving
options(kab_nosave = TRUE); .libPaths(c("/tmp/kablib", .libPaths()))
pn <- function(script) { e <- new.env(); sys.source(file.path(REV, "scripts", script), envir = e); e }
f1 <- pn("41_fig1_overview.R"); f2 <- pn("42_fig2_pitfalls.R"); f3 <- pn("43_fig3_scorecard.R"); f4 <- pn("47_fig4_power.R"); f5 <- pn("44_fig4_why_labels_fail.R")
f6 <- pn("45_fig5_signed.R"); f7 <- pn("48_fig7_tissue_source.R"); s1 <- pn("46_fig6_robustness_tool.R"); s2 <- pn("49_figS2_gtex.R")
options(kab_nosave = FALSE); FD <- file.path(REV, "benchmark", "figures", "v5_KI"); dir.create(FD, showWarnings = FALSE)
tags <- plot_annotation(tag_levels = "a")
## Figure 1 design and data; Figure 2 audit; Figure 3 scorecard (unchanged layouts)
save_fig(f1$pa / (f1$pb + f1$pc + plot_layout(widths = c(1, 1.05))) + plot_layout(heights = c(1.05, 1)) + tags, "Figure1_design", 183, 140)
save_fig(f2$pa + f2$pb + f2$pc + plot_layout(widths = c(1, .9, 1.05)) + tags, "Figure2_audit", 183, 66)
save_fig(f3$pa / f3$pb + plot_layout(heights = c(1, .14)) + tags, "Figure3_scorecard", 183, 230)
## Figure 4 the label cannot be tested set by set, and its source databases do not track kidney age
save_fig(((f4$pa / f4$pb / f4$pc + plot_layout(heights = c(1, 1, 1.05))) | s2$pc) + tags, "Figure4_power_databases", 183, 200)
## Figure 5 what tracks kidney age
V <- read.csv(file.path(REV, "results/validation/VALIDATION_SUMMARY_long.csv"))
V <- V[V$dataset %in% c("Discovery (6 GEO cohorts, 441)", "TCGA normal (129)", "CPTAC NAT RNA (75)"), ]
V$ds <- factor(c("Discovery (6 GEO cohorts, 441)" = "GEO donors (441; derivation)", "TCGA normal (129)" = "TCGA normal (129; validation)", "CPTAC NAT RNA (75)" = "CPTAC adjacent (75; validation)")[V$dataset],
  levels = c("GEO donors (441; derivation)", "TCGA normal (129; validation)", "CPTAC adjacent (75; validation)"))
plab <- c(TLS12_chemokine = "Lymphoid chemokines", PT_injury = "PT injury", IFNG_response = "IFN-γ response", SenMayo = "SenMayo (senescence)", macrophage = "Macrophage",
  PT_identity = "PT identity", OXPHOS = "Oxidative phosphorylation", FATTY_ACID = "Fatty acid metabolism")
tc <- V[V$dataset == "TCGA normal (129)" & V$model == "base", ]; V$prog <- factor(plab[V$outcome], levels = plab[tc$outcome[order(tc$rho)]])
V$model <- factor(c(base = "Unadjusted", composition = "Composition-adjusted")[V$model], levels = c("Unadjusted", "Composition-adjusted"))
dcol <- c("GEO donors (441; derivation)" = src_col[["G6"]], "TCGA normal (129; validation)" = src_col[["TCGA"]], "CPTAC adjacent (75; validation)" = src_col[["CPTAC"]])
pa5 <- ggplot(V, aes(rho, prog, colour = ds)) + geom_vline(xintercept = 0, linewidth = .3, colour = "grey40") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, linewidth = .4, orientation = "y", position = position_dodge(width = .65)) +
  geom_point(aes(shape = model), size = 1.7, stroke = .45, fill = "white", position = position_dodge(width = .65)) + facet_grid(. ~ model) +
  scale_colour_manual(values = dcol, name = NULL) + scale_shape_manual(values = c(16, 21), guide = "none") +
  labs(x = "Partial Spearman ρ with age (95% CI)", y = NULL) + theme(legend.position = "top", legend.direction = "vertical", legend.margin = margin(0, 0, 0, 0), panel.spacing.x = unit(4, "mm"))
RA <- read.csv(file.path(RR, "REVISION_raag_ckd_split.csv")); RA$set <- factor(sprintf("%s\n(%d genes)", sub("RAAG, ", "", RA$set), RA$n), levels = sprintf("%s\n(%d genes)", c("not CKD-concordant", "CKD-concordant"), c(584, 31)))
RA$an <- factor(c("All eligible datasets" = "All eligible datasets", "Without tumour-adjacent tissue" = "Without CPTAC")[RA$analysis], levels = c("All eligible datasets", "Without CPTAC"))
pb5 <- ggplot(RA, aes(ER, set, shape = an)) + geom_vline(xintercept = 1, linetype = 2, linewidth = .3, colour = "grey40") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, linewidth = .45, orientation = "y", position = position_dodge(width = .5)) +
  geom_point(size = 2.1, stroke = .4, fill = cls_col[["A"]], position = position_dodge(width = .5)) +
  scale_shape_manual(values = c(21, 24), name = NULL) + scale_x_continuous(limits = c(.95, 1.92)) +
  labs(x = "Enrichment vs matched genes (95% null interval)", y = "Renal age-associated genes") + theme(legend.position = "top", legend.direction = "vertical", legend.margin = margin(0, 0, 0, 0))
save_fig((pa5 + pb5 + plot_layout(widths = c(1.55, 1))) / (f5$pe + f6$pd + plot_layout(widths = c(.85, 1.25))) + plot_layout(heights = c(1.1, 1)) + tags, "Figure5_kidney_age_axis", 183, 170)
write.csv(V, file.path(FD, "SourceData_Figure5a.csv"), row.names = FALSE); write.csv(RA, file.path(FD, "SourceData_Figure5b.csv"), row.names = FALSE)
## Figure 6 signed signatures, composition and tissue source
save_fig((f6$pa + f6$pc + f7$pb + plot_layout(widths = c(1.3, .8, .7))) / f7$pa + plot_layout(heights = c(1.25, .8)) + tags, "Figure6_signatures_composition_source", 183, 165)
## Supplementary figures
save_fig(s1$pb + s1$pc + plot_layout(widths = c(1, 1)) + tags, "FigureS1_kidneyAgeBench", 140, 80)
save_fig(s2$pa / s2$pb + plot_layout(heights = c(1.15, 1)) + tags, "FigureS2_gtex_reference_sets", 120, 150)
save_fig((f5$pa + f5$pb + plot_layout(widths = c(1.1, 1))) / (free(f5$pd) + plot_spacer() + plot_layout(widths = c(1, .6))) + tags, "FigureS3_why_labels_look_age_associated", 183, 150)
save_fig(f6$pb + tags, "FigureS4_per_cohort", 110, 120)
save_fig(f7$pc + f4$pd + plot_layout(widths = c(1.4, 1)) + tags, "FigureS5_without_CPTAC_expected_counts", 183, 95)
cat("done\n")
