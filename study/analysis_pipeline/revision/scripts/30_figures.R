# Figures 2-3 of the benchmark manuscript.
# Fig 2 claim: most "ageing-related" gene sets used in kidney research are no more age-associated in normal human kidney
#   than expression-matched random genes; human-kidney age-derived signatures are the exception.
# Fig 3 claim: signed age-derived signatures replicate modestly across independent cohorts, roughly half of the signal
#   tracks tissue composition, and post-mortem GTEx cortex does not reproduce them.
suppressMessages({library(ggplot2); library(patchwork)})
REV <- path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030")
R <- file.path(REV, "benchmark", "results"); FD <- file.path(REV, "benchmark", "figures")
theme_set(theme_classic(base_size = 6.5, base_family = "Arial") + theme(axis.line = element_line(linewidth = .35), axis.ticks = element_line(linewidth = .35),
  legend.title = element_text(size = 6.2), legend.text = element_text(size = 5.8), strip.text = element_text(size = 6.2, face = "bold"),
  strip.background = element_blank(), plot.title = element_text(size = 7, face = "bold"), plot.tag = element_text(size = 8, face = "bold"), panel.grid = element_blank()))
save_pub_r <- function(plot, filename, width_mm = 183, height_mm = 120, dpi = 600) { w <- width_mm/25.4; h <- height_mm/25.4
  svglite::svglite(paste0(filename, ".svg"), width = w, height = h); print(plot); dev.off()
  grDevices::pdf(paste0(filename, ".pdf"), width = w, height = h, family = "Helvetica", useDingbats = FALSE); print(plot + theme(text = element_text(family = "Helvetica"))); dev.off()
  ragg::agg_tiff(paste0(filename, ".tiff"), width = w, height = h, units = "in", res = dpi, compression = "lzw"); print(plot); dev.off()
  ragg::agg_png(paste0(filename, ".png"), width = w, height = h, units = "in", res = 300); print(plot); dev.off() }
short <- c(A01_Rodwell2004="Rodwell 2004", A02_Rodwell2004_noblood="Rodwell 2004 (no blood)", A03_Melk2005_505="Melk 2005 (505)", A04_Melk2005_PAM="Melk 2005 (PAM)",
  A05_Rowland2019_37="Rowland 2019 (37)", A06_Rowland2019_rep19="Rowland 2019 (replicated 19)", A07_Kainz2007_donor="Kainz 2007 donor quality",
  A08_Perco2019_RAAG634="RAAG 2019 (634, unsigned)", A09_Perco2019_RAAG_CKD="RAAG 2019 CKD-concordant", A10_Sinning2023_MHCII="MHC-II tubuloinflammaging",
  A11_KAS94_original="KAS-94 (this study, original)", A12_KAS82_corrected="KAS-82 (this study, corrected)",
  B01_deMagalhaes2009="de Magalhães 2009", B02_Palmer2021="Palmer 2021", B03_GenAge="GenAge", B04_CellAge="CellAge", B05_SenMayo="SenMayo",
  B06_Fridman2008_up="Fridman 2008", B08_Sood2015="Sood 2015 classifier")
cls_lab <- c(A = "A  Human-kidney age-derived", B = "B  General ageing / senescence resources", C = "C  'Ageing-related' kidney disease models")
pal <- c(sig = "#2F6B9A", ns = "#BDBDBD")
## ---------------- Figure 2
F <- read.csv(file.path(R, "BENCHMARK_FINAL_all_sets.csv"))
cshort <- c(C01="DN hub genes (CCR2, VCAM1...)", C02="DKD senescence signature", C03="DN glomerular senescence (FOS, ZFP36)", C04="Transplant rejection senescence genes",
  C05="IgAN ageing hub genes", C06="Renal stone ageing genes", C07="CKD ageing hub genes", C08="Lupus nephritis senescence biomarkers", C09="ccRCC ageing prognostic genes",
  C10="ccRCC ageing/senescence model", C11="ccRCC senescence score", C12="RCC normal-tissue ageing genes", C13="Renal fibrosis ageing diagnostic genes", C14="CKD ageing hub genes (AGER...)",
  C15="CKD ferro-ageing genes", C16="DKD ageing hub genes", C17="DN ageing signature (EFEMP1...)", C18="ESRD senescence signature", C19="DKD-ageing co-differential genes",
  C20="Glomerular DN senescence genes", C21="Sepsis-AKI senescence (LCN2, TIMP1)", C22="DKD VSMC senescence (MMP14, FOS)", C24="ccRCC senescence scoring model",
  C25="DN mitochondrial/senescence hub genes", C27="ccRCC senescence risk model", C28="Lupus nephritis senescence genes", C29="Membranous nephropathy senescence genes",
  C30="Diabetic glomerulopathy senescence genes", C31="Papillary RCC senescence signature")
cid <- sub("_.*", "", F$sig_id); F$lab <- ifelse(F$sig_id %in% names(short), short[F$sig_id], paste0(cshort[cid], "  [PMID ", sub(".*PMID", "", F$sig_id), "]"))
F$sig <- ifelse(F$fdr_all < .05, "sig", "ns"); F$classf <- factor(cls_lab[F$class], levels = cls_lab)
F <- F[order(F$class, F$enrichment), ]; F$labf <- factor(F$lab, levels = unique(F$lab))
pa <- ggplot(F, aes(enrichment, labf)) + geom_vline(xintercept = 1, linewidth = .3, linetype = 2, colour = "grey40") +
  geom_segment(aes(x = 1, xend = enrichment, yend = labf), linewidth = .3, colour = "grey70") +
  geom_point(aes(fill = sig, size = n_measured), shape = 21, stroke = .25) +
  facet_grid(class ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = pal, labels = c(sig = "FDR < 0.05", ns = "FDR ≥ 0.05"), name = NULL) +
  scale_size_continuous(range = c(.9, 3.2), breaks = c(3, 30, 300), trans = "log10", name = "Genes measured") +
  scale_x_continuous(trans = "log2", breaks = c(.25, .5, 1, 1.5, 2), labels = c("0.25", "0.5", "1", "1.5", "2")) +
  labs(x = "Age-association enrichment vs matched random genes (mean |z| ratio)", y = NULL) +
  theme(strip.text.y.right = element_text(angle = 0, size = 7), axis.text.y = element_text(size = 5.2), legend.position = "top", legend.box = "vertical", legend.spacing.y = unit(0, "pt"))
S <- do.call(rbind, lapply(split(F, F$class), function(d) data.frame(class = d$class[1], med = median(d$enrichment), n = nrow(d), nsig = sum(d$fdr_all < .05), nbelow = sum(d$enrichment < 1))))
S$classf <- factor(cls_lab[S$class], levels = rev(cls_lab)); F$classr <- factor(cls_lab[F$class], levels = rev(cls_lab))
pb <- ggplot(F, aes(enrichment, classr)) + geom_vline(xintercept = 1, linewidth = .3, linetype = 2, colour = "grey40") +
  geom_point(aes(fill = sig), shape = 21, size = 1.6, stroke = .25, position = position_jitter(height = .18, seed = 1)) +
  geom_point(data = S, aes(med, classf), shape = 124, size = 6, colour = "black", inherit.aes = FALSE) +
  geom_text(data = S, aes(x = .13, y = as.numeric(classf) + .38, label = sprintf("FDR<0.05: %d/%d\nbelow random: %d/%d", nsig, n, nbelow, n)), hjust = 0, size = 1.9, lineheight = .9, inherit.aes = FALSE) +
  scale_fill_manual(values = pal, guide = "none") + scale_x_continuous(trans = "log2", breaks = c(.25, .5, 1, 2), labels = c("0.25", "0.5", "1", "2"), limits = c(.12, 2.6)) +
  scale_y_discrete(labels = function(x) sub("  ", "\n", x)) + labs(x = "Enrichment (log2); bar = class median", y = NULL) + theme(axis.text.y = element_text(size = 5.8, hjust = 0.5))
f2 <- pa + pb + plot_layout(widths = c(1.25, 1)) + plot_annotation(tag_levels = "a")
save_pub_r(f2, file.path(FD, "Figure2_benchmark_enrichment"), 183, 150)
## ---------------- Figure 3
B <- read.csv(file.path(R, "BENCHMARK_summary.csv")); N <- read.csv(file.path(R, "POSTHOC_M_signed_random_null.csv"))
G <- merge(B[is.finite(B$pooled_rho), ], N[, c("sig_id","null_95","emp_p","fdr")], by = "sig_id", suffixes = c("", "_null"))
G$lab <- short[G$sig_id]; G <- G[order(G$class, G$pooled_rho), ]; G$labf <- factor(G$lab, levels = G$lab)
G$classf <- factor(c(A = "Kidney age-derived", B = "General ageing")[G$class], levels = c("Kidney age-derived", "General ageing"))
p3a <- ggplot(G, aes(pooled_rho, labf)) + geom_vline(xintercept = 0, linewidth = .3, colour = "grey40") +
  geom_point(aes(x = null_95), shape = 124, size = 3, colour = "#C0504D") +
  geom_errorbar(aes(xmin = pooled_lo, xmax = pooled_hi), width = 0, orientation = "y", linewidth = .4, colour = "#2F6B9A") +
  geom_point(size = 1.6, colour = "#2F6B9A") + facet_grid(classf ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Pooled partial Spearman ρ with age (95% CI)", y = NULL) +
  theme(strip.text.y = element_text(angle = -90, size = 5.6), axis.text.y = element_text(size = 5.4))
L <- rbind(data.frame(sig_id = G$sig_id, labf = G$labf, classf = G$classf, ds = "GEO (6 cohorts)", rho = G$G6_rho), data.frame(sig_id = G$sig_id, labf = G$labf, classf = G$classf, ds = "TCGA normal", rho = G$TCGA_rho),
           data.frame(sig_id = G$sig_id, labf = G$labf, classf = G$classf, ds = "CPTAC NAT", rho = G$CPTAC_rho), data.frame(sig_id = G$sig_id, labf = G$labf, classf = G$classf, ds = "GTEx (post-mortem)", rho = G$GTEx_rho))
L$ds <- factor(L$ds, levels = c("GEO (6 cohorts)", "TCGA normal", "CPTAC NAT", "GTEx (post-mortem)"))
p3b <- ggplot(L, aes(ds, labf, fill = rho)) + geom_tile(colour = "white", linewidth = .4) + geom_text(aes(label = ifelse(is.na(rho), "excl.", sprintf("%.2f", rho))), size = 1.7) +
  facet_grid(classf ~ ., scales = "free_y", space = "free_y") + scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC", midpoint = 0, limits = c(-.5, .5), oob = scales::squish, na.value = "grey92", name = "ρ") + guides(fill = guide_colourbar(barwidth = unit(2, "mm"), barheight = unit(18, "mm"))) +
  labs(x = NULL, y = NULL) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank(), strip.text = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1))
P <- rbind(data.frame(labf = G$labf, classf = G$classf, model = "Unadjusted", rho = G$pooled_rho), data.frame(labf = G$labf, classf = G$classf, model = "Composition-adjusted", rho = G$pooled_rho_comp))
P$model <- factor(P$model, levels = c("Unadjusted", "Composition-adjusted")); att <- median(1 - G$pooled_rho_comp/G$pooled_rho)
p3c <- ggplot(P, aes(model, rho, group = labf)) + geom_hline(yintercept = 0, linewidth = .3, colour = "grey40") + geom_line(linewidth = .3, colour = "grey60") +
  geom_point(aes(colour = classf), size = 1.3) + scale_colour_manual(values = c("Kidney age-derived" = "#2F6B9A", "General ageing" = "#8C8C8C"), name = NULL) + scale_x_discrete(labels = c("Unadjusted", "Composition-\nadjusted")) +
  annotate("text", x = 1.5, y = -.06, label = sprintf("median attenuation %.0f%%", 100*att), size = 1.9) +
  labs(x = NULL, y = "Pooled ρ with age") + theme(legend.position = "top", legend.direction = "vertical", legend.key.size = unit(2.5, "mm"), legend.text = element_text(size = 5.2))
f3 <- (p3a | p3b | p3c) + plot_layout(widths = c(1.5, 1.05, .8)) + plot_annotation(tag_levels = "a")
save_pub_r(f3, file.path(FD, "Figure3_signed_replication"), 183, 100)
write.csv(F[, c("sig_id","class","lab","n_measured","enrichment","emp_p","fdr_all")], file.path(FD, "SourceData_Figure2.csv"), row.names = FALSE)
write.csv(G[, c("sig_id","class","lab","pooled_rho","pooled_lo","pooled_hi","pooled_p","null_95","fdr","pooled_rho_comp","G6_rho","TCGA_rho","CPTAC_rho","GTEx_rho")], file.path(FD, "SourceData_Figure3.csv"), row.names = FALSE)
cat("done\n")
