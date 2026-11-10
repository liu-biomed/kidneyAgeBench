# Shared setup for the figure set of manuscript v4 (Figures 1-7, Supplementary Figure S1). All values are read from saved benchmark outputs.
# Style rules: no panel titles (identification goes in axis labels, strips and legends); no text below 6 pt at print size;
# each colour family has one meaning across all figures:
#   class (A/B/C) = blue / teal / orange; benchmark status = grey ramp (dark = FDR < 0.05);
#   tissue source = purple / lavender / magenta / grey; warning or invalid null = dark red.
suppressMessages({library(ggplot2); library(patchwork); library(ggrepel); library(ggbeeswarm); library(ggnewscale)})
REV <- path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030")
RR <- file.path(REV, "benchmark", "results"); FD <- file.path(REV, "benchmark", "figures", "v3"); dir.create(FD, showWarnings = FALSE, recursive = TRUE)
src <- readLines(file.path(REV, "scripts", "30_figures.R"))
eval(parse(text = src[grep("^theme_set", src):grep("^pal <-", src)]))
eval(parse(text = src[grep("^cshort <-", src):(grep("^cid <-", src) - 1)]))
fs <- function(pt) pt / ggplot2::.pt   # font size in pt -> geom_text size (mm)
TX <- fs(6); TXL <- fs(6.5)
theme_update(axis.text = element_text(size = 6, colour = "black"), axis.title = element_text(size = 7), legend.title = element_text(size = 6.5), legend.text = element_text(size = 6),
  strip.text = element_text(size = 6.5, face = "bold"), plot.title = element_blank(), plot.subtitle = element_blank(), plot.tag = element_text(size = 9, face = "bold"),
  legend.key.size = unit(3, "mm"))
cls_col <- c(A = "#1F5A8A", B = "#5E9C97", C = "#C8763A")
cls_name <- c(A = "A  Kidney age-derived", B = "B  General ageing / senescence", C = "C  'Ageing-related' disease models")
src_col <- c(G6 = "#5B4A8B", TCGA = "#9C8CC4", CPTAC = "#B03A6F", GTEx = "#8C8C8C")
src_name <- c(G6 = "GEO donor cohorts", TCGA = "TCGA normal", CPTAC = "CPTAC tumour-adjacent", GTEx = "GTEx post-mortem")
warn_col <- "#B2182B"
F <- read.csv(file.path(RR, "BENCHMARK_FINAL_all_sets.csv"), colClasses = c(pmid = "character"))
S <- read.csv(file.path(RR, "SENS_gene_level_G6_TCGA_only.csv")); names(S)[names(S) %in% c("enrichment","fdr","emp_p","n")] <- c("n_sens","enrichment_sens","emp_p_sens","fdr_sens"); F <- merge(F, S[, c("sig_id","n_sens","enrichment_sens","fdr_sens")], by = "sig_id", all.x = TRUE)
N <- read.csv(file.path(RR, "POSTHOC_M_signed_random_null.csv")); names(N)[names(N) == "emp_p"] <- "emp_p_null"; F <- merge(F, N[, c("sig_id","null_95","emp_p_null")], by = "sig_id", all.x = TRUE)
cid <- sub("_.*", "", F$sig_id)
F$lab <- ifelse(F$sig_id %in% names(short), short[F$sig_id], paste0(cshort[cid], " [", sub(".*PMID", "", F$sig_id), "]"))
F$status <- ifelse(F$fdr_all < .05, "FDR < 0.05", ifelse(F$enrichment >= 1, "n.s., above random", "Below random"))
F$status_sens <- ifelse(F$fdr_sens < .05, "FDR < 0.05", ifelse(F$enrichment_sens >= 1, "n.s., above random", "Below random"))
st_col <- c("FDR < 0.05" = "#2B2B2B", "n.s., above random" = "#9A9A9A", "Below random" = "#FFFFFF")
BI <- readRDS(path.expand("~/kidney_aging_data/benchmark_intermediates.rds")); FDv <- readRDS(path.expand("~/kidney_aging_data/figure_data_v2.rds"))
NUL <- readRDS(path.expand("~/kidney_aging_data/benchmark_signed_null_full.rds"))
MS <- read.csv(file.path(RR, "M_signed_by_dataset.csv"))
# co-primary gene-level analysis (G6 + TCGA, null re-matched on the same data) and power by set size (script 36)
RV <- read.csv(file.path(RR, "REVISION_gene_level_four_configs.csv")); RV <- RV[RV$config == "noCPTAC_evalbins", c("sig_id", "enrichment", "fdr")]; names(RV)[2:3] <- c("enrichment_cop", "fdr_cop")
F <- merge(F, RV, by = "sig_id", all.x = TRUE)
F$status_cop <- ifelse(F$fdr_cop < .05, "FDR < 0.05", ifelse(F$enrichment_cop >= 1, "n.s., above random", "Below random"))
PW <- read.csv(file.path(RR, "REVISION_power_by_setsize.csv"))
pw_at <- function(n, er = 1.3, cf = "primary_locked") { s <- PW[PW$config == cf & PW$true_ER == er, ]; approx(log(s$k), s$power, log(n), rule = 2)$y }
F$power13 <- 100 * sapply(F$n_measured, pw_at)
SN <- read.csv(file.path(RR, "REVISION_signed_noCPTAC.csv")); F <- merge(F, SN[, c("sig_id", "rho_noCPTAC", "lo_noCPTAC", "hi_noCPTAC", "rho_comp_noCPTAC")], by = "sig_id", all.x = TRUE)
# The base pdf() device in save_pub_r wrote blank pages for figures with Unicode glyphs (ρ, ≥, →) and this R build has no
# cairo, so the vector PDF is converted from the SVG with rsvg-convert, which embeds Arial.
save_fig <- function(p, name, w, h) { fn <- file.path(FD, name); wi <- w / 25.4; hi <- h / 25.4
  if (isTRUE(getOption("kab_nosave"))) return(invisible(NULL))   # used when script 60 rebuilds panels
  svglite::svglite(paste0(fn, ".svg"), width = wi, height = hi); print(p); dev.off()
  stopifnot(system2(Sys.which("rsvg-convert"), c("-f", "pdf", "-o", shQuote(paste0(fn, ".pdf")), shQuote(paste0(fn, ".svg")))) == 0)
  ragg::agg_tiff(paste0(fn, ".tiff"), width = wi, height = hi, units = "in", res = 600, compression = "lzw"); print(p); dev.off()
  ragg::agg_png(paste0(fn, ".png"), width = wi, height = hi, units = "in", res = 300); print(p); dev.off() }
