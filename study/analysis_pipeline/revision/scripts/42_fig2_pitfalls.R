# Figure 2. The problem in miniature: a kidney ageing signature derived the conventional way looks strong in-sample,
# halves under nested held-out evaluation, vanishes in GTEx, and its gene-level concordance test is anticonservative.
source(file.path(path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030"), "scripts", "40_fig_common.R"))
T0 <- file.path(REV, "results/tables")
L <- rbind(cbind(set = "KAS-94\n(9 series)", read.csv(file.path(T0, "original/R1_nested_LOCO_pooled.csv"))), cbind(set = "KAS-82\n(6 eligible cohorts)", read.csv(file.path(T0, "corrected/R1_nested_LOCO_pooled.csv"))))
L <- L[L$analysis != "nested LOCO (fixed top-94)", ]; L$set <- factor(L$set, levels = unique(L$set)); L$an <- factor(ifelse(grepl("in-sample", L$analysis), "In-sample (selection data)", "Nested held-out"), levels = c("In-sample (selection data)", "Nested held-out"))
# point estimate with 95% CI rather than bars: rho is an estimate, not an amount
pa <- ggplot(L, aes(an, rho)) + geom_hline(yintercept = 0, linewidth = .3, colour = "grey60") +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = .45) + geom_point(aes(fill = an), shape = 21, size = 2.2, stroke = .4) + facet_grid(. ~ set) +
  geom_text(aes(label = sprintf("%.2f", rho)), nudge_x = .28, hjust = 0, size = TX) + scale_fill_manual(values = c("white", "grey15"), guide = "none") +
  scale_x_discrete(labels = c("In-\nsample", "Nested\nheld-out"), expand = expansion(add = c(.5, .8))) + scale_y_continuous(limits = c(0, .5), expand = c(0, 0)) +
  labs(x = NULL, y = "Pooled ρ with age (95% CI)")
G <- readRDS(path.expand("~/kidney_aging_data/harmonized/gtex_scores_original.rds")); g <- data.frame(age = G$meta$age, score = G$score$score)
st <- read.csv(file.path(T0, "original/R1_GTEx_score_level.csv"))[1, ]
pb <- ggplot(g, aes(age, score)) + geom_quasirandom(width = 2.2, size = .8, alpha = .75, colour = src_col[["GTEx"]]) +
  stat_summary(fun = median, geom = "crossbar", width = 5, linewidth = .3, colour = "black") +
  annotate("text", x = 20, y = Inf, hjust = 0, vjust = 1.2, size = TX, label = sprintf("ρ = %.2f (95%% CI %.2f to %.2f)\nn = %d", st$spearman_rho, st$boot_lo, st$boot_hi, nrow(g)), lineheight = .9) +
  scale_x_continuous(breaks = seq(25, 75, 10)) + scale_y_continuous(expand = expansion(c(.05, .2))) + labs(x = "Age, GTEx cortex (decade midpoint)", y = "KAS-94 score")
nul <- read.csv(file.path(REV, "benchmark/figures/SourceData_Figure1c_null.csv"))$permutation_concordant; obs <- 68; ng <- 93
pp <- (1 + sum(nul >= obs)) / (length(nul) + 1); pbn <- binom.test(obs, ng, alternative = "greater")$p.value
bn <- data.frame(k = 0:ng, d = dbinom(0:ng, ng, .5) * length(nul))
ex <- floor(log10(pbn)); mt <- pbn / 10^ex
pc <- ggplot(data.frame(k = nul), aes(k)) + geom_histogram(binwidth = 1, fill = "#BDBDBD", colour = NA) + geom_area(data = bn, aes(k, d), fill = warn_col, alpha = .2, colour = warn_col, linewidth = .3) +
  geom_vline(xintercept = obs, colour = "black", linewidth = .5) +
  annotate("text", x = 11, y = 900, label = "Binomial (genes independent)", colour = warn_col, size = TX, hjust = 0, vjust = 1) +
  annotate("text", x = 11, y = 835, label = sprintf("italic(P) == %.1f %%*%% 10^%d", mt, ex), parse = TRUE, colour = warn_col, size = TX, hjust = 0, vjust = 1) +
  annotate("text", x = 11, y = 700, label = "Age permutation\n(correlation kept)", colour = "grey30", size = TX, lineheight = .9, hjust = 0, vjust = 1) +
  annotate("text", x = 11, y = 570, label = sprintf("italic(P) == \"%.2f\"", pp), parse = TRUE, colour = "grey30", size = TX, hjust = 0, vjust = 1) +
  annotate("text", x = obs + 1.5, y = 600, label = sprintf("Observed\n%d/%d", obs, ng), size = TX, hjust = 0, lineheight = .9) +
  scale_y_continuous(expand = expansion(c(0, .12))) + coord_cartesian(xlim = c(10, 92)) +
  labs(x = "KAS-94 genes concordant in GTEx", y = "Null frequency")
f <- pa + pb + pc + plot_layout(widths = c(1, .9, 1.05)) + plot_annotation(tag_levels = "a")
saveRDS(list(pa = pa, pb = pb, pc = pc), file.path(FD, "panels_fig2.rds"))
save_fig(f, "Figure2_signature_audit", 183, 66)
cat("perm P", pp, "\n")
