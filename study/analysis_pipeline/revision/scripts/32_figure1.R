# Figure 1: self-audit of KAS-94. Claim: circular evaluation and independence-assuming tests overstate replication.
suppressMessages({library(ggplot2); library(patchwork)})
REV <- path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030")
Sys.setenv(KSET = "original"); source(file.path(REV, "scripts", "fun.R")); FD <- file.path(REV, "benchmark", "figures")
src <- readLines(file.path(REV, "scripts", "30_figures.R")); eval(parse(text = src[grep("^theme_set", src):grep("^pal <-", src)]))
## (a) in-sample vs nested LOCO
L <- rbind(cbind(set = "KAS-94 (9 cohorts)", read.csv(file.path(REV, "results/tables/original/R1_nested_LOCO_pooled.csv"))),
           cbind(set = "KAS-82 (6 eligible cohorts)", read.csv(file.path(REV, "results/tables/corrected/R1_nested_LOCO_pooled.csv"))))
L <- L[L$analysis != "nested LOCO (fixed top-94)", ]; L$analysis <- factor(ifelse(grepl("in-sample", L$analysis), "In-sample", "Nested held-out"), levels = c("In-sample", "Nested held-out"))
pa <- ggplot(L, aes(rho, set, colour = analysis)) + geom_vline(xintercept = 0, linewidth = .3, colour = "grey40") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, orientation = "y", linewidth = .5, position = position_dodge(.5)) + geom_point(size = 1.8, position = position_dodge(.5)) +
  scale_colour_manual(values = c("In-sample" = "#8C8C8C", "Nested held-out" = "#2F6B9A"), name = NULL) + scale_x_continuous(limits = c(0, .5)) +
  labs(x = "Pooled ρ with age (95% CI)", y = NULL) + theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0))
## (b) GTEx score vs age
G <- readRDS(file.path(D, "harmonized", "gtex_scores_original.rds")); g <- data.frame(age = G$meta$age, score = G$score$score)
st <- read.csv(file.path(REV, "results/tables/original/R1_GTEx_score_level.csv"))[1, ]
pb <- ggplot(g, aes(age, score)) + geom_point(position = position_jitter(width = 1.5, seed = 1), size = .9, alpha = .7, colour = "#2F6B9A") +
  stat_summary(fun = median, geom = "point", shape = 95, size = 6, colour = "black") +
  scale_x_continuous(breaks = seq(25, 75, 10)) + labs(x = "Age (decade midpoint, years)", y = "KAS-94 score (GTEx cortex)",
  subtitle = sprintf("ρ = %.2f (95%% CI %.2f to %.2f); n = %d", st$spearman_rho, st$boot_lo, st$boot_hi, nrow(g)))
## (c) concordance: permutation null (recomputed with the original seed and procedure)
set.seed(20261030); b <- load_bundle(); gt <- b$gtex; P <- readRDS(MP); rob <- P$rob; m <- gt$meta
ge <- compute_effects(gt$expr, m, TRUE); rownames(ge) <- ge$gene; g93 <- intersect(rob$gene, ge$gene); disc <- setNames(rob$z, rob$gene)[g93]
obs <- sum(sign(ge[g93, "r"]) == sign(disc)); Xg <- gt$expr[g93, ]; B <- 10000
nul <- replicate(B, { a <- m$age; for (sx in unique(m$sex)) { i <- which(m$sex == sx); a[i] <- a[i][sample(length(i))] }; e <- compute_effects(Xg, m, TRUE, age = a); sum(sign(e$r) == sign(disc)) })
pp <- (1 + sum(nul >= obs)) / (B + 1); pbn <- binom.test(obs, length(g93), alternative = "greater")$p.value
bn <- data.frame(k = 0:length(g93), d = dbinom(0:length(g93), length(g93), .5) * B)
pc <- ggplot(data.frame(k = nul), aes(k)) + geom_histogram(binwidth = 1, fill = "#BDBDBD", colour = NA) +
  geom_line(data = bn, aes(k, d), colour = "#C0504D", linewidth = .4) + geom_vline(xintercept = obs, colour = "#2F6B9A", linewidth = .5) +
  labs(x = "Genes concordant with discovery direction", y = "Permutations",
       subtitle = sprintf("Observed %d/%d (blue); binomial P = %.1e (red)\nage-permutation P = %.2f (grey)", obs, length(g93), pbn, pp)) + coord_cartesian(xlim = c(10, 90))
pb <- pb + theme(plot.subtitle = element_text(size = 6)); pc <- pc + theme(plot.subtitle = element_text(size = 6, lineheight = .9)); f1 <- (pa | pb | pc)
f1 <- f1 + plot_layout(widths = c(1.1, 1, 1.1)) + plot_annotation(tag_levels = "a")
save_pub_r(f1, file.path(FD, "Figure1_self_audit"), 183, 70)
write.csv(data.frame(permutation_concordant = nul), file.path(FD, "SourceData_Figure1c_null.csv"), row.names = FALSE)
cat("obs", obs, "perm P", pp, "binom P", pbn, "\n")
