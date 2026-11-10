# Reviewer 1 point 24: empirical null for the combined reversal statistic, ranking stability, and overlap genes
source(file.path(path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030"), "scripts", "fun.R"))
set.seed(20261030)
P <- readRDS(MP); rob <- P$rob; U <- unique(toupper(P$res$gene)); NU <- length(U)
up <- toupper(rob$gene[rob$direction=="up"]); dn <- toupper(rob$gene[rob$direction=="down"])
L <- readLines(file.path(D,"suppl","creeds_single_drug_v1.gmt")); ds <- list()
for (l in L) { p <- strsplit(l, "\t", fixed=TRUE)[[1]]; if (length(p) < 3) next; ds[[p[1]]] <- intersect(unique(toupper(trimws(p[-(1:2)]))), U) }
drugs <- unique(sub("-(up|dn)$", "", names(ds))); drugs <- drugs[paste0(drugs,"-up") %in% names(ds) & paste0(drugs,"-dn") %in% names(ds)]
hp <- function(a, b) phyper(length(intersect(a,b))-1, length(b), NU-length(b), length(a), lower.tail=FALSE)
stat <- function(u, d, du, dd){ r1 <- hp(u, dd); r2 <- hp(d, du); m1 <- hp(u, du); m2 <- hp(d, dd)
  c(comb=pchisq(-2*(log(max(r1,1e-300))+log(max(r2,1e-300))), 4, lower.tail=FALSE), rev=-log10(max(r1,1e-300))-log10(max(r2,1e-300)),
    net=(-log10(max(r1,1e-300))-log10(max(r2,1e-300))) - (-log10(max(m1,1e-300))-log10(max(m2,1e-300)))) }
obs <- t(sapply(drugs, function(x) stat(up, dn, ds[[paste0(x,"-up")]], ds[[paste0(x,"-dn")]])))
B <- 5000
# Label-permutation null: random 65/29 signatures from the measured universe; same drug sets
null <- mclapply(seq_len(B), function(i){ g <- sample(U, length(up)+length(dn)); u <- g[seq_along(up)]; d <- g[-seq_along(up)]
  t(sapply(drugs, function(x) stat(u, d, ds[[paste0(x,"-up")]], ds[[paste0(x,"-dn")]])[c("rev","net")])) }, mc.cores=NCORE)
NR <- sapply(null, function(m) m[, "rev"]); NN <- sapply(null, function(m) m[, "net"])
emp_rev <- (1 + rowSums(NR >= obs[, "rev"])) / (B + 1); emp_net <- (1 + rowSums(NN >= obs[, "net"])) / (B + 1)
# Max-statistic (family-wise) calibration across compounds
fw <- sapply(seq_len(B), function(i) max(NR[, i])); emp_fwer <- sapply(obs[, "rev"], function(v) (1 + sum(fw >= v)) / (B + 1))
R <- data.frame(compound=drugs, combined_reversal_p=obs[, "comb"], fisher_fdr=p.adjust(obs[, "comb"], "BH"), reversal_score=obs[, "rev"], net=obs[, "net"],
  empirical_p_reversal=emp_rev, empirical_fdr_reversal=p.adjust(emp_rev, "BH"), empirical_p_net=emp_net, empirical_fdr_net=p.adjust(emp_net, "BH"), fwer_maxT=emp_fwer)
R$rank_net <- rank(-R$net, ties.method="min"); R$rank_fisher <- rank(R$combined_reversal_p, ties.method="min"); R$rank_emp_net <- rank(R$empirical_p_net, ties.method="min")
R$overlap_ageup_drugdn <- sapply(drugs, function(x) paste(sort(intersect(up, ds[[paste0(x,"-dn")]])), collapse=";"))
R$overlap_agedn_drugup <- sapply(drugs, function(x) paste(sort(intersect(dn, ds[[paste0(x,"-up")]])), collapse=";"))
R <- R[order(R$rank_net), ]; wr(R, "R1_drug_reversal_empirical.csv")
key <- c("Calcitriol","Imiquimod","Rosiglitazone","Colchicine","Pioglitazone","Acetaminophen","Troglitazone")
print(format(R[R$compound %in% key | R$rank_net <= 10, c("compound","rank_net","rank_fisher","rank_emp_net","combined_reversal_p","fisher_fdr","empirical_p_reversal","empirical_fdr_reversal","empirical_p_net","empirical_fdr_net","fwer_maxT")], digits=2), row.names=FALSE)
cat("Spearman rank agreement net vs fisher:", round(cor(R$rank_net, R$rank_fisher, method="s"), 3), "| net vs empirical-net:", round(cor(R$rank_net, R$rank_emp_net, method="s"), 3), "\n")
cat("compounds empirical FDR<0.05 (reversal):", sum(R$empirical_fdr_reversal < 0.05), "| Fisher FDR<0.05:", sum(R$fisher_fdr < 0.05), "| FWER<0.05:", sum(R$fwer_maxT < 0.05), "\n")
