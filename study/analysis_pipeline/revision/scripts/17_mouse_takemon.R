# Addendum 3 Part A (locked): DO mouse kidney (Takemon 2021), paired RNA + protein, 188 mice, 6/12/18 months
suppressMessages({library(fgsea); library(msigdbr)})
REV <- path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030")
source(file.path(REV, "scripts", "04_locked_definitions.R")); O <- file.path(REV, "results", "validation"); set.seed(20261031)
e <- new.env(); load(path.expand("~/kidney_aging_data/validation/mouse_takemon/extracted/JAC_DO_kidney_v4.RData"), envir=e)
m <- e$dataset.mrna; p <- e$dataset.protein; s <- as.data.frame(p$annot.samples); rownames(s) <- s$Mouse.ID
Pm <- t(p$data$norm); Rm <- t(m$data$norm)               # genes x mice
ps <- toupper(p$annot.protein$symbol); rs <- toupper(m$annot.mrna$symbol)
dupP <- ps %in% ps[duplicated(ps)]; dupR <- rs %in% rs[duplicated(rs)]
Pm <- Pm[!dupP & !is.na(ps), ]; rownames(Pm) <- ps[!dupP & !is.na(ps)]; Rm <- Rm[!dupR & !is.na(rs), ]; rownames(Rm) <- rs[!dupR & !is.na(rs)]
ids <- intersect(intersect(colnames(Pm), colnames(Rm)), s$Mouse.ID); Pm <- Pm[, ids]; Rm <- Rm[, ids]; s <- s[ids, ]
cv <- data.frame(sex=factor(s$Sex), gen=factor(s$Gen), batch=factor(s$Batch), tag=factor(s$Tag)); age <- s$Age
cat("mice:", length(ids), "| ages:", paste(names(table(age)), table(age), sep=":", collapse=" "), "| proteins:", nrow(Pm), "| mRNA:", nrow(Rm), "\n")
gb <- intersect(rownames(Pm), rownames(Rm)); cat("genes with both:", length(gb), "\n")
pr <- function(y, X){ residuals(lm(rank(y, na.last="keep") ~ ., data=X)) }
res <- do.call(rbind, lapply(gb, function(g){ y <- Pm[g,]; r <- Rm[g,]; ok <- is.finite(y) & is.finite(r) & complete.cases(cv); if (sum(ok) < 150) return(NULL)
  cv_ok <- cv[ok, , drop=FALSE]; a0 <- pr(age[ok], cv_ok); p0 <- pr(y[ok], cv_ok); r0 <- pr(r[ok], cv_ok)
  X1 <- cbind(cv_ok, mrna=rank(r[ok])); a1 <- pr(age[ok], X1); p1 <- pr(y[ok], X1)
  data.frame(gene=g, prot_age=cor(a0, p0), rna_age=cor(a0, r0), ptae=cor(a1, p1)) }))
write.csv(res, file.path(O, "MOUSE_PTAE_by_gene.csv"), row.names=FALSE)
c2 <- as.data.frame(msigdbr(species="Homo sapiens", collection="C2")); core <- unique(c2$gene_symbol[c2$gs_name=="NABA_CORE_MATRISOME"])
st <- sort(setNames(res$ptae, res$gene), decreasing=TRUE)
M1 <- as.data.frame(fgsea(list(NABA_CORE_MATRISOME=core), st, minSize=10, maxSize=2000, nPermSimple=10000))
idx <- c("FGL1","OGN","CCN5","SRPX2","BGN","VWA2","MFAP5","EDIL3","COL15A1","EFEMP1","COL12A1","THBS4","DPT","FBN1","CXCL12")
ix <- res[res$gene %in% idx, ]; M2 <- binom.test(sum(ix$ptae > 0), nrow(ix), alternative="greater")
cs <- intersect(core, rownames(Pm)); M3 <- pspear(age, mz(Pm, cs), cv)
hu <- read.csv(file.path(O, "CPTAC_PTAE_by_gene.csv")); hm <- merge(hu[, c("gene","ptae")], res[, c("gene","ptae")], by="gene", suffixes=c("_human","_mouse"))
M4 <- cor.test(hm$ptae_human, hm$ptae_mouse, method="s", exact=FALSE); hmc <- hm[hm$gene %in% core, ]; M4c <- cor.test(hmc$ptae_human, hmc$ptae_mouse, method="s", exact=FALSE)
cat("\nM1 core matrisome PTAE enrichment: NES =", round(M1$NES, 2), "P =", signif(M1$pval, 3), "size", M1$size, "\n")
cat("M2 index proteins measured:", nrow(ix), "| PTAE>0:", sum(ix$ptae > 0), "| sign test P =", signif(M2$p.value, 3), "\n"); print(ix[order(-ix$ptae), ], row.names=FALSE)
cat("M3 core matrisome protein score ~ age:\n"); print(M3)
cat("M4 human vs mouse PTAE (all shared genes): rho =", round(M4$estimate, 3), "P =", signif(M4$p.value, 3), "n =", nrow(hm), "| core matrisome only: rho =", round(M4c$estimate, 3), "P =", signif(M4c$p.value, 3), "n =", nrow(hmc), "\n")
core_sum <- data.frame(median_prot_age=median(res$prot_age[res$gene %in% core]), median_rna_age=median(res$rna_age[res$gene %in% core]), median_ptae=median(res$ptae[res$gene %in% core]),
                       genome_median_ptae=median(res$ptae), pct_core_prot_up=mean(res$prot_age[res$gene %in% core] > 0), pct_core_rna_up=mean(res$rna_age[res$gene %in% core] > 0))
print(core_sum)
dec <- data.frame(hypothesis=c("M1","M2","M3","M4"), estimate=c(M1$NES, mean(ix$ptae > 0), M3$rho, M4$estimate), p=c(M1$pval, M2$p.value, M3$p, M4$p.value),
                  supported=c(M1$NES > 0 & M1$pval < .05, M2$p.value < .05, M3$rho > 0 & M3$p < .05, M4$estimate > 0 & M4$p.value < .05))
print(dec, row.names=FALSE); write.csv(dec, file.path(O, "MOUSE_decisions.csv"), row.names=FALSE); write.csv(core_sum, file.path(O, "MOUSE_core_matrisome_summary.csv"), row.names=FALSE)
