# Confirmatory validation in TCGA kidney Solid Tissue Normal, exactly per prereg/LOCKED_VALIDATION_PLAN.md
suppressMessages({library(data.table); library(metafor)})
REV <- path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030")
source(file.path(REV, "scripts", "04_locked_definitions.R")); V <- path.expand("~/kidney_aging_data/validation/tcga"); O <- file.path(REV, "results", "validation")
x <- fread(file.path(V, "tcga_normal_tpm.tsv"), data.table=FALSE); m <- read.csv(file.path(V, "tcga_normal_meta.csv"))
sx <- read.csv(file.path(V, "tcga_case_sex.csv")); m$sex <- sx$sex[match(m$patient, sx$patient)]
m$age <- ifelse(is.finite(m$age_days), m$age_days/365.25, m$age_at_index)
m <- m[order(m$sample), ]; m <- m[!duplicated(m$patient) & is.finite(m$age) & !is.na(m$sex), ]
mat <- as.matrix(x[, m$file_id]); keep <- !duplicated(x$gene_name) | TRUE
# symbol level: highest-mean Ensembl entry per symbol (protein-coding and lncRNA alike)
e <- log2(mat + 1); o <- order(rowMeans(e), decreasing=TRUE); e <- e[o,]; g <- x$gene_name[o]; e <- e[!duplicated(g),]; rownames(e) <- g[!duplicated(g)]
e <- e[rowMeans(e > 0) >= 0.2, ]
cat("TCGA normals:", ncol(e), "| genes:", nrow(e), "| age", round(min(m$age)), "-", round(max(m$age)), "| F/M", sum(m$sex=="female"), "/", sum(m$sex=="male"), "\n"); print(table(m$project))
base <- data.frame(sex=factor(m$sex), project=factor(m$project))
R <- run_locked(e, m$age, base); R$score <- "z"
# S4 robustness: within-sample rank score
rk <- apply(e, 2, function(v) rank(v)/length(v)); R2 <- run_locked(rk, m$age, base); R2$score <- "rank"
# S4: per project
PP <- do.call(rbind, lapply(split(seq_len(nrow(m)), m$project), function(i) cbind(project=m$project[i[1]], run_locked(e[, i], m$age[i], data.frame(sex=factor(m$sex[i]))))))
all <- rbind(R, R2); write.csv(all, file.path(O, "TCGA_locked_results.csv"), row.names=FALSE); write.csv(PP, file.path(O, "TCGA_locked_by_project.csv"), row.names=FALSE)
write.csv(data.frame(n=ncol(e), age_min=min(m$age), age_median=median(m$age), age_max=max(m$age), female=sum(m$sex=="female"), KIRC=sum(m$project=="TCGA-KIRC"), KIRP=sum(m$project=="TCGA-KIRP"), KICH=sum(m$project=="TCGA-KICH")), file.path(O, "TCGA_cohort.csv"), row.names=FALSE)
print(format(all, digits=2), row.names=FALSE)
## Pre-specified decisions
g <- function(o, mod, s="z") all[all$outcome==o & all$model==mod & all$score==s, ]
H1 <- g("TLS12_chemokine","base")$rho > 0 & g("TLS12_chemokine","base")$p < .05 & g("TLS12_chemokine","composition")$rho > 0 & g("TLS12_chemokine","composition")$p < .05
H2 <- g("OXPHOS","composition")$p_tost < .05 & g("FATTY_ACID","composition")$p_tost < .05
S1 <- all[all$model=="base" & all$score=="z" & all$outcome %in% c("PT_identity","PT_injury","SenMayo","IFNG_response","macrophage"), c("outcome","rho","lo","hi","p")]
S1$expected <- c(PT_identity=-1, PT_injury=1, SenMayo=1, IFNG_response=1, macrophage=1)[S1$outcome]; S1$p_holm <- p.adjust(S1$p, "holm"); S1$supported <- sign(S1$rho)==S1$expected & S1$p_holm < .05
cat("\n==== H1 (TLS-priming, composition-independent):", ifelse(H1, "CONFIRMED", "NOT CONFIRMED"), "\n==== H2 (metabolic equivalent to 0 after composition):", ifelse(H2, "CONFIRMED", "NOT CONFIRMED"), "\n"); print(S1, row.names=FALSE)
write.csv(rbind(data.frame(hypothesis=c("H1","H2"), outcome=c(H1, H2)), data.frame(hypothesis=paste0("S1_", S1$outcome), outcome=S1$supported)), file.path(O, "TCGA_decisions.csv"), row.names=FALSE)
cat("\nBy project:\n"); print(format(PP[PP$outcome %in% c("TLS12_chemokine","OXPHOS","FATTY_ACID","PT_identity"), c("project","outcome","model","n","rho","lo","hi","p")], digits=2), row.names=FALSE)
