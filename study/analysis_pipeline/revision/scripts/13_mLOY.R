# EXPLORATORY (post hoc): Y-chromosome gene expression vs age in male kidneys (mosaic loss of Y as a positive control / new angle)
suppressMessages({library(data.table)})
REV <- path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030")
source(file.path(REV, "scripts", "10_cptac_helpers.R"))
Y <- c("RPS4Y1","DDX3Y","EIF1AY","KDM5D","UTY","USP9Y","ZFY","NLGN4Y","TXLNGY")
sp <- function(a, v, C=NULL){ r <- pspear(a, v, C); r }
yscore <- function(e){ g <- intersect(Y, rownames(e)); list(n=length(g), s=mz(e, g, 3)) }
out <- list()
# CPTAC NAT protein + RNA, males
V <- path.expand("~/kidney_aging_data/validation/cptac"); mp <- cptac_rd(V, "bcm-ccrcc-mapping"); e2s <- setNames(mp$gene_name, sub("\\..*", "", mp$gene)); C <- cptac_clin(V)
P <- cptac_sym(cptac_rd(V, "bcm-ccrcc-proteomics.*Normal"), e2s); tr <- cptac_rd(V, "bcm-ccrcc-transcriptomics"); tr <- cbind(gene=tr[[1]], tr[, grep("_A$", names(tr))]); R <- cptac_sym(tr, e2s); colnames(R) <- sub("_A$","",colnames(R))
for (lay in c("protein","RNA")) { X <- if (lay=="protein") P else R; cs <- intersect(colnames(X), C$case[C$sex=="Male" & is.finite(C$age)]); m <- C[match(cs, C$case),]; y <- yscore(X[, cs])
  pt <- mz(X[, cs], LOCK$proxy$PT_content); im <- mz(X[, cs], EXT$immune)
  out[[length(out)+1]] <- cbind(dataset=paste("CPTAC NAT", lay), genes=y$n, model="age only", sp(m$age, y$s))
  out[[length(out)+1]] <- cbind(dataset=paste("CPTAC NAT", lay), genes=y$n, model="+ PT content + immune", sp(m$age, y$s, data.frame(pt=pt, im=im))) }
# TCGA normals, males
x <- fread(path.expand("~/kidney_aging_data/validation/tcga/tcga_normal_tpm.tsv"), data.table=FALSE); tm <- read.csv(path.expand("~/kidney_aging_data/validation/tcga/tcga_normal_meta.csv"))
sx <- read.csv(path.expand("~/kidney_aging_data/validation/tcga/tcga_case_sex.csv")); tm$sex <- sx$sex[match(tm$patient, sx$patient)]; tm$age <- tm$age_days/365.25
e <- log2(as.matrix(x[, tm$file_id]) + 1); o <- order(rowMeans(e), decreasing=TRUE); e <- e[o,]; g <- x$gene_name[o]; e <- e[!duplicated(g),]; rownames(e) <- g[!duplicated(g)]
mi <- which(tm$sex=="male" & is.finite(tm$age)); y <- yscore(e[, mi])
out[[length(out)+1]] <- cbind(dataset="TCGA normal RNA", genes=y$n, model="age + project", sp(tm$age[mi], y$s, data.frame(pr=factor(tm$project[mi]))))
out[[length(out)+1]] <- cbind(dataset="TCGA normal RNA", genes=y$n, model="+ PT content + immune", sp(tm$age[mi], y$s, data.frame(pr=factor(tm$project[mi]), pt=mz(e[, mi], LOCK$proxy$PT_content), im=mz(e[, mi], LOCK$proxy$immune))))
# GTEx kidney cortex males
b <- readRDS(path.expand("~/kidney_aging_data/harmonized/harmonized_bundle.rds")); gt <- b$gtex; mi <- which(gt$meta$sex=="M"); y <- yscore(gt$expr[, mi])
out[[length(out)+1]] <- cbind(dataset="GTEx kidney cortex RNA", genes=y$n, model="age only", sp(gt$meta$age[mi], y$s))
# Discovery GEO cohorts with sex, males
for (gn in c("GSE60807","GSE66272","GSE93069","GSE167093","GSE36895")) { h <- b$harm[[gn]]; mi <- which(h$meta$sex=="M"); if (length(mi) < 10) next; y <- yscore(h$expr[, mi])
  if (y$n >= 3) out[[length(out)+1]] <- cbind(dataset=paste(gn, "RNA"), genes=y$n, model="age only", sp(h$meta$age[mi], y$s)) }
Rr <- do.call(rbind, out); print(format(Rr[, c("dataset","genes","model","n","rho","lo","hi","p")], digits=2), row.names=FALSE)
write.csv(Rr, file.path(O, "EXPLORATORY_mLOY_Ygenes_vs_age.csv"), row.names=FALSE)
d <- Rr[Rr$model %in% c("age only","age + project") & !grepl("protein", Rr$dataset),]; f <- metafor::rma(yi=atanh(d$rho), sei=1/sqrt(d$n-d$q-3), method="REML")
cat(sprintf("pooled RNA (males, %d cohorts, n=%d): rho=%.3f [%.3f, %.3f] p=%.2g I2=%.0f\n", nrow(d), sum(d$n), tanh(f$b), tanh(f$ci.lb), tanh(f$ci.ub), f$pval, f$I2))
