# Addendum 3 Part B: cell of origin of the human index matrisome genes and within-fibroblast transcription vs age (Lake 2023)
suppressMessages(library(data.table)); REV <- path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030")
source(file.path(REV, "scripts", "04_locked_definitions.R")); O <- file.path(REV, "results", "validation"); L <- path.expand("~/kidney_aging_data/validation/lake")
idx <- c("FGL1","OGN","CCN5","SRPX2","BGN","VWA2","MFAP5","EDIL3","COL15A1","EFEMP1","COL12A1","THBS4","DPT","FBN1","CXCL12")
ct <- fread(file.path(L, "lake_matrisome_by_celltype.csv.gz"), data.table=FALSE); ct <- ct[ct$set=="normal" & ct$gene %in% idx, ]
top <- do.call(rbind, lapply(split(ct, ct$gene), function(d){ d <- d[order(-d$mean_cpm),]; data.frame(gene=d$gene[1], top1=d$subclass_l2[1], cpm1=round(d$mean_cpm[1],1), top2=d$subclass_l2[2], top3=d$subclass_l2[3],
  share_stroma=round(sum(d$mean_cpm[d$subclass_l1 %in% c("FIB","VSM/P")])/sum(d$mean_cpm), 2)) }))
print(top, row.names=FALSE); write.csv(top, file.path(O, "LAKE_index_gene_cell_origin.csv"), row.names=FALSE)
dd <- read.csv(file.path(L, "lake_pt_state_by_donor.csv"))
pb <- fread(file.path(L, "lake_matrisome_pseudobulk.csv.gz"), data.table=FALSE); f <- pb[pb$subclass_l1=="FIB", ]; f$lcpm <- log2(f$count/f$lib*1e6 + 1)
W <- reshape(f[, c("donor_id","gene","lcpm")], idvar="donor_id", timevar="gene", direction="wide"); rownames(W) <- W$donor_id; W <- t(as.matrix(W[,-1])); rownames(W) <- sub("^lcpm\\.", "", rownames(W))
W <- W[, colnames(W) %in% dd$donor_id[dd$n_nuclei >= 500]]; d <- dd[match(colnames(W), dd$donor_id), ]
nfib <- f$n_cells[match(paste(colnames(W), "OGN"), paste(f$donor_id, f$gene))]
sc <- mz(W, idx); core <- mz(W, rownames(W))
nd <- d$disease=="normal"
B1 <- rbind(cbind(test="B1 FIB index score ~ age, normal donors", pspear(d$age[nd], sc[nd])),
            cbind(test="B1 FIB index score ~ age, all donors | disease", pspear(d$age, sc, data.frame(dis=factor(d$disease)))),
            cbind(test="FIB matrisome-wide score ~ age, all donors | disease", pspear(d$age, core, data.frame(dis=factor(d$disease)))))
fibfrac <- read.csv(file.path(L, "lake_subclass_l2_counts.csv"), check.names=FALSE); rownames(fibfrac) <- fibfrac[[1]]; fibfrac <- fibfrac[,-1]
fs <- grep("FIB|MYOF", colnames(fibfrac), value=TRUE); fr <- rowSums(fibfrac[, fs])/rowSums(fibfrac); fr <- fr[d$donor_id]
B2 <- rbind(cbind(test="B2 fibroblast fraction ~ age, normal donors", pspear(d$age[nd], fr[nd])),
            cbind(test="B2 fibroblast fraction ~ age, all donors | disease", pspear(d$age, fr, data.frame(dis=factor(d$disease)))))
cat("FIB subclasses:", paste(fs, collapse=","), "| donors:", ncol(W), "(normal", sum(nd), ")\n"); R <- rbind(B1, B2); print(format(R, digits=2), row.names=FALSE)
write.csv(R, file.path(O, "LAKE_fibroblast_age.csv"), row.names=FALSE)
