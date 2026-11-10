# Addendum 2 R1 in an independent protein cohort: CPTAC non-ccRCC study (PDC000464, TMT11), normal adjacent kidney
source(file.path(path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030"), "scripts", "10_cptac_helpers.R"))
suppressMessages(library(jsonlite)); W <- path.expand("~/kidney_aging_data/validation/cptac_nonccrcc")
q <- fromJSON(file.path(W, "quant_PDC000464.json")); hdr <- q[1,]; X <- q[-1,]; P <- suppressWarnings(matrix(as.numeric(X[,-1]), nrow(X))); rownames(P) <- X[,1]
colnames(P) <- sub(".*:", "", hdr[-1]); b <- fromJSON(file.path(W, "bio_PDC000464.json"))
nb <- b[b$sample_type=="Solid Tissue Normal", ]; nb <- nb[nb$aliquot_submitter_id %in% colnames(P), ]; nb <- nb[!duplicated(nb$case_submitter_id), ]
P <- P[, nb$aliquot_submitter_id]; colnames(P) <- nb$case_submitter_id; P <- P[rowMeans(is.finite(P)) >= 0.5, ]
C <- cptac_clin(path.expand("~/kidney_aging_data/validation/cptac"))
cl <- suppressWarnings(data.table::fread(cmd=paste("gzcat", shQuote(list.files(path.expand("~/kidney_aging_data/validation/cptac"), pattern="mssm", full.names=TRUE))), data.table=FALSE))
m <- C[match(colnames(P), C$case), ]; m$tumour <- cl[["baseline/histologic_type"]][match(colnames(P), cl$case_id)]
cat("non-ccRCC normals:", ncol(P), "| proteins:", nrow(P), "| age", min(m$age), "-", max(m$age), "| F", sum(m$sex=="Female"), "\n"); print(table(m$tumour, useNA="a"))
FC <- full_cov(P, m$sex); cat("covariates finite per column:\n"); print(colSums(is.finite(as.matrix(FC[,-1]))))
# covariates that cannot be computed in this proteome (all NA) are dropped and reported; pspear() fails on all-NA columns
dropped <- names(FC)[sapply(FC, function(v) all(is.na(v)))]; cat("dropped (not measurable):", paste(dropped, collapse=", "), "\n"); FC <- FC[, !names(FC) %in% dropped, drop=FALSE]
sm <- LOCK$outcome$SenMayo; cat("SenMayo proteins measured:", length(intersect(sm, rownames(P))), "\n")
# n=19: the full model (10 covariates) leaves 6 df; report the locked full model and a reduced model (sex + PT content + plasma)
r <- rbind(cbind(model="R1 full model (locked)", pspear(m$age, mz(P, sm), FC)),
           cbind(model="reduced: sex + PT content + plasma", pspear(m$age, mz(P, sm), FC[, c("sex","PT_content","plasma")])),
           cbind(model="sex only", pspear(m$age, mz(P, sm), FC[, "sex", drop=FALSE])))
for (o in c("IFNG_response","PT_injury","PT_identity","OXPHOS")) r <- rbind(r, cbind(model=paste("sex only:", o), pspear(m$age, mz(P, LOCK$outcome[[o]]), FC[, "sex", drop=FALSE])))
print(format(r, digits=3), row.names=FALSE); write.csv(r, file.path(O, "R1_nonccRCC_PDC000464.csv"), row.names=FALSE)
