# Addendum 2 R1 (locked): SenMayo protein score ~ age under the full model, CPTAC ccRCC confirmatory cohort NAT (DIA),
# restricted to patients NOT in the discovery proteome. Intensities log2-transformed; >=50% completeness per protein.
source(file.path(path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030"), "scripts", "10_cptac_helpers.R"))
V <- path.expand("~/kidney_aging_data/validation/cptac"); W <- path.expand("~/kidney_aging_data/validation/cptac_conf/Supplementary_Processed_Proteomic_Data_for_PMID-36563681")
x <- fread(file.path(W, "CPTAC_ccRCC_global_proteome_DIA.tsv"), data.table=FALSE); g <- x$Genes
X <- log2(as.matrix(x[, -1])); X[!is.finite(X)] <- NA
o <- order(rowMeans(X, na.rm=TRUE), decreasing=TRUE); X <- X[o,]; g <- g[o]; keep <- !is.na(g) & g != "" & !duplicated(g); X <- X[keep,]; rownames(X) <- g[keep]
cl <- cptac_rd(V, "mssm-all"); cl <- cl[cl$tumor_code=="CCRCC",]
conf <- cl$case_id[grepl("Normal", cl[["confirmatory_study/type_of_analyzed_samples"]])]
disc_prot <- names(cptac_rd(V, "bcm-ccrcc-proteomics.*Normal"))[-1]
nc <- grep("_N$", colnames(X), value=TRUE); case <- sub("_N$", "", nc)
sel <- case %in% conf & !(case %in% disc_prot)
cat("NAT columns:", length(nc), "| confirmatory NAT:", sum(case %in% conf), "| overlap with discovery proteome:", sum(case %in% conf & case %in% disc_prot), "| used:", sum(sel), "\n")
P <- X[, nc[sel]]; colnames(P) <- case[sel]; P <- P[rowMeans(is.finite(P)) >= 0.5, ]
C <- cptac_clin(V); m <- C[match(colnames(P), C$case), ]; ok <- is.finite(m$age) & m$sex %in% c("Male","Female"); P <- P[, ok]; m <- m[ok, ]
cat("n =", ncol(P), "| proteins:", nrow(P), "| age", min(m$age), "-", max(m$age), "| F/M", sum(m$sex=="Female"), "/", sum(m$sex=="Male"), "\n")
FC <- full_cov(P, m$sex); dropped <- names(FC)[sapply(FC, function(v) all(is.na(v)))]; cat("covariates not measurable (dropped):", paste(dropped, collapse=", "), "\n"); FC <- FC[, !names(FC) %in% dropped, drop=FALSE]
sm <- LOCK$outcome$SenMayo; cat("SenMayo proteins measured:", length(intersect(sm, rownames(P))), "\n")
R1 <- cbind(test="R1: SenMayo protein ~ age | full model (LOCKED)", pspear(m$age, mz(P, sm), FC))
sens <- rbind(cbind(test="sex only", pspear(m$age, mz(P, sm), FC[, "sex", drop=FALSE])),
              cbind(test="women only, full model", pspear(m$age[m$sex=="Female"], mz(P, sm)[m$sex=="Female"], FC[m$sex=="Female", setdiff(names(FC), "sex")])))
Y <- c("RPS4Y1","DDX3Y","EIF1AY","KDM5D","UTY","USP9Y","ZFY"); ml <- m$sex=="Male"
if (length(intersect(Y, rownames(P))) >= 3) { ys <- mz(P[, ml], Y, 3); sens <- rbind(sens, cbind(test="men, full model + Y-score", pspear(m$age[ml], mz(P, sm)[ml], cbind(FC[ml, setdiff(names(FC), "sex")], y=ys))))
  cat("men: Y-score vs age rho =", round(cor(ys, m$age[ml], method="s", use="complete"), 3), "\n") }
other <- do.call(rbind, lapply(c("IFNG_response","PT_injury","PT_identity","OXPHOS","FATTY_ACID"), function(o) cbind(test=paste("secondary, full model:", o), pspear(m$age, mz(P, LOCK$outcome[[o]]), FC))))
c2 <- as.data.frame(msigdbr::msigdbr(species="Homo sapiens", collection="C2"))
for (s in c("NABA_CORE_MATRISOME","NABA_SECRETED_FACTORS")) other <- rbind(other, cbind(test=paste("secondary, full model:", s), pspear(m$age, mz(P, unique(c2$gene_symbol[c2$gs_name==s])), FC)))
for (p in c("CXCL12","OGN","BGN","IL1RN","TIMP2","FBN1")) if (p %in% rownames(P)) other <- rbind(other, cbind(test=paste("protein, full model:", p), pspear(m$age, as.numeric(scale(P[p,])), FC)))
R <- rbind(R1, sens, other); print(format(R, digits=3), row.names=FALSE)
cat("\n==== R1:", ifelse(R1$rho > 0 & R1$p < .05, "CONFIRMED", "NOT CONFIRMED"), "| R2: not testable (no NAT RNA for confirmatory patients in the downloaded resources)\n")
write.csv(R, file.path(O, "R1_confirmatory_PDC000411.csv"), row.names=FALSE)
