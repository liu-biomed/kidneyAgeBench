# POST-HOC diagnostic: discovery-cohort estimate restricted to the SenMayo proteins and covariates available in the confirmatory DIA data
source(file.path(path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030"), "scripts", "10_cptac_helpers.R"))
V <- path.expand("~/kidney_aging_data/validation/cptac"); W <- path.expand("~/kidney_aging_data/validation/cptac_conf/Supplementary_Processed_Proteomic_Data_for_PMID-36563681")
x <- fread(file.path(W, "CPTAC_ccRCC_global_proteome_DIA.tsv"), data.table=FALSE)
mp <- cptac_rd(V, "bcm-ccrcc-mapping"); e2s <- setNames(mp$gene_name, sub("\\..*", "", mp$gene)); C <- cptac_clin(V)
P <- cptac_sym(cptac_rd(V, "bcm-ccrcc-proteomics.*Normal"), e2s); P <- P[rowMeans(is.finite(P)) >= 0.5, ]
cs <- intersect(colnames(P), C$case[is.finite(C$age) & C$sex %in% c("Male","Female")]); P <- P[, cs]; m <- C[match(cs, C$case), ]
# confirmatory measured set (same filter as script 15)
X <- log2(as.matrix(x[, grep("_N$", names(x))])); X[!is.finite(X)] <- NA; conf_meas <- x$Genes[rowMeans(is.finite(X)) >= 0.5]
sm40 <- intersect(intersect(LOCK$outcome$SenMayo, conf_meas), rownames(P))
FC <- full_cov(P, m$sex); FCh <- FC[, setdiff(names(FC), c("endothelial","immunoglobulin","tumour_ccRCC"))]
R <- rbind(cbind(test="discovery, all SenMayo, full model", pspear(m$age, mz(P, LOCK$outcome$SenMayo), FC)),
           cbind(test=sprintf("discovery, %d shared SenMayo proteins, full model", length(sm40)), pspear(m$age, mz(P, sm40), FC)),
           cbind(test="discovery, shared proteins, harmonised covariates", pspear(m$age, mz(P, sm40), FCh)),
           cbind(test="discovery, core matrisome, harmonised covariates", pspear(m$age, mz(P, unique(as.data.frame(msigdbr::msigdbr(species="Homo sapiens", collection="C2"))$gene_symbol[as.data.frame(msigdbr::msigdbr(species="Homo sapiens", collection="C2"))$gs_name=="NABA_CORE_MATRISOME"])), FCh)))
print(format(R[, c("test","n","rho","lo","hi","p")], digits=3), row.names=FALSE); write.csv(R, file.path(O, "POSTHOC_R1_harmonised_discovery.csv"), row.names=FALSE)
# pooled discovery + confirmatory (random effects), for SenMayo and core matrisome
cf <- read.csv(file.path(O, "R1_confirmatory_PDC000411.csv"))
pool <- function(r1, n1, q1, r2, n2, q2){ f <- metafor::rma(yi=atanh(c(r1, r2)), sei=1/sqrt(c(n1-q1-3, n2-q2-3)), method="REML"); c(rho=tanh(f$b), lo=tanh(f$ci.lb), hi=tanh(f$ci.ub), p=f$pval, I2=f$I2) }
a <- R[3,]; b <- cf[cf$test=="R1: SenMayo protein ~ age | full model (LOCKED)",]; print(round(pool(a$rho, a$n, a$q, b$rho, b$n, b$q), 3))
a <- R[4,]; b <- cf[cf$test=="secondary, full model: NABA_CORE_MATRISOME",]; print(round(pool(a$rho, a$n, a$q, b$rho, b$n, b$q), 3))
