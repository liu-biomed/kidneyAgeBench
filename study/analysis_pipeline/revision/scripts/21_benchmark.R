# Addendum 5 (locked): benchmark of 34 kidney-ageing gene sets in leakage-free datasets
suppressMessages({library(data.table); library(metafor)})
REV <- path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030")
source(file.path(REV, "scripts", "04_locked_definitions.R")); O <- file.path(REV, "benchmark", "results"); dir.create(O, showWarnings=FALSE); set.seed(20261031)
lib <- read.csv(file.path(REV, "benchmark", "signatures", "signature_library_long.csv")); sigs <- split(lib, lib$sig_id)
## ---- datasets: list(expr, age, base covariates)
D <- list()
cb <- readRDS(path.expand("~/kidney_aging_data/harmonized/harmonized_bundle_corrected.rds"))
for (g in names(cb$harm)) { h <- cb$harm[[g]]; D[[paste0("G6_", g)]] <- list(e=h$expr, age=h$meta$age, base=if (all(is.na(h$meta$sex))) NULL else data.frame(sex=factor(h$meta$sex)), grp="G6") }
x <- fread(path.expand("~/kidney_aging_data/validation/tcga/tcga_normal_tpm.tsv"), data.table=FALSE); tm <- read.csv(path.expand("~/kidney_aging_data/validation/tcga/tcga_normal_meta.csv"))
sx <- read.csv(path.expand("~/kidney_aging_data/validation/tcga/tcga_case_sex.csv")); tm$sex <- sx$sex[match(tm$patient, sx$patient)]; tm$age <- ifelse(is.finite(tm$age_days), tm$age_days/365.25, tm$age_at_index)
tm <- tm[order(tm$sample), ]; tm <- tm[!duplicated(tm$patient) & is.finite(tm$age) & !is.na(tm$sex), ]
e <- log2(as.matrix(x[, tm$file_id]) + 1); o <- order(rowMeans(e), decreasing=TRUE); e <- e[o,]; gg <- x$gene_name[o]; e <- e[!duplicated(gg),]; rownames(e) <- gg[!duplicated(gg)]; e <- e[rowMeans(e > 0) >= 0.2, ]
D$TCGA <- list(e=e, age=tm$age, base=data.frame(sex=factor(tm$sex), project=factor(tm$project)), grp="TCGA")
source(file.path(REV, "scripts", "10_cptac_helpers.R")); O <- file.path(REV, "benchmark", "results") ; V <- path.expand("~/kidney_aging_data/validation/cptac")
mp <- cptac_rd(V, "bcm-ccrcc-mapping"); e2s <- setNames(mp$gene_name, sub("\\..*", "", mp$gene)); C <- cptac_clin(V)
tr <- cptac_rd(V, "bcm-ccrcc-transcriptomics"); tr <- cbind(gene=tr[[1]], tr[, grep("_A$", names(tr))]); R <- cptac_sym(tr, e2s); colnames(R) <- sub("_A$", "", colnames(R)); R <- R[rowMeans(R > 0) >= 0.2, ]
cs <- intersect(colnames(R), C$case[is.finite(C$age) & C$sex %in% c("Male","Female")]); mm <- C[match(cs, C$case), ]
D$CPTAC <- list(e=R[, cs], age=mm$age, base=data.frame(sex=factor(mm$sex)), grp="CPTAC")
gt <- cb$gtex; D$GTEx <- list(e=gt$expr, age=gt$meta$age, base=data.frame(sex=factor(gt$meta$sex)), grp="GTEx")
leak <- function(id) if (grepl("^A1[12]_", id)) "G6" else if (grepl("^A0[56]_Rowland", id)) "TCGA" else character()
## ---- gene-level age z per dataset (partial correlation on ranks, Fisher z)
gene_z <- function(d){ X <- if (is.null(d$base)) matrix(1, length(d$age)) else model.matrix(~., d$base); q <- ncol(X) - 1
  Hm <- X %*% solve(crossprod(X)) %*% t(X); a <- rank(d$age, na.last="keep"); a <- a - Hm %*% a
  Rk <- t(apply(d$e, 1, rank, na.last="keep")); Rk <- Rk - t(Hm %*% t(Rk)); r <- as.numeric((Rk/sqrt(rowSums(Rk^2))) %*% (a/sqrt(sum(a^2))))
  z <- atanh(pmin(pmax(r, -.9999), .9999)) * sqrt(length(d$age) - q - 3); setNames(z, rownames(d$e)) }
Z <- lapply(D, gene_z)
## ---- M-signed per dataset
score <- function(e, s){ up <- s$gene[s$dir > 0]; dn <- s$gene[s$dir < 0]; u <- mz(e, up, 3); d <- mz(e, dn, 3)
  if (all(is.na(u)) && all(is.na(d))) return(NULL); if (all(is.na(d))) u else if (all(is.na(u))) -d else u - d }
signed <- names(sigs)[sapply(sigs, function(s) any(s$dir != 0))]
MS <- do.call(rbind, lapply(signed, function(id){ s <- sigs[[id]]; do.call(rbind, lapply(names(D), function(dn){ d <- D[[dn]]; v <- score(d$e, s); if (is.null(v)) return(NULL)
  px <- as.data.frame(lapply(LOCK$proxy, function(g) mz(d$e, g))); C2 <- if (is.null(d$base)) px else cbind(d$base, px)
  rbind(cbind(sig_id=id, dataset=dn, grp=d$grp, model="base", pspear(d$age, v, d$base)), cbind(sig_id=id, dataset=dn, grp=d$grp, model="composition", pspear(d$age, v, C2))) })) }))
write.csv(MS, file.path(O, "M_signed_by_dataset.csv"), row.names=FALSE)
poolfun <- function(r, n, q) { ok <- is.finite(r); if (sum(ok) == 0) return(c(NA,NA,NA,NA)); if (sum(ok) == 1) { z <- atanh(r[ok]); se <- 1/sqrt(n[ok]-q[ok]-3); return(c(r[ok], tanh(z-1.96*se), tanh(z+1.96*se), 2*pnorm(-abs(z/se)))) }
  f <- rma(yi=atanh(r[ok]), sei=1/sqrt(n[ok]-q[ok]-3), method=if (sum(ok) > 2) "REML" else "FE"); c(tanh(f$b), tanh(f$ci.lb), tanh(f$ci.ub), f$pval) }
MSp <- do.call(rbind, lapply(split(MS, paste(MS$sig_id, MS$model)), function(d){ id <- d$sig_id[1]; lk <- leak(id)
  g6 <- d[d$grp=="G6",]; g6p <- if (nrow(g6) && !"G6" %in% lk) poolfun(g6$rho, g6$n, g6$q) else rep(NA, 4)
  parts <- rbind(if (!"G6" %in% lk && nrow(g6)) data.frame(r=g6p[1], n=NA, q=NA, z=atanh(g6p[1]), se=(atanh(g6p[3])-atanh(g6p[2]))/3.92),
                 do.call(rbind, lapply(c("TCGA","CPTAC"), function(k){ x <- d[d$grp==k,]; if (!nrow(x) || k %in% lk) return(NULL); data.frame(r=x$rho, n=x$n, q=x$q, z=atanh(x$rho), se=1/sqrt(x$n-x$q-3)) })))
  f <- rma(yi=parts$z, sei=parts$se, method=if (nrow(parts) > 2) "REML" else "FE")
  gtx <- d[d$grp=="GTEx",]
  data.frame(sig_id=id, model=d$model[1], G6_rho=g6p[1], G6_p=g6p[4], TCGA_rho=if ("TCGA" %in% lk) NA else d$rho[d$grp=="TCGA"][1], CPTAC_rho=d$rho[d$grp=="CPTAC"][1],
             pooled_rho=tanh(f$b), pooled_lo=tanh(f$ci.lb), pooled_hi=tanh(f$ci.ub), pooled_p=f$pval, k_datasets=nrow(parts), I2=if (nrow(parts) > 2) f$I2 else NA,
             GTEx_rho=gtx$rho[1], GTEx_p=gtx$p[1], leakage_excluded=paste(lk, collapse=";")) }))
## ---- M-gene (all sets) using leakage-free pooled gene z; null matched on CPTAC expression mean/variance deciles
evalsets <- c("G6","TCGA","CPTAC")
pooled_gene_z <- function(excl){ use <- names(D)[sapply(D, function(d) d$grp %in% setdiff(evalsets, excl))]
  g <- Reduce(union, lapply(Z[use], names)); M <- sapply(Z[use], function(z) z[g]); rownames(M) <- g
  k <- rowSums(is.finite(M)); s <- rowSums(M, na.rm=TRUE)/sqrt(pmax(k,1)); s[k < max(2, floor(length(use)/2))] <- NA; s }
mu <- rowMeans(D$CPTAC$e); sdv <- apply(D$CPTAC$e, 1, sd); bin <- setNames(paste(cut(rank(mu), 10, labels=FALSE), cut(rank(sdv), 10, labels=FALSE)), names(mu))
PZ <- list(none=pooled_gene_z(character()), G6=pooled_gene_z("G6"), TCGA=pooled_gene_z("TCGA"))
MG <- do.call(rbind, lapply(names(sigs), function(id){ s <- sigs[[id]]; lk <- leak(id); pz <- PZ[[if (length(lk)) lk else "none"]]
  g <- intersect(unique(s$gene), names(pz)[is.finite(pz)]); g <- intersect(g, names(bin)); if (length(g) < 2) return(NULL)
  obs <- mean(abs(pz[g])); pool_by_bin <- split(names(bin)[names(bin) %in% names(pz)[is.finite(pz)]], bin[names(bin) %in% names(pz)[is.finite(pz)]])
  nul <- replicate(2000, mean(abs(pz[vapply(g, function(x) { pl <- pool_by_bin[[bin[x]]]; pl[sample.int(length(pl), 1)] }, "")])))
  sg <- s[match(g, s$gene), "dir"]; dirok <- sg != 0
  data.frame(sig_id=id, n_measured=length(g), mean_abs_z=obs, null_mean=mean(nul), enrichment=obs/mean(nul), emp_p=(1+sum(nul >= obs))/2001,
             pct_nominal=100*mean(abs(pz[g]) > 1.96), n_signed=sum(dirok), pct_dir_correct=if (sum(dirok)) 100*mean(sign(pz[g][dirok])==sg[dirok]) else NA,
             p_dir=if (sum(dirok)) binom.test(sum(sign(pz[g][dirok])==sg[dirok]), sum(dirok), alternative="greater")$p.value else NA) }))
MG$fdr <- p.adjust(MG$emp_p, "BH")
## ---- random-set negative control for M-signed
rnd <- do.call(rbind, lapply(1:1000, function(i){ s <- data.frame(gene=sample(rownames(D$CPTAC$e), 50), dir=sample(c(-1,1), 50, TRUE)); v <- score(D$CPTAC$e, s); pspear(D$CPTAC$age, v, D$CPTAC$base)$rho }))
## ---- assemble
summ <- read.csv(file.path(REV, "benchmark", "signatures", "signature_library_summary.csv"))
B <- merge(summ[, c("sig_id","class","label","pmid","n_genes")], MG, by="sig_id", all.x=TRUE)
B <- merge(B, MSp[MSp$model=="base", setdiff(names(MSp), "model")], by="sig_id", all.x=TRUE)
cp <- MSp[MSp$model=="composition", c("sig_id","pooled_rho","pooled_p")]; names(cp)[2:3] <- c("pooled_rho_comp","pooled_p_comp"); B <- merge(B, cp, by="sig_id", all.x=TRUE)
B <- B[order(B$class, -B$enrichment), ]; write.csv(B, file.path(O, "BENCHMARK_summary.csv"), row.names=FALSE); write.csv(MG, file.path(O, "M_gene.csv"), row.names=FALSE)
print(format(B[, c("sig_id","n_measured","enrichment","fdr","pct_nominal","pct_dir_correct","pooled_rho","pooled_p","pooled_rho_comp","GTEx_rho")], digits=2), row.names=FALSE)
## ---- pre-specified expectations
Cc <- B[B$class=="C",]; A <- B[B$sig_id %in% c("A01_Rodwell2004","A02_Rodwell2004_noblood","A03_Melk2005_505","A04_Melk2005_PAM","A05_Rowland2019_37","A06_Rowland2019_rep19","A12_KAS82_corrected"),]
sg <- B[is.finite(B$pooled_rho),]; att <- 1 - sg$pooled_rho_comp/sg$pooled_rho
E <- data.frame(expectation=c("E1 C-class median enrichment < 1.2", "E1 fewer than half of C-class FDR<0.05", "E2 A-class replication rate > C-class age-informative rate", "E3 median attenuation > 30%"),
  value=c(median(Cc$enrichment), mean(Cc$fdr < .05), sprintf("%.2f vs %.2f", mean(A$pooled_rho > 0 & A$pooled_p < .05), mean(Cc$fdr < .05)), median(att[sg$pooled_rho > 0])),
  met=c(median(Cc$enrichment) < 1.2, mean(Cc$fdr < .05) < .5, mean(A$pooled_rho > 0 & A$pooled_p < .05) > mean(Cc$fdr < .05), median(att[sg$pooled_rho > 0]) > .3))
print(E, row.names=FALSE); write.csv(E, file.path(O, "EXPECTATIONS.csv"), row.names=FALSE)
cat("random 50-gene signed sets in CPTAC: rho 95% range", round(quantile(rnd, c(.025,.975)), 3), "\n")
saveRDS(list(Z=Z, PZ=PZ, MS=MS, rnd=rnd), path.expand("~/kidney_aging_data/benchmark_intermediates.rds"))
