# Shared functions for the revision analyses (effect sizes and meta-analysis follow worker-0/1.ipynb)
suppressMessages({library(metafor); library(parallel)})
D   <- path.expand("~/kidney_aging_data")
REV <- path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030")
KSET <- Sys.getenv("KSET", "original")   # "original" = 9 cohorts as submitted; "corrected" = 6 eligible cohorts
OT  <- file.path(REV, "results", "tables", KSET); OF <- file.path(REV, "results", "figures"); dir.create(OT, showWarnings=FALSE, recursive=TRUE)
PAR <- if (KSET == "original") list(mink=5, kmin=8, mink_loco=5, kmin_loco=7, mink_sex=5, kmin_sex=6, kgrid=c(7,8,9), mink_arr=4, kmin_arr=6) else
                               list(mink=4, kmin=5, mink_loco=3, kmin_loco=4, mink_sex=3, kmin_sex=4, kgrid=c(4,5,6), mink_arr=3, kmin_arr=4)
MP  <- file.path(D, "harmonized", paste0("meta_primary_", KSET, ".rds"))
wr  <- function(x, f) write.csv(x, file.path(OT, f), row.names=FALSE)
NCORE <- 8
NO_SEX <- c("GSE147451","GSE232249"); RNASEQ <- c("GSE93069","GSE217427")
load_bundle <- function() { b <- readRDS(file.path(D, "harmonized", if (KSET == "original") "harmonized_bundle.rds" else "harmonized_bundle_corrected.rds"))
  RNASEQ <<- intersect(RNASEQ, names(b$harm)); b }
first_gsm <- function(x) sub("\\+.*", "", x)

# Partial Pearson correlation of expression with age, adjusting for sex when available -> Fisher z, SE = 1/sqrt(n-k-3)
compute_effects <- function(expr, meta, adjust_sex=TRUE, age=meta$age){
  n <- length(age); has_sex <- adjust_sex && !all(is.na(meta$sex) | meta$sex=="")
  if (has_sex) { sexf <- factor(meta$sex); X <- model.matrix(~sexf); k <- ncol(X)-1
    age_res <- residuals(lm(age ~ sexf)); H <- X %*% solve(t(X) %*% X) %*% t(X); expr_res <- expr - t(H %*% t(expr))
  } else { k <- 0; age_res <- age - mean(age); expr_res <- expr - rowMeans(expr) }
  ar <- age_res/sqrt(sum(age_res^2)); er <- expr_res/sqrt(rowSums(expr_res^2))
  r <- as.numeric(er %*% ar); r[!is.finite(r)] <- NA; r <- pmax(pmin(r, 0.9999), -0.9999)
  data.frame(gene=rownames(expr), r=r, z=atanh(r), se=1/sqrt(n-k-3), n=n, k=k, stringsAsFactors=FALSE)
}
effects_all <- function(harm, sex_rule=function(g) !(g %in% NO_SEX), expr_slot="expr")
  setNames(lapply(names(harm), function(g) compute_effects(harm[[g]][[expr_slot]], harm[[g]]$meta, sex_rule(g))), names(harm))
zse <- function(eff){ cs <- names(eff); gs <- sort(unique(unlist(lapply(eff, `[[`, "gene"))))
  Z <- matrix(NA_real_, length(gs), length(cs), dimnames=list(gs, cs)); SE <- Z
  for (g in cs) { e <- eff[[g]]; Z[e$gene, g] <- e$z; SE[e$gene, g] <- e$se }; list(Z=Z, SE=SE) }

pool1 <- function(zi, sei, test="z"){
  ok <- is.finite(zi) & is.finite(sei); zi <- zi[ok]; sei <- sei[ok]
  fit <- tryCatch(rma(yi=zi, sei=sei, method="REML", test=test, control=list(maxiter=200, stepadj=0.5)),
                  error=function(e) tryCatch(rma(yi=zi, sei=sei, method="DL", test=test), error=function(e2) NULL))
  fit
}
meta_genes <- function(Z, SE, min_k=5, genes=NULL){
  k <- rowSums(is.finite(Z) & is.finite(SE)); if (is.null(genes)) genes <- names(k)[k >= min_k]
  out <- mclapply(genes, function(g){ zi <- Z[g,]; sei <- SE[g,]; f <- pool1(zi, sei); if (is.null(f)) return(NULL)
    ok <- is.finite(zi) & is.finite(sei); b <- as.numeric(f$b)
    data.frame(gene=g, z=b, se=as.numeric(f$se), pval=f$pval, k=sum(ok), I2=f$I2, tau2=f$tau2, QEp=f$QEp,
               sign_consist=mean(sign(zi[ok])==sign(b)), n_pos=sum(zi[ok]>0), n_neg=sum(zi[ok]<0)) }, mc.cores=NCORE)
  res <- do.call(rbind, out); res$r <- tanh(res$z); res$fdr <- p.adjust(res$pval, "BH"); res
}
lodo_pass <- function(gene, Z, SE, overall_sign){
  zi <- Z[gene,]; sei <- SE[gene,]; ok <- is.finite(zi) & is.finite(sei); zi <- zi[ok]; sei <- sei[ok]
  s <- sapply(seq_along(zi), function(i){ f <- tryCatch(rma(yi=zi[-i], sei=sei[-i], method="REML", control=list(maxiter=200, stepadj=0.5)),
        error=function(e) tryCatch(rma(yi=zi[-i], sei=sei[-i], method="DL"), error=function(e2) NULL))
        if (is.null(f)) c(NA, NA) else c(sign(as.numeric(f$b)), f$pval) })
  c(stable=all(s[1,]==overall_sign, na.rm=TRUE), maxp=max(s[2,], na.rm=TRUE))
}
select_robust <- function(res, Z, SE, fdr=0.05, sc=0.7, kmin=8, lodo=TRUE){
  cand <- res[res$fdr < fdr & res$sign_consist >= sc & res$k >= kmin, ]
  if (!nrow(cand)) return(cand)
  if (lodo) { L <- do.call(rbind, mclapply(cand$gene, function(g) lodo_pass(g, Z, SE, sign(cand$z[cand$gene==g])), mc.cores=NCORE))
    cand$lodo_stable <- L[,1]; cand$lodo_maxp <- L[,2]; cand <- cand[cand$lodo_stable==1 & cand$lodo_maxp < 0.05, ] }
  cand$direction <- ifelse(cand$z > 0, "up", "down"); cand
}
# Cohort-relative score (original definition): mean z of age-up genes minus mean z of age-down genes
score_z <- function(expr, up, dn){ up <- intersect(up, rownames(expr)); dn <- intersect(dn, rownames(expr))
  z <- t(scale(t(expr[c(up, dn), , drop=FALSE]))); z[!is.finite(z)] <- NA
  u <- colMeans(z[up,,drop=FALSE], na.rm=TRUE); d <- colMeans(z[dn,,drop=FALSE], na.rm=TRUE)
  data.frame(up=u, dn=d, score=u-d, n_up=length(up), n_dn=length(dn)) }
# Single-sample, rank-based score (no cohort-level normalisation): within-sample rank percentile
score_rank <- function(expr, up, dn){ up <- intersect(up, rownames(expr)); dn <- intersect(dn, rownames(expr))
  rk <- apply(expr, 2, function(v) rank(v, ties.method="average")/length(v))
  colMeans(rk[up,,drop=FALSE]) - colMeans(rk[dn,,drop=FALSE]) }
pool_rho <- function(rho, n, method="REML"){ zi <- atanh(rho); sei <- 1/sqrt(n-3); f <- rma(yi=zi, sei=sei, method=method)
  p <- predict(f); data.frame(k=length(rho), rho=tanh(f$b), lo=tanh(f$ci.lb), hi=tanh(f$ci.ub), p=f$pval, I2=f$I2, pi_lo=tanh(p$pi.lb), pi_hi=tanh(p$pi.ub)) }
spear <- function(x, y){ t <- suppressWarnings(cor.test(x, y, method="spearman", exact=FALSE)); c(rho=unname(t$estimate), p=t$p.value) }
