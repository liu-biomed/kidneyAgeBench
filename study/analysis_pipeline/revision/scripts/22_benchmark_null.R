# POST-HOC (labelled): pooled negative control for M-signed. For each signed signature, 500 random sets with the same
# numbers of up and down genes (drawn from genes measured in CPTAC), scored in exactly the same leakage-free datasets and pooled.
suppressMessages({library(metafor)}); REV <- path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030")
source(file.path(REV, "scripts", "04_locked_definitions.R")); O <- file.path(REV, "benchmark", "results"); set.seed(20261031)
src <- readLines(file.path(REV, "scripts", "21_benchmark.R")); eval(parse(text=src[1:grep("^Z <- lapply", src)-1]))   # rebuild datasets D, sigs, leak()
# null replicates use a closed-form DerSimonian-Laird pool (robust; no iterative fitting)
srma <- function(yi, sei, method) { ok <- is.finite(yi) & is.finite(sei) & sei > 0; yi <- yi[ok]; sei <- sei[ok]
  if (!length(yi)) return(list(b=NA, se=NA)); if (length(yi) == 1) return(list(b=yi, se=sei)); w <- 1/sei^2; fe <- sum(w*yi)/sum(w)
  Q <- sum(w*(yi-fe)^2); t2 <- max(0, (Q-(length(yi)-1))/(sum(w)-sum(w^2)/sum(w))); ws <- 1/(sei^2+t2); list(b=sum(ws*yi)/sum(ws), se=sqrt(1/sum(ws))) }
score <- function(e, up, dn){ u <- mz(e, up, 3); d <- mz(e, dn, 3); if (all(is.na(d))) u else if (all(is.na(u))) -d else u - d }
sp <- function(a, v, C) { if (sum(is.finite(v)) < 10) return(data.frame(rho=NA, n=NA, q=NA)); tryCatch(pspear(a, v, C), error=function(e) data.frame(rho=NA, n=NA, q=NA)) }
pooled_rho <- function(up, dn, lk){ parts <- list()
  if (!"G6" %in% lk) { g <- names(D)[grepl("^G6_", names(D))]; r <- sapply(g, function(k) sp(D[[k]]$age, score(D[[k]]$e, up, dn), D[[k]]$base)$rho)
    nq <- sapply(g, function(k) length(D[[k]]$age) - (if (is.null(D[[k]]$base)) 0 else ncol(model.matrix(~., D[[k]]$base)) - 1) - 3); ok <- is.finite(r)
    f <- srma(atanh(r[ok]), 1/sqrt(nq[ok]), "REML"); parts[["G6"]] <- c(f$b, f$se) }
  for (k in setdiff(c("TCGA","CPTAC"), lk)) { p <- sp(D[[k]]$age, score(D[[k]]$e, up, dn), D[[k]]$base); if (!is.finite(p$rho)) next; parts[[k]] <- c(atanh(p$rho), 1/sqrt(p$n - p$q - 3)) }
  P <- do.call(rbind, parts); tanh(srma(P[,1], P[,2], if (nrow(P) > 2) "REML" else "FE")$b) }
pool_genes <- rownames(D$CPTAC$e)
B <- read.csv(file.path(O, "BENCHMARK_summary.csv")); signed <- B$sig_id[is.finite(B$pooled_rho)]
N <- do.call(rbind, lapply(signed, function(id){ s <- sigs[[id]]; lk <- leak(id)
  nu <- sum(s$dir > 0 & s$gene %in% pool_genes); nd <- sum(s$dir < 0 & s$gene %in% pool_genes)
  nul <- replicate(500, { g <- sample(pool_genes, nu + nd); pooled_rho(g[seq_len(nu)], g[nu + seq_len(nd)], lk) })
  obs <- B$pooled_rho[B$sig_id==id]
  data.frame(sig_id=id, n_up=nu, n_down=nd, observed=obs, null_median=median(nul), null_95=quantile(nul, .95), emp_p=(1+sum(nul >= obs))/501) }))
N$fdr <- p.adjust(N$emp_p, "BH"); print(format(N, digits=2), row.names=FALSE); write.csv(N, file.path(O, "POSTHOC_M_signed_random_null.csv"), row.names=FALSE)
