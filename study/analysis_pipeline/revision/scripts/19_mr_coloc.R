# Addendum 4 (with recorded deviation): lead cis-pQTL Wald-ratio MR and coloc.abf (implemented directly, Wakefield ABF, default priors)
suppressMessages(library(data.table)); REV <- path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030")
O <- file.path(REV, "results", "validation"); M <- path.expand("~/kidney_aging_data/validation/mr")
# harmonised (hm_) columns are preferred: in the eGFR files raw beta refers to the unharmonised effect allele
rd <- function(f){ x <- fread(cmd=paste("gzcat", shQuote(f)), data.table=FALSE); names(x) <- tolower(names(x))
  pick <- function(...) { for (k in c(...)) if (k %in% names(x)) return(x[[k]]); rep(NA, nrow(x)) }
  data.frame(chr=as.character(pick("hm_chrom","chromosome","chr")), pos=as.numeric(pick("hm_pos","base_pair_location","pos")),
             ea=toupper(pick("hm_effect_allele","effect_allele")), oa=toupper(pick("hm_other_allele","other_allele")),
             beta=as.numeric(pick("hm_beta","beta")), se=as.numeric(pick("standard_error","se")), p=as.numeric(pick("p_value","pval")),
             eaf=as.numeric(pick("hm_effect_allele_frequency","effect_allele_frequency","eaf")), n=as.numeric(pick("n","sample_size")), region=pick("region"), stringsAsFactors=FALSE) }
fl <- list.files(M, pattern="\\.tsv\\.gz$", full.names=TRUE); key <- sub("_GCST.*", "", basename(fl)); names(fl) <- key
out <- lapply(c("eGFRcrea","eGFRcys"), function(k) rd(fl[[k]])); names(out) <- c("eGFRcrea","eGFRcys")
labf <- function(b, se, W) { r <- W/(W + se^2); z <- b/se; 0.5*(log(1 - r) + r*z^2) }   # Wakefield log ABF
lse <- function(x) { m <- max(x); m + log(sum(exp(x - m))) }
coloc_abf <- function(l1, l2, p1=1e-4, p2=1e-4, p12=1e-5){ lsum <- l1 + l2
  lH <- c(0, log(p1)+lse(l1), log(p2)+lse(l2), log(p1)+log(p2)+ (function(a,b){ x <- lse(a)+lse(b); y <- lse(a+b); if (x > y) x + log(1-exp(y-x)) else -Inf })(l1, l2), log(p12)+lse(lsum))
  pp <- exp(lH - lse(lH)); setNames(pp, paste0("PP.H", 0:4)) }
res <- list()
for (g in setdiff(key, c("eGFRcrea","eGFRcys"))) { e <- rd(fl[[g]]); if (nrow(e) < 10) { res[[length(res)+1]] <- data.frame(protein=g, outcome="both", n_snps=0, note="no cis variants in pQTL file (chrX not released)"); next }; e <- e[is.finite(e$beta) & is.finite(e$se) & e$se > 0, ]
  ne <- if (all(is.na(e$n))) c(BGN=2935, CXCL12=2935, OGN=2498, CCN5=3283)[g] else median(e$n, na.rm=TRUE); if (is.na(ne)) ne <- 3506
  for (k in names(out)) { o <- out[[k]]; o <- o[o$region==g & is.finite(o$beta) & is.finite(o$se) & o$se > 0, ]
    m <- merge(e, o, by=c("chr","pos"), suffixes=c(".e",".o"))
    m <- m[(m$ea.e==m$ea.o & m$oa.e==m$oa.o) | (m$ea.e==m$oa.o & m$oa.e==m$ea.o), ]; flip <- m$ea.e != m$ea.o; m$beta.o[flip] <- -m$beta.o[flip]
    m <- m[!duplicated(paste(m$chr, m$pos)), ]
    if (nrow(m) < 50) { res[[length(res)+1]] <- data.frame(protein=g, outcome=k, n_snps=nrow(m), note="too few shared variants"); next }
    lead <- m[which.min(m$p.e), ]; F <- (lead$beta.e/lead$se.e)^2
    wr <- lead$beta.o/lead$beta.e; wse <- abs(lead$se.o/lead$beta.e); wp <- 2*pnorm(-abs(wr/wse))
    pp <- coloc_abf(labf(m$beta.e, m$se.e, 0.15^2), labf(m$beta.o, m$se.o, 0.15^2))
    res[[length(res)+1]] <- data.frame(protein=g, outcome=k, n_snps=nrow(m), lead_pos=lead$pos, lead_p_pqtl=lead$p.e, F_stat=F, instrument_ok=lead$p.e < 5e-8,
      wald_ratio=wr, wald_se=wse, wald_p=wp, t(pp), note="") } }
R <- rbindlist(res, fill=TRUE); R <- as.data.frame(R); ok <- R$outcome=="eGFRcrea" & R$instrument_ok %in% TRUE
R$bonf_p <- NA; R$bonf_p[ok] <- pmin(1, R$wald_p[ok]*sum(ok)); R$ancestry_mismatch <- R$protein %in% c("BGN","CXCL12")
print(format(R, digits=3), row.names=FALSE); write.csv(R, file.path(O, "MR_coloc_results.csv"), row.names=FALSE)
