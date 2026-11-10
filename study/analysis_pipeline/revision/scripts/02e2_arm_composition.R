# Age association of each KAS-94 arm before/after adjustment for estimated composition (three deconvolution methods)
source(file.path(path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030"), "scripts", "fun.R"))
b <- load_bundle(); harm <- b$harm; P <- readRDS(MP); rob <- P$rob
pcor <- function(x, y, C=NULL){ if (is.null(C) || ncol(C)==0) return(cor(x, y, method="s")); cor(residuals(lm(rank(x) ~ ., data=C)), residuals(lm(rank(y) ~ ., data=C))) }
out <- list()
for (meth in c("bisque","music","nnls")) {
  dec <- read.csv(file.path(REV, "..", "tables", sprintf("deconv_proportions_%s_class.csv", meth)))
  cl <- intersect(c("epithelial","immune","stroma","endothelial"), names(dec))
  for (g in names(harm)) { h <- harm[[g]]; dd <- dec[match(first_gsm(h$meta$gsm), dec$gsm), cl]; if (any(is.na(dd))) next
    s <- score_z(h$expr, rob$gene[rob$direction=="up"], rob$gene[rob$direction=="down"]); a <- h$meta$age
    sx <- if (g %in% NO_SEX) NULL else data.frame(sex=factor(h$meta$sex)); C <- if (is.null(sx)) dd else cbind(sx, dd)
    out[[length(out)+1]] <- data.frame(method=meth, cohort=g, n=ncol(h$expr),
      up=pcor(a, s$up, sx), up_comp=pcor(a, s$up, C), dn=pcor(a, s$dn, sx), dn_comp=pcor(a, s$dn, C), score=pcor(a, s$score, sx), score_comp=pcor(a, s$score, C)) } }
O <- do.call(rbind, out); wr(O, "R1_arm_age_composition_by_cohort.csv")
pl <- do.call(rbind, lapply(split(O, O$method), function(d) do.call(rbind, lapply(c("up","up_comp","dn","dn_comp","score","score_comp"), function(v)
  cbind(method=d$method[1], measure=v, pool_rho(d[[v]], d$n - ifelse(grepl("comp", v), 5, 1)))))))
print(format(pl, digits=3), row.names=FALSE); wr(pl, "R1_arm_age_composition_pooled.csv")
