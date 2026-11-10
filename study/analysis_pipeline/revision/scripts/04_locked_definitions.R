# LOCKED definitions for the independent validation (do not edit after the plan hash is recorded)
suppressMessages(library(msigdbr))
.c2 <- as.data.frame(msigdbr(species="Homo sapiens", collection="C2")); .h <- as.data.frame(msigdbr(species="Homo sapiens", collection="H"))
LOCK <- list(
  outcome = list(
    TLS12_chemokine = c("CCL2","CCL3","CCL4","CCL5","CCL8","CCL18","CCL19","CCL21","CXCL9","CXCL10","CXCL11","CXCL13"),
    PT_identity     = c("LRP2","CUBN","SLC34A1","SLC5A12","SLC22A6","SLC22A8","HNF4A","ALDOB","PCK1","G6PC1"),
    PT_injury       = c("VCAM1","HAVCR1","LCN2","SPP1","CD44","SOX9","KRT8","KRT18","CCL2","TNC"),
    OXPHOS          = sort(unique(.h$gene_symbol[.h$gs_name=="HALLMARK_OXIDATIVE_PHOSPHORYLATION"])),
    FATTY_ACID      = sort(unique(.h$gene_symbol[.h$gs_name=="HALLMARK_FATTY_ACID_METABOLISM"])),
    IFNG_response   = sort(unique(.h$gene_symbol[.h$gs_name=="HALLMARK_INTERFERON_GAMMA_RESPONSE"])),
    SenMayo         = sort(unique(.c2$gene_symbol[.c2$gs_name=="SAUL_SEN_MAYO"])),
    macrophage      = c("CD68","CD163","C1QA","C1QB","C1QC","MS4A7","CSF1R","FCGR3A")),
  # marker-based composition proxies (applicable to any bulk RNA or proteome dataset)
  proxy = list(
    immune      = c("PTPRC","LAPTM5","CD53","CORO1A","LCP1","CD37"),
    stroma      = c("PDGFRB","PDGFRA","THY1","ACTA2","MYH11","TAGLN"),
    endothelial = c("PECAM1","VWF","CDH5","EMCN","PLVAP","KDR"),
    distal_tub  = c("PAX8","PAX2","CDH16","HNF1B","SLC12A1","SLC12A3","AQP2","UMOD"),
    PT_content  = c("LRP2","CUBN","SLC34A1","SLC5A12","SLC22A6","SLC22A8","HNF4A","ALDOB","PCK1","G6PC1")),
  # covariates used for composition adjustment of each outcome (PT outcomes are not adjusted for PT content: circular)
  adjust = list(TLS12_chemokine=c("immune","stroma","endothelial","distal_tub","PT_content"),
                PT_identity=c("immune","stroma","endothelial","distal_tub"),
                PT_injury=c("immune","stroma","endothelial","distal_tub","PT_content"),
                OXPHOS=c("immune","stroma","endothelial","distal_tub","PT_content"),
                FATTY_ACID=c("immune","stroma","endothelial","distal_tub","PT_content"),
                IFNG_response=c("immune","stroma","endothelial","distal_tub","PT_content"),
                SenMayo=c("immune","stroma","endothelial","distal_tub","PT_content"),
                macrophage=c("immune","stroma","endothelial","distal_tub","PT_content")),
  min_genes = 5, equivalence_bound = 0.20, alpha = 0.05)
mz <- function(e, g, min_genes=LOCK$min_genes){ g <- intersect(g, rownames(e)); if (length(g) < min_genes) return(rep(NA_real_, ncol(e)))
  z <- t(scale(t(e[g,,drop=FALSE]))); z[!is.finite(z)] <- NA; colMeans(z, na.rm=TRUE) }
# rank-based partial Spearman of score with age, adjusting for covariate data frame C; CI by Fisher z with n - q - 3
pspear <- function(age, y, C=NULL){ ok <- is.finite(age) & is.finite(y); if (!is.null(C)) ok <- ok & stats::complete.cases(C)
  a <- rank(age[ok]); v <- rank(y[ok]); n <- sum(ok); q <- 0
  if (!is.null(C) && ncol(C)) { Cd <- C[ok,,drop=FALSE]; Cd <- Cd[, sapply(Cd, function(x) length(unique(x)) > 1), drop=FALSE]; q <- ncol(model.matrix(~., Cd)) - 1
    if (q > 0) { a <- residuals(lm(a ~ ., data=Cd)); v <- residuals(lm(v ~ ., data=Cd)) } }
  r <- cor(a, v); se <- 1/sqrt(n - q - 3); z <- atanh(r)
  data.frame(n=n, q=q, rho=r, lo=tanh(z-1.96*se), hi=tanh(z+1.96*se), p=2*pnorm(-abs(z/se)),
             p_tost=max(pnorm((z - atanh(-LOCK$equivalence_bound))/se, lower.tail=FALSE), pnorm((z - atanh(LOCK$equivalence_bound))/se))) }
# Apply the locked analysis to one dataset. base = data.frame of non-composition covariates (sex, project, ...)
run_locked <- function(expr, age, base){
  px <- as.data.frame(lapply(LOCK$proxy, function(g) mz(expr, g)))
  do.call(rbind, lapply(names(LOCK$outcome), function(o){ y <- mz(expr, LOCK$outcome[[o]]); ng <- length(intersect(LOCK$outcome[[o]], rownames(expr)))
    if (all(is.na(y))) return(data.frame(outcome=o, genes=ng, model=c("base","composition"), n=NA, q=NA, rho=NA, lo=NA, hi=NA, p=NA, p_tost=NA))
    pc <- px[, LOCK$adjust[[o]], drop=FALSE]; pc <- pc[, colSums(is.finite(as.matrix(pc))) > 0, drop=FALSE]
    rbind(cbind(outcome=o, genes=ng, model="base", pspear(age, y, base)),
          cbind(outcome=o, genes=ng, model="composition", pspear(age, y, if (is.null(base)) pc else cbind(base, pc)))) })) }
