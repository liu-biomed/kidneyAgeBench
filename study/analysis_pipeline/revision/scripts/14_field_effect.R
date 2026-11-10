# POST-HOC: is the CPTAC protein-level ageing signal a tumour field effect? Y-loss as a field-effect marker in men;
# + VHL-loss/hypoxia programme (ccRCC hallmark) and ccRCC admixture as covariates
source(file.path(path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030"), "scripts", "10_cptac_helpers.R"))
suppressMessages(library(msigdbr)); hh <- as.data.frame(msigdbr(species="Homo sapiens", collection="H")); c2 <- as.data.frame(msigdbr(species="Homo sapiens", collection="C2"))
V <- path.expand("~/kidney_aging_data/validation/cptac"); mp <- cptac_rd(V, "bcm-ccrcc-mapping"); e2s <- setNames(mp$gene_name, sub("\\..*", "", mp$gene)); C <- cptac_clin(V)
P <- cptac_sym(cptac_rd(V, "bcm-ccrcc-proteomics.*Normal"), e2s); P <- P[rowMeans(is.finite(P)) >= 0.5, ]
cs <- intersect(colnames(P), C$case[is.finite(C$age) & C$sex %in% c("Male","Female")]); P <- P[, cs]; m <- C[match(cs, C$case), ]
FC <- full_cov(P, m$sex); Y <- c("RPS4Y1","DDX3Y","EIF1AY","KDM5D","UTY","USP9Y","ZFY","NLGN4Y","TXLNGY")
hyp <- unique(hh$gene_symbol[hh$gs_name=="HALLMARK_HYPOXIA"]); FC$hypoxia <- mz(P, hyp)
sets <- list(SenMayo=LOCK$outcome$SenMayo, core_matrisome=unique(c2$gene_symbol[c2$gs_name=="NABA_CORE_MATRISOME"]),
             secreted=unique(c2$gene_symbol[c2$gs_name=="NABA_SECRETED_FACTORS"]), IFNG=LOCK$outcome$IFNG_response, OXPHOS=LOCK$outcome$OXPHOS, PT_identity=LOCK$outcome$PT_identity)
ml <- m$sex=="Male"; ys <- rep(NA, ncol(P)); ys[ml] <- mz(P[, ml], Y, 3)
cat("Y-score vs hypoxia (men):", round(cor(ys[ml], FC$hypoxia[ml], method="s"), 3), "| Y vs ccRCC admixture:", round(cor(ys[ml], FC$tumour_ccRCC[ml], method="s"), 3), "\n")
hiY <- ml & ys >= median(ys[ml], na.rm=TRUE)
rows <- list()
for (s in names(sets)) { v <- mz(P, sets[[s]])
  rows[[length(rows)+1]] <- cbind(set=s, analysis="all, full model", pspear(m$age, v, FC[, setdiff(names(FC), "hypoxia")]))
  rows[[length(rows)+1]] <- cbind(set=s, analysis="all, full + hypoxia", pspear(m$age, v, FC))
  rows[[length(rows)+1]] <- cbind(set=s, analysis="women only, full + hypoxia", pspear(m$age[!ml], v[!ml], FC[!ml, setdiff(names(FC), "sex")]))
  rows[[length(rows)+1]] <- cbind(set=s, analysis="men, full + hypoxia + Y-score", pspear(m$age[ml], v[ml], cbind(FC[ml, setdiff(names(FC), "sex")], y=ys[ml])))
  rows[[length(rows)+1]] <- cbind(set=s, analysis="men with retained Y (>= median), full + hypoxia", pspear(m$age[hiY], v[hiY], FC[hiY, setdiff(names(FC), "sex")])) }
R <- do.call(rbind, rows); print(format(R[, c("set","analysis","n","q","rho","lo","hi","p")], digits=2), row.names=FALSE)
write.csv(R, file.path(O, "POSTHOC_field_effect_sensitivity.csv"), row.names=FALSE)
