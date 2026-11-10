# POST-HOC (not in the locked plan): CPTAC NAT protein-level robustness to blood/plasma contamination,
# residual tumour admixture, comorbidity, and dominance by single proteins
suppressMessages(library(data.table))
REV <- path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030")
source(file.path(REV, "scripts", "04_locked_definitions.R")); V <- path.expand("~/kidney_aging_data/validation/cptac"); O <- file.path(REV, "results", "validation")
rd <- function(p) suppressWarnings(fread(cmd=paste("gzcat", shQuote(list.files(V, pattern=p, full.names=TRUE)[1])), data.table=FALSE))
mp <- rd("bcm-ccrcc-mapping"); e2s <- setNames(mp$gene_name, sub("\\..*", "", mp$gene))
to_sym <- function(df){ id <- sub("\\..*", "", df[[1]]); m <- as.matrix(df[,-1]); sy <- e2s[id]; ok <- !is.na(sy) & sy != ""
  m <- m[ok,,drop=FALSE]; sy <- sy[ok]; o <- order(rowMeans(m, na.rm=TRUE), decreasing=TRUE); m <- m[o,,drop=FALSE]; sy <- sy[o]; m <- m[!duplicated(sy),,drop=FALSE]; rownames(m) <- sy[!duplicated(sy)]; m }
cl <- rd("mssm-all_cancers-clinical"); cl <- cl[cl$tumor_code=="CCRCC",]
med <- tolower(cl[["general_medical_history/medical_condition"]])
C <- data.frame(case=cl$case_id, age=suppressWarnings(as.numeric(cl[["consent/age"]])), sex=cl[["consent/sex"]],
  htn=as.integer(grepl("hypertension", med)), dm=as.integer(grepl("diabetes", med)), bmi=suppressWarnings(as.numeric(cl[["medical_history/bmi"]])),
  stage=cl[["baseline/tumor_stage_pathological"]], size=suppressWarnings(as.numeric(cl[["baseline/tumor_size_cm"]])))
P <- to_sym(rd("bcm-ccrcc-proteomics.*Normal")); P <- P[rowMeans(is.finite(P)) >= 0.5, ]
cs <- intersect(colnames(P), C$case[is.finite(C$age) & C$sex %in% c("Male","Female")]); P <- P[, cs]; m <- C[match(cs, C$case), ]
idx <- list(RBC=c("HBA1","HBA2","HBB","HBD","CA1","ALAS2","SPTA1","SLC4A1"),
            plasma=c("ALB","APOA1","APOA2","APOB","FGA","FGB","FGG","HP","HPX","TF","SERPINA1","AHSG","APOH","ITIH2","AMBP"),
            immunoglobulin=c("IGHG1","IGHG2","IGKC","IGHA1","IGHM","JCHAIN","IGLC2"),
            tumour_ccRCC=c("CA9","NDUFA4L2","ANGPTL4","EGLN3","VIM","PLIN2"))
I <- as.data.frame(lapply(idx, function(g) mz(P, g, 3))); cov_idx <- sapply(idx, function(g) paste(intersect(g, rownames(P)), collapse=";"))
cat("contamination indices measured:\n"); print(cov_idx)
# The locked proteome analysis silently dropped the immune and stromal proxies (<5 locked markers measured after the 50% completeness
# filter). Post hoc, extended protein-detectable marker panels are used for these two proxies (defined from cell-type biology, not from age data).
ext <- list(immune=c("PTPRC","LCP1","CORO1A","LAPTM5","CD53","LSP1","ARHGDIB","CD74","HLA-DRA","HLA-DRB1","LCP2","FCER1G","TYROBP","CD14","CD163","MPEG1","MNDA","CYBB","ITGB2","CD48"),
            stroma=c("PDGFRB","PDGFRA","THY1","ACTA2","MYH11","TAGLN","CNN1","COL1A1","COL1A2","COL3A1","COL6A1","COL6A2","COL6A3","DCN","LUM","POSTN","FBLN1","MFAP4"))
locked_px <- as.data.frame(lapply(LOCK$proxy, function(g) mz(P, g)))
cat("locked proxies usable in proteome:\n"); print(colSums(is.finite(as.matrix(locked_px))) > 0)
px <- locked_px; px$immune <- mz(P, ext$immune); px$stroma <- mz(P, ext$stroma)
cat("extended immune proteins:", paste(intersect(ext$immune, rownames(P)), collapse=","), "\nextended stroma proteins:", paste(intersect(ext$stroma, rownames(P)), collapse=","), "\n")
base <- data.frame(sex=factor(m$sex))
cat("index vs age (Spearman):\n"); print(round(sapply(I, function(v) cor(v, m$age, method="s", use="complete")), 3))
outs <- c("SenMayo","IFNG_response","PT_injury","PT_identity","OXPHOS","FATTY_ACID")
mods <- list(`locked proxies only (as reported)`="LOCKED", `composition, extended immune/stroma`=NULL, `+ RBC + plasma`=c("RBC","plasma"), `+ immunoglobulin`="immunoglobulin",
             `+ tumour admixture`="tumour_ccRCC", `+ all contamination`=names(idx))
rows <- list()
for (o in outs) { y <- mz(P, LOCK$outcome[[o]]); pc <- px[, LOCK$adjust[[o]], drop=FALSE]
  for (k in names(mods)) { if (identical(mods[[k]], "LOCKED")) { lp <- locked_px[, LOCK$adjust[[o]], drop=FALSE]; lp <- lp[, colSums(is.finite(as.matrix(lp))) > 0, drop=FALSE]
      rows[[length(rows)+1]] <- cbind(outcome=o, model=k, pspear(m$age, y, cbind(base, lp))); next }
    Cv <- cbind(base, pc); if (!is.null(mods[[k]])) Cv <- cbind(Cv, I[, mods[[k]], drop=FALSE])
    rows[[length(rows)+1]] <- cbind(outcome=o, model=k, pspear(m$age, y, Cv)) }
  Cv <- cbind(base, pc, I, htn=m$htn, dm=m$dm, bmi=m$bmi); rows[[length(rows)+1]] <- cbind(outcome=o, model="+ all contamination + HTN/DM/BMI", pspear(m$age, y, Cv))
  Cv <- cbind(Cv, stage=factor(m$stage), size=m$size); rows[[length(rows)+1]] <- cbind(outcome=o, model="+ ... + tumour stage/size", pspear(m$age, y, Cv)) }
R <- do.call(rbind, rows); print(format(R[, c("outcome","model","n","q","rho","lo","hi","p")], digits=2), row.names=FALSE)
write.csv(R, file.path(O, "CPTAC_protein_contamination_sensitivity.csv"), row.names=FALSE)
## SenMayo restricted to non-plasma members, and leave-one-protein-out
sm <- intersect(LOCK$outcome$SenMayo, rownames(P)); plasma_like <- unique(c(unlist(idx[c("plasma","immunoglobulin")]), "C3","CP","SERPINE1","A2M","ORM1","ORM2","SAA1","CRP","APCS","F2","PLG","KNG1","VTN","FN1"))
ce <- cor(t(P[sm, ]), I$plasma, use="pairwise")[,1]
sets <- list(all_measured=sm, minus_plasma_list=setdiff(sm, plasma_like), minus_plasma_correlated=sm[abs(ce) < 0.3 | !is.finite(ce)])
Cfull <- cbind(base, px[, LOCK$adjust$SenMayo], I)
SR <- do.call(rbind, lapply(names(sets), function(s) cbind(set=s, proteins=length(sets[[s]]), pspear(m$age, mz(P, sets[[s]]), Cfull))))
mem <- do.call(rbind, lapply(sm, function(g) data.frame(protein=g, rho_age=cor(P[g,], m$age, method="s", use="complete"), r_plasma=ce[g])))
mem <- mem[order(-mem$rho_age), ]; loo <- sapply(sm, function(g) pspear(m$age, mz(P, setdiff(sm, g)), Cfull)$rho)
SR <- rbind(SR, cbind(set="leave-one-protein-out (min rho)", proteins=length(sm)-1, pspear(m$age, mz(P, setdiff(sm, names(which.min(loo)))), Cfull)))
print(format(SR, digits=2), row.names=FALSE); cat("top SenMayo proteins by age correlation:\n"); print(head(mem, 12), row.names=FALSE)
cat("pct SenMayo proteins with positive age rho:", round(100*mean(mem$rho_age > 0), 1), "% of", nrow(mem), "\n")
write.csv(SR, file.path(O, "CPTAC_SenMayo_subset_sensitivity.csv"), row.names=FALSE); write.csv(mem, file.path(O, "CPTAC_SenMayo_members_age.csv"), row.names=FALSE)
