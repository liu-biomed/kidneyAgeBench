# CPTAC ccRCC normal adjacent tissue: H1/H2 in NAT transcriptome (second confirmatory cohort) and S2 at protein level
suppressMessages(library(data.table))
REV <- path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030")
source(file.path(REV, "scripts", "04_locked_definitions.R")); V <- path.expand("~/kidney_aging_data/validation/cptac"); O <- file.path(REV, "results", "validation")
rd <- function(p) fread(cmd=paste("gzcat", shQuote(list.files(V, pattern=p, full.names=TRUE)[1])), data.table=FALSE)
mp <- rd("bcm-ccrcc-mapping"); e2s <- setNames(mp$gene_name, sub("\\..*", "", mp$gene))
cli <- rd("mssm-all_cancers-clinical"); cli <- cli[cli$tumor_code=="CCRCC", c("case_id","consent/age","consent/sex")]; names(cli) <- c("case","age","sex")
cli$age <- suppressWarnings(as.numeric(sub(">=?", "", cli$age)))
to_sym <- function(df){ id <- sub("\\..*", "", df[[1]]); m <- as.matrix(df[,-1]); sy <- e2s[id]; ok <- !is.na(sy) & sy != ""
  m <- m[ok,,drop=FALSE]; sy <- sy[ok]; o <- order(rowMeans(m, na.rm=TRUE), decreasing=TRUE); m <- m[o,,drop=FALSE]; sy <- sy[o]; m <- m[!duplicated(sy),,drop=FALSE]; rownames(m) <- sy[!duplicated(sy)]; m }
## NAT transcriptome
tr <- rd("bcm-ccrcc-transcriptomics"); tr <- cbind(gene=tr[[1]], tr[, grep("_A$", names(tr))]); E <- to_sym(tr); colnames(E) <- sub("_A$", "", colnames(E))
E <- E[rowMeans(E > 0, na.rm=TRUE) >= 0.2, ]
## NAT proteome
pr <- rd("bcm-ccrcc-proteomics.*Normal"); Pm <- to_sym(pr); Pm <- Pm[rowMeans(is.finite(Pm)) >= 0.5, ]
run_set <- function(X, label){ cs <- intersect(colnames(X), cli$case[is.finite(cli$age) & cli$sex %in% c("Male","Female")]); m <- cli[match(cs, cli$case),]
  r <- run_locked(X[, cs], m$age, data.frame(sex=factor(m$sex))); cbind(dataset=label, r, age_min=min(m$age), age_max=max(m$age), female=sum(m$sex=="Female")) }
R <- rbind(run_set(E, "CPTAC NAT transcriptome"), run_set(Pm, "CPTAC NAT proteome"))
print(format(R, digits=2), row.names=FALSE); write.csv(R, file.path(O, "CPTAC_locked_results.csv"), row.names=FALSE)
g <- function(d, o, mod) R[R$dataset==d & R$outcome==o & R$model==mod, ]
tx <- "CPTAC NAT transcriptome"; px <- "CPTAC NAT proteome"
H1 <- with(list(a=g(tx,"TLS12_chemokine","base"), b=g(tx,"TLS12_chemokine","composition")), a$rho>0 & a$p<.05 & b$rho>0 & b$p<.05)
H2 <- g(tx,"OXPHOS","composition")$p_tost < .05 & g(tx,"FATTY_ACID","composition")$p_tost < .05
S2 <- c(protein_metabolic_equivalence = g(px,"OXPHOS","composition")$p_tost < .05 & g(px,"FATTY_ACID","composition")$p_tost < .05,
        protein_IFNG_up = g(px,"IFNG_response","base")$rho > 0 & g(px,"IFNG_response","base")$p < .05,
        protein_PTidentity_down = g(px,"PT_identity","base")$rho < 0 & g(px,"PT_identity","base")$p < .05)
dec <- data.frame(hypothesis=c("H1_CPTAC_transcriptome","H2_CPTAC_transcriptome", paste0("S2_", names(S2))), outcome=c(H1, H2, S2))
print(dec, row.names=FALSE); write.csv(dec, file.path(O, "CPTAC_decisions.csv"), row.names=FALSE)
cat("NAT transcriptome n:", sum(colnames(E) %in% cli$case), "| proteome n:", sum(colnames(Pm) %in% cli$case), "| proteome genes:", nrow(Pm), "\n")
