# Loaders and the fixed "full model" covariates for CPTAC NAT proteome analyses (marker lists identical to script 09)
suppressMessages(library(data.table))
REV <- path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030")
source(file.path(REV, "scripts", "04_locked_definitions.R")); O <- file.path(REV, "results", "validation")
CONTAM <- list(RBC=c("HBA1","HBA2","HBB","HBD","CA1","ALAS2","SPTA1","SLC4A1"),
  plasma=c("ALB","APOA1","APOA2","APOB","FGA","FGB","FGG","HP","HPX","TF","SERPINA1","AHSG","APOH","ITIH2","AMBP"),
  immunoglobulin=c("IGHG1","IGHG2","IGKC","IGHA1","IGHM","JCHAIN","IGLC2"), tumour_ccRCC=c("CA9","NDUFA4L2","ANGPTL4","EGLN3","VIM","PLIN2"))
EXT <- list(immune=c("PTPRC","LCP1","CORO1A","LAPTM5","CD53","LSP1","ARHGDIB","CD74","HLA-DRA","HLA-DRB1","LCP2","FCER1G","TYROBP","CD14","CD163","MPEG1","MNDA","CYBB","ITGB2","CD48"),
  stroma=c("PDGFRB","PDGFRA","THY1","ACTA2","MYH11","TAGLN","CNN1","COL1A1","COL1A2","COL3A1","COL6A1","COL6A2","COL6A3","DCN","LUM","POSTN","FBLN1","MFAP4"))
full_cov <- function(P, sex){ px <- as.data.frame(lapply(LOCK$proxy[c("endothelial","distal_tub","PT_content")], function(g) mz(P, g)))
  px$immune <- mz(P, EXT$immune); px$stroma <- mz(P, EXT$stroma); I <- as.data.frame(lapply(CONTAM, function(g) mz(P, g, 3)))
  cbind(data.frame(sex=factor(sex)), px, I) }
cptac_rd <- function(V, p) suppressWarnings(fread(cmd=paste("gzcat", shQuote(list.files(V, pattern=p, full.names=TRUE)[1])), data.table=FALSE))
cptac_sym <- function(df, e2s){ id <- sub("\\..*", "", df[[1]]); m <- as.matrix(df[,-1]); sy <- e2s[id]; ok <- !is.na(sy) & sy != ""
  m <- m[ok,,drop=FALSE]; sy <- sy[ok]; o <- order(rowMeans(m, na.rm=TRUE), decreasing=TRUE); m <- m[o,,drop=FALSE]; sy <- sy[o]; m <- m[!duplicated(sy),,drop=FALSE]; rownames(m) <- sy[!duplicated(sy)]; m }
cptac_clin <- function(V){ cl <- cptac_rd(V, "mssm-all"); cl <- cl[cl$tumor_code=="CCRCC",]
  data.frame(case=cl$case_id, age=suppressWarnings(as.numeric(cl[["consent/age"]])), sex=cl[["consent/sex"]]) }
