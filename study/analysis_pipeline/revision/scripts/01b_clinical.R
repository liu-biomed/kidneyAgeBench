# Rebuild the three clinical cohorts as in worker-2.ipynb cells 137-153
suppressMessages({library(data.table); library(edgeR); library(AnnotationDbi); library(org.Hs.eg.db)})
D <- path.expand("~/kidney_aging_data"); OUT <- file.path(D, "harmonized")
read_sm <- function(f){ con <- gzfile(f); L <- readLines(con, warn=FALSE); close(con)
  b <- grep("^!series_matrix_table_begin", L); e <- grep("^!series_matrix_table_end", L)
  list(hdr=L[seq_len(b-1)], tab=if (e-b > 2) fread(text=L[(b+1):(e-1)], data.table=FALSE) else NULL) }
hfield <- function(hdr, key) { l <- hdr[startsWith(hdr, key)]; lapply(l, function(x) gsub('^"|"$', "", strsplit(x, "\t")[[1]][-1])) }
parse_field <- function(hdr, key){ ch <- hfield(hdr, "!Sample_characteristics_ch1"); out <- rep(NA_character_, length(ch[[1]]))
  for (v in ch) { hit <- grepl(key, v, ignore.case=TRUE); out[hit] <- sub(paste0("^.*", key, "[:\\s]*"), "", v[hit], ignore.case=TRUE, perl=TRUE) }; trimws(out) }
clin <- list()
## GSE175759 (RNA-seq counts, logCPM)
sm <- read_sm(file.path(D,"geo","GSE175759_series_matrix.txt.gz"))
meta <- data.frame(gsm=hfield(sm$hdr,"!Sample_geo_accession")[[1]], title=hfield(sm$hdr,"!Sample_title")[[1]],
  egfr=suppressWarnings(as.numeric(parse_field(sm$hdr, "estimated gfr \\(ckd-epi\\)"))), diagnosis=parse_field(sm$hdr,"diagnosis"),
  outlier=parse_field(sm$hdr,"technical outlier"), stringsAsFactors=FALSE)
td <- file.path(D,"suppl","raw175"); dir.create(td, showWarnings=FALSE); untar(file.path(D,"suppl","GSE175759_RAW.tar"), exdir=td)
fs <- list.files(td, pattern="\\.txt\\.gz$", full.names=TRUE)
lst <- lapply(fs, function(f){ dt <- fread(f, header=TRUE); setnames(dt, c("gene", sub("_.*","",basename(f)))); dt })
mat <- Reduce(function(a,b) merge(a,b,by="gene",all=TRUE), lst); cm <- as.matrix(mat[,-1]); rownames(cm) <- mat$gene
sym <- suppressMessages(mapIds(org.Hs.eg.db, keys=sub("\\..*","",rownames(cm)), keytype="ENSEMBL", column="SYMBOL", multiVals="first"))
agg <- rowsum(cm[!is.na(sym),], sym[!is.na(sym)]); meta <- meta[match(colnames(agg), meta$gsm),]
keep <- meta$outlier=="No"; agg <- agg[,keep]; meta <- meta[keep,]
dge <- DGEList(agg); dge <- dge[suppressWarnings(filterByExpr(dge)),, keep.lib.sizes=FALSE]; dge <- calcNormFactors(dge)
meta$grp <- ifelse(meta$diagnosis=="Control","Reference","Disease")
clin$GSE175759 <- list(expr=cpm(dge, log=TRUE, prior.count=1), meta=meta, counts=agg)
## GSE104954 (two Affymetrix ENTREZ-CDF platforms)
e2s <- function(m){ s <- suppressMessages(mapIds(org.Hs.eg.db, keys=sub("_at$","",rownames(m)), keytype="ENTREZID", column="SYMBOL", multiVals="first"))
  k <- !is.na(s); rs <- rowsum(m[k,,drop=FALSE], s[k]); rs / as.numeric(table(s[k])[rownames(rs)]) }
subs <- lapply(c("GPL22945","GPL24120"), function(p){ sm <- read_sm(file.path(D,"geo",paste0("GSE104954-",p,"_series_matrix.txt.gz")))
  ex <- as.matrix(sm$tab[,-1]); rownames(ex) <- sm$tab[[1]]
  list(expr=e2s(ex), meta=data.frame(gsm=colnames(ex), title=hfield(sm$hdr,"!Sample_title")[[1]], diagnosis=parse_field(sm$hdr,"diagnosis"), platform=p, stringsAsFactors=FALSE)) })
cg <- Reduce(intersect, lapply(subs, function(s) rownames(s$expr)))
ex <- do.call(cbind, lapply(subs, function(s) s$expr[cg,])); meta <- do.call(rbind, lapply(subs, `[[`, "meta"))
meta$group <- ifelse(grepl("-LD", meta$title), "Reference_LD", ifelse(!is.na(meta$diagnosis) & meta$diagnosis=="Tumor nephrectomy", "Reference_TN", "Disease"))
meta$diagnosis[is.na(meta$diagnosis) & grepl("-LD", meta$title)] <- "Living donor"
meta$grp <- ifelse(grepl("Reference", meta$group), "Reference", "Disease")
clin$GSE104954 <- list(expr=ex, meta=meta)
## GSE66494 (Agilent GPL6480, linear -> log2, max-variance probe)
sm <- read_sm(file.path(D,"geo","GSE66494_series_matrix.txt.gz")); ex <- as.matrix(sm$tab[,-1]); rownames(ex) <- sm$tab[[1]]
con <- gzfile(file.path(D,"gpl","GPL6480.annot.gz")); L <- readLines(con, warn=FALSE); close(con)
b <- grep("^!platform_table_begin", L); e <- grep("^!platform_table_end", L); an <- fread(text=L[(b+1):(e-1)], data.table=FALSE, quote="", sep="\t", fill=TRUE)
sym <- setNames(an[["Gene symbol"]], an$ID)[rownames(ex)]; k <- !is.na(sym) & sym!=""; ex <- ex[k,]; sym <- sym[k]
ex[ex < 1] <- 1; ex <- log2(ex); o <- order(apply(ex,1,var), decreasing=TRUE); ex <- ex[o,]; sym <- sym[o]
ex <- ex[!duplicated(sym),]; rownames(ex) <- sym[!duplicated(sym)]
st <- parse_field(sm$hdr, "disease status")
meta <- data.frame(gsm=colnames(ex), status=st, grp=ifelse(grepl("normal", st, ignore.case=TRUE), "Reference", "Disease"),
   source=hfield(sm$hdr,"!Sample_source_name_ch1")[[1]], title=hfield(sm$hdr,"!Sample_title")[[1]], stringsAsFactors=FALSE)
clin$GSE66494 <- list(expr=ex, meta=meta, header_series=grep("^!Series_(summary|overall_design)", sm$hdr, value=TRUE),
                      sample_chars=do.call(cbind, hfield(sm$hdr, "!Sample_characteristics_ch1")))
for (n in names(clin)) cat(n, dim(clin[[n]]$expr), table(clin[[n]]$meta$grp), "\n")
saveRDS(clin, file.path(OUT, "clinical_bundle.rds"))
