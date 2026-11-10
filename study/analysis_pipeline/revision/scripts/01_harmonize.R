# Rebuild the nine discovery cohorts and GTEx Kidney Cortex exactly as in the original
# execution notebooks (worker-0.ipynb, cells 17-79): same sample filters, probe-to-symbol
# rules (max-mean collapse), transformations and age/sex parsing.
suppressMessages({library(data.table); library(DESeq2); library(AnnotationDbi); library(org.Hs.eg.db); library(readxl)})
D <- path.expand("~/kidney_aging_data"); OUT <- file.path(D, "harmonized"); dir.create(OUT, showWarnings=FALSE)

read_sm <- function(f){            # series matrix -> list(header lines, expression table)
  con <- gzfile(f); L <- readLines(con, warn=FALSE); close(con)
  b <- grep("^!series_matrix_table_begin", L); e <- grep("^!series_matrix_table_end", L)
  hdr <- L[seq_len(b-1)]
  tab <- if (e-b > 2) fread(text=L[(b+1):(e-1)], data.table=FALSE) else NULL
  list(hdr=hdr, tab=tab)
}
sample_records <- function(hdr){   # mimic notebook full_meta(): gsm, title, source, characteristics
  sp <- function(key) { l <- hdr[startsWith(hdr, key)]; lapply(l, function(x) gsub('^"|"$', "", strsplit(x, "\t")[[1]][-1])) }
  gsm <- sp("!Sample_geo_accession")[[1]]; title <- sp("!Sample_title")[[1]]; src <- sp("!Sample_source_name_ch1")[[1]]
  ch <- sp("!Sample_characteristics_ch1")
  lapply(seq_along(gsm), function(i) c(gsm=gsm[i], title=title[i], source=src[i], sapply(ch, `[`, i)))
}
scan_token <- function(rec, key, kind){
  kr <- paste0("(?:", key, ")")
  for (v in rec) { if (is.na(v) || !grepl(kr, v, ignore.case=TRUE, perl=TRUE)) next
    if (kind=="age") { m <- regmatches(v, regexec(paste0(kr, "\\s*[:=]?\\s*(?:age\\s*)?(\\d{1,3}(?:\\.\\d+)?)"), v, ignore.case=TRUE, perl=TRUE))[[1]]
      if (length(m)) { a <- as.numeric(m[2]); if (a>=0 && a<=110) return(a) } }
    else { m <- regmatches(v, regexec(paste0(kr, "\\s*[:=]?\\s*([mf]|male|female)\\b"), v, ignore.case=TRUE, perl=TRUE))[[1]]
      if (length(m)) return(ifelse(startsWith(tolower(m[2]), "f"), "F", "M")) } }
  NA
}
build_meta <- function(hdr, filt, age_key="\\bage\\b"){
  recs <- sample_records(hdr); recs <- Filter(function(r) filt(tolower(paste(r, collapse=" "))), recs)
  m <- data.frame(gsm=sapply(recs, `[[`, "gsm"), title=sapply(recs, `[[`, "title"),
                  age=sapply(recs, scan_token, key=age_key, kind="age"),
                  sex=sapply(recs, scan_token, key="sex|gender", kind="sex"), stringsAsFactors=FALSE)
  m[!is.na(m$age), ]
}
collapse_maxmean <- function(mat, symbols){
  symbols <- trimws(sub("///.*$", "", symbols))
  keep <- !is.na(symbols) & symbols!="" & symbols!="---"; mat <- mat[keep,,drop=FALSE]; symbols <- symbols[keep]
  o <- order(rowMeans(mat, na.rm=TRUE), decreasing=TRUE); mat <- mat[o,,drop=FALSE]; symbols <- symbols[o]
  d <- duplicated(symbols); mat <- mat[!d,,drop=FALSE]; rownames(mat) <- symbols[!d]; mat
}
read_gpl <- function(f, symcol){
  con <- gzfile(f); L <- readLines(con, warn=FALSE); close(con)
  b <- grep("^!platform_table_begin", L); e <- grep("^!platform_table_end", L)
  if (!length(b)) { b <- max(grep("^[#!^]", L)); e <- length(L) + 1 }
  t <- fread(text=L[(b+1):(e-1)], data.table=FALSE, quote="", sep="\t", fill=TRUE)
  setNames(as.character(t[[symcol]]), as.character(t$ID))
}

filters <- list(
  GSE36895 = list(f=function(s) grepl("normal cortex", s, fixed=TRUE), gpl="GPL570.annot.gz", sym="Gene symbol"),
  GSE60807 = list(f=function(s) TRUE, gpl="GPL6480.annot.gz", sym="Gene symbol"),
  GSE21785 = list(f=function(s) TRUE, gpl="GPL96.annot.gz", sym="Gene symbol"),
  GSE147451= list(f=function(s) TRUE, gpl="GPL571.annot.gz", sym="Gene symbol", age="donor age", adult=TRUE),
  GSE66272 = list(f=function(s) grepl("adj.normal", s, fixed=TRUE) | grepl("adj. normal", s, fixed=TRUE), gpl="GPL570.annot.gz", sym="Gene symbol"),
  GSE232249= list(f=function(s) grepl("kidney", s, fixed=TRUE), gpl="GPL25759_table.txt", sym=NA, age="age"))
harm <- list()
for (g in names(filters)) {
  p <- filters[[g]]; sm <- read_sm(file.path(D, "geo", paste0(g, "_series_matrix.txt.gz")))
  meta <- build_meta(sm$hdr, p$f, if (is.null(p$age)) "\\bage\\b" else p$age)
  if (isTRUE(p$adult)) meta <- meta[meta$age >= 18, ]
  if (g %in% c("GSE147451","GSE232249")) meta$sex <- NA   # sex unavailable in source metadata (as originally analysed)
  ex <- as.matrix(sm$tab[,-1]); rownames(ex) <- as.character(sm$tab[[1]])
  symcol <- p$sym
  if (is.na(symcol)) { # submitter GPL; locate symbol column
    con <- gzfile(file.path(D,"gpl",p$gpl)); L <- readLines(con, warn=FALSE); close(con)
    b <- grep("^!platform_table_begin", L); if (!length(b)) b <- max(grep("^[#!^]", L)); h <- strsplit(L[b+1], "\t")[[1]]
    symcol <- h[grep("^GENE_?SYMBOL$|^GENE SYMBOL$", toupper(h))[1]] }
  map <- read_gpl(file.path(D, "gpl", p$gpl), symcol)
  gm <- collapse_maxmean(ex, map[rownames(ex)])
  common <- intersect(meta$gsm, colnames(gm)); gm <- gm[, common, drop=FALSE]; meta <- meta[match(common, meta$gsm), ]
  sy <- trimws(sub("///.*$", "", map[rownames(ex)])); ok <- !is.na(sy) & sy!="" & sy!="---"
  alt <- rowsum(ex[ok, common, drop=FALSE], sy[ok]) / as.vector(table(sy[ok])[rownames(rowsum(ex[ok, common, drop=FALSE], sy[ok]))])
  harm[[g]] <- list(expr=gm, meta=meta, type="array", expr_meanprobe=alt)
  cat(sprintf("%-10s genes=%5d n=%3d age[%g-%g] sexNA=%d\n", g, nrow(gm), ncol(gm), min(meta$age), max(meta$age), sum(is.na(meta$sex))))
}
ensg2sym <- function(x) mapIds(org.Hs.eg.db, keys=sub("\\..*$","",x), column="SYMBOL", keytype="ENSEMBL", multiVals="first")
# GSE93069: counts -> VST
sm <- read_sm(file.path(D,"geo","GSE93069_series_matrix.txt.gz")); meta <- build_meta(sm$hdr, function(s) TRUE)
cnt <- fread(file.path(D,"suppl","GSE93069_count_matrix.txt.gz"), data.table=FALSE); rownames(cnt) <- cnt[[1]]; cnt <- cnt[,-1]
ct <- intersect(meta$title, colnames(cnt)); cnt <- cnt[, ct]; meta <- meta[match(ct, meta$title), ]; colnames(cnt) <- meta$gsm
cm <- round(as.matrix(cnt)); cm[cm<0] <- 0; mode(cm) <- "integer"; cm <- cm[rowSums(cm>=10) >= 0.2*ncol(cm), ]
vm <- assay(vst(DESeqDataSetFromMatrix(cm, data.frame(row.names=colnames(cm), s=factor(rep("a",ncol(cm)))), ~1), blind=TRUE))
harm$GSE93069 <- list(expr=collapse_maxmean(vm, suppressMessages(ensg2sym(rownames(vm)))), meta=meta, type="rnaseq")
# GSE217427: cortex TPM -> log2(TPM+1)
sm <- read_sm(file.path(D,"geo","GSE217427_series_matrix.txt.gz")); meta <- build_meta(sm$hdr, function(s) grepl("cortex", s, fixed=TRUE))
gct <- fread(file.path(D,"suppl","GSE217427_expression_matrix.TPM.gct.gz"), skip=2, data.table=FALSE)
ct <- intersect(meta$title, colnames(gct)); meta <- meta[match(ct, meta$title), ]
gm <- collapse_maxmean(log2(as.matrix(gct[, ct]) + 1), gct$Description); colnames(gm) <- meta$gsm
harm$GSE217427 <- list(expr=gm, meta=meta, type="rnaseq")
# GSE167093: non-tumour columns of the pairs-normalised matrix, matched by patient number
sm <- read_sm(file.path(D,"geo","GSE167093_series_matrix.txt.gz")); meta <- build_meta(sm$hdr, function(s) grepl("type: normal", s, fixed=TRUE))
x <- as.data.frame(read_excel(file.path(D,"suppl","GSE167093_Pairs_normalization_matrix.xlsx"), sheet=1))
ex <- as.matrix(x[,-1]); rownames(ex) <- x[[1]]
nl <- grep("NL$", colnames(ex), value=TRUE); idx <- match(as.numeric(sub("patient\\s+(\\d+).*","\\1", meta$title)), as.numeric(gsub("2020_|NL","",nl)))
stopifnot(!anyNA(idx)); ex <- ex[, nl[idx]]; colnames(ex) <- meta$gsm
map <- read_gpl(file.path(D,"gpl","GPL10558_table.txt"), "Symbol"); sym <- map[rownames(ex)]
sym[grepl("phage|:low|:high|ERCC|^NA$", sym)] <- NA
okp <- !is.na(sym) & sym!=""
harm$GSE167093 <- list(expr=collapse_maxmean(ex, sym), meta=meta, type="array",
  expr_meanprobe=rowsum(ex[okp,], sym[okp]) / as.vector(table(sym[okp])[rownames(rowsum(ex[okp,], sym[okp]))]))
for (g in c("GSE93069","GSE217427","GSE167093")) { h <- harm[[g]]; cat(sprintf("%-10s genes=%5d n=%3d age[%g-%g] sexNA=%d\n", g, nrow(h$expr), ncol(h$expr), min(h$meta$age), max(h$meta$age), sum(is.na(h$meta$sex)))) }
harm <- harm[c("GSE36895","GSE60807","GSE21785","GSE147451","GSE66272","GSE232249","GSE93069","GSE217427","GSE167093")]
# GTEx v8 Kidney Cortex (RNASEQ freeze)
sa <- fread(file.path(D,"gtex","SampleAttributesDS.txt"), data.table=FALSE, quote="")
sj <- fread(file.path(D,"gtex","SubjectPhenotypesDS.txt"), data.table=FALSE)
kc <- sa$SAMPID[sa$SMTSD=="Kidney - Cortex" & sa$SMAFRZE=="RNASEQ"]
g <- fread(file.path(D,"gtex","gene_tpm_v8_kidney_cortex.gct.gz"), skip=2, data.table=FALSE)
kc <- intersect(kc, colnames(g))
gm <- collapse_maxmean(log2(as.matrix(g[, kc]) + 1), g$Description)
subj <- sub("^(GTEX-[^-]+).*","\\1", kc); mt <- sj[match(subj, sj$SUBJID), ]
mid <- c("20-29"=25,"30-39"=35,"40-49"=45,"50-59"=55,"60-69"=65,"70-79"=75)
gtex <- list(expr=gm, meta=data.frame(gsm=kc, SUBJID=subj, AGE=mt$AGE, age=unname(mid[mt$AGE]), sex=ifelse(mt$SEX==1,"M","F"), DTHHRDY=mt$DTHHRDY), type="rnaseq")
# GTEx sample-level technical covariates for sensitivity analyses
gtex$meta <- cbind(gtex$meta, sa[match(kc, sa$SAMPID), c("SMRIN","SMTSISCH","SMCENTER","SMNABTCH")])
cat(sprintf("GTEx genes=%d n=%d age[%g-%g]\n", nrow(gm), ncol(gm), min(gtex$meta$age), max(gtex$meta$age)))
saveRDS(list(harm=harm, gtex=gtex), file.path(OUT, "harmonized_bundle.rds"))
cat("saved\n")
