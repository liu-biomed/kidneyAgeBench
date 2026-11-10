# Corrected discovery set: only cohorts meeting the study's own eligibility rule
# (>=12 independent adult normal/reference donors with genome-wide total-RNA expression).
#  - GSE232249 excluded: 6 kidneys, each split into m6A-IP and supernatant fractions (two-colour arrays), not expression profiles
#  - GSE21785  excluded: 6 living donors, each microdissected into tubulus and glomerulus
#  - GSE217427 excluded: 22 cortex samples from 15 donors (left/right kidneys), 8 donors with AKI; 7 non-AKI donors
#  - GSE147451 collapsed to donor level: 177 kidneys from 166 deceased donors (11 donors contributed both kidneys)
D <- path.expand("~/kidney_aging_data"); b <- readRDS(file.path(D, "harmonized", "harmonized_bundle.rds")); h <- b$harm
con <- gzfile(file.path(D, "geo", "GSE147451_series_matrix.txt.gz")); L <- readLines(con, n=300); close(con)
gsm <- gsub('"', "", strsplit(grep("^!Sample_geo_accession", L, value=TRUE), "\t")[[1]][-1])
ch <- grep("^!Sample_characteristics_ch1", L, value=TRUE); C <- t(sapply(ch, function(x) gsub('"', "", strsplit(x, "\t")[[1]][-1])))
rec <- apply(C, 2, paste, collapse=" | ")
get <- function(k) sapply(rec, function(r){ m <- regmatches(r, regexec(paste0(k, ": ([^|]*)"), r))[[1]]; if (length(m)) trimws(m[2]) else NA })
ann <- data.frame(gsm=gsm, bmi=suppressWarnings(as.numeric(get("donor bmi"))), race=get("donor race"), kdpi=get("kdpi"),
                  gfr_group=get("gfr group"), batch=suppressWarnings(as.numeric(get("batch effect"))), stringsAsFactors=FALSE)
g <- h$GSE147451; m <- merge(g$meta, ann, by="gsm", sort=FALSE); m <- m[match(colnames(g$expr), m$gsm), ]
m$donor <- paste(m$age, round(m$bmi, 6), m$race, m$kdpi)
ex <- sapply(split(seq_len(nrow(m)), m$donor), function(i) rowMeans(g$expr[, i, drop=FALSE]))
mp <- if (!is.null(g$expr_meanprobe)) sapply(split(seq_len(nrow(m)), m$donor), function(i) rowMeans(g$expr_meanprobe[, i, drop=FALSE])) else NULL
md <- do.call(rbind, lapply(split(m, m$donor), function(d) data.frame(gsm=paste(d$gsm, collapse="+"), title=d$title[1], age=d$age[1], sex=NA,
        bmi=d$bmi[1], race=d$race[1], gfr_group=d$gfr_group[1], batch=mean(d$batch, na.rm=TRUE), n_kidneys=nrow(d))))
md <- md[match(colnames(ex), rownames(md)), ]; colnames(ex) <- md$gsm; if (!is.null(mp)) colnames(mp) <- md$gsm
h$GSE147451 <- list(expr=ex, meta=md, type="array", expr_meanprobe=mp, kidney_level=list(expr=g$expr, meta=m))
keep <- c("GSE36895","GSE60807","GSE147451","GSE66272","GSE93069","GSE167093")
cat("GSE147451 donors:", ncol(ex), "| kidneys:", nrow(m), "\n")
for (k in keep) cat(k, ncol(h[[k]]$expr), "\n"); cat("total donors:", sum(sapply(h[keep], function(x) ncol(x$expr))), "\n")
saveRDS(list(harm=h[keep], gtex=b$gtex, excluded=b$harm[c("GSE232249","GSE21785","GSE217427")]), file.path(D, "harmonized", "harmonized_bundle_corrected.rds"))
