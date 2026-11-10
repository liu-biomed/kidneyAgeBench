# Summary across discovery and independent validation layers (descriptive synthesis; not a pre-specified test)
suppressMessages({library(data.table); library(metafor)})
REV <- path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030")
source(file.path(REV, "scripts", "04_locked_definitions.R")); O <- file.path(REV, "results", "validation"); V <- path.expand("~/kidney_aging_data/validation/cptac")
mp <- fread(cmd=paste("gzcat", shQuote(list.files(V, pattern="bcm-ccrcc-mapping", full.names=TRUE))), data.table=FALSE)
pr <- fread(cmd=paste("gzcat", shQuote(list.files(V, pattern="proteomics.*Normal", full.names=TRUE))), data.table=FALSE)
psym <- unique(mp$gene_name[match(sub("\\..*","",pr[[1]]), sub("\\..*","",mp$gene))])
cov <- data.frame(proxy=names(LOCK$proxy), markers=sapply(LOCK$proxy, length), measured_in_proteome=sapply(LOCK$proxy, function(g) sum(g %in% psym)),
                  measured_genes=sapply(LOCK$proxy, function(g) paste(intersect(g, psym), collapse=";")))
print(cov, row.names=FALSE); write.csv(cov, file.path(O, "CPTAC_proteome_proxy_coverage.csv"), row.names=FALSE)
disc <- read.csv(file.path(REV, "prereg", "discovery_estimates_pooled.csv")); disc$dataset <- "Discovery (6 GEO cohorts, 441)"
tc <- read.csv(file.path(O, "TCGA_locked_results.csv")); tc <- tc[tc$score=="z",]; tc$dataset <- "TCGA normal (129)"
cp <- read.csv(file.path(O, "CPTAC_locked_results.csv")); cp$dataset <- ifelse(grepl("proteome", cp$dataset), "CPTAC NAT protein (80)", "CPTAC NAT RNA (75)")
cols <- c("dataset","outcome","model","rho","lo","hi","p")
S <- rbind(disc[, cols], tc[, cols], cp[, cols]); S <- S[is.finite(S$rho),]
# random-effects pooling of the two independent RNA validation cohorts
val <- rbind(tc[, c("outcome","model","rho","n","q")], cp[cp$dataset=="CPTAC NAT RNA (75)", c("outcome","model","rho","n","q")])
vp <- do.call(rbind, lapply(split(val, paste(val$outcome, val$model)), function(d){ d <- d[is.finite(d$rho),]; f <- rma(yi=atanh(d$rho), sei=1/sqrt(d$n-d$q-3), method="REML")
  data.frame(dataset="Validation RNA pooled (204)", outcome=d$outcome[1], model=d$model[1], rho=tanh(f$b), lo=tanh(f$ci.lb), hi=tanh(f$ci.ub), p=f$pval) }))
S <- rbind(S, vp); w <- reshape(S[, c("dataset","outcome","model","rho","p")], idvar=c("outcome","model"), timevar="dataset", direction="wide")
write.csv(S, file.path(O, "VALIDATION_SUMMARY_long.csv"), row.names=FALSE)
fmt <- function(r, p) sprintf("%+.2f%s", r, ifelse(is.na(p), "", ifelse(p<.001,"***",ifelse(p<.01,"**",ifelse(p<.05,"*","")))))
tab <- data.frame(outcome=w$outcome, model=w$model)
for (d in c("Discovery (6 GEO cohorts, 441)","TCGA normal (129)","CPTAC NAT RNA (75)","Validation RNA pooled (204)","CPTAC NAT protein (80)")) tab[[d]] <- fmt(w[[paste0("rho.",d)]], w[[paste0("p.",d)]])
tab <- tab[order(tab$outcome, tab$model),]; print(tab, row.names=FALSE); write.csv(tab, file.path(O, "VALIDATION_SUMMARY_table.csv"), row.names=FALSE)
