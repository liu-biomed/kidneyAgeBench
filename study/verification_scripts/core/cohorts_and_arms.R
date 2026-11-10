# Materialize the identical reference datasets by reading (not sourcing) script 51's loader.
suppressPackageStartupMessages({library(data.table); library(Matrix)})
revroot <- Sys.getenv("KIDNEY_REVISION_ROOT",unset="/Users/zheranliu/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030")
orig <- function(file) file.path(revroot,"scripts",file)
assignments <- function(file,targets) {
 for(ex in parse(file))if(is.call(ex)&&identical(ex[[1]],as.name("<-"))&&is.symbol(ex[[2]])&&as.character(ex[[2]])%in%targets)eval(ex,.GlobalEnv)
}
assignments(orig("10_cptac_helpers.R"),c("cptac_rd","cptac_sym","cptac_clin"))
assignments(orig("04_locked_definitions.R"),"mz")
ss <- readLines(orig("51_genomemed_upgrades.R"))
eval(parse(text=ss[which(ss=="D <- list()"):which(startsWith(ss,"rm(x, e, R, tr)"))]),.GlobalEnv)
load(file.path(revroot,"kidneyAgeBench/data/kab_library.rda"))
load(file.path(revroot,"kidneyAgeBench/data/kab_bins.rda"))
sigs <- split(kab_library,kab_library$sig_id)
signed <- sigs[vapply(sigs,function(s)any(s$dir!=0),logical(1))]
mats <- lapply(D,function(x)x$e)
mats$GTEx <- cb$gtex$expr
cl <- readRDS(path.expand("~/kidney_aging_data/harmonized/clinical_bundle.rds"))
mats$ERCB <- cl$GSE104954$expr[,cl$GSE104954$meta$diagnosis!="Tumor nephrectomy"]
mats$GSE175759 <- cl$GSE175759$expr
armcheck <- list()
for(ds in names(mats))for(id in names(signed)){
 s <- signed[[id]];e <- mats[[ds]]
 arm <- function(m) {
    u<-mz(e,s$gene[s$dir>0],m);d<-mz(e,s$gene[s$dir<0],m)
    if(all(is.na(u))) -d else if(all(is.na(d))) u else u-d
 }
 x<-arm(2);y<-arm(3);ok<-is.finite(x)&is.finite(y)
 armcheck[[length(armcheck)+1]]<-data.frame(dataset=ds,sig_id=id,up=sum(s$gene[s$dir>0]%in%rownames(e)),down=sum(s$gene[s$dir<0]%in%rownames(e)),n_profiles=ncol(e),usable_min2=sum(is.finite(x)),usable_min3=sum(is.finite(y)),max_abs_difference=if(any(ok))max(abs(x[ok]-y[ok])) else NA)
}
write.csv(do.call(rbind,armcheck),"calibration/thresholds_2_vs_3.csv",row.names=FALSE)
rm(mats,cl)

