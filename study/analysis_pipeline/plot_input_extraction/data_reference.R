# Run after data_clinical.R. Recreate plotting tables directly from archived inputs.
root <- ".."; D <- path.expand(Sys.getenv("KIDNEY_DATA","~/kidney_aging_data"))
br <- function(fn)read.csv(file.path(root,"benchmark/results",fn))
out <- function(x,f)write.csv(x,file.path("source_data",paste0(f,".csv")),row.names=FALSE,na="")
assign_rhs <- function(file,name){
 expr <- as.list(parse(file.path(root,"scripts",file)))
 hit <- Filter(function(x)is.call(x) && identical(x[[1]],as.name("<-")) && identical(x[[2]],as.name(name)),expr)
 stopifnot(length(hit)==1);hit[[1]][[3]]
}
ages <- readRDS(file.path(D,"figure_data_v2.rds"))$ages
out(as.data.frame(table(ages$grp)),"F1_reference_n")
folds <- pool <- list()
for(v in c("original","corrected")){
 nm <- if(v=="original")"KAS-94" else "KAS-82"
 z <- readRDS(file.path(D,"harmonized",paste0("loco_",v,".rds")))
 a <- z$folds[,c("heldout","n","rho","sig_rule","sig_size")]
 a$in_sample <- unname(z$insample[paste0(a$heldout,".rho")])
 stopifnot(all(is.finite(a$rho)),all(is.finite(a$in_sample)))
 folds[[v]] <- cbind(signature=nm,a)
 f <- read.csv(file.path(root,"results/tables",v,"R1_nested_LOCO_pooled.csv"))
 pool[[v]] <- cbind(signature=nm,f[!grepl("fixed top",f$analysis),])
}
out(do.call(rbind,folds),"F1_folds");out(do.call(rbind,pool),"F1_pooled_KAS")
gt <- readRDS(file.path(D,"harmonized/gtex_scores_original.rds"))
stopifnot(identical(rownames(gt$score),as.character(gt$meta$gsm)))
out(data.frame(sample=gt$meta$gsm,age=gt$meta$age,age_bin=gt$meta$AGE,score=gt$score$score),"F1_GTEx")
gs <- read.csv(file.path(root,"results/tables/original/R1_GTEx_score_level.csv"))
out(subset(gs,axis=="score"),"F1_GTEx_estimate")
null <- read.csv(file.path(root,"benchmark/figures/SourceData_Figure1c_null.csv"))
stopifnot(nrow(null)==10000)
out(null,"F1_sex_perm")
observed <- read.csv(file.path(root,"results/tables/original/R1_GTEx_concordance_nulls.csv"))
out(observed[,c("n_genes","concordant","p_perm_concordance")],"F1_concordance_observed")
pools <- br("REVISION_gene_pools_four_configs.csv")
out(subset(pools,config=="primary_locked"),"F1_background")
f <- br("BENCHMARK_FINAL_all_sets.csv")
k <- read.csv("source_data/sets_key.csv",check.names=FALSE)
m <- match(k$sig_id,f$sig_id);stopifnot(length(m)==48,!anyNA(m))
out(cbind(k,f[m,c("enrichment","emp_p","fdr_all")]),"F2_48_sets")
v <- subset(br("REVISION_gene_level_composition_adjusted.csv"),config=="base_locked")
out(v[,c("class","sig_id","fdr","fdr_vif","camera_VIF")],"F2_VIF")
out(subset(br("REVISION_power_by_setsize.csv"),config=="primary_locked" & true_ER %in% c(1,1.2,1.3,1.5)),"F2_power")
alt <- subset(br("REVISION_gene_level_four_configs.csv"),config=="noCPTAC_evalbins")
q <- merge(f[,c("sig_id","class","enrichment")],alt[,c("sig_id","enrichment")],by="sig_id",suffixes=c("_primary","_noCPTAC"),all.x=TRUE)
stopifnot(sum(is.finite(q$enrichment_noCPTAC))==46);out(q,"F2_without_CPTAC")
MS <- br("M_signed_by_dataset.csv")
coh <- unique(MS$dataset); ids <- read.csv("source_data/signed_key.csv")$sig_id
raw <- subset(MS,model=="base")
cells <- merge(expand.grid(sig_id=ids,dataset=coh,stringsAsFactors=FALSE),raw,by=c("sig_id","dataset"),all.x=TRUE)
leak <- eval(assign_rhs("21_benchmark.R","leak"))
cells$use <- mapply(function(id,ds,r){
 class <- if(grepl("^G6_",ds))"G6" else ds
 if(class %in% leak(id))"Derivation" else if(!is.finite(r))"Unavailable" else "Observed"
},cells$sig_id,cells$dataset,cells$rho)
stopifnot(identical(as.integer(table(factor(cells$use,levels=c("Observed","Derivation","Unavailable")))),c(120L,14L,1L)))
out(cells,"F2_signed_cohorts")
# Exact two-stage pooling from the archived study, including the original
# signature-specific exclusions and the choice of REML versus fixed effects.
suppressPackageStartupMessages(library(metafor))
poolfun <- eval(assign_rhs("21_benchmark.R","poolfun"))
MSp <- eval(assign_rhs("21_benchmark.R","MSp"))
names(MSp)[names(MSp)=="pooled_rho"] <- "rho"
names(MSp)[names(MSp)=="pooled_lo"] <- "lo";names(MSp)[names(MSp)=="pooled_hi"] <- "hi"
for(model in c("base","composition")){
 z <- MSp[MSp$model==model,];g <- f[match(z$sig_id,f$sig_id),]
 expected <- if(model=="base") g$pooled_rho else g$pooled_rho_comp
 stopifnot(max(abs(z$rho-expected))<1e-10)
 if(model=="base")stopifnot(max(abs(z$lo-g$pooled_lo))<1e-10,max(abs(z$hi-g$pooled_hi))<1e-10)
}
out(MSp,"F2_F3_signed_intervals")
px <- read.csv(file.path(root,"benchmark/figures/v3/SourceData_Figure6d.csv"))
out(px,"F3_age_proxies")
out(subset(br("REVISION_gene_pools_comp_and_ckd.csv"),analysis %in% c("Age (unadjusted, locked)","Age (composition-adjusted, locked)")),"F3_pools")
# Largest retained GEO cohort: actual ages and marker scores from donor-level profiles.
cb <- readRDS(file.path(D,"harmonized/harmonized_bundle_corrected.rds"))
cohort <- names(which.max(vapply(cb$harm,function(g)ncol(g$expr),integer(1))))
stopifnot(cohort=="GSE147451")
h <- cb$harm[[cohort]];stopifnot(identical(as.character(h$meta$gsm),colnames(h$expr)))
spec <- assign_rhs("04_locked_definitions.R","LOCK")
LOCK <- list(min_genes=5,proxy=eval(spec[["proxy"]]),equivalence_bound=.20)
mz <- eval(assign_rhs("04_locked_definitions.R","mz"));score <- eval(assign_rhs("21_benchmark.R","score"))
load(file.path(root,"kidneyAgeBench/data/kab_library.rda")); sel <- ids[c(1,6,14)]
s <- data.frame(cohort,sample=h$meta$gsm,age=h$meta$age,as.data.frame(lapply(LOCK$proxy,function(g)mz(h$expr,g))))
for(id in sel)s[[id]] <- score(h$expr,kab_library[kab_library$sig_id==id,])
out(s,"F3_donor_samples")


# Descriptive correlations use the same paired donor scores as the 3E scatter.
# No cutoffs, covariate adjustment or additional hypothesis tests.
f3e <- do.call(rbind,lapply(sel,function(id){
  ok <- is.finite(s$immune) & is.finite(s[[id]])
  stopifnot(sum(ok)==166,!anyDuplicated(s$sample[ok]))
  data.frame(cohort=cohort,sig_id=id,proxy="immune",n=sum(ok),
             rho=cor(s$immune[ok],s[[id]][ok],method="spearman"))
}))
out(f3e,"F3E_descriptive_rho")
