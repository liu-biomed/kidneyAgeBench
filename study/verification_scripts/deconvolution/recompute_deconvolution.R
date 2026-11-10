# Local-only correction of the published-working deconvolution sensitivity calculation.
# Arguments: [read-only package root] [output directory]. No sampling/remote lookups.
args <- commandArgs(trailingOnly=TRUE)
input <- if (length(args)) args[1] else "/Users/zheranliu/Kidney_Manuscript_Packages/Genome_Medicine_终稿资料包_20261106"
out <- if (length(args) > 1) args[2] else "/Users/zheranliu/Kidney_Manuscript_Packages/submission_revision/deconvolution"
dir.create(out, recursive=TRUE, showWarnings=FALSE)
input <- normalizePath(input); out <- normalizePath(out)
stopifnot(!startsWith(out, paste0(input,"/")), input != out)
data_root <- file.path(input,"04_Source_data/c_Datasets")
code_root <- file.path(input,"04_Source_data/b_Original_code")
deconv_root <- file.path(data_root,"02_Derived/upstream/tables")
hpath <- file.path(data_root,"02_Derived/kidney_aging_data/harmonized/harmonized_bundle_corrected.rds")
original_hpath <- file.path(dirname(hpath),"harmonized_bundle.rds")
series <- file.path(data_root,"01_Study_inputs/geo/GSE147451_series_matrix.txt.gz")
library_file <- file.path(code_root,"kidneyAgeBench/data/kab_library.rda")
lock_file <- file.path(code_root,"revision/scripts/04_locked_definitions.R")
old_code <- file.path(code_root,"revision/scripts/51_genomemed_upgrades.R")
pair_code <- file.path(code_root,"revision/scripts/01c_corrected_set.R")
before_file <- file.path(data_root,"02_Derived/revision/benchmark/results/REVISION_deconv_method_concordance.csv")
fig_link <- file.path(data_root,"03_Figure_source_tables/F5_donor_linkage.csv")
methods <- list(MuSiC="music_class",Bisque5="bisque_subclass",Bisque_class="bisque_class",NNLS="nnls_class")
covars <- list(MuSiC=c("immune","stroma","endothelial"), Bisque5=c("PT","DCT","IMM","FIB","EC"),
               Bisque_class=c("immune","stroma","endothelial"), NNLS=c("immune","stroma","endothelial"))
files <- vapply(methods,function(x)file.path(deconv_root,paste0("deconv_proportions_",x,".csv")),character(1))
archives <- c(hpath,original_hpath,series,library_file,lock_file,pair_code,old_code,before_file,fig_link,files)
stopifnot(all(file.exists(archives)))
csv <- function(x,name) write.csv(x,file.path(out,name),row.names=FALSE,na="")
csv(data.frame(path=substring(archives,nchar(input)+2),sha256=vapply(archives,function(f)
  digest::digest(file=f,algo="sha256"),character(1))),"input_manifest.csv")

# Resolve the frozen proxy list and functions from the original file, without initializing
# unrelated msigdbr pathways or invoking any external API.
exprs <- parse(lock_file)
defs <- new.env(parent=globalenv())
for(expr in exprs) {
  if(!is.call(expr) || as.character(expr[[1]])[1] != "<-") next
  nm <- as.character(expr[[2]])[1]
  if(nm == "LOCK") {
    fields <- as.list(expr[[3]])[-1]
    defs$LOCK <- list(proxy=eval(fields$proxy),min_genes=eval(fields$min_genes),equivalence_bound=.20)
  }
  if(nm %in% c("mz","pspear")) eval(expr, envir=defs)
}
stopifnot(length(defs$LOCK$proxy)==5, is.function(defs$mz),is.function(defs$pspear))
lib <- new.env(); load(library_file,envir=lib)
sigs <- split(lib$kab_library,lib$kab_library$sig_id)
ids <- setdiff(names(sigs)[vapply(sigs,function(s)any(s$dir != 0),logical(1))], c("A11_KAS94_original","A12_KAS82_corrected"))
stopifnot(length(ids)==13)
h <- readRDS(hpath)$harm
stopifnot(length(h)==6, all(vapply(h,function(v) identical(colnames(v$expr),as.character(v$meta$gsm)) &&
       !anyDuplicated(v$meta$gsm),logical(1))))
# Reconstruct component membership from the original GEO fields and the original 01c
# adult-kidney set, independently of delimiters in the donor-level expression IDs.
con <- gzfile(series,"rt"); hdr <- readLines(con,n=300); close(con)
row <- function(line) sub('^"(.*)"$','\\1',strsplit(line,"\t",fixed=TRUE)[[1]][-1])
gsm <- row(grep("^!Sample_geo_accession\t",hdr,value=TRUE))
chars <- lapply(grep("^!Sample_characteristics_ch1\t",hdr,value=TRUE),row)
field <- function(name) vapply(seq_along(gsm),function(i) {
  x <- vapply(chars,"[",character(1),i); k <- paste0(name,": ")
  v <- x[startsWith(x,k)]; stopifnot(length(v)<=1)
  if(!length(v)) return(NA_character_)
  substring(v,nchar(k)+1)
},character(1))
ann <- data.frame(gsm=gsm,age=as.numeric(field("donor age")),bmi=suppressWarnings(as.numeric(field("donor bmi"))),
                  race=field("donor race"),kdpi=field("kdpi"),gfr=field("gfr group"),title=row(grep("^!Sample_title\t",hdr,value=TRUE)))
stopifnot(!anyDuplicated(ann$gsm))
gp <- h$GSE147451; kid <- gp$kidney_level; km <- kid$meta
original_kid <- readRDS(original_hpath)$harm$GSE147451
stopifnot(identical(original_kid$expr,kid$expr),identical(as.character(original_kid$meta$gsm),as.character(km$gsm)))
rm(original_kid); invisible(gc())
stopifnot(identical(colnames(kid$expr),as.character(km$gsm)),!anyDuplicated(km$gsm),
          all(km$gsm %in% ann$gsm),identical(rownames(kid$expr),rownames(gp$expr)))
ann <- ann[match(km$gsm,ann$gsm),]
ann$inferred_key <- paste(ann$age,round(ann$bmi,6),ann$race,ann$kdpi)
stopifnot(identical(ann$inferred_key,as.character(km$donor)),
          isTRUE(all.equal(ann$age,as.numeric(km$age),check.attributes=FALSE)),
          identical(ann$title,as.character(km$title)))
groups <- split(ann$gsm,ann$inferred_key)
links <- read.csv(fig_link,stringsAsFactors=FALSE)
stopifnot(length(groups)==166,nrow(ann)==177,max(lengths(groups))==2,sum(lengths(groups)==2)==11)
fraction_data <- lapply(files,function(f){
 x <- read.csv(f,stringsAsFactors=FALSE)
 stopifnot(!anyDuplicated(paste(x$cohort,x$gsm)))
 x
})
eligible <- names(h)[vapply(names(h),function(g)any(vapply(fraction_data,function(t)g %in% t$cohort,logical(1))),logical(1))]
stopifnot(length(eligible)==5,!('GSE66272' %in% eligible))
# Join the independently reconstructed kidney groups to corrected donor metadata using
# the original (age, BMI, race, KDPI) grouping key. Then check the complete expression
# vector, original source IDs and clinical linkage table for every donor.
stopifnot(setequal(rownames(gp$meta),names(groups)))
pair_evidence <- do.call(rbind,lapply(seq_len(nrow(gp$meta)),function(i){
 dm <- gp$meta[i,,drop=FALSE]; key <- rownames(dm); components <- groups[[key]]
 id <- as.character(dm$gsm)
 contained <- strsplit(id,"+",fixed=TRUE)[[1]]
 stopifnot(identical(components,contained),as.integer(dm$n_kidneys)==length(components),
           all(ann$age[match(components,ann$gsm)]==dm$age))
 delta <- max(abs(gp$expr[,i]-rowMeans(kid$expr[,components,drop=FALSE])))
 stopifnot(is.finite(delta),delta<1e-10)
 z <- links[match(id,links$gsm),]
 stopifnot(nrow(z)==1,!is.na(z$kidneys),as.integer(z$kidneys)==length(components))
 fractions <- vapply(names(fraction_data),function(m){
   t <- fraction_data[[m]];d <- t[t$cohort=="GSE147451",]
   j <- match(components,d$gsm)
   stopifnot(!anyNA(j), all(d$age[j]==dm$age))
   sum(complete.cases(d[j,covars[[m]],drop=FALSE]))
 },integer(1))
 a <- ann[match(components,ann$gsm),]
 data.frame(donor_record=id,inferred_key=key,kidney_profiles=length(components),geo_components=paste(components,collapse=";"),
            age=dm$age,bmi=a$bmi[1],race=a$race[1],kdpi=a$kdpi[1],
            max_abs_expression_mean_difference=delta,genes_checked=nrow(gp$expr),
            gfr_labels=paste(a$gfr,collapse=";"),MuSiC_inputs=fractions["MuSiC"],
            Bisque5_inputs=fractions["Bisque5"],Bisque_class_inputs=fractions["Bisque_class"],NNLS_inputs=fractions["NNLS"])
}))
stopifnot(length(unique(unlist(groups)))==177)
csv(pair_evidence,"GSE147451_donor_reconciliation.csv")

aggregate_method <- function(m,g,ids_g){
 t <- fraction_data[[m]]
 t <- t[t$cohort==g,,drop=FALSE]; vars <- covars[[m]]
 stopifnot(all(vars %in% names(t)))
 mat <- matrix(NA_real_,length(ids_g),length(vars),dimnames=list(ids_g,vars))
 for(i in seq_along(ids_g)){
  components <- if(g=="GSE147451")groups[[rownames(h[[g]]$meta)[i]]] else ids_g[i]
  j <- match(components,t$gsm)
  if(anyNA(j)) next    # no borrowing/partial pair aggregation
  d <- as.matrix(t[j,vars,drop=FALSE])
  if(any(!is.finite(d))) next
  stopifnot(all(abs(t$age[j]-h[[g]]$meta$age[i])<1e-8))
  mat[i,] <- colMeans(d)
 }
 as.data.frame(mat)
}
score <- function(e,s){
 up <- defs$mz(e,s$gene[s$dir>0],2); dn <- defs$mz(e,s$gene[s$dir<0],2)
 if(all(is.na(up))&&all(is.na(dn))) return(rep(NA_real_,ncol(e)))
 if(all(is.na(dn))) return(up)
 if(all(is.na(up))) return(-dn)
 up-dn
}
allmodels <- c("base","marker5",names(methods))
stopifnot(!anyDuplicated(allmodels))
cohort_data <- list(); linkage <- list(); proxies <- list()
for(g in names(h)){
 m <- h[[g]]$meta; ep <- h[[g]]$expr
 a <- if(all(is.na(m$sex)))data.frame(row.names=seq_len(nrow(m))) else data.frame(sex=factor(m$sex))
 px <- as.data.frame(lapply(defs$LOCK$proxy,function(gs)defs$mz(ep,gs)))
 agg <- setNames(lapply(names(methods),aggregate_method,g=g,ids_g=as.character(m$gsm)), names(methods))
 allcov <- c(list(base=a,marker5=cbind(a,px)),lapply(agg,function(x)cbind(a,x)))
 # For legacy replay: use the original semicolon-cut as written, without fixing identifiers.
 legacy <- lapply(names(methods),function(method) {
  t <- fraction_data[[method]]; t <- t[t$cohort==g,,drop=FALSE]
  v <- t[match(sub(";.*","",m$gsm),t$gsm),covars[[method]],drop=FALSE]
  rownames(v) <- seq_len(nrow(m)); cbind(a,v)
 });names(legacy) <- names(methods)
 for(method in names(methods)){
  pp <- agg[[method]];rownames(pp)<-NULL
  pp <- data.frame(cohort=g,donor_record=m$gsm,method=method,pp)
  proxies[[paste(g,method,sep="_")]] <- pp
 }
 scores <- sapply(ids,function(id)score(ep,sigs[[id]]))
 cohort_data[[g]] <- list(m=m,cov=allcov,legacy=c(list(base=a,marker5=cbind(a,px)),legacy),scores=scores)
 oldable <- do.call(cbind,lapply(legacy,complete.cases))
 shared <- Reduce('&', lapply(allcov,complete.cases))
 linkage[[g]] <- data.frame(cohort=g,donor_records=nrow(m),historical_deconv_match=sum(apply(oldable,1,all)),
              verified_all_methods=sum(shared),has_fraction_tables=g %in% eligible)
}
csv(do.call(rbind,linkage),"cohort_linkage.csv")
csv(data.table::rbindlist(proxies,fill=TRUE),"donor_covariates.csv")

# Same partial Spearman estimand as 04_locked_definitions: rank age and signature
# within the rows fitted, then project both ranks off the raw covariates. Use QR rank
# for effective q; report regression, correlation and Fisher degrees separately.
partial <- function(age,v,C,selected=rep(TRUE,length(age))){
 ok <- selected & is.finite(age)&is.finite(v)&complete.cases(C)
 n <- sum(ok)
 if(n<4)return(c(n=n,q=NA,reg_df=NA,cor_df=NA,fisher_info=NA,rho=NA,lo=NA,hi=NA))
 cd <- C[ok,,drop=FALSE]
 use <- vapply(cd,function(x)length(unique(x))>1,logical(1))
 X <- if(any(use))model.matrix(~.,droplevels(cd[,use,drop=FALSE])) else matrix(1,nrow=n,ncol=1)
 fit <- lm.fit(X,cbind(rank(age[ok]),rank(v[ok])))
 q <- fit$rank-1
 d <- n-q-3
 rr <- if(d>0)cor(fit$residuals[,1],fit$residuals[,2]) else NA_real_
 ci <- tanh(atanh(rr)+c(-1,1)*1.96/sqrt(d))
 c(n=n,q=q,reg_df=n-q-1,cor_df=n-q-2,fisher_info=d,rho=rr,lo=ci[1],hi=ci[2])
}
audits <- list(); use_ids <- list()
scopes <- c("historical_replay","historical_common_correct_weights","common_five","max_five","max_six")
for(scope in scopes){
 gs <- if(scope=="max_six")names(h) else eligible
 for(g in gs){
  dat <- cohort_data[[g]]; cov <- if(scope=="historical_replay")dat$legacy else dat$cov
  n0 <- nrow(dat$m)
  for(id in ids){
   v <- dat$scores[,id]
   sel <- switch(scope,common_five=Reduce('&',lapply(dat$cov,complete.cases)),
                         historical_common_correct_weights=Reduce('&',lapply(dat$legacy,complete.cases)),rep(TRUE,n0))
   for(method in allmodels){
    result <- partial(dat$m$age,v,cov[[method]],sel)
    ok <- sel & is.finite(dat$m$age) & is.finite(v) & complete.cases(cov[[method]])
    # Archived code skips an unevaluable signature/cohort before any pooling. Weighting
    # wrongly kept cohort N and nominal q even when the actual fitted rows were fewer.
    oldq <- as.integer(!all(is.na(dat$m$sex)))+switch(method,base=0,marker5=5,Bisque5=5,3)
    eff_n <- if(scope=="historical_replay")n0 else unname(result['n'])
    eff_q <- if(scope=="historical_replay")oldq else unname(result['q'])
    key <- paste(scope,g,id,method,sep="/")
    audits[[key]] <- data.frame(scope=scope,cohort=g,sig_id=id,method=method,
                       as.data.frame(as.list(result)),n_for_weight=eff_n,q_for_weight=eff_q,
                       fisher_variance=if(is.finite(result["rho"])) 1/(eff_n-eff_q-3) else NA_real_)
    use_ids[[key]] <- data.frame(scope=scope,cohort=g,sig_id=id,method=method,
      n_used=sum(ok),donor_set_sha256=digest::digest(paste(sort(dat$m$gsm[ok]),collapse=";"),algo="sha256",serialize=FALSE))
   }
  }
 }
}
a <- do.call(rbind,audits); donors <- do.call(rbind,use_ids)
csv(a,"cohort_effects.csv");csv(donors,"row_set_checks.csv")
same <- subset(donors,scope=="common_five")
stopifnot(all(vapply(split(same$donor_set_sha256,paste(same$cohort,same$sig_id)),function(z)length(unique(z))==1,logical(1))))
pool <- list();weights <- list(); warning_checks <- list()
for(k in unique(paste(a$scope,a$sig_id,a$method,sep="|"))){
 key <- strsplit(k,"|",fixed=TRUE)[[1]];v <- a[a$scope==key[1]&a$sig_id==key[2]&a$method==key[3]&is.finite(a$rho)&is.finite(a$fisher_variance),]
 stopifnot(nrow(v)>=2,all(abs(v$rho)<1),all(v$fisher_variance>0))
 z <- atanh(v$rho)
 warnings <- character()
 model <- withCallingHandlers(tryCatch(metafor::rma(yi=z,vi=v$fisher_variance,method="REML"),
                   error=function(e)metafor::rma(yi=z,vi=v$fisher_variance,method="DL")),
          warning=function(w){warnings <<- c(warnings,conditionMessage(w));invokeRestart("muffleWarning")})
 if(length(warnings)){
   # Check the boundary against a dense log grid and refined restricted-likelihood optimum.
   # Intercept-only REML profile, ignoring constants.
   reml_ll <- function(tau2){
     w <- 1/(v$fisher_variance+tau2)
     -.5*(sum(log(v$fisher_variance+tau2))+log(sum(w))+sum(w*(z-weighted.mean(z,w))^2))
   }
   upper <- max(1,var(z)*10)
   grid <- c(0,exp(seq(log(1e-12),log(upper),length.out=2001)))
   vals <- vapply(grid,reml_ll,numeric(1)); imax <- which.max(vals)
   opt <- optimize(reml_ll,c(grid[max(1,imax-1)],grid[min(length(grid),imax+1)]),maximum=TRUE,tol=1e-13)
   candidates <- c(0,model$tau2,opt$maximum,grid[imax],upper)
   ll <- vapply(candidates,reml_ll,numeric(1)); independent_tau2 <- candidates[which.max(ll)]
   gap <- max(ll)-reml_ll(model$tau2)
   warning_checks[[k]] <- data.frame(scope=key[1],sig_id=key[2],method=key[3],message=paste(warnings,collapse=" | "),
             fitted_tau2=model$tau2,grid_refined_tau2=independent_tau2,log_likelihood_gap=gap)
   stopifnot(gap<1e-6)
 }
 ww <- 1/(v$fisher_variance+model$tau2)
 stopifnot(abs(drop(model$b)-weighted.mean(z,ww))<1e-8)
 pool[[k]] <- data.frame(scope=key[1],sig_id=key[2],method=key[3],cohorts=nrow(v),
                          donors_actually_used=sum(v$n),donors_in_weight=sum(v$n_for_weight),pool_method=model$method,
                          tau2=model$tau2,rho=tanh(as.numeric(model$b)),lo=tanh(model$ci.lb),hi=tanh(model$ci.ub),p=model$pval)
 weights[[k]] <- data.frame(scope=key[1],sig_id=key[2],method=key[3],cohort=v$cohort,n=v$n,q=v$q,
             n_for_weight=v$n_for_weight,q_for_weight=v$q_for_weight,vi=v$fisher_variance,tau2=model$tau2,
             inverse_variance_weight=ww,weight_percent=100*ww/sum(ww))
}
p <- do.call(rbind,pool);csv(p,"pooled_effects.csv");csv(do.call(rbind,weights),"meta_weights.csv")
if(length(warning_checks)) csv(do.call(rbind,warning_checks),"REML_boundary_checks.csv")
orig <- read.csv(before_file)
file.copy(before_file,file.path(out,"archived_REVISION_deconv_method_concordance.csv"),overwrite=TRUE)
cn <- c(base="G6_unadjusted",marker5="G6_marker5",MuSiC="G6_KPMP_MuSiC",Bisque5="G6_KPMP_Bisque_5sub",
        Bisque_class="G6_KPMP_Bisque_class",NNLS="G6_KPMP_NNLS")
s10 <- data.frame(sig_id=ids)
for(m in allmodels)s10[[cn[[m]]]] <- vapply(ids,function(id)p$rho[p$scope=="common_five" & p$method==m & p$sig_id==id],numeric(1))
s10$n_common <- vapply(ids,function(id)p$donors_actually_used[p$scope=="common_five" & p$method=="marker5" & p$sig_id==id],numeric(1))
s10$cohorts_common <- vapply(ids,function(id)p$cohorts[p$scope=="common_five" & p$method=="marker5" & p$sig_id==id],numeric(1))
csv(s10,"S10_common_five_source.csv")
rhod <- function(scope,id,method)p$rho[p$scope==scope & p$sig_id==id & p$method==method]
change <- do.call(rbind,lapply(ids,function(id)do.call(rbind,lapply(allmodels,function(m) {
  old <- orig[match(id,orig$sig_id),cn[m]]
  data.frame(sig_id=id,method=m,archived_rho=old,historical_replay=rhod("historical_replay",id,m),
             historical_common_correct_weights=rhod("historical_common_correct_weights",id,m),
             corrected_common_five=rhod("common_five",id,m),corrected_max_five=rhod("max_five",id,m),
             corrected_max_six=rhod("max_six",id,m),change_common_minus_archive=rhod("common_five",id,m)-old)
}))))
csv(change,"before_after.csv")
error <- max(abs(change$archived_rho-change$historical_replay))
stopifnot(is.finite(error),error<1e-7)
concord <- do.call(rbind,lapply(scopes,function(scope)do.call(rbind,lapply(names(methods),function(m){
  cmp1 <- vapply(ids,function(id)rhod(scope,id,"marker5"),numeric(1));cmp2 <- vapply(ids,function(id)rhod(scope,id,m),numeric(1))
  data.frame(scope=scope,comparison=paste0("marker5_vs_",m),signatures=length(ids),spearman=cor(cmp1,cmp2,method="spearman"),
             identical_support_by_signature=all(vapply(ids,function(id){
               z<-donors[donors$scope==scope & donors$sig_id==id,,drop=FALSE];a1<-z[z$method=="marker5",];a2<-z[z$method==m,]
               identical(a1$donor_set_sha256,a2$donor_set_sha256)},logical(1))))
}))))
stopifnot(all(concord$identical_support_by_signature[concord$scope %in% c("common_five","historical_common_correct_weights","max_five")]),
          !any(concord$identical_support_by_signature[concord$scope %in% c("historical_replay","max_six")]))
csv(concord,"cross_signature_concordance.csv")
primary <- subset(p,scope=="common_five")
for(m in allmodels){
 v <- primary[primary$method==m,]
 cat(sprintf("%-14s donor N %s; studies %s; pooled rho %.4f to %.4f\n",m,paste(sort(unique(v$donors_actually_used)),collapse="/"),
             paste(sort(unique(v$cohorts)),collapse="/"),min(v$rho),max(v$rho)))
}
cat("Archived replay maximum absolute error:",error,"\n")
print(do.call(rbind,linkage),row.names=FALSE)
print(concord,row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
