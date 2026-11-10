# All outputs relative to the output working directory. Never source the original I/O scripts.
rev <- Sys.getenv("PROJECT_REV", file.path(path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta"),"aging_cell_revision_20261030"))
dat <- path.expand(Sys.getenv("KIDNEY_DATA","~/kidney_aging_data"))
rhs <- function(file,nm) {
  terms <- as.list(parse(file.path(rev,"scripts",file)))
  hit <- Filter(function(i)is.call(i)&&identical(i[[1]],as.name("<-"))&&identical(i[[2]],as.name(nm)),terms)
  stopifnot(length(hit)==1); hit[[1]][[3]]
}
lock <- rhs("04_locked_definitions.R","LOCK")
LOCK <- list(min_genes=5,equivalence_bound=.2,proxy=eval(lock[["proxy"]]))
injury <- eval(lock[["outcome"]][["PT_injury"]])
markers <- c(LOCK$proxy,list(PT_injury=injury))
mz <- eval(rhs("04_locked_definitions.R","mz"))
pspear <- eval(rhs("04_locked_definitions.R","pspear"))
load(file.path(rev,"kidneyAgeBench/data/kab_library.rda"))
lib <- split(kab_library,kab_library$sig_id)
signed <- names(lib)[vapply(lib,function(x)any(x$dir!=0),logical(1))]
stopifnot(length(lib)==48, length(signed)==15)
out <- function(x,fn) write.csv(x,fn,row.names=FALSE,na="")
archive <- read.csv(file.path(rev,"benchmark/results/REVISION_clinical_CKD_and_transplant.csv"))
cl <- readRDS(file.path(dat,"harmonized/clinical_bundle.rds"))
tx <- readRDS(file.path(dat,"harmonized/harmonized_bundle_corrected.rds"))$harm$GSE147451
pairs <- do.call(rbind,lapply(split(tx$kidney_level$meta,tx$kidney_level$meta$donor),function(dd) {
  data.frame(gsm=paste(dd$gsm,collapse="+"),kidneys=nrow(dd),titles=paste(dd$title,collapse=" / "),
    matching_key=dd$donor[1], label_sequence=paste(ifelse(is.na(dd$gfr_group),"missing",dd$gfr_group),collapse="/"),
    first_outcome=dd$gfr_group[1], discordant=length(unique(na.omit(dd$gfr_group)))>1, missing_BMI=any(is.na(dd$bmi)))
}))
meta_tx <- tx$meta
pairs <- pairs[match(meta_tx$gsm,pairs$gsm),]
stopifnot(identical(pairs$gsm,as.character(meta_tx$gsm)),identical(as.character(pairs$first_outcome),as.character(meta_tx$gfr_group)),
 nrow(pairs)==166,sum(pairs$kidneys)==177,sum(pairs$kidneys==2)==11,sum(pairs$discordant)==1)
out(pairs,"transplant_linkage_audit.csv")
# Verify the stored arithmetic-mean expression and predictor grouping using the underlying kidneys.
dev <- vapply(seq_len(nrow(pairs)),function(i){
  ix<-strsplit(pairs$gsm[i],"+",fixed=TRUE)[[1]]
  max(abs(rowMeans(tx$kidney_level$expr[,ix,drop=FALSE])-tx$expr[,i]),na.rm=TRUE)
},numeric(1));stopifnot(max(dev)<1e-10)
r1<-cl$GSE104954; keep<-r1$meta$diagnosis!="Tumor nephrectomy";r1$meta<-r1$meta[keep,];r1$expr<-r1$expr[,keep,drop=FALSE]
cohorts<-list(ERCB=r1,GSE175759=cl$GSE175759,GSE147451=tx)
overlaps<-arms<-results<-list()
join_genes <- function(x)paste(sort(unique(x)),collapse="|")
rank_test <- function(y,x,C=NULL) {
  if(!all(is.finite(x))||length(unique(x))<2) return(c(n=length(y),q=NA_real_,estimate=NA,lo=NA,hi=NA,p=NA))
  z<-pspear(y,x,C); c(n=z$n,q=z$q,estimate=z$rho,lo=z$lo,hi=z$hi,p=z$p)
}
ols_test <- function(y,disease,platform,adj) {
  dd<-data.frame(score=y,disease=disease)
  if(!is.null(platform))dd$platform<-factor(platform)
  if(!is.null(adj)) dd<-cbind(dd,adj)
  if(!all(is.finite(y)))return(c(n=nrow(dd),q=ncol(dd)-1,estimate=NA,lo=NA,hi=NA,p=NA))
  fit<-lm(score~., data=dd); co<-coef(summary(fit))["disease",]; ci<-confint(fit,"disease")[1,]
  c(n=nobs(fit),q=fit$rank-1,estimate=unname(co[1]),lo=unname(ci[1]),hi=unname(ci[2]),p=unname(co[4]))
}
for(cohort in names(cohorts)) {
  h<-cohorts[[cohort]];e<-h$expr;m<-h$meta
  stopifnot(identical(colnames(e),as.character(m$gsm)))
  xpx<-as.data.frame(lapply(LOCK$proxy,function(g)mz(e,g)))
  inj<-mz(e,injury);adj<-cbind(xpx,PT_injury=inj)
  ref<-if(cohort!="GSE147451")m$grp!="Disease" else rep(FALSE,nrow(m))
  # A membership listing for every set/proxy/cohort including zero overlaps.
  for(id in names(lib))for(pn in names(markers)) {
   sg<-lib[[id]]; genes<-intersect(sg$gene,markers[[pn]])
   overlaps[[length(overlaps)+1L]]<-data.frame(cohort,sig_id=id,class=sg$class[1],proxy=pn,
    panel_genes=join_genes(markers[[pn]]),all_shared=join_genes(genes),measured_shared=join_genes(intersect(genes,rownames(e))),
    n_shared=length(genes),n_measured_shared=length(intersect(genes,rownames(e))))
  }
  drops<-if(cohort=="GSE147451") unique(unlist(LOCK$proxy)) else unique(unlist(markers))
  for(id in signed) for(trim in c(FALSE,TRUE)) for(cutoff in 2:3) {
   sg<-lib[[id]];eligible<- if(trim) sg[!sg$gene %in% drops,] else sg
   up<-intersect(unique(eligible$gene[eligible$dir>0]),rownames(e)); dn<-intersect(unique(eligible$gene[eligible$dir<0]),rownames(e))
   u<-mz(e,up,cutoff);d<-mz(e,dn,cutoff)
   use_u<-any(is.finite(u));use_d<-any(is.finite(d))
   v<-if(use_u&&use_d)u-d else if(use_u)u else if(use_d)-d else rep(NA_real_,ncol(e))
   arm<-if(use_u&&use_d)"up-minus-down" else if(use_u)"up only" else if(use_d)"down only" else "unscorable"
   variant<-if(trim)"drop_proxy_union" else "intact"
   arms[[length(arms)+1]]<-data.frame(cohort,sig_id=id,variant,min_arm=cutoff,n_up=length(up),n_down=length(dn),
      use_up=use_u,use_down=use_d,arm,removed=join_genes(setdiff(sg$gene,eligible$gene)),measured_up=join_genes(up),measured_down=join_genes(dn))
   emit<-function(endpoint,population,model,stat) {
     results[[length(results)+1L]] <<- data.frame(cohort,sig_id=id,variant,min_arm=cutoff,arm,endpoint,population,model,as.list(stat),check.names=FALSE)
   }
   if(cohort=="GSE147451"){
    pops<-list(primary=!is.na(m$gfr_group),omit_discordant=!is.na(m$gfr_group)&!pairs$discordant,omit_all_bilateral=!is.na(m$gfr_group)&(pairs$kidneys==1))
    for(pop in names(pops)){
     k<-pops[[pop]];y<-as.numeric(m$gfr_group[k]=="poor");age<-data.frame(age=m$age[k]); full<-cbind(age,xpx[k,,drop=FALSE])
     for(mod in c("raw","age","age_composition"))emit("poor_function_group",pop,mod,rank_test(y,v[k],switch(mod,raw=NULL,age=age,age_composition=full)))
    }
   }else {
    vs<-(v-mean(v[ref]))/sd(v[ref]); plat<-if(cohort=="ERCB")m$platform else NULL
    emit("CKD_reference_difference","all","base",ols_test(vs,as.numeric(!ref),plat,NULL))
    emit("CKD_reference_difference","all","composition_injury",ols_test(vs,as.numeric(!ref),plat,adj))
    for(pop in c("all","CKD_only")){
     k<-if(pop=="all")rep(TRUE,length(ref)) else !ref
     for(pn in names(adj))emit(paste0("score_vs_",pn),pop,"raw",rank_test(adj[[pn]][k],v[k]))
     if(cohort=="GSE175759"){
      emit("eGFR",pop,"raw",rank_test(m$egfr[k],v[k]))
      emit("eGFR",pop,"composition_injury",rank_test(m$egfr[k],v[k],adj[k,,drop=FALSE]))
     }
    }
   }
  }
}
overlaps<-do.call(rbind,overlaps);out(overlaps,"overlap_all48_by_cohort.csv");arms<-do.call(rbind,arms);out(arms,"arm_availability.csv")
r<-do.call(rbind,results)
fam<-interaction(r$cohort,r$variant,r$min_arm,r$population,r$model,r$endpoint,drop=TRUE)
stopifnot(all(table(fam)==15))
r$q_BH15<-ave(r$p,fam,FUN=function(v)p.adjust(v,method="BH",n=15))
out(r,"all_clinical_effects.csv")
orig<-subset(r,variant=="intact" & min_arm==2)
out(subset(orig,cohort=="GSE147451"),"transplant_sensitivities.csv")
out(subset(orig,cohort!="GSE147451"),"CKD_population_comparisons.csv")
# Guarantee baseline effects agree with archived output, not just summaries.
verify<-function(rows,col,val="estimate"){
  rows<-rows[match(archive$sig_id,rows$sig_id),]
  stopifnot(identical(rows$sig_id,archive$sig_id), max(abs(rows[[val]]-archive[[col]])) <1e-9)
}
for(co in c("ERCB","GSE175759")) {
 for(mod in c("base","composition_injury")){
 x<-subset(orig,cohort==co&endpoint=="CKD_reference_difference"&model==mod)
 name<-paste0(if(co=="ERCB")"ERCB" else "GSE175","_beta_",if(mod=="base")"raw" else "adj")
 verify(x,name); verify(x,sub("beta","p",name),"p")
 }
}
for(mod in c("raw","composition_injury")){
 x<-subset(orig,cohort=="GSE175759"&endpoint=="eGFR"&model==mod&population=="all")
 name<-paste0("GSE175_eGFR_rho_",if(mod=="raw")"raw" else "adj")
 verify(x,name);verify(x,sub("rho","p",name),"p")
}
verify(subset(orig,cohort=="ERCB"&endpoint=="score_vs_PT_injury"&population=="all"),"ERCB_cor_PT_injury")
for(mod in c("raw","age","age_composition")){
 x<-subset(orig,cohort=="GSE147451"&model==mod&population=="primary")
 name<-paste0("TX_poorGFR_rho_",switch(mod,raw="raw",age="adj_age",age_composition="adj_age_comp"))
 verify(x,name);verify(x,sub("rho","p",name),"p")
}
cat("archived clinical estimates/nominal P reproduced. Rows",nrow(r)," / arm configurations",nrow(arms)," / overlap rows",nrow(overlaps),"\n")
writeLines(capture.output(sessionInfo()),"R_session.txt")
# Compact freeze-ready extract for figures and their source tables.
main <- subset(orig, (cohort=="ERCB" & endpoint %in% c("CKD_reference_difference","score_vs_PT_injury") & population=="all") |
                     (cohort=="GSE175759" & endpoint %in% c("CKD_reference_difference","eGFR") & population=="all") |
                     (cohort=="GSE147451" & population=="primary"))
main$panel <- ifelse(main$cohort=="ERCB",ifelse(main$endpoint=="CKD_reference_difference","4D","4B"),
              ifelse(main$cohort=="GSE175759",ifelse(main$endpoint=="CKD_reference_difference","4E","4C/4F"),
                     ifelse(main$model=="raw","5D",ifelse(main$model=="age","5E","5F"))))
main$unit <- ifelse(main$cohort=="GSE147451","inferred donor with first kidney recipient outcome","biopsy")
main$measure <- ifelse(main$endpoint=="CKD_reference_difference","linear coefficient, reference-group SD","Spearman rho or model-adjusted partial rank rho")
main$inputs <- ifelse(main$cohort=="GSE147451","harmonized_bundle_corrected.rds; original clinical model scripts/51_genomemed_upgrades.R","clinical_bundle.rds; original clinical model scripts/51_genomemed_upgrades.R")
out(main,"clinical_verified_figures.csv")
