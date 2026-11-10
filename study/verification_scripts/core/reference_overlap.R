# Sensitivity only; keep the five proxies and injury definition fixed.
source("code/cohorts_and_arms.R")
library(metafor)
defs<-parse(orig("04_locked_definitions.R"))
# Read only literal marker vectors; no external lookup or global script execution.
for(ex in defs)if(is.call(ex)&&identical(ex[[1]],as.name("<-"))&&identical(ex[[2]],as.name("LOCK"))){
  proxy<-eval(ex[[3]][["proxy"]]); injury<-eval(ex[[3]][["outcome"]][["PT_injury"]])
  LOCK<-list(equivalence_bound=eval(ex[[3]][["equivalence_bound"]]))
}
assignments(orig("04_locked_definitions.R"),"pspear")
assignments(orig("21_benchmark.R"),"poolfun")
stopifnot(length(proxy)==5,length(injury)==10)
D$GTEx<-list(e=cb$gtex$expr,age=cb$gtex$meta$age,base=data.frame(sex=factor(cb$gtex$meta$sex)),grp="GTEx")
out<-list()
for(ds in names(D)){
 d<-D[[ds]];px<-as.data.frame(lapply(proxy,function(g)mz(d$e,g,5)))
 for(id in names(signed))for(rule in c("original","remove_shared")){
   s<-signed[[id]];if(rule=="remove_shared")s<-s[!(s$gene%in%c(unlist(proxy),injury)),]
   u<-mz(d$e,s$gene[s$dir>0],3);v<-mz(d$e,s$gene[s$dir<0],3)
   score<-if(all(is.na(u)))-v else if(all(is.na(v)))u else u-v
   for(model in c("base","composition")){
     C<-if(model=="base") d$base else if(is.null(d$base))px else cbind(d$base,px)
     est<-if(any(is.finite(score)))pspear(d$age,score,C) else data.frame(n=0,q=NA,rho=NA,lo=NA,hi=NA,p=NA,p_tost=NA)
     out[[length(out)+1]]<-cbind(data.frame(ds,id,rule,model,G6=startsWith(ds,"G6_"),immune_score_rho=if(any(is.finite(score)))cor(px$immune,score,method="spearman") else NA),est)
   }
 }
}
df<-do.call(rbind,out)
write.csv(df,"calibration/reference_overlap_by_cohort.csv",row.names=FALSE)
pool<-list()
for(id in names(signed))for(rule in unique(df$rule))for(model in unique(df$model)){
 # explicit indexing (avoid data-mask self-comparisons)
 a<-df[df$id==id & df$rule==rule & df$model==model & is.finite(df$rho),]
 exc<-if(grepl("^A1[12]_",id))"G6" else if(grepl("^A0[56]_Rowland",id))"TCGA" else ""
 parts<-list();six<-a[a$G6,]
 if(exc!="G6" && nrow(six)){
  t<-poolfun(six$rho,six$n,six$q);parts[[1]]<-c(z=atanh(t[1]),se=(atanh(t[3])-atanh(t[2]))/3.92)
 }
 for(k in setdiff(c("TCGA","CPTAC"),exc)){
  r<-a[a$ds==k,];if(nrow(r))parts[[length(parts)+1]]<-c(z=atanh(r$rho),se=1/sqrt(r$n-r$q-3))
 }
 if(!length(parts))next
 P<-do.call(rbind,parts);m<-rma(yi=P[,"z"],sei=P[,"se"],method=if(nrow(P)>2)"REML" else "FE")
 pool[[length(pool)+1]]<-data.frame(sig_id=id,rule,model,rho=tanh(coef(m)),lo=tanh(m$ci.lb),hi=tanh(m$ci.ub),p=m$pval)
}
pool<-do.call(rbind,pool);write.csv(pool,"calibration/reference_overlap_pooled.csv",row.names=FALSE)
old<-read.csv(file.path(revroot,"genome_medicine_polished/source_data/F2_F3_signed_intervals.csv"))
cmp<-merge(subset(pool,rule=="original"),old,by=c("sig_id","model"))
stopifnot(nrow(cmp)==30,max(abs(cmp$rho.x-cmp$rho.y))<1e-6)

