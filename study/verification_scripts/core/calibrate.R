# Fixed competitive-null calibration. These simulations are never used as patient data.
source("code/cohorts_and_arms.R")
v<-subset(read.csv(file.path(revroot,"benchmark/results/REVISION_gene_level_composition_adjusted.csv")),config=="base_locked")
stopifnot(nrow(v)==48)
set.seed(20261106);B<-1000L;K<-2000L
Z<-readRDS(path.expand("~/kidney_aging_data/benchmark_intermediates.rds"))
domains<-lapply(Z$PZ,function(z)intersect(names(z)[is.finite(z)],names(kab_bins)))
excl<-function(id)if(grepl("^A1[12]_",id))"G6" else if(grepl("^A0[56]_Rowland",id))"TCGA" else "none"
G<-Reduce(union,domains)
res<-list();rhoRows<-list()
for(ds in names(D)){
 d<-D[[ds]];e<-d$e;g<-intersect(G,rownames(e));e<-e[g,,drop=FALSE]
 ranks<-t(apply(e,1,rank,na.last="keep"))
 X<-if(is.null(d$base))matrix(1,ncol(e),1) else model.matrix(~.,d$base)
 X<-cbind(X,rank(d$age))
 z<-t(qr.resid(qr(X),t(ranks)));norm<-sqrt(rowSums(z^2))
 z<-z[norm>1e-7 & is.finite(norm),,drop=FALSE];z<-z/sqrt(rowSums(z^2));res[[ds]]<-z
 for(id in names(sigs)){
    g<-intersect(sigs[[id]]$gene,rownames(z));k<-length(g)
    # Sum of off-diagonal inner products without forming k-by-k matrices.
    rho<-if(k>1)(sum(colSums(z[g,,drop=FALSE])^2)-k)/(k*(k-1)) else NA
    rhoRows[[length(rhoRows)+1]]<-data.frame(sig_id=id,cohort=ds,n=ncol(z),genes=k,residual_rho=rho)
 }
}
write.csv(do.call(rbind,rhoRows),"calibration/residual_correlations.csv",row.names=FALSE)
sets<-lapply(names(sigs),function(id){
 dom<-domains[[excl(id)]];g<-intersect(unique(sigs[[id]]$gene),dom)
 stopifnot(length(g)==v$n[match(id,v$sig_id)])
 pool<-split(dom,kab_bins[dom])
 draws<-vapply(g,function(x){p<-pool[[as.character(kab_bins[x])]];p[sample.int(length(p),K,replace=TRUE)]},character(K))
 sparseMatrix(i=c(rep(seq_len(K),length(g)),rep(K+1,length(g))),j=match(c(draws,g),G),x=1/length(g),dims=c(K+1,length(G)))
})
names(sets)<-names(sigs)
uses<-lapply(names(domains),function(key)names(D)[!vapply(D,function(d)d$grp==key,logical(1))]);names(uses)<-names(domains)
ps<-ph<-list()
for(sc in c("no_signal","equal_nonzero_effect_distribution"))
 ps[[sc]]<-ph[[sc]]<-matrix(NA_real_,length(sets),B,dimnames=list(names(sets),NULL))
for(start in seq(1,B,by=25)){
 idx<-start:min(start+24,B);n<-length(idx)
 noise<-lapply(res,function(z)z%*%matrix(rnorm(ncol(z)*n),ncol(z),n))
 # Pooled-score true effects have the same distribution for target and matched background.
 signal<-matrix(rnorm(length(G)*n,sd=1.5)*rbinom(length(G)*n,1,.25),length(G),n)
 for(key in names(domains)){
  pooled<-matrix(0,length(G),n);ng<-integer(length(G))
  for(ds in uses[[key]]){
   j<-match(rownames(noise[[ds]]),G);pooled[j,]<-pooled[j,,drop=FALSE]+noise[[ds]];ng[j]<-ng[j]+1L
  }
  stopifnot(all(ng[match(domains[[key]],G)]>=max(2,floor(length(uses[[key]])/2))))
  pooled<-pooled/sqrt(pmax(ng,1))
  for(sc in names(ps)){
   magnitude<-abs(if(sc=="no_signal")pooled else pooled+signal)
   for(id in names(sets)[vapply(names(sets),excl,character(1))==key]){
    vals<-as.matrix(sets[[id]]%*%magnitude);nul<-vals[seq_len(K),,drop=FALSE];obs<-vals[K+1,]
    ps[[sc]][id,idx]<-pnorm((obs-colMeans(nul))/(apply(nul,2,sd)*sqrt(v$camera_VIF[match(id,v$sig_id)])),lower.tail=FALSE)
    ph[[sc]][id,idx]<-(1+colSums(sweep(nul,2,obs,">=")))/(K+1)
   }
  }
 }
}
rates<-family<-list()
for(sc in names(ps))for(test in c("archived_VIF","empirical_matched")){
 p<-if(test=="archived_VIF")ps[[sc]] else ph[[sc]]
 rates[[length(rates)+1]]<-data.frame(scenario=sc,test=test,sig_id=rownames(p),reject=rowSums(p<.05),replicates=B,rate=rowMeans(p<.05))
 q<-apply(p,2,p.adjust,"BH")
 family[[length(family)+1]]<-data.frame(scenario=sc,test=test,B=B,mean_marginal_rate=mean(p<.05),min_rate=min(rowMeans(p<.05)),max_rate=max(rowMeans(p<.05)),fraction_with_any_BH=mean(colSums(q<.05)>0))
}
write.csv(do.call(rbind,rates),"calibration/false_positive_by_set.csv",row.names=FALSE)
write.csv(do.call(rbind,family),"calibration/false_positive_family.csv",row.names=FALSE)

