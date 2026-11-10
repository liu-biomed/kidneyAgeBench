# Execute with Rscript --vanilla from any location. This writes only to submission_revision/raag/.
suppressPackageStartupMessages({library(readxl);library(xml2);library(metafor)})
project<-file.path(path.expand("~"),"Kidney_Manuscript_Packages")
out<-file.path(project,"submission_revision","raag")
package<-Sys.glob(file.path(project,"Genome_Medicine_*"))[1]
source_root<-file.path(package,"04_Source_data")
derived<-file.path(source_root,"c_Datasets/02_Derived")
data_root<-file.path(derived,"kidney_aging_data")
inputs<-file.path(source_root,"c_Datasets/01_Study_inputs")
scripts<-file.path(source_root,"b_Original_code/revision/scripts")
br<-function(n)read.csv(file.path(derived,"revision/benchmark/results",paste0(n,".csv")),check.names=FALSE)
figdata<-function(n)read.csv(file.path(source_root,"c_Datasets/03_Figure_source_tables",paste0(n,".csv")),check.names=FALSE)
emit<-function(x,fn)write.csv(x,file.path(out,fn),row.names=FALSE,na="")
s1<-file.path(package,"02_Supplementary/TableS1.xlsx")
s1genes<-as.data.frame(read_excel(s1,"S1b_gene_lists"))
lib<-read.csv(file.path(derived,"revision/benchmark/signatures/signature_library_long.csv"))
load(file.path(source_root,"b_Original_code/kidneyAgeBench/data/kab_library.rda"))
id<-"A09_Perco2019_RAAG_CKD";allid<-"A08_Perco2019_RAAG634"
# Primary Table 3, which supplies CKD fold changes and is sorted by them.
xm<-read_xml(file.path(out,"PMC6611921.xml"))
tr<-xml_find_all(xm,"//table-wrap[label='Table 3']//tbody/tr")
t3<-as.data.frame(t(vapply(tr,function(row)xml_text(xml_find_all(row,"./td")),character(8))),stringsAsFactors=FALSE)
names(t3)<-c("table3_gene","description","reported_aging_source","glomerular_protein","tubular_protein","ckd_fc","ckd_p","ckd_fdr_band")
t3$ckd_fc<-as.numeric(gsub("−","-",t3$ckd_fc));t3$ckd_p<-as.numeric(t3$ckd_p)
g<-s1genes[s1genes$sig_id==id,];stopifnot(nrow(g)==31, identical(t3$table3_gene,g$gene_orig),all(g$gene_orig==g$gene),all(sign(t3$ckd_fc)==g$dir))
for(other in list(lib,kab_library)){v<-other[other$sig_id==id,]; ii<-match(g$gene,v$gene);stopifnot(nrow(v)==31,!anyNA(ii),all(v$dir[ii]==g$dir))}
# Public figure 3 labels were inspected visually. Read each label's red/green/grey
# colour relative to its Table-3 row; do not use the colour of the bars.
im<-jpeg::readJPEG(file.path(out,"Koppelstaetter_Fig3.jpg"))
stopifnot(identical(dim(im),c(921L,788L,3L)))
pix<-t(vapply(seq(19,829,length.out=31),function(row){
 a<-im[round(row)+(-9:9),1:111,];r<-a[,,1];g<-a[,,2];b<-a[,,3]
 c(red=sum(r-g>.20 & r-b>.20),green=sum(g-r>.17 & g-b>.14))},numeric(2)))
color<-ifelse(pix[,"red"]>5,"red",ifelse(pix[,"green"]>5,"green","grey"))
stopifnot(sum(color=="red")==19,sum(color=="green")==4,sum(color=="grey")==8)
stopifnot(all(t3$ckd_fc[color=="red"]>0),all(t3$ckd_fc[color=="green"]<0))
sup<-as.data.frame(read_excel(file.path(inputs,"published/supp/raag/mmc1.xlsx")))
i<-match(t3$table3_gene,sup$Symbol);stopifnot(!anyNA(i),nrow(sup)==634)
parent<-s1genes[s1genes$sig_id==allid,];stopifnot(nrow(parent)==634,setequal(parent$gene_orig,sup$Symbol),length(unique(parent$gene))==634,all(parent$dir==0))
I<-readRDS(file.path(data_root,"benchmark_intermediates.rds"))
cb<-readRDS(file.path(data_root,"harmonized/harmonized_bundle_corrected.rds"))
cl<-readRDS(file.path(data_root,"harmonized/clinical_bundle.rds"))
clinical_matrix<-list(ERCB_GSE104954=cl$GSE104954$expr[,cl$GSE104954$meta$diagnosis!="Tumor nephrectomy"],GSE175759=cl$GSE175759$expr,TX_GSE147451=cb$harm$GSE147451$expr[,!is.na(cb$harm$GSE147451$meta$gfr_group)])
addcoverage<-function(df){
 for(nm in names(I$Z)){
  v<-I$Z[[nm]]; df[[paste0(nm,"_present")]]<-df$current_gene%in%names(v)
  df[[paste0(nm,"_finite_gene_age_z")]]<-is.finite(v[df$current_gene])
 }
 for(nm in names(clinical_matrix)){
  e<-clinical_matrix[[nm]]; df[[paste0(nm,"_present")]]<-df$current_gene%in%rownames(e)
  vf<-apply(e[intersect(df$current_gene,rownames(e)),,drop=FALSE],1,sd)
  df[[paste0(nm,"_variable")]]<-is.finite(vf[df$current_gene])&vf[df$current_gene]>0
 }
 df
}
result<-data.frame(table3_order=1:31,source_symbol=t3$table3_gene,current_gene=g$gene,ensembl_id=sup[i,"Ensembl GeneID"],table3_aging_evidence=t3$reported_aging_source,source_ckd_fc=t3$ckd_fc,source_ckd_progression_direction=sign(t3$ckd_fc),source_ckd_p=t3$ckd_p,source_ckd_FDR_band=t3$ckd_fdr_band,figure3_label_color=color,explicit_age_direction_in_table3=NA_integer_,age_direction_supported_by_named_concordant_class=ifelse(color=="grey",NA_integer_,sign(t3$ckd_fc)),age_direction_implied_by_figure3_category=ifelse(color=="grey",-sign(t3$ckd_fc),sign(t3$ckd_fc)),figure3_age_ckd_relationship=ifelse(color=="grey","discordant (author figure 3 caption)","concordant (author figure 3 caption)"),is_source_concordant23=color!="grey",source_age_direction_basis="Inferred from published Figure 3 concordance label plus signed Table 3 CKD FC; no per-gene numeric age FC in Table 3/S1",current_signed_score_direction=g$dir,current_signed_direction_source="CKD progression (progressive vs stable); S1b dir and script 20 rc",stringsAsFactors=FALSE)
result<-addcoverage(result);emit(result,"RAAG_provenance_audit.csv")
full<-data.frame(source_symbol=sup$Symbol,current_gene=parent$gene[match(sup$Symbol,parent$gene_orig)],ensembl_id=sup[["Ensembl GeneID"]],source_count=sup$Sum,subset31=sup$Symbol%in%t3$table3_gene,subset23=sup$Symbol%in%t3$table3_gene[color!="grey"])
full<-addcoverage(full);emit(full,"full_RAAG_membership.csv")
dfcov<-do.call(rbind,lapply(grep("_present$",names(result),value=TRUE),function(v){
 prefix<-sub("_present$","",v);sv<-if(paste0(prefix,"_finite_gene_age_z")%in%names(result))paste0(prefix,"_finite_gene_age_z") else paste0(prefix,"_variable")
 data.frame(cohort=prefix,present_A09=sum(result[[v]]),usable_A09=sum(result[[sv]]),usable_CKD_up=sum(result[[sv]]&g$dir==1),usable_CKD_down=sum(result[[sv]]&g$dir== -1),usable_A08=sum(full[[sv]]))}))
emit(dfcov,"RAAG_cohort_coverage.csv")
# Named functions only; no writer/network setup is sourced.
rhs<-function(fname,name){
 e<-as.list(parse(file.path(scripts,fname)));hits<-Filter(function(x)is.call(x)&&identical(x[[1]],as.name("<-"))&&identical(x[[2]],as.name(name)),e);stopifnot(length(hits)==1);hits[[1]][[3]]
}
lock<-rhs("04_locked_definitions.R","LOCK")
LOCK<-list(min_genes=5,equivalence_bound=.20,proxy=eval(lock[["proxy"]]))
inj<-eval(lock[["outcome"]][["PT_injury"]]);mz<-eval(rhs("04_locked_definitions.R","mz"));pspear<-eval(rhs("04_locked_definitions.R","pspear"));score<-eval(rhs("51_genomemed_upgrades.R","score_sig"));s<-kab_library[kab_library$sig_id==id,]
test<-function(a,b)stopifnot(all(is.finite(a)),all(is.finite(b)),length(a)==length(b),max(abs(a-b))<1e-10)
keys<-data.frame();res<-list()
add<-function(panel,endpoint,model,est,lo=NA_real_,hi=NA_real_,p=NA_real_,unit,nn=NA_integer_,nunit="",interval="",source,cohort="",fdr=NA_real_,sig_id=id,rank_in_family=NA_integer_){
 res[[length(res)+1]]<<-data.frame(panel,endpoint,model,cohort,sig_id,est,lo,hi,p,fdr,unit,n=nn,n_unit=nunit,interval,source,rank_in_family)
}
bn<-br("BENCHMARK_FINAL_all_sets");bn<-bn[order(-bn$enrichment),]
for(kid in c(allid,id)){z<-bn[bn$sig_id==kid,]; add("F2A","absolute gene-level age association","original",z$enrichment,p=z$emp_p,unit="enrichment ratio",nn=z$n_measured,nunit="measured genes",source="BENCHMARK_FINAL_all_sets.csv: enrichment,emp_p,fdr_all,n_measured",sig_id=kid,rank_in_family=match(kid,bn$sig_id),fdr=z$fdr_all)}
four<-br("REVISION_gene_level_four_configs");vif<-br("REVISION_gene_level_composition_adjusted")
for(kid in c(allid,id)){
 for(j in which(four$sig_id==kid)){z<-four[j,];add(if(z$config%in%c("primary_locked","noCPTAC_evalbins"))"F2F/sensitivity" else "S2","unsigned gene-level enrichment",z$config,z$enrichment,p=z$emp_p,fdr=z$fdr,unit="enrichment ratio",nn=z$n,nunit="measured genes",source="REVISION_gene_level_four_configs.csv: enrichment,emp_p,fdr",sig_id=kid)}
 for(j in which(vif$sig_id==kid)){z<-vif[j,]; add("F2D/S2","co-expression (approximate sensitivity)",z$config,z$ER,p=z$p_vif,fdr=z$fdr_vif,unit="enrichment ratio (unchanged by VIF)",nn=z$n,nunit="measured genes",source="REVISION_gene_level_composition_adjusted.csv: ER,p_vif,fdr_vif",sig_id=kid)}
}
# Re-pool the original archived per-cohort correlations with the unchanged two-stage procedure.
MS<-br("M_signed_by_dataset");MS<-MS[MS$sig_id==id,];leak<-eval(rhs("21_benchmark.R","leak"));poolfun<-eval(rhs("21_benchmark.R","poolfun"))
MSp<-eval(rhs("21_benchmark.R","MSp"));f<-figdata("F2_F3_signed_intervals")
for(j in seq_len(nrow(MSp))){z<-MSp[j,];v<-f[f$sig_id==id&f$model==z$model,]
 test(c(z$pooled_rho,z$pooled_lo,z$pooled_hi,z$pooled_p),c(v$rho,v$lo,v$hi,v$pooled_p))
 add(if(z$model=="base") "F2B / F3C" else "F3C","signed score vs age",z$model,z$pooled_rho,z$pooled_lo,z$pooled_hi,z$pooled_p,"pooled rank correlation",nn=645,nunit="reference profiles available in contributing datasets",interval="pooled Fisher scale 95% CI",source="M_signed_by_dataset.csv; script21 MSp; F2_F3_signed_intervals.csv",rank_in_family=1L+sum(f$rho[f$model==z$model]>z$pooled_rho))}
for(j in seq_len(nrow(MS))){z<-MS[j,];add(if(z$model=="base") "F2C" else "S1d","cohort signed score vs age",z$model,z$rho,z$lo,z$hi,z$p,"rank correlation",nn=z$n,nunit="profiles",interval="Fisher 95% CI",source="M_signed_by_dataset.csv",cohort=z$dataset)}
pools<-br("REVISION_gene_pools_comp_and_ckd")
for(j in grep("^RAAG,",pools$pool)){z<-pools[j,];co<-if(grepl("concordant",z$pool))"CKD progression-associated RAAGs" else "Other RAAGs"; add("F3D / S10","gene pool enrichment",z$analysis,z$ER,z$lo,z$hi,z$emp_p,"enrichment ratio",nn=z$n,nunit="measured genes",interval="matched-null range, NOT effect CI",source="REVISION_gene_pools_comp_and_ckd.csv: ER,lo,hi,emp_p",cohort=co,sig_id=if(co=="Other RAAGs")"A08 minus A09" else id)}
# Refit only unchanged CKD/eGFR/Tx models over exact original units. No sensitivity analysis.
tab<-br("REVISION_clinical_CKD_and_transplant");ref<-tab[tab$sig_id==id,]
f4<-figdata("F4_ckd_intervals");eg<-figdata("F4_egfr_intervals");tx<-figdata("F5_outcome_intervals")
for(cg in c("GSE104954","GSE175759")){
 d<-cl[[cg]];m<-d$meta;if(cg=="GSE104954")m<-m[m$diagnosis!="Tumor nephrectomy",]; e<-d$expr[,m$gsm]
 dx<-as.numeric(m$grp=="Disease");v<-score(e,s);z<-(v-mean(v[dx==0]))/sd(v[dx==0]);fd<-figdata(if(cg=="GSE104954")"F4_ercb_biopsies" else "F4_gse175759_biopsies"); test(z,fd[[id]][match(m$gsm,fd$gsm)]);px<-as.data.frame(lapply(LOCK$proxy,function(g)mz(e,g)));injury<-mz(e,inj)
 fits<-if(cg=="GSE104954")list(lm(z~dx+factor(m$platform)),lm(z~dx+factor(m$platform)+injury+.,data=px)) else list(lm(z~dx),lm(z~dx+injury+.,data=px))
 for(j in 1:2){name<-c("raw","adj")[j];fit<-fits[[j]];cf<-coef(summary(fit))["dx",];CI<-confint(fit,"dx");cgout<-if(cg=="GSE104954")"ERCB" else cg
 v0<-f4[f4$sig_id==id&f4$cohort==cgout&f4$model==name,]
 test(c(cf[1],CI[1],CI[2],cf[4]),c(v0$beta,v0$lo,v0$hi,v0$p))
 add(if(cg=="GSE104954")"F4D" else "F4E","CKD vs reference",name,cf[1],CI[1],CI[2],cf[4],"reference-group SD",nobs(fit),"biopsies","OLS t 95% CI","51:296–303; F4_ckd_intervals.csv",cgout)}
 if(cg=="GSE104954"){
  rho<-cor(v,injury,method="spearman");test(rho,ref$ERCB_cor_PT_injury);add("S3","signed score vs tubular injury","unadjusted",rho,unit="Spearman correlation",nn=nrow(m),nunit="biopsies",source="51:299; REVISION_clinical_CKD_and_transplant.csv",cohort="ERCB")
 }else for(j in 1:2){mod<-c("raw","adj")[j];q<-pspear(m$egfr,v,if(j==1) NULL else cbind(px,injury));v0<-eg[eg$sig_id==id&eg$model==mod,];test(unlist(q[,c("rho","lo","hi","p")]),unlist(v0[,c("rho","lo","hi","p")]))
 add("F4F","score vs eGFR",mod,q$rho,q$lo,q$hi,q$p,"rank correlation",q$n,"biopsies","Fisher 95% CI","51:304–305; F4_egfr_intervals.csv",cg)}
}
h<-cb$harm$GSE147451;m<-h$meta;ok<-!is.na(m$gfr_group); v<-score(h$expr,s)[ok];outcome<-as.numeric(m$gfr_group[ok]=="poor");age<-data.frame(age=m$age[ok]);px<-as.data.frame(lapply(LOCK$proxy,function(g)mz(h$expr,g)))[ok,,drop=FALSE]
f5d<-figdata("F5_donors");test(v,f5d[[id]][match(m$gsm[ok],f5d$gsm)])
for(j in 1:3){mod<-c("raw","age","full")[j];q<-pspear(outcome,v,list(NULL,age,cbind(age,px))[[j]]);v0<-tx[tx$sig_id==id&tx$model==mod,];test(unlist(q[,c("rho","lo","hi","p")]),unlist(v0[,c("rho","lo","hi","p")]))
 add(c("F5D","F5E","F5F")[j],"score vs recorded 24-month graft group",mod,q$rho,q$lo,q$hi,q$p,"rank correlation",q$n,"inferred donors","Fisher 95% CI","51:307–310; F5_outcome_intervals.csv","GSE147451")}
emit(do.call(rbind,res),"RAAG_verified_results.csv")
cat("31 symbols and CKD signs agree with original Table 3, project Table S1, signature CSV and analysis RDA.\n")
print(table(figure3=color,CKDsign=sign(t3$ckd_fc)));print(dfcov,row.names=FALSE)
cat("A08: 634 source rows, 634 unique current symbols, ",sum(parent$gene!=parent$gene_orig)," changed aliases.\n")
print(data.frame(original_code=id,original_rank_48=match(id,bn$sig_id),original_rank_class_A=match(id,bn$sig_id[bn$class=="A"])))
