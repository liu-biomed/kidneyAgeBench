# Direct accession/subject crosswalk audit. Never equate two different local ID namespaces.
suppressPackageStartupMessages(library(readxl))
project<-file.path(path.expand("~"),"Kidney_Manuscript_Packages");out<-file.path(project,"submission_revision/raag")
pack<-Sys.glob(file.path(project,"Genome_Medicine_*"))[1];inputs<-file.path(pack,"04_Source_data/c_Datasets/01_Study_inputs")
base<-file.path(pack,"04_Source_data/c_Datasets/02_Derived/kidney_aging_data")
cb<-readRDS(file.path(base,"harmonized/harmonized_bundle_corrected.rds"));cl<-readRDS(file.path(base,"harmonized/clinical_bundle.rds"))
records_soft<-function(path){all<-readLines(path,warn=FALSE); ix<-which(startsWith(all,"^SAMPLE = ")); end<-c(ix[-1]-1L,length(all))
 do.call(rbind,lapply(seq_along(ix),function(k){ss<-all[ix[k]:end[k]]; ti<-grep("^!Sample_title = ",ss,value=TRUE)
   item<-ss[grepl("^!Sample_(title|description|source_name_ch1|characteristics_ch1|relation) = ",ss)]
   data.frame(gsm=sub("^\\^SAMPLE = ","",ss[1]),title=sub("^!Sample_title = ","",ti[1]),metadata=paste(item,collapse="; "))}))}
rec<-records_soft(file.path(out,"GSE60861_geo_header.soft"))
pt<-as.data.frame(read_excel(file.path(inputs,"published/supp/raag/mmc2.xlsx")))
src<-rec[match(pt[[1]],rec$title),];stopifnot(nrow(src)==63,!anyNA(src$gsm),!anyDuplicated(src$title),!anyDuplicated(src$gsm),all(rec$title %in% rec$title))
counts<-table(pt[[2]]);stopifnot(unname(counts["stable"])==39,unname(counts["progressive"])==24)
write.csv(data.frame(source_nc=src$title,source_gsm=src$gsm,source_group=pt[[2]]),file.path(out,"GSE60861_63_source_keys.csv"),row.names=FALSE)
records_matrix<-function(f){
 z<-gzfile(f,open="rt"); on.exit(close(z)); ss <- character();repeat{line<-readLines(z,n=1,warn=FALSE);if(!length(line)||line=="!series_matrix_table_begin")break;if(startsWith(line,"!Sample_"))ss<-c(ss,line)}
 pick<-function(k){idx<-startsWith(ss,paste0("!Sample_",k,"\t"));lapply(strsplit(ss[idx],"\t",fixed=TRUE),function(i)gsub('^"|"$',"",i[-1]))}
 gsm<-pick("geo_accession")[[1]];ti<-pick("title")[[1]]
 vv<-c(pick("title"),pick("description"),pick("source_name_ch1"),pick("characteristics_ch1"),pick("relation"))
 data.frame(gsm=gsm,title=ti,metadata=apply(do.call(cbind,vv),1,paste,collapse="; "))}
geos<-names(cb$harm); ids<-lapply(cb$harm,function(h)unique(unlist(strsplit(h$meta$gsm,"+",fixed=TRUE))))
ids[["GSE104954"]]<-with(cl$GSE104954$meta,gsm[diagnosis!="Tumor nephrectomy"]);ids[["GSE175759"]]<-cl$GSE175759$meta$gsm
rows<-list()
audit<-function(nm,ss,note){
 co<-intersect(src$gsm,ss$gsm);nc<-src$title[vapply(src$title,function(id)any(grepl(paste0("(?<![A-Za-z0-9])",id,"(?![A-Za-z0-9])"),ss$metadata,perl=TRUE)),TRUE)]
 data.frame(source="63 identified CKD progression patients (GSE60861)",tested=nm,evaluation_records=nrow(ss),source_gsm_matches=length(co),source_nc_token_matches=length(nc),disclosure=note)}
for(nm in names(ids)){
 files<-Sys.glob(file.path(inputs,"geo",paste0(nm,"*series_matrix.txt.gz")))
 s<-do.call(rbind,lapply(files,records_matrix));used<-s[s$gsm%in%ids[[nm]],]
 stopifnot(setequal(used$gsm,ids[[nm]]),!anyDuplicated(used$gsm))
 note<-switch(nm,GSE60807="Same Innsbruck institution, deceased-donor time-zero biopsies; different clinical setting, no common person key.",GSE147451="177 adult kidney profiles subsequently collapsed to 166 inferred donors; transplant subset uses these profiles.",GSE104954="ERCB metadata document reuse of GSE47184 and GSE32591; checked their headers separately.",GSE175759="Seoul National University Hospital kidney biopsies/reference; no shared cross-institution patient key.", "No common person-level identifier links study-specific IDs.")
 rows[[nm]]<-audit(nm,used,note)
}
for(nm in c("GSE47184","GSE32591"))rows[[nm]]<-audit(paste(nm,"ERCB predecessor"),records_soft(file.path(out,paste0(nm,"_geo_header.soft"))),"Earlier GEO deposit linked by ERCB; check all released sample records.")
# TCGA case identifiers (129) and CPTAC case IDs (75) as evaluated; use original metadata.
q<-read.csv(file.path(base,"validation/tcga/tcga_normal_meta.csv"));sx<-read.csv(file.path(base,"validation/tcga/tcga_case_sex.csv"));q$sex<-sx$sex[match(q$patient,sx$patient)]
q$age<-ifelse(is.finite(q$age_days),q$age_days/365.25,q$age_at_index);q<-q[order(q$sample),];q<-q[!duplicated(q$patient)&is.finite(q$age)&!is.na(q$sex),]
rows[["TCGA"]]<-audit("TCGA reference",data.frame(gsm=q$sample,metadata=apply(q[,c("sample","patient")],1,paste,collapse="; ")),"TCGA-barcode metadata; no shared person crosswalk.");stopifnot(nrow(q)==129)
# These original helper functions read the package's public CPTAC matrices without sourcing its global I/O.
suppressMessages(library(data.table))
script<-file.path(pack,"04_Source_data/b_Original_code/revision/scripts/10_cptac_helpers.R");h<-as.list(parse(script))
for(n in c("cptac_rd","cptac_clin","cptac_sym")){expr<-Filter(function(x)is.call(x)&&identical(x[[1]],as.name("<-"))&&identical(x[[2]],as.name(n)),h);eval(expr[[1]])}
V<-file.path(inputs,"validation/cptac"); if(!dir.exists(V))V<-file.path(base,"validation/cptac")
mp<-cptac_rd(V,"bcm-ccrcc-mapping");tr<-cptac_rd(V,"bcm-ccrcc-transcriptomics");cc<-cptac_clin(V)
cs<-intersect(sub("_A$","",grep("_A$",names(tr),value=TRUE)),cc$case[is.finite(cc$age)&cc$sex%in%c("Male","Female")])
stopifnot(length(cs)==75);rows[["CPTAC"]]<-audit("CPTAC reference",data.frame(gsm=cs,metadata=cs),"CPTAC case keys; no shared person crosswalk.")
gt<-cb$gtex$meta;rows[["GTEx"]]<-audit("GTEx reference",data.frame(gsm=gt$gsm,metadata=paste(gt$gsm,gt$SUBJID)),"GTEx keys do not match progression samples. HPA expression filter is a separate possible indirect reuse; see HPA sample-key audit.")
all<-do.call(rbind,rows);write.csv(all,file.path(out,"overlap_keys_by_cohort.csv"),row.names=FALSE)
# HPA v18 archived kidney tissue data: distinguish access to expression-only reference from age/outcome derivation.
html<-paste(readLines(file.path(out,"HPA_v18_IL7R_kidney.html"),warn=FALSE),collapse=" ")
hs<-unique(regmatches(html,gregexpr("GTEX-[A-Z0-9]+-[A-Z0-9]+-SM-[A-Z0-9]+",html))[[1]])
su<-sub("(^GTEX-[^-]+).*","\\1",hs);maps<-data.frame(hpa_sample=hs,hpa_subject=su,match_sample=hs%in%gt$gsm,match_subject=su%in%gt$SUBJID)
write.csv(maps,file.path(out,"HPA_GTEx_expression_filter_overlap.csv"),row.names=FALSE)
stopifnot(nrow(maps)==32,sum(maps$match_sample)==28,length(unique(maps$hpa_subject[maps$match_subject]))==29)
print(all[,c("tested","evaluation_records","source_gsm_matches","source_nc_token_matches")],row.names=FALSE)
cat("Source patients: ",nrow(src),", clinical groups: stable 39 / progressive 24.\n")
cat("Public HPA-v18 IL7R kidney records ",length(hs),"; shared exact sample IDs ",sum(maps$match_sample),"; shared donor IDs ",length(unique(maps$hpa_subject[maps$match_subject])),".\n")
