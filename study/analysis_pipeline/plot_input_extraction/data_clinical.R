# From redesign_v2. Read original model definitions without sourcing scripts that write or fetch.
original <- ".."
inputs <- path.expand(Sys.getenv("KIDNEY_DATA","~/kidney_aging_data"))
read_result <- function(f) read.csv(file.path(original,"benchmark/results",f),check.names=FALSE)
write_data <- function(x,n) write.csv(x,file.path("source_data",paste0(n,".csv")),row.names=FALSE,na="")
definition <- function(file,name){
 e <- as.list(parse(file.path(original,"scripts",file)))
 hit <- Filter(function(x) is.call(x) && identical(x[[1]],as.name("<-")) && identical(x[[2]],as.name(name)),e)
 stopifnot(length(hit)==1);hit[[1]][[3]]
}
spec <- definition("04_locked_definitions.R","LOCK")
LOCK <- list(min_genes=5, equivalence_bound=.20, proxy=eval(spec[["proxy"]]))
injgenes <- eval(spec[["outcome"]][["PT_injury"]])
mz <- eval(definition("04_locked_definitions.R","mz"))
pspear <- eval(definition("04_locked_definitions.R","pspear"))
signed_score <- eval(definition("51_genomemed_upgrades.R","score_sig"))
load(file.path(original,"kidneyAgeBench/data/kab_library.rda"))
library_sets <- split(kab_library,kab_library$sig_id)
ids <- names(Filter(function(s)any(s$dir!=0),library_sets))
label_key <- read.csv("source_data/sets_key.csv")
sigmap <- data.frame(sig_id=ids,label=label_key$short_label[match(ids,label_key$sig_id)])
stopifnot(!anyNA(sigmap$label))
stopifnot(identical(ids,sigmap$sig_id))
write_data(sigmap,"signed_key")
cl <- readRDS(file.path(inputs,"harmonized/clinical_bundle.rds"))
obs <- read_result("REVISION_clinical_CKD_and_transplant.csv")
stopifnot(identical(ids,obs$sig_id))
diagnoses <- c("Living donor"="LD","Thin membrande disease"="TMD","Minimal change disease"="MCD",
 "Focal segmental glomerulosclerosis and minimal change disease"="FSGS/MCD",
 "Focal segmental glomerulosclerosis"="FSGS","Membranous glomerulonephropathy"="MN","IgA nephropathy"="IgAN",
 "Hypertensive nephropathy"="HTN","Systemic lupus erythematosus"="LN","Diabetic nephropathy"="DKD",
 "ANCA-associated vasculitis"="ANCA")
refz <- function(x,r) (x-mean(x[r]))/sd(x[r])
same <- function(a,b)stopifnot(all(is.finite(a)),all(is.finite(b)),max(abs(a-b))<1e-9)
allfit <- list(); renal <- list()
for(cohort in c("ERCB","GSE175759")){
 x <- cl[[if(cohort=="ERCB")"GSE104954" else "GSE175759"]]
 keep <- if(cohort=="ERCB") x$meta$diagnosis!="Tumor nephrectomy" else rep(TRUE,nrow(x$meta))
 m <- x$meta[keep,]; e <- x$expr[,keep]
 stopifnot(identical(as.character(m$gsm),colnames(e)))
 disease <- as.numeric(m$grp=="Disease"); ref <- disease==0
 injury <- mz(e,injgenes); proxies <- as.data.frame(lapply(LOCK$proxy,function(g)mz(e,g)))
 scores <- sapply(ids,function(id) signed_score(e,library_sets[[id]]))
 display <- cbind(score_injury=injury,proxies,scores)
 samples <- data.frame(gsm=m$gsm,diagnosis=m$diagnosis,group=ifelse(ref,"Reference","CKD"),
                       egfr=if(cohort=="ERCB")NA_real_ else m$egfr,
                       apply(display,2,refz,r=ref),check.names=FALSE)
 if(cohort=="ERCB"){
   samples$dx <- unname(diagnoses[samples$diagnosis])
   stopifnot(!anyNA(samples$dx))
   samples <- samples[order(match(samples$dx,diagnoses),samples$score_injury,samples$gsm),]
 }
 write_data(samples,if(cohort=="ERCB")"F4_ercb_biopsies" else "F4_gse175759_biopsies")
 for(id in ids){
  i <- match(id,ids); z <- refz(scores[,id],ref)
  fits <- if(cohort=="ERCB")list(lm(z~disease+factor(m$platform)),lm(z~disease+factor(m$platform)+injury+.,data=proxies)) else
    list(lm(z~disease),lm(z~disease+injury+.,data=proxies))
  for(k in 1:2){
   mod <- fits[[k]]; ci <- confint(mod,"disease");co <- coef(summary(mod))["disease",];type <- c("raw","adj")[k]
   stem <- if(cohort=="ERCB")"ERCB" else "GSE175"
   same(co[1],obs[[paste(stem,"beta",type,sep="_")]][i]);same(co[4],obs[[paste(stem,"p",type,sep="_")]][i])
   allfit[[length(allfit)+1]] <- data.frame(cohort,sig_id=id,model=type,n=nobs(mod),df=df.residual(mod),
    beta=unname(co[1]),se=unname(co[2]),lo=ci[1],hi=ci[2],p=unname(co[4]))
  }
  if(cohort=="GSE175759"){
   es <- rbind(pspear(m$egfr,z),pspear(m$egfr,z,cbind(proxies,inj=injury)))
   same(es$rho,unlist(obs[i,c("GSE175_eGFR_rho_raw","GSE175_eGFR_rho_adj")]))
   renal[[id]] <- data.frame(sig_id=id,model=c("raw","adj"),es)
  }else same(cor(z,injury,method="spearman"),obs$ERCB_cor_PT_injury[i])
 }
}
write_data(do.call(rbind,allfit),"F4_ckd_intervals")
write_data(do.call(rbind,renal),"F4_egfr_intervals")
# Donor grouping, expression and endpoint are inherited rather than reassigned.
bundle <- readRDS(file.path(inputs,"harmonized/harmonized_bundle_corrected.rds"))
donor <- bundle$harm$GSE147451;m <- donor$meta; e <- donor$expr;ok <- !is.na(m$gfr_group)
stopifnot(identical(as.character(m$gsm),colnames(e)),sum(ok)==151,sum(ok & m$gfr_group=="poor",na.rm=TRUE)==62)
kid <- donor$kidney_level$meta
pairs <- do.call(rbind,lapply(split(kid,kid$donor),function(g)data.frame(gsm=paste(g$gsm,collapse="+"),
   kidneys=nrow(g),kidney_labels=paste(ifelse(is.na(g$gfr_group),"missing",g$gfr_group),collapse="/"),
   retained_label=as.character(g$gfr_group[1]))))
stopifnot(nrow(pairs)==166,sum(pairs$kidneys)==177,sum(pairs$kidneys==2)==11,
 sum(pairs$kidney_labels %in% c("good/poor","poor/good"))==1,
 identical(as.character(pairs$retained_label[match(m$gsm,pairs$gsm)]),as.character(m$gfr_group)))
write_data(pairs,"F5_donor_linkage")
v <- sapply(ids,function(id)signed_score(e,library_sets[[id]]))
px <- as.data.frame(lapply(LOCK$proxy,function(g)mz(e,g)))[ok,,drop=FALSE]
tx <- data.frame(gsm=m$gsm[ok],age=m$age[ok],outcome=m$gfr_group[ok],v[ok,,drop=FALSE],check.names=FALSE)
write_data(tx,"F5_donors")
txfit <- do.call(rbind,lapply(ids,function(id){
 y <- as.numeric(tx$outcome=="poor");x <- tx[[id]]
 rows <- rbind(pspear(y,x),pspear(y,x,data.frame(age=tx$age)),pspear(y,x,cbind(age=tx$age,px)))
 same(rows$rho,unlist(obs[match(id,ids),c("TX_poorGFR_rho_raw","TX_poorGFR_rho_adj_age","TX_poorGFR_rho_adj_age_comp")]))
 data.frame(sig_id=id,model=c("raw","age","full"),rows)
}))
write_data(txfit,"F5_outcome_intervals")

