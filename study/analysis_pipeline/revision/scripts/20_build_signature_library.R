# Benchmark signature library: every gene list with its source, class, and (where defined) direction with age.
# Symbols are updated to current HGNC via limma::alias2SymbolTable (org.Hs.eg.db).
suppressMessages({library(limma); library(org.Hs.eg.db); library(msigdbr)})
REV <- path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030")
P <- path.expand("~/kidney_aging_data/published"); OUTD <- file.path(REV, "benchmark", "signatures")
L <- list(); add <- function(id, cls, label, pmid, genes, dir=0, note="") { if (!length(genes)) { cat("EMPTY:", id, "\n"); return(invisible()) }
  L[[id]] <<- data.frame(sig_id=id, class=cls, label=label, pmid=pmid, gene_orig=genes, dir=rep_len(dir, length(genes)), note=note, stringsAsFactors=FALSE) }
c2 <- as.data.frame(msigdbr(species="Homo sapiens", collection="C2")); g2 <- function(n) unique(c2$gene_symbol[c2$gs_name==n])
## A: kidney, human, age-derived (or age-confounded kidney) signatures
add("A01_Rodwell2004", "A", "Rodwell 2004 human kidney ageing (cortex+medulla, 74 kidneys)", "15562319", c(g2("RODWELL_AGING_KIDNEY_UP"), g2("RODWELL_AGING_KIDNEY_DN")), c(rep(1, length(g2("RODWELL_AGING_KIDNEY_UP"))), rep(-1, length(g2("RODWELL_AGING_KIDNEY_DN")))))
add("A02_Rodwell2004_noblood", "A", "Rodwell 2004, blood-cell genes removed", "15562319", c(g2("RODWELL_AGING_KIDNEY_NO_BLOOD_UP"), g2("RODWELL_AGING_KIDNEY_NO_BLOOD_DN")), c(rep(1, length(g2("RODWELL_AGING_KIDNEY_NO_BLOOD_UP"))), rep(-1, length(g2("RODWELL_AGING_KIDNEY_NO_BLOOD_DN")))))
mk <- read.csv(file.path(REV, "benchmark", "signatures", "Melk2005_supp2_parsed.csv"))
add("A03_Melk2005_505", "A", "Melk 2005 renal-ageing set (SAM multiclass, old vs adult direction)", "16316342", mk$symbol, sign(mk$old_vs_adult), "Supp Table 2 recovered from web.archive.org 2010-07-21")
t2 <- c(VCAN=1, IGLJ3=1, CD52=1, TNC=1, JCHAIN=1, CDC25C=1, ESPL1=1, FN1=1, SLC31A1=-1, GDA=-1, PLGLB1=-1, CALML3=-1, CDR1=-1, BPHL=-1, ASS1=-1, AK3=-1, CYP8B1=-1, BHMT2=-1,
        AIFM1=-1, WNK1=-1, ALDH6A1=-1, SLC13A3=-1, GPD1=-1, FABP1=-1, PLG=-1, HCN3=-1, DDC=-1, APOE=-1, PSAT1=-1, POLR3H=-1, DIO1=-1, PCYT1B=-1, SLC7A7=-1, SLC37A4=-1, HNF4A=-1, SEPTIN2=-1)
add("A04_Melk2005_PAM", "A", "Melk 2005 Table 2 age-flagged PAM genes (high-pathology old vs adult)", "16316342", names(t2), unname(t2), "Table 2 genes with Age=Yes; legacy symbols updated manually")
rw <- c(ZNF518B=-1, PTCHD4=1, NCS1=1, LYG1=1, `LHFPL3-AS2`=-1, EHD3=-1, PDZK1IP1=1, APLNR=-1, EGF=-1, FAM150B=1, SLC30A8=-1, TMEM54=1, `RP11-352M15.2`=-1, HCG11=-1, SLC16A5=1,
        ATP1B2=-1, EIF4G1=1, VIM=1, SNX24=-1, GUCY1A2=-1, MAP4=1, PPP1R3C=1, ADAMTS5=-1, ARPC1A=1, NREP=-1, SARS2=-1, PHLDA3=1, NUDT3=1, USH1C=1, ERP27=1, SLPI=1, SNHG19=-1, TRIM26=-1, TSPYL5=-1, ARFGEF3=-1, LTF=1)
add("A05_Rowland2019_37", "A", "Rowland 2019 TRANSLATE+TCGA 37 age genes (direction from TRANSLATE beta; WDR55 excluded, beta~0)", "30784661", names(rw), unname(rw))
rr <- c(EGF=-1, SLC16A5=1, VIM=1, TSPYL5=-1, LTF=1, PHLDA3=1, MAP4=1, ATP1B2=-1, SLPI=1, EHD3=-1, TMEM54=1, LYG1=1, ZNF518B=-1, ARFGEF3=-1, ERP27=1, TRIM26=-1, SNX24=-1, PPP1R3C=1, NCS1=-1)
add("A06_Rowland2019_rep19", "A", "Rowland 2019 genes replicated in Nephroseq (FDR<0.05, either compartment)", "30784661", names(rr), unname(rr), "replication partly used Rodwell data")
kz_up <- c("ITGAV","REG1A","REG1B","IGKC","MALT1","CR2","CFB","DCDC2","RDH5","TKT","PSMF1","ADAM30","GRK5","STAT1","MS4A7","SHB","TNRC4","KMT2C","TF")
kz_dn <- c("PBRM1","OGN","PROK1","JCHAIN","IGLC2","MT1H","TEX27","PLA2G12B","REG4","UPP2","POLA1","STATIP1","RAB12","CLDN2","MEF2C","E2F1","ASCL1","VAMP1","SLC6A13","DAO","TMED6","DNAJC9","ZFAND2B")
add("A07_Kainz2007_donor", "A", "Kainz 2007 donor biopsy poor vs good 1-yr graft function (age-confounded donor quality)", "17452894", c(kz_up, kz_dn), c(rep(1, length(kz_up)), rep(-1, length(kz_dn))), "direction = higher in poor-function group; not age-derived")
ra <- read.csv(file.path(P, "supp", "raag_S1.csv"))
add("A08_Perco2019_RAAG634", "A", "Renal age-associated genes (6 sources incl. Rodwell, Melk, GenAge, literature)", "31316728", ra$Symbol, 0, "unsigned; partly non-transcriptomic sources")
rc <- c(C3=1,LTF=1,IL7R=1,NNMT=1,CD52=1,UBD=1,APBB1IP=1,CFB=1,MMP7=1,VCAN=1,CLDN1=1,ITGB2=1,MICALL2=1,INPP5D=1,CLU=1,ITPR3=1,TMPRSS4=1,WFDC2=1,CDH6=1,TSPAN1=1,ANXA1=1,GLS=1,TNFRSF11B=1,CXCL6=1,TPBG=1,NOS3=-1,CGNL1=-1,STARD8=-1,CDKN2A=-1,LPL=-1,EGF=-1)
add("A09_Perco2019_RAAG_CKD", "A", "RAAGs associated with CKD progression (direction = CKD progression)", "31316728", names(rc), unname(rc))
add("A10_Sinning2023_MHCII", "A", "Tubuloinflammaging MHC-II genes (human scRNA)", "37675124", c("CD74","HLA-DQB1","HLA-DRB1"), 1)
k94 <- read.csv(file.path(REV, "..", "tables", "robust_aging_signature.csv"))
add("A11_KAS94_original", "A", "This study, original 94-gene signature (9 cohorts)", "this study", k94$gene, ifelse(k94$direction=="up", 1, -1), "in-sample for discovery cohorts")
kc <- readRDS(path.expand("~/kidney_aging_data/harmonized/meta_primary_corrected.rds"))$rob
add("A12_KAS82_corrected", "A", "This study, corrected 82-gene signature (6 cohorts)", "this study", kc$gene, ifelse(kc$direction=="up", 1, -1), "in-sample for discovery cohorts")
## B: cross-tissue / ageing-biology references
add("B01_deMagalhaes2009", "B", "de Magalhaes 2009 cross-species ageing meta-analysis", "19189975", c(g2("DEMAGALHAES_AGING_UP"), g2("DEMAGALHAES_AGING_DN")), c(rep(1, length(g2("DEMAGALHAES_AGING_UP"))), rep(-1, length(g2("DEMAGALHAES_AGING_DN")))))
pm <- list(up=read.csv(file.path(P, "hagr", "palmer_over.csv")), dn=read.csv(file.path(P, "hagr", "palmer_under.csv")))
add("B02_Palmer2021", "B", "Palmer 2021 mammalian multi-tissue ageing signature (q<0.05)", "33568573", c(pm$up$Gene, pm$dn$Gene), c(rep(1, nrow(pm$up)), rep(-1, nrow(pm$dn))))
ga <- read.csv(file.path(P, "hagr", "genage_human.csv")); add("B03_GenAge", "B", "GenAge human ageing genes (build 21)", "29121237", ga$symbol, 0, "unsigned")
ca <- read.delim(file.path(P, "hagr", "cellage3.tsv")); add("B04_CellAge", "B", "CellAge senescence genes", "29121237", ca$Gene.symbol, 0, "unsigned")
add("B05_SenMayo", "B", "SenMayo senescence/SASP (Saul 2022)", "35974106", g2("SAUL_SEN_MAYO"), 1)
add("B06_Fridman2008_up", "B", "Fridman 2008 senescence-up", "18332867", g2("FRIDMAN_SENESCENCE_UP"), 1)
mo <- read.csv(file.path(REV, "results", "validation", "MOUSE_PTAE_by_gene.csv")); n <- 188; q_mo <- 6; z <- atanh(mo$rna_age)*sqrt(n-3-q_mo); fq <- p.adjust(2*pnorm(-abs(z)), "BH") # q_mo=6 matches locked Addendum 5 df specification for DO mouse covariates
add("B07_Takemon2021_mouseRNA", "B", "DO mouse kidney RNA age genes (6-18 mo, FDR<0.05, orthologues)", "33687326", mo$gene[fq < .05], sign(mo$rna_age[fq < .05]))
## C: 'ageing/senescence-related' gene models from kidney disease bioinformatics studies (unsigned w.r.t. age)
C <- list(C01=list("37760894","DN ageing-related hub genes", c("CCR2","VCAM1","CSF1R","ITGAM")),
  C02=list("37396184","DKD cellular-senescence signature", c("LIMA1","ZFP36","FOS","IGFBP6","CKB")),
  C03=list("37957539","DN glomerular senescence genes", c("FOS","ZFP36")),
  C04=list("39602449","Kidney transplant rejection senescence genes", c("CYBB","ACKR1","IFNG","IL2RB","C1QA","TNFSF13B","BCL2A1","IRF1","MAPK1","TNFRSF1B","CCL4","STAT1","IL15")),
  C05=list("40671964","IgAN ageing-related hub genes", c("AGT","ALB","CD36","EGF","KDR","LPL","MYC","PPARGC1A")),
  C06=list("40593989","Renal stone ageing-related genes", c("CNR1","KIT","HTR2A","DES","IL33","UCP2","PPT1")),
  C07=list("40319112","CKD ageing-related hub genes", c("SOD2","FGF21","FOS","RELA","DDIT4","BMI1","DUSP6","LGALS3","CXCR2","CEBPB")),
  C08=list("40290492","Lupus nephritis senescence biomarkers", c("ALOX5","PTGER2","PRKCB")),
  C09=list("42152390","ccRCC ageing-related prognostic genes", c("PCK1","TOP2A","TFAP2A","CCNA2","FOXM1","CDKN2A","PLAU")),
  C10=list("35646056","ccRCC ageing/senescence prognostic model", c("FOXG1","FOXM1","GNRH1","HAMP","IGFBP2","IL10","MPEG1","VASH1")),
  C11=list("35720278","ccRCC senescence score", c("P3H1","PROX1","HJURP","HK3","CDKN1A","AR","VENTX","MAGOHB","MAP2K6")),
  C12=list("34207247","RCC normal-tissue ageing prognostic genes", c("DUSP22","MAPK14","MAPKAPK3","STAT1","VCP")),
  C13=list("42102994","Renal fibrosis ageing-related diagnostic genes", c("MMP7","AGR2","RASSF5","AHR","HIF1A","AXL","PTTG1")),
  C14=list("42227369","CKD ageing-related hub genes", c("AGER","PDGFRB","SOD2")),
  C15=list("42792788","CKD ferro-ageing hub genes", c("IRF7","SOCS2")),
  C16=list("42464365","DKD ageing-related hub genes", c("CLU","EGF","SLC16A7","MYC","RPA1","RB1","APOC3","SYK","NR3C1")))
for (k in names(C)) add(paste0(k, "_PMID", C[[k]][[1]]), "C", C[[k]][[2]], C[[k]][[1]], C[[k]][[3]], 0, "unsigned; gene list from abstract/main text")
lib <- do.call(rbind, L); rownames(lib) <- NULL
lib$gene_orig <- trimws(lib$gene_orig); lib <- lib[!is.na(lib$gene_orig) & lib$gene_orig != "" & !grepl("^EST|^Hs\\.|^[A-Z]{1,2}[0-9]{5,6}$", lib$gene_orig), ]
map <- suppressWarnings(alias2SymbolTable(lib$gene_orig, species="Hs")); lib$gene <- ifelse(is.na(map), lib$gene_orig, map)
lib <- lib[!duplicated(paste(lib$sig_id, lib$gene)), ]
s <- aggregate(gene ~ sig_id + class + label + pmid, lib, length); s$n_up <- tapply(lib$dir > 0, lib$sig_id, sum)[s$sig_id]; s$n_down <- tapply(lib$dir < 0, lib$sig_id, sum)[s$sig_id]
names(s)[names(s)=="gene"] <- "n_genes"; s <- s[order(s$sig_id), ]
write.csv(lib, file.path(OUTD, "signature_library_long.csv"), row.names=FALSE); write.csv(s, file.path(OUTD, "signature_library_summary.csv"), row.names=FALSE)
print(s, row.names=FALSE); cat("signatures:", nrow(s), "| unmapped symbols kept as-is:", sum(is.na(map)), "\n")
