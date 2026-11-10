# Reviewer 1 point 25: direct comparison with previously published ageing signatures. Because these signatures were
# derived without our data, their age association in our cohorts is free of selection bias and is the fair benchmark
# for the nested (held-out) performance of the present signature.
source(file.path(path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030"), "scripts", "fun.R"))
suppressMessages(library(msigdbr))
b <- load_bundle(); harm <- b$harm; gt <- b$gtex; P <- readRDS(MP); rob <- P$rob; res <- P$res
Lc <- read.csv(file.path(OT, "R1_nested_LOCO_folds.csv"))
M <- as.data.frame(msigdbr(species="Homo sapiens", collection="C2")); gs <- function(n) unique(M$gene_symbol[M$gs_name==n])
pub <- list(`Rodwell 2004 kidney`=list(up=gs("RODWELL_AGING_KIDNEY_UP"), dn=gs("RODWELL_AGING_KIDNEY_DN")),
            `Rodwell 2004 kidney (blood-filtered)`=list(up=gs("RODWELL_AGING_KIDNEY_NO_BLOOD_UP"), dn=gs("RODWELL_AGING_KIDNEY_NO_BLOOD_DN")),
            `de Magalhaes 2009 cross-tissue`=list(up=gs("DEMAGALHAES_AGING_UP"), dn=gs("DEMAGALHAES_AGING_DN")))
rf <- file.path(D, "published", "rowland2019_37genes.csv")
if (file.exists(rf)) { rw <- read.csv(rf); pub$`Rowland 2019 kidney (37 genes)` <- list(up=rw$gene[rw$direction=="up"], dn=rw$gene[rw$direction=="down"]) }
pub$`This study (in-sample)` <- list(up=rob$gene[rob$direction=="up"], dn=rob$gene[rob$direction=="down"])
rows <- list()
for (nm in names(pub)) { s <- pub[[nm]]
  rr <- sapply(names(harm), function(g){ h <- harm[[g]]; v <- score_z(h$expr, s$up, s$dn); c(cor(v$score, h$meta$age, method="s"), ncol(h$expr), v$n_up[1]+v$n_dn[1]) })
  pr <- pool_rho(rr[1,], rr[2,]); gv <- score_z(gt$expr, s$up, s$dn); gs_ <- spear(gv$score, gt$meta$age)
  m <- res[res$gene %in% c(s$up, s$dn), ]; exp_dir <- ifelse(m$gene %in% s$up, 1, -1)
  rows[[nm]] <- data.frame(signature=nm, n_up=length(s$up), n_down=length(s$dn), measured_median=median(rr[3,]), pooled_rho=pr$rho, lo=pr$lo, hi=pr$hi, p=pr$p, I2=pr$I2,
    gtex_rho=gs_[1], gtex_p=gs_[2], overlap_with_signature=length(intersect(c(s$up, s$dn), rob$gene)),
    direction_agreement_in_meta_pct=100*mean(sign(m$z)==exp_dir), n_FDR05_in_meta=sum(m$fdr < 0.05)) }
ns <- pool_rho(Lc$rho, Lc$n)
rows$nested <- data.frame(signature="This study (nested held-out)", n_up=NA, n_down=NA, measured_median=NA, pooled_rho=ns$rho, lo=ns$lo, hi=ns$hi, p=ns$p, I2=ns$I2,
  gtex_rho=median(Lc$gtex_rho_fold_signature), gtex_p=NA, overlap_with_signature=NA, direction_agreement_in_meta_pct=NA, n_FDR05_in_meta=NA)
R <- do.call(rbind, rows); print(format(R, digits=3), row.names=FALSE); wr(R, "R1_published_signatures.csv")
