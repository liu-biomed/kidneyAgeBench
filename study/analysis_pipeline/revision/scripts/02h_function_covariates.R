# Age-adjusted functional associations inside discovery cohorts with per-donor age, using the nested (held-out)
# fold signature so that the tested cohort never contributed to gene selection; plus GSE147451 covariate/batch sensitivity
source(file.path(path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030"), "scripts", "fun.R"))
b <- load_bundle(); harm <- b$harm; P <- readRDS(MP); Lc <- readRDS(file.path(D, "harmonized", paste0("loco_", KSET, ".rds")))
hdr <- function(g){ con <- gzfile(file.path(D,"geo",paste0(g,"_series_matrix.txt.gz"))); L <- readLines(con, n=300); close(con)
  gsm <- gsub('"',"",strsplit(grep("^!Sample_geo_accession",L,value=TRUE),"\t")[[1]][-1]); ch <- grep("^!Sample_characteristics_ch1",L,value=TRUE)
  rec <- apply(t(sapply(ch, function(x) gsub('"',"",strsplit(x,"\t")[[1]][-1]))), 2, paste, collapse=" | ")
  list(gsm=gsm, get=function(k) setNames(sapply(rec, function(r){ m <- regmatches(r, regexec(paste0(k,": ([^|]*)"), r))[[1]]; if (length(m)) trimws(m[2]) else NA }), gsm)) }
out <- list()
## GSE60807: last donor creatinine (zero-hour biopsies), adjusted for age, sex, cause of death
h <- harm$GSE60807; a <- hdr("GSE60807"); cr <- as.numeric(a$get("last donor creatinine")[h$meta$gsm]); dth <- a$get("donor death")[h$meta$gsm]
for (sig in c("heldout","primary")) { s <- if (sig=="heldout") Lc$sigs$GSE60807 else list(up=P$rob$gene[P$rob$direction=="up"], dn=P$rob$gene[P$rob$direction=="down"])
  sc <- score_z(h$expr, s$up, s$dn)$score; f <- lm(scale(log(cr)) ~ scale(sc) + h$meta$age + h$meta$sex + dth)
  out[[length(out)+1]] <- data.frame(cohort="GSE60807", outcome="log last donor creatinine", signature=sig, n=sum(is.finite(cr)),
    rho_unadj=cor(sc, cr, method="s"), beta_adj_age_sex_death=coef(f)[2], p_adj=summary(f)$coef[2,4]) }
## GSE147451: post-transplant graft function group (good vs poor), adjusted for donor age, BMI, race (KDPI omitted: contains age)
h <- harm$GSE147451; km <- h$meta
for (sig in c("heldout","primary")) { s <- if (sig=="heldout") Lc$sigs$GSE147451 else list(up=P$rob$gene[P$rob$direction=="up"], dn=P$rob$gene[P$rob$direction=="down"])
  sc <- score_z(h$expr, s$up, s$dn)$score; y <- as.integer(km$gfr_group=="poor"); ok <- !is.na(y) & !is.na(km$bmi)
  f <- glm(y[ok] ~ scale(sc[ok]) + km$age[ok] + km$bmi[ok] + factor(km$race[ok]), family=binomial)
  out[[length(out)+1]] <- data.frame(cohort="GSE147451", outcome="poor post-transplant GFR group (OR per SD)", signature=sig, n=sum(ok),
    rho_unadj=cor(sc[ok], y[ok], method="s"), beta_adj_age_sex_death=exp(coef(f)[2]), p_adj=summary(f)$coef[2,4]) }
F <- do.call(rbind, out); print(F, row.names=FALSE); wr(F, "R1_function_age_adjusted.csv")
## GSE147451 age effects with additional adjustment for batch score, BMI and race (donor level), for the signature genes
if (!is.null(km$batch)) { ok <- is.finite(km$batch) & is.finite(km$bmi)
  e0 <- compute_effects(h$expr[, ok], km[ok, ], FALSE)
  X <- model.matrix(~ batch + bmi + factor(race), km[ok, ]); ar <- residuals(lm(km$age[ok] ~ X - 1)); Hm <- X %*% solve(crossprod(X)) %*% t(X)
  er <- h$expr[, ok] - t(Hm %*% t(h$expr[, ok])); r1 <- as.numeric((er/sqrt(rowSums(er^2))) %*% (ar/sqrt(sum(ar^2))))
  g <- intersect(P$rob$gene, rownames(h$expr)); i <- match(g, rownames(h$expr))
  cv <- data.frame(n=sum(ok), genes=length(g), same_sign=sum(sign(e0$r[i])==sign(r1[i])), cor_r=cor(e0$r[i], r1[i]), median_abs_ratio=median(abs(r1[i])/abs(e0$r[i])))
  print(cv); wr(cv, "R1_GSE147451_covariate_adjustment.csv") }
