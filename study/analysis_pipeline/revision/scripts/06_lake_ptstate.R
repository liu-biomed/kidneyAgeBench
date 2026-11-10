# S3: PT altered-state fraction vs age (Lake 2023 snRNA-seq), per locked plan; extra: cortex/PT-count restriction (reported as non-prespecified)
REV <- path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030")
d <- read.csv(path.expand("~/kidney_aging_data/validation/lake/lake_pt_state_by_donor.csv"))
sp <- function(x, y){ t <- suppressWarnings(cor.test(x, y, method="spearman", exact=FALSE)); n <- sum(is.finite(x) & is.finite(y)); z <- atanh(t$estimate); se <- 1/sqrt(n-3)
  data.frame(n=n, rho=unname(t$estimate), lo=tanh(z-1.96*se), hi=tanh(z+1.96*se), p=t$p.value) }
pr <- function(y, a, C){ ok <- complete.cases(C) & is.finite(y); n <- sum(ok); q <- ncol(model.matrix(~., C[ok,,drop=FALSE]))-1
  r <- cor(residuals(lm(rank(a[ok]) ~ ., C[ok,,drop=FALSE])), residuals(lm(rank(y[ok]) ~ ., C[ok,,drop=FALSE]))); z <- atanh(r); se <- 1/sqrt(n-q-3)
  data.frame(n=n, rho=r, lo=tanh(z-1.96*se), hi=tanh(z+1.96*se), p=2*pnorm(-abs(z/se))) }
nd <- d[d$disease=="normal",]; ok20 <- d$n_PT >= 100
out <- rbind(
  cbind(analysis="S3 primary: normal donors, altered PT fraction", prespecified=TRUE, sp(nd$age, nd$PT_altered_frac)),
  cbind(analysis="S3 sensitivity: all donors, disease-adjusted", prespecified=TRUE, pr(d$PT_altered_frac, d$age, data.frame(dis=factor(d$disease)))),
  cbind(analysis="normal donors, >=100 PT nuclei", prespecified=FALSE, sp(nd$age[nd$n_PT>=100], nd$PT_altered_frac[nd$n_PT>=100])),
  cbind(analysis="all donors >=100 PT nuclei, disease + cortex fraction adjusted", prespecified=FALSE, pr(d$PT_altered_frac[ok20], d$age[ok20], data.frame(dis=factor(d$disease[ok20]), ctx=d$cortex_frac[ok20]))),
  cbind(analysis="normal donors, aPT fraction", prespecified=FALSE, sp(nd$age, nd$PT_aPT_frac)),
  cbind(analysis="normal donors, dPT fraction", prespecified=FALSE, sp(nd$age, nd$PT_dPT_frac)),
  cbind(analysis="all donors, aPT fraction, disease-adjusted", prespecified=FALSE, pr(d$PT_aPT_frac, d$age, data.frame(dis=factor(d$disease)))))
print(format(out, digits=2), row.names=FALSE); write.csv(out, file.path(REV, "results", "validation", "Lake_PT_state_vs_age.csv"), row.names=FALSE)
