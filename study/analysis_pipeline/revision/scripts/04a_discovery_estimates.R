# Discovery-set estimates for the locked definitions (corrected 6-cohort set), pooled by random effects
source(file.path(path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030"), "scripts", "fun.R"))
source(file.path(REV, "scripts", "04_locked_definitions.R"))
b <- load_bundle(); harm <- b$harm
R <- do.call(rbind, lapply(names(harm), function(g){ h <- harm[[g]]; base <- if (g %in% NO_SEX) NULL else data.frame(sex=factor(h$meta$sex)); cbind(cohort=g, run_locked(h$expr, h$meta$age, base)) }))
pl <- do.call(rbind, lapply(split(R, paste(R$outcome, R$model)), function(d){ d <- d[is.finite(d$rho),]
  f <- rma(yi=atanh(d$rho), sei=1/sqrt(d$n - d$q - 3), method="REML"); data.frame(outcome=d$outcome[1], model=d$model[1], k=nrow(d), n=sum(d$n), rho=tanh(f$b), lo=tanh(f$ci.lb), hi=tanh(f$ci.ub), p=f$pval, I2=f$I2) }))
pl <- pl[order(pl$outcome, pl$model),]; print(format(pl, digits=2), row.names=FALSE)
write.csv(R, file.path(REV, "prereg", "discovery_estimates_by_cohort.csv"), row.names=FALSE); write.csv(pl, file.path(REV, "prereg", "discovery_estimates_pooled.csv"), row.names=FALSE)
