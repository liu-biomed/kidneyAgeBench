# Post hoc revision analysis (2026-11-01): signed-score replication and composition attenuation re-pooled without
# CPTAC tumour-adjacent tissue (G6 random-effects estimate + TCGA), from the saved per-dataset estimates of script 21.
suppressMessages(library(metafor))
REV <- path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030")
O <- file.path(REV, "benchmark", "results"); MS <- read.csv(file.path(O, "M_signed_by_dataset.csv"))
leak <- function(id) if (grepl("^A1[12]_", id)) "G6" else if (grepl("^A0[56]_Rowland", id)) "TCGA" else character()
pool <- function(z, se) { if (length(z) == 1) return(c(z, z - 1.96 * se, z + 1.96 * se)); f <- rma(yi = z, sei = se, method = if (length(z) > 2) "REML" else "FE"); c(f$b, f$ci.lb, f$ci.ub) }
est <- function(d, lk) { parts <- NULL; g6 <- d[d$grp == "G6", ]
  if (!"G6" %in% lk && nrow(g6)) { p <- pool(atanh(g6$rho), 1 / sqrt(g6$n - g6$q - 3)); parts <- rbind(parts, data.frame(z = p[1], se = (p[3] - p[2]) / 3.92)) }
  tc <- d[d$grp == "TCGA", ]; if (!"TCGA" %in% lk && nrow(tc)) parts <- rbind(parts, data.frame(z = atanh(tc$rho), se = 1 / sqrt(tc$n - tc$q - 3)))
  p <- pool(parts$z, parts$se); tanh(p) }
R <- do.call(rbind, lapply(split(MS, MS$sig_id), function(d) { id <- d$sig_id[1]; lk <- leak(id)
  b <- est(d[d$model == "base", ], lk); cpn <- est(d[d$model == "composition", ], lk)
  data.frame(sig_id = id, rho_noCPTAC = b[1], lo_noCPTAC = b[2], hi_noCPTAC = b[3], rho_comp_noCPTAC = cpn[1], datasets = paste(setdiff(c("G6", "TCGA"), lk), collapse = "+")) }))
B <- read.csv(file.path(O, "BENCHMARK_summary.csv")); R <- merge(R, B[, c("sig_id", "class", "pooled_rho", "pooled_rho_comp")], by = "sig_id")
R$att_noCPTAC <- 1 - R$rho_comp_noCPTAC / R$rho_noCPTAC; R$att_primary <- 1 - R$pooled_rho_comp / R$pooled_rho
write.csv(R, file.path(O, "REVISION_signed_noCPTAC.csv"), row.names = FALSE); print(R[order(R$class, -R$rho_noCPTAC), ], digits = 2, row.names = FALSE)
for (k in c("A", "B")) { d <- R[R$class == k & R$sig_id != "A10_Sinning2023_MHCII", ]; cat(k, "range rho noCPTAC:", round(range(d$rho_noCPTAC), 2), "; CI excludes 0:", sum(d$lo_noCPTAC > 0), "of", nrow(d), "\n") }
cat("median attenuation without CPTAC:", round(median(R$att_noCPTAC), 2), "IQR", round(quantile(R$att_noCPTAC, c(.25, .75)), 2), "; primary:", round(median(R$att_primary), 2), "\n")
