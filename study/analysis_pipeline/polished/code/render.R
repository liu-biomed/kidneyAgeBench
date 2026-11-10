# The five layouts share a 170 mm width and physical font/marker sizes.
# Vertical anchors are panel or chart boundaries. Numeric y scales remain
# linear within each measured scatterplot, histogram or power curve.
source("code/draw.R")
bands <- list(
  # Remove two blank inter-row bands in F1; keep the height available to data.
  list(old=c(0,4,50,57,95,105,144,148),
       new=c(0,3.8,55.4,62.4,105.8,112.8,156.5,160)),
  # F2A gets two explicit inter-class breaks and a common two-column baseline.
  list(old=c(0,3.6,13,82.7,88,94,133,143,147,151,165,171.7,175),
       new=c(0,4.5,15.2,98.4,107.5,114.8,161.3,172.4,180,185,204.2,212.5,216.2)),
  list(old=c(0,3.5,13,28,35,41,57,64.2,71.5,77.5,117,127,133,140,155,162,165),
       new=c(0,4.1,15.1,32.7,40.5,47.6,68.4,76.4,83.4,90.9,139.8,151.6,158.4,167.2,187.8,195.8,199.6)),
  list(old=c(0,3.5,16.9,69.4,73.7,81.6,102,108.5,114,122.5,160.3,169.8,172.5),
       new=c(0,4.5,19.1,85.5,90.8,99.9,125.9,133.8,140.5,150.1,198.6,209.2,213)),
  list(old=c(0,3.5,22.5,29.5,36,44.5,65.5,72.6,83.2,87.7,131.8,142,148,151),
       new=c(0,4.1,25.2,33,40.3,50.3,76.3,84.2,92.4,99.2,151.5,164.2,171,175))
)

figs <- as.integer(commandArgs(trailingOnly=TRUE))
if(!length(figs))figs <- 1:5
stopifnot(all(figs %in% 1:5))
for(i in figs){
  layout <- bands[[i]]
  stopifnot(length(layout$old)==length(layout$new), all(diff(layout$old)>0),
            all(diff(layout$new)>0),tail(layout$new,1)<=225)
  page_y <- approxfun(layout$old,layout$new)
  source(sprintf("code/figure%d.R",i))
  export(sprintf("figures/Figure%d",i),tail(layout$new,1),get(sprintf("f%d",i)))
}
