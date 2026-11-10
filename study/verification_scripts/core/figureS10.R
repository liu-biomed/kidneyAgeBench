source("code/layout_supplement.R")
s10 <- function(){
 headfig("S10","SENSITIVITY ANALYSES","Gene pools and deconvolution checks")
 tag("A",5,25,"Enrichment relative to expression-matched genes")
 tx(5,32,"Entries show the observed ratio and its matched-null range. Age results use the primary configuration.",8.2,col=subink)
 q <- read.csv("source_data/S10_pools.csv")
 nm <- c("Kidney age-derived","CKD-prog. RAAGs","RAAG, other","Disease-derived","CellAge","GenAge")
 cols <- c(58,100,149)
 for(i in 1:3)tx(cols[i],42,c("Age, baseline","Age, composition-adjusted","CKD versus reference")[i],8,TRUE,just="centre")
 analyses <- c("Age (unadjusted, locked)","Age (composition-adjusted, locked)","CKD vs reference (GSE104954 + GSE175759)")
 for(j in 1:6){
  y <- 50+(j-1)*7; if(j%%2==1)box(5,y-3.3,170,6.6,faint);tx(5,y,nm[j],8)
  for(i in 1:3){z <- q[q$analysis==analyses[i],][j,];tx(cols[i],y,sprintf("%.2f (%.2f–%.2f)",z$ER,z$lo,z$hi),8.1,just="centre")}
 }
 tx(5,94,"Ranges divide the observed statistic by the 97.5th and 2.5th matched-null percentiles.",8.1,col=subink)
 tx(5,100,"CKD–reference comparisons pool ERCB and GSE175759. Pool sizes are specific to each comparison.",8.1,col=subink)
 seg(5,105,175,105)
 tag("B",5,112,"Marker adjustment compared with reference-based deconvolution")
 q <- read.csv("source_data/S10_deconvolution.csv");dom <- c(-.06,.30);ploty <- c(170,130)
 for(i in 1:2){
  r <- c(25,75)+82*(i-1); v <- q[[c("G6_KPMP_MuSiC","G6_KPMP_Bisque_5sub")[i]]]
  tx(mean(r),123,sprintf("%s  ·  ρ = %.2f",c("MuSiC","Bisque, 5 subclasses")[i],cor(q$G6_marker5,v,method="spearman")),8.4,TRUE,just="centre")
  seg(r[1],ploty[1],r[2],ploty[2],subink,.6,"dashed")
  seg(r[1],ploty[1],r[1],ploty[2],subink,.7)
  for(t in 0:3/10){y <- scale_to(t,dom,ploty);tx(r[1]-1.8,y,sprintf("%.1f",t),7.6,just="right")}
  axis_x(r,ploty[1],dom,0:3/10,function(x)sprintf("%.1f",x))
  for(j in seq_len(nrow(q)))dot(scale_to(q$G6_marker5[j],dom,r),scale_to(v[j],dom,ploty),1,2)
 }
 tx(6,150,"Deconvolution-adjusted age ρ",7.8,rot=90,just="centre")
 tx(90,181,"Marker-adjusted age ρ",8.3,just="centre")
 tx(5,190,"Each point is a signature (n = 13). The five-cohort comparison uses the same profiles for each method.",8.2)
 tx(5,196,"414 profiles contributed for 12 signatures; 261 for tubular MHC II. Weights use actual n and model rank.",8.2)
}
expfig("Supplementary/FigureS10",203,s10)

