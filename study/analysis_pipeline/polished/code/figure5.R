f5 <- function(){
 d <- getdata("F5_donors");p <- getdata("F5_donor_linkage")
 panel("A",3,3.5)
 # Horizontal path represents the source unit conversions, not a clinical mechanism.
 xs <- c(10,36,66)
 words <- c(sprintf("%d kidney\nprofiles",sum(p$kidneys)),sprintf("%d inferred\ndonors",nrow(p)),sprintf("%d with\n24-mo. group",nrow(d)))
 for(i in 1:3){
  txt(xs[i],6,words[i],7.7,bold=TRUE)
  if(i<3){x <- c(28,57)[i];line(x,6,x+5,6,muted,.85);line(x+3.7,5,x+5,6,muted,.85);line(x+3.7,7,x+5,6,muted,.85)}
 }
 txt(10,14.5,"11 inferred pairs · 1 discordant pair",7.5)
 txt(10,19.5,sprintf("%d without a 24-month group",nrow(p)-nrow(d)),7.5)
 panel("B",89,3.5)
 xx <- c(107,176);dom <- c(16,76)
 # Stacks are within exact age and outcome; coordinates never jitter along the age axis.
 for(j in 1:2){
  g <- c("good","poor")[j];v <- d$age[d$outcome==g];y <- c(6.5,16.3)[j]
  txt(91,y,paste0(if(j==1)"Higher" else "Lower"," (",length(v),")"),7.4)
  for(age in unique(v)){
   n <- sum(v==age);ys <- y+((1:n)-(n+1)/2)*1.12
   point(rep(map(age,dom,xx),length(ys)),ys,if(j==1)21 else 24,.85,if(j==1)"white" else dxcol,if(j==1)base else dxcol,.65)
  }
 }
 axisx(xx,22.5,dom,c(20,40,60));txt(mean(xx),29.5,"Donor age (years)",7.5,"centre")
 panel("C",3,36)
 point(12,36,21,1.7,"white",base);txt(15,36,"Higher function",7.5)
 point(42,36,24,1.7,dxcol,dxcol);txt(45,36,"Lower function",7.5)
 for(i in 1:3){
  xr <- c(18,59)+(i-1)*57;yr <- c(65.5,44.5); v <- d[[reps[i]]]
  txt(mean(xr),40.5,keys$label[match(reps[i],keys$sig_id)],7.7,"centre",TRUE)
  axisx(xr,yr[1],dom,c(20,40,60));axisy(xr[1],yr,c(-1.7,2.9),c(-1,0,1,2))
  for(j in 1:2){k <- d$outcome==c("good","poor")[j]
   point(map(d$age[k],dom,xr),map(v[k],c(-1.7,2.9),yr),c(21,24)[j],1.08,
         if(j==1)"white" else "#89503956",if(j==1)"#293D4C85" else "#89503985",.65)}
 }
 txt(4.5,55,"Signed score",7.3,"centre",rot=90);txt(91,72.6,"Donor age (years)",7.7,"centre")
 s <- getdata("F5_outcome_intervals")
 ys <- seq(87.7,131.8,length.out=15)
 starts <- c(3,86,134); areas <- list(c(38,80),c(86,128),c(134,176))
 model_labels <- c("Unadjusted","Age-adjusted","Age + composition")
 for(j in 1:3){
   panel(c("D","E","F")[j],starts[j],83.2)
   legendpoint(starts[j]+8,83.2,model_labels[j],c(1,3,2)[j])
   md <- c("raw","age","full")[j]
   forests(s[s$model==md,],keys$sig_id,areas[[j]],c(-.26,.46),ys,"rho",c(-.2,0,.2,.4),models=md,styles=c(1,3,2)[j])
 }
 forestlabels(ys,ifelse(grepl("^KAS",keys$label),paste0(keys$label,"*"),keys$label))
 txt(89,142,"Correlation with lower 24-month graft-function group, ρ (95% CI)",7.8,"centre")
 txt(3,148,"* GSE147451 derivation overlap",7.3,col=muted)
}

