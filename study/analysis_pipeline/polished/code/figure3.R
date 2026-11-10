f3 <- function(){
 panel("A",3,3.5)
 bar_key(25.2,3,10.8,c(-.5,.5),"Marker–age\nρ")
 p <- getdata("F3_age_proxies");coh <- c("GSE36895","GSE60807","GSE147451","GSE66272","GSE93069","GSE167093","TCGA","CPTAC","GTEx")
 xs <- seq(43,168,length.out=9)
 for(j in 1:9){txt(xs[j],3.5,sub("GSE","",coh[j]),7.5,"centre",TRUE);txt(xs[j],8,paste0("n=",unique(p$n[p$ds==coh[j]])),7.1,"centre")}
 labs <- c("Immune","Stroma","Endothelial","Distal tubule","Proximal tubule")
 for(i in 1:5){y <- 13.0+(i-1)*3.75;txt(3,y,labs[i],7.7)
  for(j in 1:9){
   v <- p$rho[p$ds==coh[j]&p$proxy==labs[i]]
   rect(xs[j]-7.45,y-1.76,14.9,3.52,heatcols(v,c(-.5,.5)))
   txt(xs[j],y,sprintf("%.2f",v),7.4,"centre",col=if(v < -.38)"white" else "#15272D")
  }}
 s <- getdata("F3_donor_samples")
 panel("B",3,35);txt(11,35,"GSE147451 donors (n = 166)",8.4,bold=TRUE)
 for(i in 1:2){
  var <- c("immune","PT_content")[i];xx <- c(18.5,84)+88*(i-1);ya <- c(57,41)
  val <- s[[var]];xdom <- c(16,76);ydom <- if(i==1)c(-1.7,4.5) else c(-3.4,1.4)
  axisx(xx,ya[1],xdom,c(20,40,60));axisy(xx[1],ya,ydom,if(i==1)c(-1,1,3) else c(-3,-1,1))
  txt(mean(xx),38.2,sprintf("%s  (ρ = %.2f)",c("Immune proxy","Proximal tubular proxy")[i],cor(s$age,val,method="spearman")),7.7,"centre",TRUE)
  point(map(s$age,xdom,xx),map(val,ydom,ya),16,1.45,"#4F788578","#4F788578")
  
 }
 txt(90,64.2,"Donor age (years)",7.5,"centre")
 txt(5,49,"Marker score",7.4,"centre",rot=90)
 panel("C",3,71.5)
 panel("D",91,71.5)
 legendpoint(13,71.5,"Baseline",1);legendpoint(41,71.5,"+ Composition",2)
 txt(99,71.5,"Whiskers: matched-null range",7.5,col=muted)
 z <- getdata("F2_F3_signed_intervals");ys <- seq(77.5,117,length.out=15)
 forestlabels(ys);forests(z,keys$sig_id,c(38,84),c(-.16,.44),ys,"rho",c(-.1,0,.1,.2,.3,.4),models=c("base","composition"))
 txt(61,127,"Score–age ρ (95% CI)",7.5,"centre")
 q <- getdata("F3_pools");u <- unique(q$pool)[c(1,4,5,6,2,3)]
 short <- c("Kidney age-derived","Disease-derived","CellAge","GenAge",keys$label[keys$sig_id=="A09_Perco2019_RAAG_CKD"],"Other RAAGs")
 rr <- c(126,175);dom <- c(.83,2.22);sy <- seq(78,113,length.out=6)
 line(map(1,dom,rr),75,map(1,dom,rr),119,gridcol,.85,2)
 for(i in 1:6){
  a <- q[q$pool==u[i],]
  txt(92,sy[i],short[i],7.45);txt(92,sy[i]+3.5,paste0("n = ",a$n[1]," genes"),7,col=muted)
  for(k in 1:2){
   y <- sy[i]+c(-1,1)[k];color <- c(base,adj)[k]
   line(map(a$lo[k],dom,rr),y,map(a$hi[k],dom,rr),y,color,.9)
   point(map(a$ER[k],dom,rr),y,c(21,23)[k],1.75,if(k==1)"white" else color,color)
  }
 }
 axisx(rr,120,dom,c(1,1.4,1.8,2.2));txt(mean(rr),127,"Enrichment ratio",7.5,"centre")
 panel("E",3,133)
 rhos <- getdata("F3E_descriptive_rho")
 for(i in 1:3){
  xr <- c(15,44)+(i-1)*34.6;yr <- c(155,140)
  xd <- c(-1.6,4.3);yd <- c(-1.6,3.0);v <- s[[reps[i]]]
  txt(mean(xr),133,c("Rodwell","Rowland replicated","SenMayo")[i],7.3,"centre",TRUE)
  txt(mean(xr),137,sprintf("ρ = %.2f",rhos$rho[match(reps[i],rhos$sig_id)]),7.3,"centre")
  axisx(xr,yr[1],xd,c(-1,1,3));axisy(xr[1],yr,yd,c(-1,1,3))
  point(map(s$immune,xd,xr),map(v,yd,yr),16,1.25,"#4F788578","#4F788578")
 }
 txt(4,148,"Signed score",7.2,"centre",rot=90);txt(64.1,162,"Immune marker score (GSE147451)",7.5,"centre")
 panel("F",116,133)
 # Six illustrated cells per field; type is encoded by outline/shape.
 # Only the intracellular dots vary between the reference and state examples.
 tubular <- function(x,y,marked=FALSE){
   grid.roundrect(hx(x),yy(y),mm(4.3),mm(3.8),r=mm(.55),
     gp=gpar(fill="#CEE1DD",col="#4A7773",lwd=.7))
   if(marked){
     for(k in 1:3){
       dx <- c(-.95,.85,-.05)[k];dy <- c(.62,.6,-.62)[k]
       grid.circle(hx(x+dx),yy(y)-mm(dy),mm(.40),gp=gpar(fill="#9B4C28",col=NA))
     }
   }else grid.circle(hx(x),yy(y),mm(.4),gp=gpar(fill="#9B4C28",col=NA))
 }
 interstitial <- function(x,y){
   grid.polygon(hx(x+c(-2.4,0,2.4,0)),yy(y)-mm(c(0,-1.45,0,1.45)),
       gp=gpar(fill="#CBB9AA",col="#8C7464",lwd=.7))
 }
 immune <- function(x,y){
   grid.circle(hx(x),yy(y),mm(1.48),gp=gpar(fill="#8493A3",col="#52657A",lwd=.7))
 }
 examples <- list(c("T","T","T","T","S","I"),c("T","I","T","I","S","I"),c("T","T","T","T","S","I"))
 cx <- c(126,146,166)
 for(j in 1:3){
   for(k in 1:6){
     x <- cx[j]+c(-5,0,5)[(k-1)%%3+1]; y <- if(k<=3)137 else 141.5
     switch(examples[[j]][k],T=tubular(x,y,marked=j==3),S=interstitial(x,y),I=immune(x,y))
   }
   txt(cx[j],145.8,c("Reference","Cell mix","Cell state")[j],7.3,"centre")
 }
 tubular(119,151.5);txt(123,151.5,"Tubular epithelium",7.2)
 interstitial(150,151.5);txt(154,151.5,"Interstitial cell",7.2)
 immune(119,156.9);txt(123,156.9,"Immune cell",7.2)
 grid.circle(hx(150),yy(156.9),mm(.46),gp=gpar(fill="#9B4C28",col=NA))
 txt(154,156.9,"State marker",7.2)
}
