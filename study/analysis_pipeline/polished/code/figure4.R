f4 <- function(){
 a <- getdata("F4_ercb_biopsies");b <- getdata("F4_gse175759_biopsies")
 panel("A",3,3.5);txt(10,3.5,"ERCB biopsies (n = 190)",8.6,bold=TRUE)
 bar_key(147,1.3,27,c(-4,4),"Living-donor SD",TRUE)
 groups <- unique(a$dx);ng <- as.integer(table(factor(a$dx,levels=groups)))
 x0 <- 38; w <- 138/nrow(a);x <- x0+w*(seq_len(nrow(a))-1)
 start <- c(0,cumsum(ng)[-length(ng)])
 # Direct labels sit over their true block centres. Two narrow groups use
 # short leaders to an upper tier; every number is the group sample count.
 for(k in seq_along(groups)){
   mid <- x0+w*(start[k]+ng[k]/2)
   raised <- groups[k] %in% c("TMD","FSGS/MCD")
   if(raised)line(mid,8.2,mid,11.8,muted,.55)
   txt(mid,if(raised)6.7 else 10.8,groups[k],7.2,"centre",TRUE)
   rect(x0+w*start[k],12.5,w*ng[k],2.9,if(k%%2)"#DDE6E6" else "#F3F5F4")
   txt(mid,13.95,ng[k],7.2,"centre")
 }
 txt(3,10.8,"Diagnosis",7.8);txt(3,13.95,"Biopsies",7.2)
 drawrow <- function(y,name,values,h){
   rect(x,y-h/2,w+.018,h,heatcols(values,c(-4,4)))
   txt(3,y,name,7.2)
 }
 tr <- c("score_injury","immune","stroma","endothelial","distal_tub","PT_content")
 labs <- c("Proximal tubular injury","Immune","Stromal","Endothelial","Distal tubular","Proximal tubular")
 for(i in 1:6)drawrow(16.9+(i-1)*2.52,labs[i],a[[tr[i]]],2.42)
 for(i in 1:15)drawrow(32.9+(i-1)*2.52,keys$label[i],a[[keys$sig_id[i]]],2.44)
 for(k in start[-1])line(x0+w*k,15.7,x0+w*k,69.4,"white",.72)
 # Each dot is one measured biopsy; no jitter or fitted trend is added.
 # B and C share one paper row, with the same three signatures in each cohort.
 triplets <- function(letter,data,xfield,xlim,xticks,xlabel,ylim,yticks){
   ercb <- letter=="B"
   origin <- if(ercb)3 else 93
   panel(letter,origin,73.7)
   txt(origin+7,73.7,if(ercb)"ERCB (n = 190)" else "GSE175759 (n = 87)",8.2,bold=TRUE)
   if(ercb){
     point(46,73.7,21,1.65,"white",base);txt(49,73.7,"Reference",7.4)
     point(72,73.7,16,1.65,dxcol,dxcol);txt(75,73.7,"CKD",7.4)
   }
   ya <- c(102,81.6)
   lefts <- if(ercb)c(13.5,39.5,65.5) else c(103.5,129.5,155.5)
   if(ercb)txt(2.9,mean(ya),"Score (ref. SD)",7.5,"centre",rot=90)
   for(i in seq_along(reps)){
     id <- reps[i];left <- lefts[i];ra <- c(left,left+20.2)
     val <- data[[id]];xx <- data[[xfield]]
     txt(mean(ra),77.1,keys$label[match(id,keys$sig_id)],7.6,"centre",TRUE)
     axisx(ra,ya[1],xlim,xticks)
     if(i==1)axisy(left,ya,ylim,yticks) else line(left,ya[1],left,ya[2],muted)
     for(g in c("CKD","Reference")){
       j <- data$group==g
       co <- if(g=="CKD")"#89503982" else "#233A46"
       point(map(xx[j],xlim,ra),map(val[j],ylim,ya),if(g=="CKD")16 else 21,1.1,if(g=="CKD")co else "white",co,.65)
     }
     txt(mean(ra),80.15,sprintf("ρ = %.2f",cor(xx,val,method="spearman")),7.2,"centre")
   }
   txt(mean(c(lefts[1],lefts[3]+20.2)),108.5,xlabel,7.6,"centre")
 }
 triplets("B",a,"score_injury",c(-2.2,7.1),c(-2,2,6),"Proximal tubular injury (living-donor SD)",c(-2.8,10.5),c(0,5,10))
 triplets("C",b,"egfr",c(37,144),c(40,80,120),"eGFR (mL/min/1.73 m²)",c(-2.2,4.65),c(-2,0,2,4))
 fit <- getdata("F4_ckd_intervals");er <- getdata("F4_egfr_intervals")
 ys <- seq(122.5,160.3,length.out=nrow(keys))
 panel("D",3,114);txt(10,114,"ERCB",8.3,bold=TRUE)
 panel("E",86,114);txt(93,114,"GSE175759",8.3,bold=TRUE)
 panel("F",134,114);txt(141,114,"eGFR · GSE175759",8,bold=TRUE)
 legendpoint(6,118.6,"Baseline",1);legendpoint(32,118.6,"+ Composition and injury",2)
 line(86,118.6,92,118.6,base,.95);txt(94,118.6,"Pointwise 95% CI (D–F)",7.5,col=muted)
 forestlabels(ys)
 forests(fit[fit$cohort=="ERCB",],keys$sig_id,c(38,80),c(-1.55,4.3),ys,"beta",c(-1,0,1,2,3,4))
 forests(fit[fit$cohort=="GSE175759",],keys$sig_id,c(86,128),c(-1.55,4.3),ys,"beta",c(-1,0,1,2,3,4))
 forests(er,keys$sig_id,c(134,176),c(-.7,.4),ys,"rho",c(-.6,-.3,0,.3))
 txt(83,169.8,"CKD − reference (reference SD)",7.8,"centre");txt(155,169.8,"Score–eGFR ρ",7.8,"centre")
}

