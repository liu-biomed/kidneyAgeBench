f2 <- function(){
 pal <- c(A="#315D85",B="#347667",C="#A8532A"); pchs <- c(A=21,B=23,C=24)
 panel("A",3,3.6)
 d <- getdata("F2_48_sets")
 # Each class has a visible start. Class C continues in ID order in the right column.
 txt(100,3.6,"Filled: BH FDR < 0.05",7.3);point(97,3.6,21,1.8,ink,ink)
 txt(139,3.6,"Open: BH FDR ≥ 0.05",7.3);point(136,3.6,21,1.8,"white",ink)
 a <- d[d$class=="A",];b <- d[d$class=="B",];c_sets <- d[d$class=="C",]
 chunks <- list(list(a,b,head(c_sets,4)),list(tail(c_sets,-4)))
 class_names <- c(A="Kidney age-derived",B="General ageing",C="Disease-derived")
 for(half in seq_along(chunks)){
   at <- (half-1)*89; rr <- at+c(47.5,85.1);dom <- c(0,2.35)
   txt(at+3,8.6,"ID",7.3,bold=TRUE);txt(at+11,8.6,"Set",7.3,bold=TRUE)
   txt(at+41,8.4,"Measured\ngenes",7.2,"centre",TRUE)
   txt(mean(rr),8.4,"Matched enrichment\nratio",7.2,"centre",TRUE)
   line(map(1,dom,rr),12,map(1,dom,rr),81.0,muted,.72,2)
   y <- 13
   for(j in seq_along(chunks[[half]])){
     chunk <- chunks[[half]][[j]]
     if(j>1) y <- y+1.20
     cl <- chunk$class[1]
     header <- if(half==2) "Disease-derived (continued)" else class_names[[cl]]
     txt(at+3,y,header,7.4,bold=TRUE,col=pal[cl]);y <- y+2.65
     for(i in seq_len(nrow(chunk))){
       txt(at+3,y,chunk$code[i],7.2,col=pal[cl]);txt(at+11,y,chunk$short_label[i],7.2)
       txt(at+41,y,chunk$n_measured[i],7.2,"centre",col=muted)
       fill <- if(chunk$fdr_all[i]<.05)pal[cl] else "white"
       point(map(chunk$enrichment[i],dom,rr),y,pchs[cl],1.85,fill,pal[cl])
       y <- y+if(half==1)2.59 else 2.69
     }
   }
   axisx(rr,82.7,dom,c(0,1,2),labelgap=2.1)
 }

 panel("B",3,88)
 panel("C",84.9,87.4)
 b <- getdata("F2_F3_signed_intervals");b <- b[b$model=="base",]
 rr <- c(39,81);dom <- c(-.03,.44);ys <- seq(94,133,length.out=15)
 forestlabels(ys)
 forests(b,keys$sig_id,rr,dom,ys,"rho",c(0,.1,.2,.3,.4),models="base",styles=1)
 txt(mean(rr),143,"Score–age ρ (95% CI)",7.6,"centre")
 pc <- getdata("F2_signed_cohorts");ds <- c("G6_GSE36895","G6_GSE60807","G6_GSE147451","G6_GSE66272","G6_GSE93069","G6_GSE167093","TCGA","CPTAC","GTEx")
 xs <- seq(93.3,171.7,length.out=9)
 for(j in 1:9){
  txt(xs[j],89,sub("G6_GSE","",ds[j]),7.0,"centre",rot=30)
  for(i in 1:15){
   row <- pc[pc$sig_id==keys$sig_id[i] & pc$dataset==ds[j],]
   valid <- row$use=="Observed"
   rect(xs[j]-4.7,ys[i]-1.30,9.4,2.60,if(valid)heatcols(row$rho,c(-.6,.6)) else "#EDF0F0")
   if(!valid)txt(xs[j],ys[i],if(row$use=="Derivation")"×" else "–",8.1,"centre")
  }
 }
 bar_key(105,137,23,c(-.6,.6),"Score–age ρ");txt(132,138.3,"×  Derivation",7.1);txt(152,138.3,"–  Unavailable",7.1)
 panel("D",3,147)
 panel("E",57,147)
 panel("F",122,147)
 v <- getdata("F2_VIF");sb <- tapply(v$fdr<.05,v$class,sum);sa <- tapply(v$fdr_vif<.05,v$class,sum)
 txt(11,147,"FDR < 0.05",7.4);txt(34,149,"Sensitivity\nbaseline",7.2,"centre",TRUE);txt(48,149,"VIF",7.2,"centre",TRUE)
 for(i in 1:3){y <- 155+(i-1)*4.5;txt(4,y,c("Kidney","General","Disease")[i],7.4)
  txt(34,y,sprintf("%d/%d",sb[i],c(12,7,29)[i]),7.3,"centre")
  txt(48,y,sprintf("%d/%d",sa[i],c(12,7,29)[i]),7.3,"centre")}
 txt(4,171,"Retained: KAS-94",7.3,col=muted)
 pow <- getdata("F2_power");px <- c(67,113);py <- c(165,152.7)
 xd <- log10(c(2,100));axisx(px,py[1],xd,log10(c(2,5,13,30,100)),lab=c("2","5","13","30","100"));axisy(px[1],py,c(0,1),c(0,.5,1))
 txt(59,mean(py),"Power",7.2,"centre",rot=90)
 txt(65,147,"True enrichment ratio",7.2)
 ers <- c(1,1.2,1.3,1.5);cols <- c("#7D888B","#587493","#00796F","#A55328")
 for(i in 1:4){row <- pow[pow$true_ER==ers[i],];row <- row[order(row$k),]
   x <- map(log10(row$k),xd,px); y <- map(row$power,c(0,1),py)
   line(head(x,-1),head(y,-1),tail(x,-1),tail(y,-1),cols[i],1.15,c(3,2,1,4)[i])
   xx <- 61+(i-1)*13.9
   line(xx,150.2,xx+3.5,150.2,cols[i],1,c(3,2,1,4)[i]);txt(xx+4.3,150.2,format(ers[i],trim=TRUE),7.05)
 }
 txt(mean(px),171.7,"Measured genes (log scale)",7.3,"centre")
 m <- getdata("F2_without_CPTAC");m <- subset(m,is.finite(enrichment_noCPTAC));sx <- c(137,175);sy <- c(165,151.4);dm <- c(0,2.6)
 axisx(sx,sy[1],dm,0:2);axisy(sx[1],sy,dm,0:2)
 line(sx[1],sy[1],sx[2],sy[2],muted,.7,2)
 for(cl in c("C","A","B")){sel <- m$class==cl;point(map(m$enrichment_primary[sel],dm,sx),map(m$enrichment_noCPTAC[sel],dm,sy),pchs[cl],1.3,grDevices::adjustcolor(pal[cl],alpha.f=.25),pal[cl],.65)}
 txt(127.2,mean(sy),"Enrichment\nwithout CPTAC",7.1,"centre",rot=90)
 txt(mean(sx),171.7,"Primary enrichment",7.3,"centre")
 txt(130,147,sprintf("46 sets · ρ = %.2f",cor(m$enrichment_primary,m$enrichment_noCPTAC,method="spearman")),7.3)
}

