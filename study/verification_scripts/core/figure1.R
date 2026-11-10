f1 <- function(){
 # A common gene-set library is examined in three settings; no totals across settings.
 panel("A",3,4)
 f <- getdata("F1_reference_n");nref <- setNames(f$Freq,f$Var1)
 genes <- getdata("F2_48_sets");library_n <- table(genes$class)
 txt(11,4,paste(nrow(genes),"gene sets"),8.4,bold=TRUE)
 txt(3,10,sprintf("%d kidney age-derived · %d general ageing · %d disease-derived",
     library_n["A"],library_n["B"],library_n["C"]),7.4)
 er <- getdata("F4_ercb_biopsies");ck <- getdata("F4_gse175759_biopsies");dd <- getdata("F5_donors")
 branch_y <- c(17,33,44)
 line(6,12.5,6,tail(branch_y,1),muted,.8)
 for(y in branch_y){
   line(6,y,12,y,muted,.8);line(10.9,y-.9,12,y,muted,.8);line(10.9,y+.9,12,y,muted,.8)
 }
 txt(14,17,"Reference kidney",7.9,bold=TRUE)
 txt(14,22,sprintf("GEO %d · TCGA %d · CPTAC %d  (n = %d)",nref["G6"],nref["TCGA"],nref["CPTAC"],
     sum(nref[c("G6","TCGA","CPTAC")])),7.4)
 txt(14,26.2,sprintf("GTEx cortex %d¹",nref["GTEx"]),7.4)
 txt(14,33,"CKD and reference biopsies",7.9,bold=TRUE)
 txt(14,37.6,sprintf("ERCB %d · GSE175759 %d",nrow(er),nrow(ck)),7.4)
 txt(14,44,"Transplant donors",7.9,bold=TRUE)
 txt(14,48.6,sprintf("GSE147451 · %d donor records²",nrow(dd)),7.4)
 panel("B",88,4)
 z <- getdata("F1_folds");co <- z$heldout[z$signature=="KAS-94"]
 # Eligible cohorts first; excluded discovery series last.
 keep <- z$heldout[z$signature=="KAS-82"];co <- c(keep,setdiff(co,keep))
 xx <- c(112.5,175);dom <- c(-.25,.80)
 legendpoint(102,4,"KAS-94",1);point(136,4,22,1.9,ink,ink);txt(139,4,"KAS-82",7.5)
 yb <- seq(9.7,40.5,length.out=length(co));line(map(0,dom,xx),8,map(0,dom,xx),42,gridcol,.75,2)
 cut <- mean(yb[c(length(keep),length(keep)+1)])
 line(89,cut,175,cut,gridcol,.5,3)
 for(i in seq_along(co)){
  v <- z[z$heldout==co[i],];txt(89,yb[i],co[i],7.6)
  for(k in seq_len(nrow(v))){
   first <- v$signature[k]=="KAS-94"
   point(map(v$rho[k],dom,xx),yb[i]+if(first)-.6 else .6,if(first)21 else 22,1.75,if(first)"white" else ink,base)
  }
 }
 axisx(xx,43,dom,c(-.2,0,.2,.4,.6,.8));txt(mean(xx),50,"Held-out score–age ρ",7.8,"centre")
 panel("C",3,57)
 q <- getdata("F1_pooled_KAS");xr <- c(19.5,83);dom <- c(0,.5)
 legendpoint(13,57,"In sample",1);legendpoint(41,57,"Nested held out",2)
 for(i in 1:2){
  nm <- c("KAS-94","KAS-82")[i];y <- 68+(i-1)*12
  txt(3,y,nm,7.8)
  for(k in 1:2){
   v <- q[q$signature==nm & grepl(c("in-sample","nested")[k],q$analysis),]
   yk <- y+c(-1.65,1.65)[k];cl <- c(base,adj)[k]
   line(map(v$lo,dom,xr),yk,map(v$hi,dom,xr),yk,cl,.95)
   point(map(v$rho,dom,xr),yk,c(21,23)[k],2,if(k==1)"white" else cl,cl)
  }
 }
 axisx(xr,88,dom,seq(0,.5,.1));txt(mean(xr),95,"Score–age ρ (95% CI)",7.8,"centre")
 panel("D",88,57);txt(96,57,"GTEx cortex (n = 85)",8.5,bold=TRUE)
 g <- getdata("F1_GTEx");e <- getdata("F1_GTEx_estimate")
 txt(96,62,sprintf("ρ = %.2f (95%% CI %.2f to %.2f)",e$spearman_rho,e$boot_lo,e$boot_hi),7.5)
 xx <- c(103,175);ya <- c(88,65);xd <- c(19,81);yd <- c(-2,1.6)
 axisx(xx,ya[1],xd,seq(25,75,10));axisy(xx[1],ya,yd,c(-2,-1,0,1))
 point(map(g$age,xd,xx),map(g$score,yd,ya),16,1.65,"#546E8890","#546E8890")
 for(i in unique(g$age)){x <- map(i,xd,xx);y <- map(median(g$score[g$age==i]),yd,ya);line(x-1.8,y,x+1.8,y,ink,1.1)}
 txt(91,77,"KAS-94 score",7.4,"centre",rot=90)
 txt(mean(xx),95,"Age-bin midpoint (years)",7.8,"centre")
 panel("E",3,105)
 perm <- getdata("F1_sex_perm")[[1]];counts <- tabulate(perm+1,94)
 ra <- c(18,83);ya <- c(137,108);xd <- c(0,94);yd <- c(0,260)
 axisx(ra,ya[1],xd,c(0,20,40,60,80));axisy(ra[1],ya,yd,c(0,100,200))
 rect(map(0:93-.48,xd,ra),map(counts,yd,ya),diff(ra)/94,ya[1]-map(counts,yd,ya),"#8CA3AE")
 observed <- getdata("F1_concordance_observed");xo <- map(observed$concordant,xd,ra);line(xo,ya[2],xo,ya[1],agecol,1.2)
 txt(xo-1.1,106.5,sprintf("Observed %d/%d · P = %.3f",observed$concordant,observed$n_genes,(1+sum(perm>=observed$concordant))/(length(perm)+1)),7.2,"right",col=agecol)
 txt(4,122,"Permutations",7.5,"centre",rot=90);txt(mean(ra),144,"Concordant genes",7.8,"centre")
 panel("F",88,105)
 p <- getdata("F1_background");p <- p[match(c("Class A, external sets","Class C","GenAge","CellAge"),p$pool),]
 xx <- c(127.5,175);dom <- c(.8,1.58);ys <- seq(113,132,length.out=4)
 legendpoint(98,105,"Genome-wide",1);legendpoint(135,105,"Expression-matched",2)
 line(map(1,dom,xx),110,map(1,dom,xx),134,gridcol,.8,2)
 for(i in 1:4){
  txt(89,ys[i],c("Kidney age-derived","Disease-derived","GenAge","CellAge")[i],7.5)
  point(map(p$naive_ratio[i],dom,xx),ys[i]-1,21,1.8,"white",base)
  y <- ys[i]+1;line(map(p$lo[i],dom,xx),y,map(p$hi[i],dom,xx),y,adj,1)
  point(map(p$matched_ratio[i],dom,xx),y,23,1.8,adj,adj)
 }
 axisx(xx,137,dom,c(.8,1,1.2,1.4));txt(mean(xx),144,"Enrichment ratio",7.8,"centre")
}

