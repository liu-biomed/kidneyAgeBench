# Vector primitives in physical units; the page is the actual journal-scale plotting area.
library(grid)
ink <- "#15272D"; muted <- "#43585F"; gridcol <- "#D3DEDF"
base <- "#52616B"; adj <- "#00796F"; agecol <- "#B16030"; dxcol <- "#895039"
keys <- read.csv("source_data/signed_key.csv")
reps <- keys$sig_id[c(1,6,14)]
getdata <- function(name) read.csv(file.path("source_data",paste0(name,".csv")),check.names=FALSE)
mm <- function(x)unit(x,"mm")
# Horizontal anchor locations use the available two-column width.
# Marker diameters, text sizes and strokes remain in absolute paper units.
hx <- function(x)unit(x*170/180,"mm")
yy <- function(y)unit(1,"npc")-mm(page_y(y))
txt <- function(x,y,s,fs=8,just="left",bold=FALSE,col=ink,rot=0){
 grid.text(s,hx(x),yy(y),just=c(just,"centre"),rot=rot,gp=gpar(fontfamily="Arial",fontsize=fs,
          fontface=if(bold)"bold" else "plain",col=col,lineheight=.99))
}
line <- function(x,y,xe,ye,col=gridcol,lw=.65,lty=1){
 grid.segments(hx(x),yy(y),hx(xe),yy(ye),gp=gpar(col=col,lwd=lw,lty=lty,lineend="butt"))
}
rect <- function(x,y,w,h,fill,stroke=NA,lw=.45){
 top <- page_y(y); bottom <- page_y(y+h)
 grid.rect(hx(x+w/2),unit(1,"npc")-mm((top+bottom)/2),hx(w),mm(bottom-top),
           gp=gpar(fill=fill,col=stroke,lwd=lw))
}
point <- function(x,y,shape=21,diam=2,fill="white",col=base,lw=.75){
 grid.points(hx(x),yy(y),pch=shape,size=mm(diam*1.38),gp=gpar(fill=fill,col=col,lwd=lw))
}
panel <- function(letter,x,y)txt(x,y,letter,11,bold=TRUE)
map <- function(x,domain,area)area[1]+(x-domain[1])/diff(domain)*diff(area)
axisx <- function(area,y,domain,ticks,fs=7.3,lab=NULL,labelgap=3.1){
 line(area[1],y,area[2],y,muted)
 if(is.null(lab))lab <- format(ticks,trim=TRUE)
 for(i in seq_along(ticks)){x <- map(ticks[i],domain,area);line(x,y,x,y+.9,muted);txt(x,y+labelgap,lab[i],fs,"centre",col=muted)}
}
axisy <- function(x,area,domain,ticks,fs=7.3){
 line(x,area[1],x,area[2],muted)
 for(t in ticks){y <- map(t,domain,area);line(x-.85,y,x,y,muted);txt(x-1.6,y,format(t,trim=TRUE),fs,"right",col=muted)}
}
dvg <- grDevices::colorRamp(c("#396F99","#E0EBF1","#FBFAF6","#EBC9A6","#A94E29"),space="Lab")
heatcols <- function(v,lim){
 vals <- pmax(0,pmin(1,(v-lim[1])/diff(lim)))
 rgb(dvg(vals)/255)
}
bar_key <- function(x,y,w,lim,label,end=FALSE){
 n <- 100; vals <- seq(lim[1],lim[2],length.out=n)
 rect(x+(seq_len(n)-1)*w/n,y,w/n+.06,2.3,heatcols(vals,lim))
 txt(x-2,y+1,label,7.6,"right")
 breaks <- c(lim[1],0,lim[2]); labs <- format(breaks,trim=TRUE)
 if(end){labs[1] <- paste0("≤",lim[1]);labs[3] <- paste0("≥",lim[2])}
 for(i in 1:3)txt(map(breaks[i],lim,c(x,x+w)),y+5,labs[i],7,"centre")
}
forests <- function(d,ids,rangeX,domain,ys,metric,ticks,offset=.72,models=c("raw","adj"),styles=c(1,2),axisgap=3){
 xx <- function(a)map(a,domain,rangeX)
 stopifnot(min(d$lo)>=domain[1] && max(d$hi)<=domain[2])
 line(xx(0),ys[1]-2.4,xx(0),tail(ys,1)+2.4,"#ABBBBD",.65,2)
 for(k in seq_along(models)){
  from <- d[d$model==models[k],];set <- from[match(ids,from$sig_id),]
  y <- ys + if(length(models)==1)0 else c(-offset,offset)[k]
  st <- styles[k];color <- c(base,adj,agecol)[st]
  line(xx(set$lo),y,xx(set$hi),y,color,.9)
  point(xx(set[[metric]]),y,c(21,23,24)[st],if(length(models)>1)1.55 else 2.0,
        if(st==1)"white" else color,color,.7)
 }
 axisx(rangeX,tail(ys,1)+axisgap,domain,ticks)
}
legendpoint <- function(x,y,label,style){co <- c(base,adj,agecol)[style];point(x,y,c(21,23,24)[style],1.9,if(style==1)"white" else co,co);txt(x+2.7,y,label,7.5)}
forestlabels <- function(y,nm=keys$label,x=3){
 for(i in seq_along(y))txt(x,y[i],nm[i],7.6)
}
export <- function(path,height,fun){
 dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
 svglite::svglite(paste0(path,".svg"),width=170/25.4,height=height/25.4,bg="white",system_fonts=list(sans="Arial"))
 grid.newpage();fun();dev.off()
 stopifnot(system2("rsvg-convert",c("-f","pdf","-o",shQuote(paste0(path,".pdf")),shQuote(paste0(path,".svg"))))==0)
 stopifnot(system2("pdftoppm",c("-singlefile","-png","-r","300",shQuote(paste0(path,".pdf")),shQuote(path)))==0)
}

