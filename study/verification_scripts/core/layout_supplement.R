# Editable page geometry in millimetres. Coordinates run from the top-left.
library(grid)
ink <- "#192C33"; subink <- "#455E65"; pale <- "#D5DEDE"; faint <- "#F2F6F5"
basecol <- "#485C69"; adjusted <- "#006E69"; accent <- "#A54D24"
classcol <- c(A="#325A80",B="#397669",C="#A5542C")
repnames <- c("Rodwell","Rowland replicated","SenMayo")
sdfile <- function(id) read.csv(file.path("source_data",paste0(id,".csv")),check.names=FALSE)
X <- function(x) unit(x,"mm")
Y <- function(y) unit(1,"npc")-unit(y,"mm")
tx <- function(x,y,label,size=8.5,bold=FALSE,col=ink,just="left",rot=0){
 grid.text(label,X(x),Y(y),just=c(just,"centre"),rot=rot,
           gp=gpar(fontfamily="Arial",fontsize=size,fontface=if(bold)"bold" else "plain",col=col,lineheight=1.03))
}
seg <- function(x,y,x2,y2,col=pale,lwd=.75,lty=1){
 grid.segments(X(x),Y(y),X(x2),Y(y2),gp=gpar(col=col,lwd=lwd,lty=lty,lineend="round"))
}
box <- function(x,y,w,h,fill="white",stroke=NA,lwd=.7){
 grid.rect(X(x+w/2),Y(y+h/2),X(w),X(h),gp=gpar(fill=fill,col=stroke,lwd=lwd))
}
dot <- function(x,y,style=1,diam=2.15,col=NULL){
 shape <- c(21,23,24)[style]; colour <- if(is.null(col)) c(basecol,adjusted,accent)[style] else col
 grid.points(X(x),Y(y),pch=shape,size=X(diam*1.46),gp=gpar(fill=if(style==1)"white" else colour,col=colour,lwd=.95))
}
tag <- function(l,x,y,s){tx(x,y,l,12,TRUE);tx(x+6,y,s,10,TRUE)}
headfig <- function(n,topic,title){
 tx(5,5.6,paste("FIGURE",n),8,TRUE,adjusted);tx(26,5.6,topic,8,TRUE,subink)
 tx(5,14.0,title,15.0,TRUE);seg(5,19,175,19,ink,.9)
}
scale_to <- function(x,domain,range) range[1]+diff(range)*(x-domain[1])/diff(domain)
axis_x <- function(range,y,domain,ticks,fmt=function(n) format(n,trim=TRUE),size=8){
 seg(range[1],y,range[2],y,subink,.65)
 for(v in ticks){x <- scale_to(v,domain,range);seg(x,y,x,y+1.1,subink,.65);tx(x,y+4.0,fmt(v),size,col=subink,just="centre")}
}
zero <- function(range,domain,top,bot,ref=0){x <- scale_to(ref,domain,range);seg(x,top,x,bot,"#7D9095",.75,"dashed")}
marklegend <- function(x,y,style,label){dot(x,y,style,2.2);tx(x+3.2,y,label,8,col=subink)}
rowlabels <- function(labels,ys,x=5){
 for(i in seq_along(labels))tx(x,ys[i],labels[i],8.2,labels[i] %in% repnames)
}
paired <- function(df,a,b,range,domain,ys){
 x1 <- scale_to(df[[a]],domain,range); x2 <- scale_to(df[[b]],domain,range)
 stopifnot(all(x1 >= range[1] & x1 <= range[2]),all(x2 >= range[1] & x2 <= range[2]))
 for(i in seq_along(ys)){seg(x1[i],ys[i],x2[i],ys[i],"#8A9B9C",1);dot(x1[i],ys[i],1,2.2);dot(x2[i],ys[i],2,2.0)}
}
legible_text <- function(fill){
 chan <- col2rgb(fill)/255
 lin <- ifelse(chan<=.04045,chan/12.92,((chan+.055)/1.055)^2.4)
 L <- sum(lin*c(.2126,.7152,.0722))
 if((L+.05)/.05 >= 1.05/(L+.05)) "#000000" else "#FFFFFF"
}
heat_colour <- function(v,domain,cols){
 ramp <- grDevices::colorRamp(cols, space="Lab")
 q <- max(0,min(1,scale_to(v,domain,c(0,1))))
 rgb(ramp(q),maxColorValue=255)
}
# Tubule cross-section and position markers are illustrative and carry no data.
tubule <- function(x,y,r=5,variant=FALSE,surrounding=variant){
 inner <- r*.43
 for(i in 1:8){
  a <- seq((i-1)*pi/4+.035,i*pi/4-.035,length.out=12)
  xx <- c(x+r*cos(a),x+inner*cos(rev(a))); yy <- c(y+r*sin(a),y+inner*sin(rev(a)))
  fill <- if(variant && i %in% c(2,3,4,7)) "#D9A781" else "#BBD6D1"
  grid.polygon(X(xx),Y(yy),gp=gpar(fill=fill,col=if(variant && i %in% c(2,3,4,7))"#996044" else "#557C77",lwd=.65))
 }
 if(surrounding){
  for(offset in list(c(-1,-1),c(.8,-1.32),c(1.2,.1)))grid.circle(X(x+r*offset[1]),Y(y+r*offset[2]),X(.65),gp=gpar(fill="#627080",col=NA))
 }
}
expfig <- function(prefix,height,draw){
 dir.create(dirname(prefix),showWarnings=FALSE,recursive=TRUE)
 svglite::svglite(paste0(prefix,".svg"),width=180/25.4,height=height/25.4,bg="white",
                  system_fonts=list(sans="Arial"))
 grid.newpage();draw();dev.off()
 convert <- Sys.which("rsvg-convert")
 stopifnot(nzchar(convert),system2(convert,c("-f","pdf","-o",shQuote(paste0(prefix,".pdf")),shQuote(paste0(prefix,".svg"))))==0)
 stopifnot(system2("pdftoppm",c("-singlefile","-png","-r","300",shQuote(paste0(prefix,".pdf")),shQuote(prefix)))==0)
}

