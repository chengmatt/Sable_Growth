sbcomps_pred<-function(data.o,data.p,data.p2,len,surv,title,nr,nc,new=T,maxy=NULL){

#Setting up for new graph parameters, typically if new dataset, then new=T
if(new){
par.old<-par(mfcol=c(nr,nc),mai=c(0.5, 0.25, 0.25, 0.25), omi=c(0.25,0.25,0.25,0))
}

#Assigning row and column lengths of dataset, and specifying number of graphsheets
n.s<-length(data.o[1,])
n.y<-length(data.o[,1])
cells <- nr*nc
if(n.y - trunc(n.y/cells)*cells!=0)
  {npages <- trunc(n.y/(nr*nc))+1 }
else
  {npages <- trunc(n.y/(nr*nc))}

#Assigning x-axis labels
if(len){
  labx<-"Length(cm)"}
 else
  { labx<-"Age"}

#Assigning max for y-axis based on observed and predicted data
if(is.null(maxy))
    {maxy<-max(as.numeric(unlist(data.o)),
         as.numeric(unlist(data.p)))}

#Assigning size bins
size.bins<-as.numeric(dimnames(data.o)[[2]])

for(j in 1:npages){
low<-(j-1)*(nr*nc)+1
hi<-min(n.y,j*(nr*nc)) 
for(i in low:hi){
barplot(as.numeric(unlist(data.o[i,])),names.arg=size.bins,space=0,
ylim=c(0,maxy),xlab=labx,ylab="")
ii<-as.numeric(dimnames(data.o)[[1]][i])
lines(as.numeric(as.factor(size.bins))-0.5,
as.numeric(unlist(data.p[as.numeric(dimnames(data.p)[[1]])==ii])),
pch=19,type="o",lty=1,lwd=1.5)
lines(as.numeric(as.factor(size.bins))-0.5,
as.numeric(unlist(data.p2[as.numeric(dimnames(data.p2)[[1]])==ii])),
pch=22,type="o",col=4,lty=2,lwd=1.5)
mtext(as.character(dimnames(data.o)[[1]][i]),side=4,cex=1,las=1,adj=3,padj=-3)
}
mtext(title, side=3, outer=T, line=-1, cex=0.9)
}
#if(new)on.exit(par(par.old))
}
