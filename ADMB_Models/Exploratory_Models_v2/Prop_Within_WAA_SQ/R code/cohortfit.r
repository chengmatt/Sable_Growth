#Function to use color bars to follow cohorts on age composition plots
#Requires matrix of observed (data.o) and predicted (data.p) composition data with year as rownames and age/size as colnames
cohortfit <- function(data.o,data.p,recage,a,max,est,plotcolor) {
   library(fields)
   subtle.color <- plotcolor
   pred.data = data.p
   obs.data  = data.o
   years     = as.numeric(rownames(data.o))
   #print(years)
   nyears <- length(years)
   
   #=========================================================================================  
   #Use top for age data, use bottom (after else statement) for length data
   #=========================================================================================
   
   if(a=="a"){
   ages <- c(1,length(obs.data[1,] ))+recage-1 #age range
   ages.list <- ages[1]:ages[2]
   #print(ages.list)
   nages <- length(ages.list)

   #setup page for plot, sort of depends on number of years
   #mfcol <- c(ceiling(nyears/3),3)
   #mfcol=c(nyears,1)
   if(nyears<=max)
   {mfcol=c(nyears,1)}
   if(nyears>max)  
   {mfcol=c(max,1)}
   
   par(mfcol=mfcol,oma=c(4,5,3.5,1),mar=c(0,0,0,0))
   # cohort.color <- rainbow(mfcol[1]+2)[-c(1:2)]   #use hideous rainbow colors because they loop more gracefully than rich.colors
   # test colors with pie(rep(1,n),col=rainbow(n)), and define n for the number of colors to evaluate
   #cohort.color <- rainbow(ages[2]+2)[-c(1:2)]   
   cohort.color <- tim.colors(nages)   
   ncolors <- length(cohort.color)

   ylim <- c(0,1.05*max(obs.data,pred.data))
   for (yr in 1:nyears) {
      names.arg <- rep("",nages) #Allows for no names on bottom of plot at first
      if(yr<nyears){
      x <- barplot(obs.data[yr,],space=0.2,ylim=ylim,las=1,names.arg=names.arg, cex.names=1, xaxs="i",yaxs="i",border=subtle.color,
                  col=cohort.color[1:nages],axes=F,ylab="",xlab="")
      #Establishes next color palette for following year based on cohort of current year
      #Takes difference of current year and next year to loop color palette for following yr
      cohort.color <- c(cohort.color[(ncolors-(years[yr+1]-years[yr]-1)):ncolors],
                        cohort.color[1:(ncolors-(years[yr+1]-years[yr]))])
      }
      
      if(yr==nyears){
      x <- barplot(obs.data[yr,],space=0.2,ylim=ylim,las=1,names.arg=names.arg, cex.names=1, xaxs="i",yaxs="i",border=subtle.color,
                  col=cohort.color[1:nages],axes=F,ylab="",xlab="")
      }
                                     
      #Statements to set the x-axis labels only at the bottom of the set of barplots
      #Use Modulo (%%) to define based on a remainder of a division e.g. 5 %% 3 has a remainder while 6 %% 3 does not and == 0
      if (yr %% mfcol[1] == 0) {
         axis(side=1,at=x,lab=ages.list, line=-0.1,col.axis=subtle.color, col="black",lwd=0,lwd.ticks=0,cex.axis=1.5)  #just use for the labels, to allow more control than names.arg
      }
      #Statement to set top plot y-axis label complete
      if (yr == 1 | yr == max+1) {
        axis(2,las=1,at=round(c(ylim[1],ylim[2]/2,ylim[2]),2),col=subtle.color,col.axis=subtle.color,lwd=0.5,cex.axis=1.5)
      }
      #Statement to set y-axis label at left side of plot window
      if (yr <= mfcol[1] & yr > 1) {
        axis(2,las=1,at=round(ylim/2,2),col=subtle.color,col.axis=subtle.color,lwd=0.5,cex.axis=1.5)
      }
      if (yr > max+1) {
        axis(2,las=1,at=round(ylim/2,2),col=subtle.color,col.axis=subtle.color,lwd=0.5,cex.axis=1.5)
      }
      
      #Statements to put year text in each plot and box around whole plot window
      par(new=T)
      par(xpd=NA)
      if(est){
      plot(x=x,y=pred.data[yr,],ylim=ylim, xlim=par("usr")[1:2],type="b",las=1,xaxs="i",yaxs="i",col=subtle.color,bg="white",pch=19,cex=1.5,axes=F,ylab="",xlab="")}
      box(col=subtle.color,lwd=0.5)
      x.pos <- par("usr")[1] + 0.85*diff(par("usr")[1:2])   #par("usr") spits out the current coordinates of the plot window
      y.pos <- par("usr")[3] + 0.75*diff(par("usr")[3:4])   #par("usr") spits out the current coordinates of the plot window
      text(x=x.pos,y=y.pos,years[yr],cex=1.5, col=subtle.color)
      par(xpd=T)
   }
   mtext(side=1,outer=T,"Age",line=2.5)
   mtext(side=2,outer=T,"Proportion",line=4)
   #mtext(side=3,outer=T,line=1.2,paste(dat$Index_names[f],"index age composition data"))
   #mtext(side=3,outer=T,line=0.2,paste("(",case_label,")",sep=""),cex=0.6)
   #detach(dat)
   par(mfcol=c(1,1), mar=c(5.1,4.1,4.1,2.1))
   }
   
   #=========================================================================================
   #Use following for length comp data
   #=========================================================================================
   if(a!="a") {
   sizes.list <- as.numeric(unlist(dimnames(obs.data)[2])) #size range
   #print(sizes.list)
   nsizes <- length(sizes.list)

   #setup page for plot, sort of depends on number of years
   #mfcol <- c(ceiling(nyears/3),3)
   if(nyears<=max)
     {mfcol=c(nyears,1)}
   if(nyears>max)  
   {mfcol=c(max,1)}
   
   par(mfcol=mfcol,oma=c(4,5,3.5,1),mar=c(0,0,0,0))
   # cohort.color <- rainbow(mfcol[1]+2)[-c(1:2)]   #use hideous rainbow colors because they loop more gracefully than rich.colors
   # test colors with pie(rep(1,n),col=rainbow(n)), and define n for the number of colors to evaluate
   #cohort.color <- rainbow(sizes[2]+2)[-c(1:2)]   
   cohort.color<-rep("dark green",nsizes)
   #cohort.color <- gray(0:(nsizes-1)/(nsizes-1))  
   ncolors <- length(cohort.color)
   ylim <- c(0,1.05*max(obs.data,pred.data))
   for (yr in 1:nyears) {
      names.arg <- rep("",nsizes) #Allows for no names on bottom of plot at first
      if(yr<nyears){
      x <- barplot(obs.data[yr,],space=0.2,ylim=ylim,las=1,names.arg=names.arg, cex.names=0.5, xaxs="i",yaxs="i",border=subtle.color,
                  col=cohort.color[1:nsizes],axes=F,ylab="",xlab="")
      #Establishes next color palette for following year based on cohort of current year
      #Takes difference of current year and next year to loop color palette for following yr
      cohort.color <- c(cohort.color[(ncolors-(years[yr+1]-years[yr]-1)):ncolors],
                        cohort.color[1:(ncolors-(years[yr+1]-years[yr]))])
      }
      
      if(yr==nyears){
      x <- barplot(obs.data[yr,],space=0.2,ylim=ylim,las=1,names.arg=names.arg, cex.names=0.5, xaxs="i",yaxs="i",border=subtle.color,
                  col=cohort.color[1:nsizes],axes=F,ylab="",xlab="")
      }
                                     
      #Statements to set the x-axis labels only at the bottom of the set of barplots
      #Use Modulo (%%) to define based on a remainder of a division e.g. 5 %% 3 has a remainder while 6 %% 3 does not and == 0
      if (yr %% mfcol[1] == 0) {
         axis(side=1,at=x,lab=sizes.list, line=-0.1,col.axis=subtle.color, col="black",lwd=0,lwd.ticks=0,cex.axis=1.5)  #just use for the labels, to allow more control than names.arg
      }
      #Statement to set top plot y-axis label complete
      if (yr == 1 | yr == max+1) {
        axis(2,las=1,at=round(c(ylim[1],ylim[2]/2,ylim[2]),2),col=subtle.color,col.axis=subtle.color,lwd=0.5,cex.axis=1.5)
      }
      #Statement to set y-axis label at left side of plot window
      if (yr <= mfcol[1] & yr > 1) {
        axis(2,las=1,at=round(ylim/2,2),col=subtle.color,col.axis=subtle.color,lwd=0.5,cex.axis=1.5)
      }
      if (yr > max+1) {
        axis(2,las=1,at=round(ylim/2,2),col=subtle.color,col.axis=subtle.color,lwd=0.5,cex.axis=1.5)
      }
      #Statements to put year text in each plot and box around whole plot window
      par(new=T)
      par(xpd=NA)
      if(est){
      plot(x=x,y=pred.data[yr,],ylim=ylim, xlim=par("usr")[1:2],type="b",las=1,xaxs="i",yaxs="i",col=subtle.color,bg="white",pch=19,cex=1.5,axes=F,ylab="",xlab="")}
      box(col=subtle.color,lwd=0.5)
      x.pos <- par("usr")[1] + 0.85*diff(par("usr")[1:2])   #par("usr") spits out the current coordinates of the plot window
      y.pos <- par("usr")[3] + 0.75*diff(par("usr")[3:4])   #par("usr") spits out the current coordinates of the plot window
      text(x=x.pos,y=y.pos,years[yr],cex=1.5, col=subtle.color)
      par(xpd=T)
   }
   mtext(side=1,outer=T,"Size",line=2.5)
   mtext(side=2,outer=T,"Proportion",line=4)
   #mtext(side=3,outer=T,line=1.2,paste(dat$Index_names[f],"index age composition data"))
   #mtext(side=3,outer=T,line=0.2,paste("(",case_label,")",sep=""),cex=0.6)
   #detach(dat)
   par(mfcol=c(1,1), mar=c(5.1,4.1,4.1,2.1)) 
   }
}
