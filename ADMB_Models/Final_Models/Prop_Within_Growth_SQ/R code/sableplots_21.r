sableplots<-function(x,x2,x3,output="screen"){
    
    #---------------------------------------------------------------------------------
    # SOURCE THE FILES WITH THE FUNCTION
    # (From E.Williams samples)
    #---------------------------------------------------------------------------------
    
    library(MASS)# Venerables & Ripley MASS library
    library(Hmisc)
    library(gplots)
  

    #---------------------------------------------------------------------------------
    # SET UP THE GRAPHICS DEVICE AND OPTIONS
    # (From E.Williams samples)
    #---------------------------------------------------------------------------------
    bs<-16
    if(dev.cur()>1) dev.off()                     # Close any current graphics device
    
 #   if(output=="screen") {
  #   op <- options("graphicsRecord"=T)             # Turn on plot & save old options
  #  windows(width=12,height=12,record=T)           # Open a graphics device
  #  assign(".SavedPlots", NULL, "package:base")   # Erase any existing plot history
  #  old.par <- par() # Save old par settings
  #  }
  #   else {
       png(filename = "sable%02d.png", width=7, height=9,units="in",res=400,bg="transparent")
       
         #---------------------------------------------------------------------------------
    # INPUT DATA AND VARIABLE ASSIGNMENTS
    # (From D. Hanselman plots)
    #---------------------------------------------------------------------------------
    
    fishacu.o <- x$oac.fish1# Observed fishery age comps
    fishacu.p <- x$eac.fish1# Estimated fishery age comps
    surv1acu.o <- x$oac.srv1# Observed trawl survey age comps
    surv1acu.p <- x$eac.srv1# Estimated trawl survey age comps
    fishlc.o <- x$olc.fish1.f# Observed fishery length comps
    fishlc.p <- x$elc.fish1.f# Estimated fishery length comps
    surv2lc.o <- x$olc.srv1.f# Observed longline survey length comps
    surv2lc.p <- x$elc.srv1.f# Estimated longline survey length comps
    surv1lc.o <- x$olc.srv2.f# Observed trawl survey length comps
    surv1lc.p <- x$elc.srv2.f# Estimated trawl survey length comps
    parm<-as.numeric(unlist(x$parameters))# Relevant parameter estimates vector
    mcmc.df<-mcmc$mcmc# MCMC trials
    
    styear<-as.numeric(min(row.names(x$t.series)))
    endyear<-as.numeric(max(row.names(x$t.series)))
    ys<-length(as.numeric(row.names(x$t.series)))# number of years of assessment
    
    print(styear)
    
    # Defining upper and lower CI for mcmc trials <?>
    uci <- length(mcmc.df)
    lci <- length(mcmc.df)
    for (i in 1:length(mcmc.df)){
      uci[i]<-quantile(as.numeric(unlist(mcmc.df[[i]])) ,0.975)
      lci[i]<-quantile(as.numeric(unlist(mcmc.df[[i]])),0.025)}
    
    
    
    #---------------------------------------------------------------------------------
    # Likelihood Components
    # 
    #---------------------------------------------------------------------------------
    
    par(omi=c(0,0,0,0),mgp=c(2.5,1,0),mai=c(0.8,0.7,0.1,0.1))
    layout(matrix(c(1,1,1,1,1,1,2,2,3,3), 10, 1, byrow = TRUE))
    term_SSB_tem<-as.numeric(unlist(x$t.series["spbiom"]))
    term_SSB<-term_SSB_tem[ys]
    B_B40<-term_SSB/B40
    term_F_tem<-as.numeric(unlist(x$t.series["fmort"]))
    term_F<-term_F_tem[ys]
    F_F40<-term_F/F40
    likes<-as.numeric(unlist(x$likecomp))
    total_lik<-likes[23]
    total_lik<-round(total_lik,digits=2)
    like<-likes[-c(5,8,22,23)]
    like_names<-c("LL Fish Age","LL Srv Age", "LL Fish Size_F","LL Fish Size_M","TRWL Fish Size_F","TRWL Fish Size_M",
                  "LL Srv Size_F","LL Srv Size_M","Coop Srv Size_F","Coop Srv Size_M", "TRWL Srv Size_F","TRWL Srv Size_M",
                  "LL Srv RPW","Coop Srv RPW","LL Srv RPN","Coop Srv RPN","LL CPUE RPN", "JPN CPUE RPN", "Trawl Survey RPW","Catch",
                  "Recruit_Pen","F_Pen","M_Prior")           
    barplot(like,names.arg=like_names,
            col=c(rep("red",times=2),rep("darkolivegreen",times=10),rep("blue",times=4),rep("deeppink",times=3),"cyan",rep("blueviolet", times=3)),cex.names=.95,las=2)
    legend("topright",ncol=2, c("Age Comps", "Length Comps","Survey Indices","CPUE","Catch","Penalties"), 
           col = c("red","darkolivegreen","blue","deeppink","cyan","blueviolet"),lwd=2,cex=1)
    title("Component Contributions to Objective Function",line=1)
    
    summary<-matrix(c(as.character(model_conv), max_grad,total_lik,num_par),ncol=4,byrow=TRUE) 
    colnames(summary)<-c("Converged?","Max Grad","nLL","# Pars")
    summary<-as.table(summary)
    textplot(summary,show.rownames=FALSE, valign="top",halign="center",cex=1.75)
    
    summary1<-matrix(round(c( term_SSB,B40, B_B40,mean_recruit,term_F, F40,F_F40, F.ABC, ABC),digits=2),ncol=9,byrow=TRUE) 
    colnames(summary1)<-c(paste0(endyear," SSB (kt)"),
                         "SSB_40 (kt)",paste0(endyear," SSB/SSB_40"), "Mean_Recruit",paste0(endyear," F"),"F_40", paste0(endyear," F/F_40"),"F_ABC",paste0(endyear+1," ABC (kt)"))
    summary1<-as.table(summary1)
    textplot(summary1,show.rownames=FALSE, valign="top",halign="center",cex=0.9)
    
    
    #---------------------------------------------------------------------------------
    # FISHERY CATCH PLOTS LONG AND RECENT HISTORY
    # (From D. Hanselman plots)
    #---------------------------------------------------------------------------------
    
    par(omi=c(0,0,0,0),mfrow=c(2,1),mgp=c(2.5,1,0),mai=c(0.8,0.7,0.1,0.1))
    
    #par(oma=c(5,5,3.75,3.75),mfrow=c(2,1))
    
    #Complete HAL Catch predicted
    plot(as.numeric(row.names(x$t.series)),
         as.numeric(unlist(x$t.series["Catch_HAL"])+unlist(x$t.series["Catch_TWL"])),type="l",
         lty=1,lwd=3,col="black",xlab="Year",ylab="Catch (kt)",ylim=c(0,60))
    
    lines(as.numeric(row.names(x$t.series)),
          as.numeric(unlist(x$t.series["Catch_TWL"])),type="l",
          lty=2,lwd=3,col="blue",xlab="Year",ylab="Catch (kt)")
    
    lines(as.numeric(row.names(x$t.series)),
          as.numeric(unlist(x$t.series["Catch_HAL"])),
          lty=3,lwd=3,col="red",xlab="Year",ylab="Catch (kt)")
    legend(1990,50,legend=c("Total","Trawl Gear","Fixed Gear"),col=c("black","blue","red"),lwd=c(3,2.5,2.5),lty=c(1,2,3))
    #Complete Trawl catch
    
    plot(as.numeric(row.names(x$t.series))[(ys/2):ys+1],
         as.numeric(unlist(x$t.series["Catch_HAL"])+unlist(x$t.series["Catch_TWL"]))[(ys/2):ys+1],type="l",
         lty=1,lwd=3,col="black",xlab="Year",ylab="Catch (kt)",ylim=c(0,40))
    
    lines(as.numeric(row.names(x$t.series))[(ys/2):ys+1],
          as.numeric(unlist(x$t.series["Catch_TWL"]))[(ys/2):ys+1],type="l",
          lty=2,lwd=3,col="blue",xlab="Year",ylab="Catch (kt)")
    
    lines(as.numeric(row.names(x$t.series))[(ys/2):ys+1],
          as.numeric(unlist(x$t.series["Catch_HAL"]))[(ys/2):ys+1],
          lty=3,lwd=3,col="red",xlab="Year",ylab="Catch (kt)")
    legend(2000,35,legend=c("Total","Trawl Gear","Fixed Gear"),col=c("black","blue","red"),lwd=c(3,2.5,2.5),lty=c(1,2,3))
    
    
    #Catch versus Survey Biomass
    #plot(as.numeric(row.names(x$t.series)),
    #  as.numeric(unlist(x$t.series["catchbio"])),
    #type="l",lty=1,lwd=1.5,xlab="Year",ylab="Catch (mt)")
    
    
    #---------------------------------------------------------------------------------
    # SURVEY BIOMASS PLOTS WITH MODEL FIT
    # (From D. Hanselman plots)
    #---------------------------------------------------------------------------------
    
    par(omi=c(0,0,0,0), mfrow=c(2,2),mgp=c(2.5,1,0))
    
    #----- SURVEY ESTIMATES ----------------------------------------------------------
    
    #Plot with CI intervals as points
    maxy<-1.05*max(as.numeric(unlist(x$obssrv1["predsrv1"])),as.numeric(unlist(x$obssrv1["obssrv1.uci"])))
    plot(as.numeric(row.names(x$obssrv1)), as.numeric(unlist(x$obssrv1["obssrv1.uci"])), 
         ylim=range(c(0,maxy)),type='n',xlab="Year", ylab="Domestic LL Relative Population Weight")
    points(as.numeric(row.names(x$obssrv1)), as.numeric(unlist(x$obssrv1["obssrv1"])),cex=1.5,col="blue")
    arrows(as.numeric(row.names(x$obssrv1)),as.numeric(unlist(x$obssrv1["obssrv1"])),
           as.numeric(row.names(x$obssrv1)),as.numeric(unlist(x$obssrv1["obssrv1.uci"])),
           length=0.04, angle=90)
    arrows(as.numeric(row.names(x$obssrv1)),as.numeric(unlist(x$obssrv1["obssrv1"])),
           as.numeric(row.names(x$obssrv1)),as.numeric(unlist(x$obssrv1["obssrv1.lci"])),
           length=0.04, angle=90)
    lines(as.numeric(row.names(x$obssrv1)),as.numeric(unlist(x$obssrv1["predsrv1"])), 
          type="l", lwd=3,lty=1,col="red")
    
    #Plot with CI intervals as points
    maxy<-1.05*max(as.numeric(unlist(x$obssrv2["predsrv2"])),as.numeric(unlist(x$obssrv2["obssrv2.uci"])))
    plot(as.numeric(row.names(x$obssrv2)), as.numeric(unlist(x$obssrv2["obssrv2.uci"])), 
         ylim=range(c(0,maxy)),type='n',xlab="Year", ylab="Japanese LL Relative Population Weight")
    points(as.numeric(row.names(x$obssrv2)), as.numeric(unlist(x$obssrv2["obssrv2"])),cex=1.5,col="blue")
    arrows(as.numeric(row.names(x$obssrv2)),as.numeric(unlist(x$obssrv2["obssrv2"])),
           as.numeric(row.names(x$obssrv2)),as.numeric(unlist(x$obssrv2["obssrv2.uci"])),
           length=0.04, angle=90)
    arrows(as.numeric(row.names(x$obssrv2)),as.numeric(unlist(x$obssrv2["obssrv2"])),
           as.numeric(row.names(x$obssrv2)),as.numeric(unlist(x$obssrv2["obssrv2.lci"])),
           length=0.04, angle=90)
    lines(as.numeric(row.names(x$obssrv2)),as.numeric(unlist(x$obssrv2["predsrv2"])), 
          type="l", lwd=3,lty=1,col="red")
    
    #Plot with CI intervals as points
    maxy<-1.05*max(as.numeric(unlist(x$obssrv3["predsrv3"])),as.numeric(unlist(x$obssrv3["obssrv3.uci"])))
    plot(as.numeric(row.names(x$obssrv3)), as.numeric(unlist(x$obssrv3["obssrv3.uci"])), 
         ylim=range(c(0,maxy)),type='n',xlab="Year", ylab="Domestic LL Relative Population Number")
    points(as.numeric(row.names(x$obssrv3)), as.numeric(unlist(x$obssrv3["obssrv3"])),cex=1.5,col="blue")
    arrows(as.numeric(row.names(x$obssrv3)),as.numeric(unlist(x$obssrv3["obssrv3"])),
           as.numeric(row.names(x$obssrv3)),as.numeric(unlist(x$obssrv3["obssrv3.uci"])),
           length=0.04, angle=90)
    arrows(as.numeric(row.names(x$obssrv3)),as.numeric(unlist(x$obssrv3["obssrv3"])),
           as.numeric(row.names(x$obssrv3)),as.numeric(unlist(x$obssrv3["obssrv3.lci"])),
           length=0.04, angle=90)
    lines(as.numeric(row.names(x$obssrv3)),as.numeric(unlist(x$obssrv3["predsrv3"])), 
          type="l", lwd=3,lty=1,col="red")
    
    #Plot with CI intervals as points
    maxy<-1.05*max(as.numeric(unlist(x$obssrv4["predsrv4"])),as.numeric(unlist(x$obssrv4["obssrv4.uci"])))
    plot(as.numeric(row.names(x$obssrv4)), as.numeric(unlist(x$obssrv4["obssrv4.uci"])), 
         ylim=range(c(0,maxy)),type='n',xlab="Year", ylab="Japanese LL Survey Relative Population Number")
    points(as.numeric(row.names(x$obssrv4)), as.numeric(unlist(x$obssrv4["obssrv4"])),cex=1.5,col="blue")
    arrows(as.numeric(row.names(x$obssrv4)),as.numeric(unlist(x$obssrv4["obssrv4"])),
           as.numeric(row.names(x$obssrv4)),as.numeric(unlist(x$obssrv4["obssrv4.uci"])),
           length=0.04, angle=90)
    arrows(as.numeric(row.names(x$obssrv4)),as.numeric(unlist(x$obssrv4["obssrv4"])),
           as.numeric(row.names(x$obssrv4)),as.numeric(unlist(x$obssrv4["obssrv4.lci"])),
           length=0.04, angle=90)
    lines(as.numeric(row.names(x$obssrv4)),as.numeric(unlist(x$obssrv4["predsrv4"])), 
          type="l", lwd=3,lty=1,col="red")
    
    #Plot with CI intervals as points
    maxy<-1.05*max(as.numeric(unlist(x$obssrv5["predsrv5"])),as.numeric(unlist(x$obssrv5["obssrv5.uci"])))
        plot(as.numeric(row.names(x$obssrv5)), as.numeric(unlist(x$obssrv5["obssrv5.uci"])), 
         ylim=range(c(0,maxy)),type='n',xlab="Year", ylab="Domestic Fishery CPUE Index")
    points(as.numeric(row.names(x$obssrv5)), as.numeric(unlist(x$obssrv5["obssrv5"])),cex=1.5,col="blue")
    arrows(as.numeric(row.names(x$obssrv5)),as.numeric(unlist(x$obssrv5["obssrv5"])),
           as.numeric(row.names(x$obssrv5)),as.numeric(unlist(x$obssrv5["obssrv5.uci"])),
           length=0.04, angle=90)
    arrows(as.numeric(row.names(x$obssrv5)),as.numeric(unlist(x$obssrv5["obssrv5"])),
           as.numeric(row.names(x$obssrv5)),as.numeric(unlist(x$obssrv5["obssrv5.lci"])),
           length=0.04, angle=90)
    lines(as.numeric(row.names(x$obssrv5)),as.numeric(unlist(x$obssrv5["predsrv5"])), 
          type="l", lwd=3,lty=1,col="red")
    
    #Plot with CI intervals as points
    maxy<-1.05*max(as.numeric(unlist(x$obssrv6["predsrv6"])),as.numeric(unlist(x$obssrv6["obssrv6.uci"])))
    plot(as.numeric(row.names(x$obssrv6)), as.numeric(unlist(x$obssrv6["obssrv6.uci"])), 
         ylim=range(c(0,maxy)),type='n',xlab="Year", ylab="Japanese Fishery CPUE Index")
    points(as.numeric(row.names(x$obssrv6)), as.numeric(unlist(x$obssrv6["obssrv6"])),cex=1.5,col="blue")
    arrows(as.numeric(row.names(x$obssrv6)),as.numeric(unlist(x$obssrv6["obssrv6"])),
           as.numeric(row.names(x$obssrv6)),as.numeric(unlist(x$obssrv6["obssrv6.uci"])),
           length=0.04, angle=90)
    arrows(as.numeric(row.names(x$obssrv6)),as.numeric(unlist(x$obssrv6["obssrv6"])),
           as.numeric(row.names(x$obssrv6)),as.numeric(unlist(x$obssrv6["obssrv6.lci"])),
           length=0.04, angle=90)
    lines(as.numeric(row.names(x$obssrv6)),as.numeric(unlist(x$obssrv6["predsrv6"])), 
          type="l", lwd=3,lty=1,col="red")
    
    #Plot with CI intervals as points
    maxy<-1.05*max(as.numeric(unlist(x$obssrv7["predsrv7"])),as.numeric(unlist(x$obssrv7["obssrv7.uci"])))
    plot(as.numeric(row.names(x$obssrv7)), as.numeric(unlist(x$obssrv7["obssrv7.uci"])), 
         ylim=range(c(0,maxy)),type='n',xlab="Year", ylab="GOA Trawl Survey Biomass (kt)")
    points(as.numeric(row.names(x$obssrv7)), as.numeric(unlist(x$obssrv7["obssrv7"])),cex=1.5,col="blue")
    arrows(as.numeric(row.names(x$obssrv7)),as.numeric(unlist(x$obssrv7["obssrv7"])),
           as.numeric(row.names(x$obssrv7)),as.numeric(unlist(x$obssrv7["obssrv7.uci"])),
           length=0.04, angle=90)
    arrows(as.numeric(row.names(x$obssrv7)),as.numeric(unlist(x$obssrv7["obssrv7"])),
           as.numeric(row.names(x$obssrv7)),as.numeric(unlist(x$obssrv7["obssrv7.lci"])),
           length=0.04, angle=90)
    lines(as.numeric(row.names(x$obssrv7)),as.numeric(unlist(x$obssrv7["predsrv7"])), 
          type="l", lwd=3,lty=1,col="red")
    
    
    #---------------------------------------------------------------------------------
    # PREDICTED BIOMASS PLOT OVER TIME WITH MCMC CIs
    # (From D. Hanselman plots)
    #---------------------------------------------------------------------------------
    
    #axis(side=4)
    #mtext(side=4, "Catch(t)", outer=T, line=1.0)
    #par(xaxs="r")
    
    #par(new=T, xaxs="d")
    par(omi=c(0,0,0,0),mfrow=c(2,1),mgp=c(2.5,1,0))
    print(ys)
    plot(as.numeric(row.names(x$t.series)),
         as.numeric(unlist(x$t.series["totbiom"])),type="l",col="blue",lwd=2,
         lty=1,xlab ="Year",ylab="Predicted Biomass (kt)",yaxs="r",ylim=c(0,(1.2*max(as.numeric(unlist(x$t.series["totbiom"]))))))
    lines(as.numeric(row.names(x$t.series)),
          as.numeric(uci[15:(15+ys-1)]), type="l", lty=3,lwd=2)
    lines(as.numeric(row.names(x$t.series)),
          as.numeric(lci[15:(15+ys-1)]), type="l", lty=3,lwd=2)
    
    #mtext("Figure 7-15.  Predicted total biomass for GOA Pacific ocean perch.",
    #side=1, outer=F, line=5, adj =0, cex=0.9)
    #mtext("Dashed lines are 95% confidence intervals from 5,000,000 MCMC runs.",
    #side=1, outer=F, line=6, adj =0, cex=0.9)
    print(9+2*ys)
    print(length(mcmc.df)-1)
    # Now plot spawning biomass w/ MCMC CIs
    plot(as.numeric(row.names(x$t.series)),
         as.numeric(unlist(x$t.series["spbiom"])),type="l",col="blue",lwd=2,
         lty=1,xlab ="Year",ylab="Spawning Biomass (kt)",yaxs="r",ylim=c(0,(1.2*max(as.numeric(unlist(x$t.series["spbiom"]))))))
    lines(as.numeric(row.names(x$t.series)),
          as.numeric(uci[(15+ys):(15+2*ys-1)]), type="l", lty=3,lwd=2)
    lines(as.numeric(row.names(x$t.series)),
          as.numeric(lci[(15+ys):(15+2*ys-1)]), type="l", lty=3,lwd=2)
    
    #---------------------------------------------------------------------------------
    # SELECTIVITY CURVES
    #---------------------------------------------------------------------------------
    par(omi=c(0,0,0,0),mfrow=c(3,2),mgp=c(2.5,1,0))
    
    #par(omi=c(0.25,0.25,0.25,0.25), mfrow=c(2,2))
    sellabs<-c("Derby fishery female","Derby fishery male","Trawl fishery female" ,"Trawl fishery male" ,"IFQ fishery female" ,"IFQ fishery male", "IFQ Recent fishery female" ,
               "IFQ Recent fishery male", "Domestic LL survey female","Domestic LL survey male", "Domestic LL Recent survey female","Domestic LL Recent survey male","Cooperative LL survey female" ,
               "Cooperative LL survey male" ,"GOA trawl survey female","GOA trawl survey male" )
    
    for(i in 1:length(x$agesel[1,])){
      plot(as.numeric(dimnames(x$agesel)[[1]]),as.numeric(x$agesel[,i]),
           pch=19,type="o",ylim=c(0.0,1.1),xlab="Age",ylab="Proportion selected")
      text(20,0.8,sellabs[i],cex=1.3,col="blue")
    }
    
        #---------------------------------------------------------------------------------
    # STOCK-RECRUIT DATA PLOT, MANAGEMENT PLOT (F/F40 VS B/B40)
    # (From D. Hanselman plots)
    #---------------------------------------------------------------------------------
    par(omi=c(0.25,0.25,0,0),mfrow=c(1,1))
    
    #par(omi=c(0,0,0,0), mfrow=c(1,1))
    #plot(Brat,as.numeric(unlist(x$t.series["Frat"]))[1:ys],type="l",lty=1,lwd=0.5,xlab="B/B40",ylab="F/F40")
    #plot(as.numeric(unlist(x$t.series["Brat"]))[1:ys],as.numeric(unlist(x$t.series["Frat"]))[1:ys],type="l",lty=1,lwd=0.5,xlab="B/B40",ylab="F/F40")
    #text(as.numeric(unlist(x$t.series["Brat"]))[1:ys], 
    #text(Brat, 
    #as.numeric(unlist(x$t.series["Frat"]))[1:ys], 
    #as.numeric(unlist(x$t.series["year"]))[1:ys],cex=0.5)
    #mtext("Figure 10A-8.  Time series of estimated fishing mortality over F40 ", 
    #side=1, outer=F, line=5, adj =0, cex=0.9)
    #mtext("versus estimated spawning biomass over B40", 
    #side=1, outer=F, line=6, adj =0, cex=0.9)
    
    
    maxx<-max(as.numeric(unlist(x$t.series["spbiom"])))
    maxy<-max(as.numeric(unlist(x$t.series["Recr"])))
    par(omi=c(0,0,0,0), mfrow=c(2,1))
    plot(as.numeric(unlist(x$t.series["spbiom"]))[1:(ys-2)], 
         as.numeric(unlist(x$t.series["Recr"]))[3:ys], type = "p", pch=" ", 
         xlab="SSB(kt)", ylab="Age 2 Recruits (millions)", xlim=c(0,maxx),ylim=c(0,maxy))
    text(as.numeric(unlist(x$t.series["spbiom"]))[1:(ys-2)], 
         as.numeric(unlist(x$t.series["Recr"]))[3:ys], 
         as.numeric(unlist(x$t.series["year"]))[1:(ys-2)],cex=0.75,col="dark blue")
 
    #---------------------------------------------------------------------------------
    # AGE AND LENGTH COMPOSITION PLOTS WITH MODEL FIT
    # (Uses recomps.r package, from D. Hanselman plots)
    #---------------------------------------------------------------------------------
  
    
    pres<-"black"
    
    cohortfit(x$oac.fish1,x$eac.fish1,2,"a",7,T,plotcolor=pres)
    cohortfit(x$oac.srv1,x$eac.srv1,2,"a",6,T,plotcolor=pres)
    cohortfit(x$oac.srv2,x$eac.srv2,2,"a",8,T,plotcolor=pres)
    cohortfit(x$olc.fish1.f,x$elc.fish1.f,2,"l",15,T,plotcolor=pres)
    cohortfit(x$olc.fish1.m,x$elc.fish1.m,2,"l",15,T,plotcolor=pres)
    cohortfit(x$olc.fish3.f,x$elc.fish3.f,2,"l",17,T,plotcolor=pres)
    cohortfit(x$olc.fish3.m,x$elc.fish3.m,2,"l",17,T,plotcolor=pres)
    cohortfit(x$olc.srv1.f,x$elc.srv1.f,2,"l",15,T,plotcolor=pres)
    cohortfit(x$olc.srv1.m,x$elc.srv1.m,2,"l",15,T,plotcolor=pres)
    cohortfit(x$olc.srv2.f,x$elc.srv2.f,2,"l",16,T,plotcolor=pres)
    cohortfit(x$olc.srv2.m,x$elc.srv2.m,2,"l",16,T,plotcolor=pres)
    cohortfit(x$olc.srv7.f,x$elc.srv7.f,2,"l",15,T,plotcolor=pres)
    cohortfit(x$olc.srv7.m,x$elc.srv7.m,2,"l",15,T,plotcolor=pres)
    
    
    dev.off()
    png(filename = "sable_end2%02d.png", width=8, height=5,units="in",res=400,bg="transparent")
    
    
    ### MCMC projection
    newswath(mcmc,styear,endyear)
     
    #### base version 
    d<-mcmc.df[,(3*ys+15):(3*ys+29)]
    belowb40<-d/B40
    belowb35<-d/B35
    belowmmst<-d/(0.5*B35)
    pbelowb35<-seq(1,14)
    pbelowb40<-seq(1,14)
    pbelowmsst<-seq(1,14)
    for(i in 1:14) {
      pbelowb40[i]<- length(which(belowb40[,i]<1))/length(belowb40[,i])
      pbelowb35[i]<- length(which(belowb35[,i]<1))/length(belowb40[,i])
      pbelowmsst[i]<- length(which(belowmmst[,i]<1))/length(belowb40[,i]) }
    
    plot(seq(endyear+1,endyear+14),pbelowmsst,ylim=c(0,1),col=4,xlab="Year",ylab="Probability",pch="")
    lines(seq(endyear+1,endyear+14),pbelowmsst,lty=2,col=4,lwd=3)
    points(seq(endyear+1,endyear+14),pbelowb35,col=3)
    lines(seq(endyear+1,endyear+14),pbelowb35,col=3,lty=3,lwd=3)
    points(seq(endyear+1,endyear+14),pbelowb40,col=9)
    lines(seq(endyear+1,endyear+14),pbelowb40,col=9,lty=1,lwd=3)
    legend(endyear+2,0.48,c("P<B40%","P<B35%","P<B17.5%"),lty=c(1,3,2),col=c(1,3,4),cex=0.9,lwd=c(3,3,3))
    
    #### More MCMC probabilities
    
    truehist(mcmc.df[,(3*ys+15)],xlim=c(100,200),xlab="Female Spawning Biomass (kt)",ylab="Density",nbins=25,col="royal blue")
    lines(x=c(B40,B40),y=c(0,0.12),lwd=3.5,lty=2,col="red")
    lines(c(B35,B35),y=c(0,0.12),lwd=3.5,lty=2,col="black")
    text(B40-7,0.021,expression(paste(italic(B)["40%"])),col="red",cex=1.4)
    text(B35-7,0.021,expression(paste(italic(B)["35%"])),col="black",cex=1.4)

    dev.off()
    png(filename = "sable_end%02d.png", width=6, height=8,units="in",res=400,bg="transparent")
    
   ### Plot correlation of MCMC matrix
    ## panel.cor function for making loess and numbers on graph
    panel.cor <- function(x, y, digits=2, prefix="", cex.cor) 
    {
      usr <- par("usr"); on.exit(par(usr)) 
      par(usr = c(0, 1, 0, 1)) 
      r <- abs(cor(x, y)) 
      txt <- format(c(r, 0.123456789), digits=digits)[1] 
      txt <- paste(prefix, txt, sep="") 
      if(missing(cex.cor)) cex <- 0.8/strwidth(txt) 
      test <- cor.test(x,y) 
      # borrowed from printCoefmat
      Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                       cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                       symbols = c("***", "**", "*", ".", " ")) 
      text(0.5, 0.5, txt, cex = cex * r) 
      #text(.8, .8, Signif, cex=cex, col=2) # for adding stars for significance
    }
    
    
    b<-cbind(mcmc.df[,4],mcmc.df[,5],mcmc.df[,8],mcmc.df[,9],mcmc.df[,(2*ys+14)],mcmc.df[,(3*ys+30)],mcmc.df[,(3*ys+4)],mcmc.df[,12])   
    colnames(b)<-c("Dom.LL.q","Jap.LL.q","Dom.Fish.q","Jap.Fish.q","Ending SSB","ABC","2000 YC","F40%")
    pairs(b,lower.panel=panel.smooth, upper.panel=panel.cor)
    
     
    ####Management phase plane plot
    ####
    # Set up plot panel and colors
    par(mfrow=c(2,1),mai=c(1.1,1.1,0.1,0.1))
    FB<-data.frame(cbind(as.numeric(rownames(x$t.series)),x$t.series[,5],x$t.series[,6]))
    colnames(FB)<-c("Year","B","F")
    nyrs<-endyear-styear+1
    FB_proj<-data.frame(matrix(nrow=2,ncol=3))
    colnames(FB_proj)<-c("Year","B","F")
    FB_proj[,1]<-(endyear+1):(endyear+2)
    FB_proj[,2]<-c(mean(mcmc.df[,(3*ys+14)]),mean(mcmc.df[,(3*ys+15)]))
    FB_proj[,3]<-FB[length(FB[,1]),3]
    FB<-rbind(FB,FB_proj)
    
    st_yr<-min(FB$Year)
    end_yr<-max(FB$Year)
    F35<-0.40/0.35*x$parameters[11]
    F_B<-FB
    F_B[,2]<-FB[,2]/B35
    F_B[,3]<-FB[,3]/F35
    #phasecolor<- tim.colors(6)
    #Use colors() to see color names list that R knows about
    phasecolor<- c("black","red","darkorange","forestgreen","deepskyblue","blue")
    #dec<-ceiling((max_yr-min_yr)/10)
    plotcolor<-"black"
    # Plot full chart
    yl<-c(0,ceiling(1.05*max(F_B$F)))
    xl<-c(0,ceiling(1.05*max(F_B$B)))
    
    plot(F_B$B[1],F_B$F[1],ylim=c(0,ceiling(max(F_B$F))),xlim=c(0,ceiling(max(F_B$B))),pch=NA,type="n",lwd=2,las=2,xaxt="n",ylab="",xlab="",
         col=plotcolor, col.axis=plotcolor,col.lab=plotcolor, fg=plotcolor)
    
    for(i in 1:dim(F_B)[1]){
      
      if((F_B$Year[i]) %% 100 >= 60 & (F_B$Year[i]) %% 100 < 70)
      {text(F_B$B[i],F_B$F[i],labels=F_B$Year[i]-1900,cex=0.75,col=phasecolor[1],font=2)}  
      if((F_B$Year[i]) %% 100 >= 70 & (F_B$Year[i]) %% 100 < 80)
      {text(F_B$B[i],F_B$F[i],labels=F_B$Year[i]-1900,cex=0.75,col=phasecolor[2],font=2)}  
      if((F_B$Year[i]) %% 100 >= 80 & (F_B$Year[i]) %% 100 < 90)
      {text(F_B$B[i],F_B$F[i],labels=F_B$Year[i]-1900,cex=0.75,col=phasecolor[3],font=2)}  
      if((F_B$Year[i]) %% 100 >= 90)
      {text(F_B$B[i],F_B$F[i],labels=F_B$Year[i]-1900,cex=0.75,col=phasecolor[4],font=2)}  
      if((F_B$Year[i]) %% 100 >= 0 & (F_B$Year[i]) %% 100 < 10)
      {text(F_B$B[i],F_B$F[i],labels=F_B$Year[i]-2000,cex=0.75,col=phasecolor[5],font=2)}  
      if((F_B$Year[i]) %% 100 >= 10 & (F_B$Year[i]) %% 100 < 20)
      {text(F_B$B[i],F_B$F[i],labels=F_B$Year[i]-2000,cex=0.75,col=phasecolor[6],font=2)}  
      if((F_B$Year[i]) %% 100 >= 20 & (F_B$Year[i]) %% 100 < 30)
      {text(F_B$B[i],F_B$F[i],labels=F_B$Year[i]-2000,cex=0.75,col=phasecolor[6],font=2)}  
      
    }
    
    points(F_B$B[1:dim(F_B)[1]],F_B$F[1:dim(F_B)[1]],type="b",col=plotcolor,lwd=2,pch=NA)
    par(fg=plotcolor)
    axis(1,0,0, col.axis=plotcolor)
    axis(1,0.5,0.5,col.axis=plotcolor)
    axis(1,1,1,col.axis=plotcolor)
    axis(1,1.5,1.5,col.axis=plotcolor)
    axis(1,2,2,col.axis=plotcolor)
    axis(1,2.5,2.5,col.axis=plotcolor)
    axis(1,3,3,col.axis=plotcolor)
    
    # Redone by D. Hanselman to be linked to control rule
    Fabc<-0.35/0.4
    
    segments(0.05,0,0.4/0.35,1,col="red",lwd=2.5)
    segments(0.4/0.35,1,4,1,col="blue",lwd=2.5)
    
    segments(0.05,0,0.4/0.35,Fabc,col="red",lwd=2.5,lty=2)
    segments(0.4/0.35,Fabc,4,Fabc,col="blue",lwd=2.5,lty=2)
    par(fg=pres)
    legend("topleft",legend=c(expression(F[OFL]),expression(F[ABC])),lwd=c(2.5,2.5),col=c("blue","blue"),lty=c(1,2),
           box.col=plotcolor, title.col=plotcolor)
    #legend("topleft",legend=c("1961-69","1970-79","1980-89","1990-99","2000-09","2010-11",expression(F[OFL]),expression(F[ABC])),lwd=c(2,2,2,2,2,2,2.5,2.5),col=c(phasecolor[1:6],"blue","blue"),lty=c(1,1,1,1,1,1,1,2))
    
    # Plot bottom panel zooooomed in on the most recent year
    # Only need this plot for species that have low catch for whole time series
    
    #  par(mar=c(4,4,0.5,0.5))
    segments(xl[1]-0.07,1,xl[2],1,lty=2,lwd=0.7)
    segments(1,yl[1]-0.025,1,yl[2]+0.05,lty=2,lwd=0.7)
    #  abline(h=1,lty=2,lwd=0.7)
    # abline(v=1,lty=2,lwd=0.7)
    yl<-c(0,1)
    xl<-c(0,2.25)
    plot(F_B$B[1],F_B$F[1],ylim=yl,xlim=xl,pch=NA,type="n",lwd=2,las=2,xaxt="n",ylab="",xlab="",
         col=plotcolor, col.axis=plotcolor,col.lab=plotcolor, fg=plotcolor)
    
    for(i in 1:dim(F_B)[1]){
      
      if((F_B$Year[i]) %% 100 >= 60 & (F_B$Year[i]) %% 100 < 70)
      {text(F_B$B[i],F_B$F[i],labels=F_B$Year[i]-1900,cex=0.75,col=phasecolor[1],font=2)}  
      if((F_B$Year[i]) %% 100 >= 70 & (F_B$Year[i]) %% 100 < 80)
      {text(F_B$B[i],F_B$F[i],labels=F_B$Year[i]-1900,cex=0.75,col=phasecolor[2],font=2)}  
      if((F_B$Year[i]) %% 100 >= 80 & (F_B$Year[i]) %% 100 < 90)
      {text(F_B$B[i],F_B$F[i],labels=F_B$Year[i]-1900,cex=0.75,col=phasecolor[3],font=2)}  
      if((F_B$Year[i]) %% 100 >= 90)
      {text(F_B$B[i],F_B$F[i],labels=F_B$Year[i]-1900,cex=0.75,col=phasecolor[4],font=2)}  
      if((F_B$Year[i]) %% 100 >= 0 & (F_B$Year[i]) %% 100 < 10)
      {text(F_B$B[i],F_B$F[i],labels=F_B$Year[i]-2000,cex=0.75,col=phasecolor[5],font=2)}  
      if((F_B$Year[i]) %% 100 >= 10 & (F_B$Year[i]) %% 100 < 20 & (F_B$Year[i])<=end_yr-1)
      {text(F_B$B[i],F_B$F[i],labels=F_B$Year[i]-2000,cex=0.75,col=phasecolor[6],font=2)}  
      if((F_B$Year[i]) %% 100 >= 20 & (F_B$Year[i]) %% 100 < 30 & (F_B$Year[i])<=end_yr-1)
      {text(F_B$B[i],F_B$F[i],labels=F_B$Year[i]-2000,cex=0.75,col=phasecolor[6],font=2)}  
      # For RE, made last year position slightly above the location using pos, since label was on top of others
      if((F_B$Year[i]) >= end_yr)
      {text(F_B$B[i],F_B$F[i],labels=F_B$Year[i],cex=1,pos=1,col=phasecolor[6],font=2)}  
      
    }
    
    points(F_B$B[1:dim(F_B)[1]],F_B$F[1:dim(F_B)[1]],type="b",col=plotcolor,lwd=2,pch=NA)
    par(fg=pres)
    axis(1,0,0,col.axis=plotcolor)
    axis(1,0.5,0.5,col.axis=plotcolor)
    axis(1,1,1,col.axis=plotcolor)
    axis(1,1.5,1.5,col.axis=plotcolor)
    axis(1,2,2,col.axis=plotcolor)
    axis(1,2.5,2.5,col.axis=plotcolor)
    axis(1,3,3,col.axis=plotcolor)
    
    # Linked to control rule
    Fabc<-0.35/0.4
    
    segments(0.05,0,0.4/0.35,1,col="red",lwd=2.5)
    segments(0.4/0.35,1,4,1,col="blue",lwd=2.5)
    
    segments(0.05,0,0.4/0.35,Fabc,col="red",lwd=2.5,lty=2)
    segments(0.4/0.35,Fabc,4,Fabc,col="blue",lwd=2.5,lty=2)
    
    par(fg=pres)
    # Write margin text
    mtext(expression(italic(F/F["35%"])),at=1.2,side=2,line=3,cex=1.5,col=plotcolor)
    mtext(expression(italic(SSB/B["35%"])),side=1,line=3,cex=1.5,col=plotcolor)
    
    segments(xl[1]-0.07,yl[2],xl[2],yl[2],lty=2,lwd=0.7)
    segments(1,yl[1]-0.025,1,yl[2]+0.025,lty=2,lwd=0.7)
    #abline(h=1,lty=2,lwd=0.7)
    #abline(v=1,lty=2,lwd=0.7)
    
    ####
    
    plot(rnorm(10),rnorm(10),type="n")  # DUMMY PLOT so last plot gets recorded
    plot(rnorm(10),rnorm(10),type="n")  # DUMMY PLOT so last plot gets recorded
    
   
    dev.off()
 
    #### GG version of below
    ProjectedSSB=mcmc.df[,c((3*ys+15),(3*ys+16),(3*ys+17))]
    names(ProjectedSSB)<-c(endyear+1,endyear+2,endyear+3)
    PSSB<-melt(ProjectedSSB)
    names(PSSB)<-c("Year","Proj.SSB")
   p<- ggplot(PSSB,aes(x=Proj.SSB)) + geom_density(aes(group=Year,fill=Year), alpha=0.3)+theme_bw(base_size=bs)+
      ylab("Density")+geom_vline(xintercept=B40,linetype=5,color="green",size=1.8,alpha=0.7)+
      geom_vline(xintercept=B35,linetype=5,color="red",size=1.8,alpha=0.7)+
      xlab("Projected Female Spawning Biomass")+
      geom_text(aes(x=B40+15,y=0.035),label="B40%",size=6,color="dark green")+
      geom_text(aes(x=B35-15,y=0.035),label="B35%",size=6,color="dark red")+
      xlim(80,300)
    ggsave("Projected_SSB_Density.png",dpi=400,width=10,height=6,units="in",plot=p)
    
    #---------------------------------------------------------------------------------
    # Natural MORTALITY RATE
    # (From D. Hanselman plots)
    #---------------------------------------------------------------------------------
    
    par(omi=c(0,0.1,0,0.25),mgp=c(2.5,1,0))
    maxy<-max(as.numeric(unlist(x$M)))
    M_yr<-as.numeric(row.names(x$t.series))
    M_age<-as.numeric(colnames(x$M))
    
    
    M<- as.numeric(unlist(x$M))
    M<- matrix(M, byrow = FALSE, nrow = length(M_yr), ncol = length(M_age))
    M_melt<-data.frame(cbind(M,M_yr),check.names=FALSE)

    year_M<-as.numeric(unlist(M[,15]))
    mean_M_yr<-mean(year_M)
    M_melt_yr<-data.frame(cbind(year_M,M_yr))    
    
    if(as.numeric(unlist(x$phases["ph.Mdevs.age"]))>0)
    {
      keep.years <- c(seq(1960, endyear, 20), endyear)
      M<- matrix(M, byrow = FALSE, nrow = length(M_yr), ncol = length(M_age), 
                 dimnames= list(c(M_yr), c(M_age)))
      
      M_new <- M %>% as.data.frame() %>%
        rownames_to_column("Year") %>%
        pivot_longer(-Year, names_to = "Age", values_to= "Values") %>%   # combo of transpose and melt
        mutate(Plot_years= ifelse(Year %in% keep.years, 1, 0),           # keep 
               Age= as.numeric(Age))
      M_sub <- M_new %>% filter(Plot_years == 1)
      p<-ggplot(M_sub,aes(x=Age,y=Values, type=Year, color=Year))+
        geom_line(lwd=1.5)+
        ylim(0,1.1*max(M_sub$Values))+
        ylab("M")+  scale_x_continuous(breaks = seq(2, 31, by = 4))+
        xlab("Age")
      ggsave("Age_Natural_Mortality.png",width=8,height=5,dpi=300,plot=p)
     }
    
    
   p<-ggplot(M_melt_yr,aes(x=M_yr,y=year_M))+
            geom_bar(stat="identity", position=position_dodge(),fill="yellow",colour="dark blue") +
            scale_fill_brewer(palette="Paired") + theme_bw(base_size = 16)+
            geom_hline(yintercept=mean_M_yr,colour="red",linetype=1,size=1.)+
            ylab("M at Age-16")+  scale_x_continuous(breaks = seq(min(M_melt_yr$M_yr), max(M_melt_yr$M_yr), by = 4))+
            xlab("Year")
    ggsave("Yearly_Natural_Mortality.png",width=8,height=5,dpi=300,plot=p)
    
    
    #---------------------------------------------------------------------------------
    # Fully Selected Fishing MORTALITY RATE
    # (From D. Hanselman plots)
    #---------------------------------------------------------------------------------
    
    par(omi=c(0,0.1,0,0.25),mfrow=c(2,1),mgp=c(2.5,1,0))
    maxy<-max(as.numeric(unlist(x$t.series["fmort"])))
    #par(omi=c(0.25,0.25,0,0), mfrow=c(1,1))
    F<- as.numeric(unlist(x$t.series["fmort"]))
    F_yr<-as.numeric(row.names(x$t.series))
    F_melt<-data.frame(cbind(F,F_yr))
    mean_F<-mean(F)
    p<-ggplot(F_melt,aes(x=F_yr,y=F))+
      geom_bar(stat="identity", position=position_dodge(),fill="yellow",colour="dark blue") +
      scale_fill_brewer(palette="Paired") + theme_bw(base_size = 16)+
      geom_hline(yintercept=mean_F,colour="red",linetype=1,size=1.)+
      ylab("Fully-selected F")+  scale_x_continuous(breaks = seq(min(F_melt$F_yr), max(F_melt$F_yr), by = 4))+
      xlab("Year")
    ggsave("Fully_Selected_F.png",width=8,height=5,dpi=300,plot=p)
    
    
    # RECRUITMENT PLOTS
    # (From D. Hanselman plots)
    #---------------------------------------------------------------------------------
    
    years <- as.character(row.names(x$t.series))
    #years[seq(2,length(years),4)] <-" "
    #years[seq(3,length(years),4)] <-" "
    #years[seq(4,length(years),4)] <-" "
    
    recLY<-as.numeric(unlist(x3$t.series["Recr"]))[1:ys-1]
    rec <- as.numeric(unlist(x$t.series["Recr"]))[1:ys-1]
    #barplot(rec, names=years[1:ys-1], width=0.8,space=0.25,
    #   xlab=" ", ylab="Age 3 recruits (millions)", ylim=c(0,500))
    
    rec_comp<-data.frame(cbind(recLY,rec))
    names(rec_comp)<-c(paste0(as.numeric(tail(years,1))-1," Model"),paste0(as.numeric(tail(years,1))," Model"))
    rec_comp2<-melt(rec_comp)
    rec_comp2$Year<-c(seq(1958,as.numeric(tail(years,1))-3))
    rec_comp2<-rec_comp2[rec_comp2$Year>1976,]
    names(rec_comp2)<-c("Model.Year","Recruitment","Year")
    ### gather data for second plot
    recg<-data.frame(cbind(as.numeric(years[-length(years)]),as.numeric(rec)))
    names(recg)<-c("Year","Recruitment")
    recg$lci<-as.numeric(lci[(15+2*ys):(15+3*ys-2)])
    recg$uci<-as.numeric(uci[(15+2*ys):(15+3*ys-2)])
    
    ###### why is uci reset to max recruitment????
   # recg[recg$uci>max(recg$Recruitment),"uci"]<-max(recg$Recruitment)
    ######
    
    recg$Year<-recg$Year-2 
    
    rec_mean<-mean(recg$Recruitment[recg$Year>1976&recg$Year<recg$Year[nrow(recg)-2]])
    rec_mean_all<-mean(recg$Recruitment)
    ### 1977 on compared to prior year
    p<-ggplot(rec_comp2) + 
      geom_bar(aes(x=Year, y=Recruitment,fill=Model.Year),stat="identity", position=position_dodge(),colour="dark blue") +
      scale_fill_brewer(palette="Paired") + theme_bw(base_size = 16)+
      geom_hline(yintercept=rec_mean,colour="black",linetype=1,size=1.)+ylim(c(0,1.05*max(rec_comp2$Recruitment)))+
      ylab("Recruitment (Millions)")+  scale_x_continuous(breaks = seq(min(recg$Year), max(recg$Year), by = 4))+
      xlab("Year Class")+scale_fill_manual(values=c("dark blue","yellow"))
    ggsave("rec_plot_comparison.png",width=8,height=5,dpi=300,plot=p)
    ### all years
    p<-ggplot(recg, aes(x=Year, y=Recruitment)) + 
      geom_bar(stat="identity", position=position_dodge(),fill="yellow",colour="dark blue") +
      geom_errorbar(aes(ymin=lci, ymax=uci,  width=.8),colour="dark blue",alpha=0.42)+
      scale_fill_brewer(palette="Paired") + theme_bw(base_size = 16)+
      geom_hline(yintercept=rec_mean_all,colour="red",linetype=1,size=1.)+
      geom_hline(yintercept=rec_mean,colour="black",linetype=1,size=1.)+ #ylim(c(0,1.05*max(rec_comp2$Recruitment)))+
      ylab("Recruitment (Millions)")+  scale_x_continuous(breaks = seq(min(recg$Year), max(recg$Year), by = 4))+
      xlab("Year Class") + ylim(0,(max(recg$uci)+.01*max(recg$uci)))
    ggsave("rec_plot_all.png",width=8,height=5,dpi=300,plot=p)
    
    
    
}
