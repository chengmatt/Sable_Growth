##########################################################
# Code for running Alaska Sablefish Graphics for Stock Assessment
# Created by Dana Hanselman
# Commented by Dan Goethel
##############################################################
rm(list=(ls()))

############### INPUTS TO BE ALTERED ######################################
dir_master<-"/Users/matthewcheng/Desktop/Side projects/2023-Sablefish-SAFE/Model/2. Base (23.5)_final model"      ## location of SA Outputs
dir_R<-paste0(dir_master,'/R code',sep='')  ## location of R graphics files
SA_curr_YR<-2022 #### enter terminal year for previous stock assessment model
SA_prev_YR<-2021 #### enter terminal year for previous stock assessment model
#####################################################################################################

#########################################################################
###### Required Data Files #############################################
########################################################################
  # tem.rdat                     // current year SA report file output as rdat file
  # tem_2020.rdat                // previous year SA report file output as rdat file, "2019" should be replaced with year of prev. SA
  # sable.rep                    // ADMB output report file of current year SA
  # evalout.sdat                 // MCMC ADMB output from current year SA, need to run Bayesian version of assessment to obtain
  # IPHC_RPN_Sablefish.csv       // IPHC survey RPN from most recent year (not used in assessment)

##################################################################
###### Required R Scripts  ############################################
#####################################################################

 # sableplots.r             // main plotting function, year indicates last time updated (does not need to be updated yearly)
 # newswath.r               // plot quantiles of timeseries (from MCMC), year indicates last time updated (does not need to be updated yearly)
 # cohortfit.r              // composition plots with different colors for cohorts
 # sbcomps_pred.r           // graphing composition data

#####################################################################
##############################################
######### START CODE #####################
##########################################################
##########################################################

library(reshape2)
library(tidyverse)
library(ggsci)
library(miceadds)
library(viridis)

dir.create(paste0(dir_master,"\\plots",sep="")) # create plot folder in master directory
dir_plots<-paste0(dir_master,"\\plots",sep="")

####### load graphing functions ############################
setwd(dir_R)
source("sableplots_21.r") 
source("newswath.r")
source("cohortfit.r")
source("sbcomps_pred.r")

######## Load Data for Plotting #######################################
setwd(dir_master)
sab_curr<-dget("tem.rdat") 
sab_prev<-dget(paste0("tem_",SA_prev_YR,".rdat",sep=""))
mcmc<-dget("evalout.sdat")
write.csv(mcmc,paste0(dir_plots,'\\mcmc.csv',sep=''))

temp<-readLines("sable.rep")
B40<-as.numeric(unlist(strsplit(temp[grep("B_40",temp)+1],split=" "))) #Get B40 from report file
mean_recruit<-as.numeric(unlist(strsplit(temp[grep("Mean_Recruitment",temp)+1],split=" "))) #Get mean recruitment used for B40 calc from report file

B35<-B40*0.35/0.40
num_par<-as.numeric(unlist(strsplit(temp[grep("Number parameters estimated",temp)+1],split=" "))) #Get B40 from report file
F40<-as.numeric(unlist(strsplit(temp[grep("F_40",temp)+1],split=" "))) 
F40<-F40[1]
F.ABC<-as.numeric(unlist(strsplit(temp[grep("F_ABC",temp)+1],split=" "))) 
F.ABC<-F.ABC[1]
ABC<-as.numeric(unlist(strsplit(temp[grep("ABC for",temp)+1],split=" "))) 
ABC<-ABC[3] # first 2 are the F_ABC values, this grabs ABC for first projection year
temp.par<-readLines("tem.par")
max_grad<-unlist(strsplit(temp.par[grep("Maximum gradient",temp.par)],split=" "))
max_grad<-as.numeric(max_grad[length(max_grad)]) # grab the max grad value which is the last value of first line of par file, and turn into numeric value
model_conv<-file.exists("tem.std")

######### Graph Main Figures for SA Doc #######################################
setwd(dir_plots)
sableplots(sab_curr,sab_curr,sab_prev,"disk") # call to function created in sableplots2018.r




bs<-22
setwd(dir_master)
#x<-readLines("\\\\nmfs.local/AKC-ABL/Users2/Dana.Hanselman/My Documents/AD2017/Sablefish/16.5_new/sable.rep")
x<-readLines("sable.rep")
pos<-c(15,55,55,55,55,15,15,55,55,55,55,55,55)
labels<-c('seq(2,31)',rep("seq(41,99,by=2)",4),"seq(2,31)","seq(2,31)",rep("seq(41,99,by=2)",6))
labels2<-c('Age',rep("Length(cm)",4),"Age","Age",rep("Length(cm)",6))
labels3<-c("IFQ fishery age","Fixed gear lengths female","Fixed gear lengths male","Trawl fishery lengths female","Trawl fishery lengths male","Domestic survey ages","Cooperative survey ages","Domestic survey lengths females","Domestic survey lengths males","Cooperative survey length females","Cooperative survey lengths males","GOA trawl survey lengths females","GOA trawl survey lengths males")

p<-function(expr) eval(parse(text=expr))
spots<-grep("Pred_P_",x)
spots2<-grep("Obs_P",x)
setwd(dir_plots)

for(i in 1:length(spots)){
  obs<-x[(spots2[i]+1):(spots[i]-2)]
  obs2<-(data.frame((strsplit(obs," ")[seq(1,length(obs),by=2)]))) 
  
  
  obs3<-(obs2[3:32,])
  obs4<-obs3
  
  pred<-x[(spots[i]+1):(spots[i]+length(obs2))]
  pred2<-(data.frame((strsplit(pred," "))))
  pred3<-(pred2[3:32,])
  pred4<-pred3 
  for(j in 1:dim(obs3)[2]) {
    obs4[,j]<-as.numeric(as.character(obs3[,j]))
    pred4[,j]<-as.numeric(as.character(pred3[,j])) }
  obs5<-rowMeans(obs4);names(obs5)<-p(labels[i])
  pred5<-rowMeans(pred4); names(pred5)<-p(labels[i]) 
  b<-apply(obs4, 1, quantile, probs = c(0.05, 0.95),  na.rm = TRUE)
  d<-rbind(obs5,b)
  rownames(d)<-c("Observed","LCI","UCI")
  #melt(t(d))
  e<-data.frame(t(d),stringsAsFactors=FALSE)
  e$Age<-as.numeric(rownames(e))
  e<-(cbind(e,matrix(NA,nrow=30,ncol=1)))
  e[,5]<-pred5 
  colnames(e)[5]<-c("M_21.12")
  cols2<-rainbow(9,start=0.7)
  cols <- c("Observed"="dark green","M_21.12"="dark blue")
  ggplot(e)+geom_line(aes(x=Age,y=Observed,colour="Observed"),size=1.5)+geom_point(aes(x=Age,y=Observed),size=5,colour="dark green")+
    geom_ribbon(aes(x=Age,ymin=LCI, ymax=UCI),fill="green",alpha=0.2)+geom_path(aes(x=Age,y=M_21.12,colour="M_21.12"),size=1.4,alpha=0.5)+
    theme_bw(base_size=bs)+xlab(labels2[i])+geom_text(mapping=aes(x=pos[i],y=0.9*max(b),label=labels3[i]),size=8,colour="red")+
    ylab("Proportion")+
    scale_colour_manual(name="Source",values=cols) + theme(legend.position = c(0.9, 0.7))+
    ggtitle(paste("Aggregated observed compositions and predictions"))
  
  
  ggsave(paste("comps_2021_",i,".png"),dpi=300,width=12.6,height=6.3,units="in") }




####size age plots

clr.pos = rgb(255, 100, 100, maxColorValue=255)    # color for positive residuals
clr.neg = rgb(200, 100, 200, maxColorValue=255)    # color for positive residuals
clr.neg2 = rgb(155, 100, 255, maxColorValue=255)    # color for positive residuals
#
png(file="Fig. 3.12b. Age-length conversion matrix.png",res=500,width=8,height=8.,units="in")
par(omi=c(0,0,0,0),mfrow=c(2,2),mgp=c(2.5,1,0),mai=c(0.8,0.7,.5,0.1))
#par(mar=c(4,4,2,1),mfrow=c(2,2),oma=c(5,3,5,3), cex=1.1, cex.main=1.0, cex.lab=1.1)
z<-as.matrix(sab_curr$sizeage.f.block1)*10
x1<-as.numeric(dimnames(sab_curr$olc.fish1.f)[[2]])
y1<-as.numeric(dimnames(sab_curr$sizeage.f.block1)[[2]])
#y<-as.numeric(dimnames(z)[[2]])
x2 <- rep(y1,ncol(z))
y2 <- sort(rep(x1,nrow(z)))
plot(x2,y2,xlab="Age (yrs)",ylab="Length (cm)",type="n", las=T)
title(main="Size-Age Transition Block 1, Female")
for(p in 1:length(x2))
{ if(z[p]>0.5*max(z)) {  clr=clr.pos}
  else {if(z[p]>0.25*max(z)) {clr=clr.neg} else {clr=clr.neg2}}
  points(x2[p],y2[p],cex=z[p],col=clr,pch=16)
  points(x2[p],y2[p],cex=z[p],col=1,pch=1) }
z<-as.matrix(sab_curr$sizeage.f.block2)*10
x1<-as.numeric(dimnames(sab_curr$olc.fish1.f)[[2]])
y1<-as.numeric(dimnames(sab_curr$sizeage.f.block2)[[2]])
#y<-as.numeric(dimnames(z)[[2]])
x2 <- rep(y1,ncol(z))
y2 <- sort(rep(x1,nrow(z)))
plot(x2,y2,xlab="Age (yrs)",ylab="Length (cm)",type="n", las=T, main="Size-Age Transition Block 2, Female")
for(p in 1:length(x2))
{ if(z[p]>0.5*max(z)) {  clr=clr.pos}
  else {if(z[p]>0.25*max(z)) {clr=clr.neg} else {clr=clr.neg2}}
  points(x2[p],y2[p],cex=z[p],col=clr,pch=16)
  points(x2[p],y2[p],cex=z[p],col=1,pch=1) }
#
#par(mar=c(5,5,2,1),mfrow=c(2,1),oma=c(5,3,5,3), cex=1.1, cex.main=1.0, cex.lab=1.1)
z<-as.matrix(sab_curr$sizeage.m.block1)*10
x1<-as.numeric(dimnames(sab_curr$olc.fish1.m)[[2]])
y1<-as.numeric(dimnames(sab_curr$sizeage.m.block1)[[2]])
#y<-as.numeric(dimnames(z)[[2]])
x2 <- rep(y1,ncol(z))
y2 <- sort(rep(x1,nrow(z)))
plot(x2,y2,xlab="Age (yrs)",ylab="Length (cm)",type="n", las=T, main="Size-Age Transition Block 1, Male")
for(p in 1:length(x2))
{ if(z[p]>0.5*max(z)) {  clr=clr.pos}
  else {if(z[p]>0.25*max(z)) {clr=clr.neg} else {clr=clr.neg2}}
  points(x2[p],y2[p],cex=z[p],col=clr,pch=16)
  points(x2[p],y2[p],cex=z[p],col=1,pch=1) }
z<-as.matrix(sab_curr$sizeage.m.block2)*10
x1<-as.numeric(dimnames(sab_curr$olc.fish1.m)[[2]])
y1<-as.numeric(dimnames(sab_curr$sizeage.m.block2)[[2]])
#y<-as.numeric(dimnames(z)[[2]])
x2 <- rep(y1,ncol(z))
y2 <- sort(rep(x1,nrow(z)))
plot(x2,y2,xlab="Age (yrs)",ylab="Length (cm)",type="n", las=T, main="Size-Age Transition Block 2, Male")
for(p in 1:length(x2))
{ if(z[p]>0.5*max(z)) {  clr=clr.pos}
  else {if(z[p]>0.25*max(z)) {clr=clr.neg} else {clr=clr.neg2}}
  points(x2[p],y2[p],cex=z[p],col=clr,pch=16)
  points(x2[p],y2[p],cex=z[p],col=1,pch=1) }

dev.off()




############### getting mcmc results for table 3.14

uci <- length(mcmc$mcmc)
lci <- length(mcmc$mcmc)
for (i in 1:length(mcmc$mcmc)){
  uci[i]<-quantile(as.numeric(unlist(mcmc$mcmc[[i]])) ,0.975)
  lci[i]<-quantile(as.numeric(unlist(mcmc$mcmc[[i]])),0.025)}
write.csv(uci,"uci.csv")
write.csv(lci,"lci.csv")


### 99th confidence interval
uci <- length(mcmc$mcmc)
lci <- length(mcmc$mcmc)
for (i in 1:length(mcmc$mcmc)){
  uci[i]<-quantile(as.numeric(unlist(mcmc$mcmc[[i]])) ,0.975)
  lci[i]<-quantile(as.numeric(unlist(mcmc$mcmc[[i]])),0.001)}



#############################################################################################################################################

# Get the remaining plots (mostly data plots)

#############################################################################################################################################

load("C:\\Users\\daniel.goethel\\Desktop\\2022 SAFE\\Data\\Sablefish_Data_2022.RData")

#catch_area_hist <- read.csv("catch_by_area_1989.csv")    #should be read in during data pull but can also do it manually

# catch by gear (3.2)

catch_area_plot <- raw_catch %>%
  dplyr::mutate(Gear = dplyr::case_when(fmp_gear == "TRW" ~ "Trawl",
                                        fmp_gear == "BTR" ~ "Trawl",
                                        fmp_gear == "PTR" ~ "Trawl",
                                        fmp_gear == "NPT" ~ "Trawl",
                                        fmp_gear == "HAL" ~ "HAL",
                                        fmp_gear == "POT" ~ "Pot")) %>%
  dplyr::group_by(year, area, Gear) %>%
  dplyr::summarise(catch = sum(weight_posted) / 1000) %>%
  dplyr::bind_rows(tibble(year = 1990, fixed90_complete, Gear ="HAL")) %>%
  dplyr::bind_rows(tibble(year = 1990, trawl90_complete, Gear ="Trawl")) %>%
  dplyr::arrange(year, Gear) %>%
  dplyr::filter(!is.na(Gear))

sabl_fixed_abund_subset <- sabl_fixed_abundance %>%
  dplyr::filter(fleet == "domestic", variable == "catch") %>%
  dplyr::select(year, catch= value, gear) %>%
  dplyr::mutate(Gear = recode(gear, 'llf' = 'HAL', 'tf' = 'Trawl')) 

catch_gear_historical <- catch_area_plot %>%
  dplyr::group_by(year, Gear) %>%
  dplyr::summarize(catch = sum(catch, na.rm=T)) %>%
  dplyr::bind_rows(sabl_fixed_abund_subset) %>%
  dplyr::arrange(year, Gear) %>%
  dplyr::rename(Year = year)


catch_area_hist_filter <- catch_area_hist %>%
  dplyr::select(!Grand.total) %>%
  dplyr::rename(year = Year,BS = Bering.Sea, AI = Aleu.tians, WG = Western , CG = Central, EG = Eastern) 


catch_area_plot_final <- catch_area_plot %>%
  dplyr::group_by(year, area) %>%
  dplyr::summarize(catch = 1000*sum(catch, na.rm=T)) %>%
  tidyr::pivot_wider(names_from= area, values_from= catch, names_expand = TRUE, values_fill = 0) %>%
  dplyr::select(year,BS, AI,WG,CG,EG) %>%
  dplyr::bind_rows(catch_area_hist_filter) %>%
  dplyr::arrange(year) %>%
  tidyr::pivot_longer(cols=!year, values_to= "catch", names_to ="area") %>%
  dplyr::mutate(catch = catch/1000) %>%
  dplyr::rename(Year = year, Area = area)
  

par(omi=c(0,0.1,0,0.25),mgp=c(2.5,1,0))

catch_gear_plot <- catch_gear_historical %>% 
  ggplot(aes(x=Year,y=catch,fill=Gear)) +
  geom_bar(stat='identity',position="dodge") + 
  #facet_wrap(~gear) + 
  ggtitle("Catch by Gear Type")+
  ylab("Catch (kilotons)")+
  scale_fill_jco()

ggsave(paste0(dir_plots , "//Fig. 3.2. Catch by gear.png"),plot=catch_gear_plot) 



# catch by area (3.3)


par(omi=c(0,0.1,0,0.25),mgp=c(2.5,1,0))

catch_area_graph <- catch_area_plot_final %>% 
  ggplot(aes(x=Year,y=catch,fill=Area)) +
  geom_bar(stat='identity') + 
  #facet_wrap(~gear) + 
  ggtitle("Catch by NPFMC Area")+
  ylab("Catch (kilotons)")+
  scale_fill_jco()

ggsave(paste0(dir_plots , "//Fig. 3.3. Catch by area.png"),plot=catch_area_graph) 



# survey comparison plot (3.4)

png(filename = "Fig. 3.4. survey_trends.png", width=8, height=5,units="in",res=400)


#Plot with CI intervals as points
maxy<-1.05*max(as.numeric(unlist(sab_curr$obssrv3["obssrv3.uci"]))/mean(as.numeric(unlist(sab_curr$obssrv3["obssrv3"]))),
               as.numeric(unlist(sab_curr$obssrv5["obssrv5.uci"]))/mean(as.numeric(unlist(sab_curr$obssrv5["obssrv5"]))),
               as.numeric(unlist(sab_curr$obssrv7["obssrv7.uci"]))/mean(as.numeric(unlist(sab_curr$obssrv7["obssrv7"]))))



plot(as.numeric(row.names(sab_curr$obssrv7)), as.numeric(unlist(sab_curr$obssrv7["obssrv7.uci"]))/mean(as.numeric(unlist(sab_curr$obssrv7["obssrv7.uci"]))), 
     ylim=range(c(0,maxy)),type='n',xlab="Year", ylab="Relative Index")
lines(as.numeric(row.names(sab_curr$obssrv3)), as.numeric(unlist(sab_curr$obssrv3["obssrv3"]))/mean(as.numeric(unlist(sab_curr$obssrv3["obssrv3"]))),cex=1.5,col="blue",type='l',lwd=2,pch=16)
lines(as.numeric(row.names(sab_curr$obssrv3)), as.numeric(unlist(sab_curr$obssrv3["obssrv3"]))/mean(as.numeric(unlist(sab_curr$obssrv3["obssrv3"]))),cex=1.5,col="blue",type='p',lwd=2,pch=16)
points(as.numeric(row.names(sab_curr$obssrv5)), as.numeric(unlist(sab_curr$obssrv5["obssrv5"]))/mean(as.numeric(unlist(sab_curr$obssrv5["obssrv5"]))),cex=1.5,col="darkred",type='l',lwd=2,pch=17)
points(as.numeric(row.names(sab_curr$obssrv5)), as.numeric(unlist(sab_curr$obssrv5["obssrv5"]))/mean(as.numeric(unlist(sab_curr$obssrv5["obssrv5"]))),cex=1.5,col="darkred",type='p',lwd=2,pch=17)
points(as.numeric(row.names(sab_curr$obssrv7)), as.numeric(unlist(sab_curr$obssrv7["obssrv7"]))/mean(as.numeric(unlist(sab_curr$obssrv7["obssrv7"]))),cex=1.5,col="darkgreen",type='l',lwd=2,pch=15)
points(as.numeric(row.names(sab_curr$obssrv7)), as.numeric(unlist(sab_curr$obssrv7["obssrv7"]))/mean(as.numeric(unlist(sab_curr$obssrv7["obssrv7"]))),cex=1.5,col="darkgreen",type='p',lwd=2,pch=15)
abline(h=1,lty=2,col='black')
legend('top',legend=c("LL Survey RPNs","Fishery CPUE RPW","Trawl Survey RPW"),col=c("blue","darkred","darkgreen"),lwd=c(2,2,2),lty=c(1,1,1),pch=c(16,17,15),cex=1.25)

dev.off()


# lls by area (3.6)

par(omi=c(0,0.1,0,0.25),mgp=c(2.5,1,0))

lls_area_fig <- lls_area_rpw_apport %>%
  dplyr::select(year,area,rpn) %>%
  dplyr::rename(Year = year, Area = area) %>%
  dplyr::group_by(Year,Area) %>%
  dplyr::summarize(RPN = rpn/1000) %>%
  dplyr::mutate(Area = recode_factor(Area, 'Aleutians' = "AI" , 'Bering Sea' = 'BS', 'Western Gulf of Alaska' = 'WG',
                                     'Central Gulf of Alaska' = 'CG', 'West Yakutat' = 'WY', 'East Yakutat/Southeast' = 'EY/SE')) %>%
  ggplot(aes(x=Year,y=RPN,fill=Area)) +
  geom_bar(stat='identity') + 
  #facet_wrap(~gear) + 
  ggtitle("AFSC Longline Survey Relative Population Numbers (RPNs) by NPFMC Area")+
  ylab("RPNs (1000s)")+
  scale_fill_jco()

ggsave(paste0(dir_plots , "//Fig. 3.6. LLS by area.png"),plot=lls_area_fig) 



# lls depredation (3.7)

par(omi=c(0,0.1,0,0.25),mgp=c(2.5,1,0))


lls_dep_fig <- dplyr::full_join(lls_rpn_dep_alt,
                        lls_rpn_no_dep) %>%
  select(year,RPN,RPN_no_dep) %>%
  tidyr::pivot_longer(cols=!year, values_to= "rpn", names_to ="type") %>%
  dplyr::mutate(type = recode(type, 'RPN' = 'Corrected', 'RPN_no_dep' = 'Uncorrected')) %>%
  ggplot(aes(x=year,y=rpn, color=type))+
  geom_line(lwd= 1.4)+
  labs(y = "RPNs (1000s)", x = "Year", color = 'Whale Correction')+
  ggtitle("Impact of Whale Depredation Corrections on Longline Survey RPNs")+
  scale_color_jco()

ggsave(paste0(dir_plots , "//Fig. 3.7. LLS depredation correction.png"),plot=lls_dep_fig) 


# fishery depredation (3.8)

par(omi=c(0,0.1,0,0.25),mgp=c(2.5,1,0))

fish_whale_dep_plot <- final_dep_area %>% 
  ggplot(aes(year, pred_dep)) +
  geom_col() + 
  facet_wrap(~factor(fmp_area,levels=c('AI','BS','WG','CG','WY','SE'))) + 
  ggtitle("Fishery Whale Depredation")+
  labs(y = "Depredation (tons)", x = "Year")+
  scale_fill_jco()
  
ggsave(paste0(dir_plots , "//Fig. 3.8. Fishery Whale Dep by area.png"),plot=fish_whale_dep_plot) 


# lls abund old fish (3.51)

par(omi=c(0,0.1,0,0.25),mgp=c(2.5,1,0))

lls_age_dist_fig <- lls_age_filter %>%
  dplyr::select(year,area,Age,obs) %>%
  dplyr::group_by(year,Age) %>%
  dplyr::summarize(Freq = sum(obs)) %>%
  dplyr::mutate(Bin = dplyr::case_when(Age < 12 ~ '2-11',
                                       Age >= 12 & Age < 22 ~ '12-21',
                                       Age >= 22 ~ '22+') ,
                Bin = factor(Bin, levels = c('2-11', '12-21', '22+'), ordered = TRUE)) %>%
  dplyr::group_by(year, Bin) %>%
  dplyr::summarize(Freq = sum(Freq, na.rm=T)) %>%
  ggplot(aes(x = year, y = Freq, color = Bin)) +
  geom_line(lwd= 1.5)+
  ggtitle("Trends in Abundance by Age Group on the Domestic Longline Survey")+
  labs(y = "Number of Otoliths Aged", x = "Year", color = 'Age Group')+
  scale_color_jco()


ggsave(paste0(dir_plots , "//Fig. 3.51. LLS freq of old fish.png"),plot=lls_age_dist_fig) 


#### Calculate spawners per year class (Figure 3.50)
p_mature<-sab_curr$growthmat[,5]
wt_m<-sab_curr$growthmat[,2]
wt_f<-sab_curr$growthmat[,1]

setwd(dir_master)
library(scales)
rep_out<-readLines("sable.rep")

setwd(dir_plots)
ages21<-as.numeric(unlist(strsplit(rep_out[grep("N_proj_f",rep_out)+1],split=" "))[2:31]) # adjust number in strsplit by+3
mature<-ages21*wt_f*p_mature
years<-seq(endyr-1,endyr-30,by=-1)
labels<-mature/sum(mature)
bs<-16
spages<-data.frame(cbind(years,ages,ages21,p_mature,mature,labels))
names(spages)<-c("Year_Class","Age","2023 Abundance","Proportion Mature","SSB","Proportion SSB")
labs2<-spages$Year_Class
labs2[which(spages$Spawning_Biomass<6)]<-" "
write.csv(spages, "percent_contribution_SSB.csv")

yc_ssb_graph <- spages %>% pivot_longer(cols = 3:6, names_to = "Source", values_to = "Value")


par(omi=c(0,0.1,0,0.25),mgp=c(2.5,1,0))

ssb_cont_plot<-ggplot(yc_ssb_graph,aes(x=Year_Class,y= Value,color= Source)) +
  geom_point(cex=1.9,pch=19)+
  facet_grid(factor(Source,levels=c('Proportion Mature','2023 Abundance','SSB','Proportion SSB'))~., scales="free")+
  geom_text(data=yc_ssb_graph %>% filter(Year_Class=='2014' | Year_Class=='2016' | Year_Class=='2017' |  Year_Class=='2019'), 
            aes(label= Year_Class), cex=2., nudge_x=1.2)+
  geom_segment( aes(xend=Year_Class,yend=0))+
  ylab("")+ 
  #ylim(0,1.01)+
  scale_x_continuous(breaks = seq(min(yc_ssb_graph$Year_Class),max(yc_ssb_graph$Year_Class), by = 2))+
  xlab("Year Class")+
  ggtitle("Contribution to 2023 SSB by Year Class")+ 
  theme(legend.position ="none",axis.line=element_line(),strip.text.y = element_text(size = 7),axis.text = element_text(size = 7))+ 
  scale_fill_jco()+ scale_color_jco()

ggsave(paste0(dir_plots , "//Fig. 3.50. Percent Contr to SSB by YC.png"),plot=ssb_cont_plot) #width=8,height=8,dpi=600,




########### Get Projection Outputs for PT Final SUmmary Table (last in text table with age-4+ bio)  ###########################################

n.f<-sab_curr$natage.female
n.m<-sab_curr$natage.male
growth.mat<-sab_curr$growthmat
wt_f<-growth.mat$wt.f.block1  
wt_m<-growth.mat$wt.m.block1

#age4bio<-n.f*wt_f+n.m*wt_m
#age4bio<-age4bio[,-c(1,2)] # only age-4+
#age4bio<-rowSums(age4bio)


n.f.proj.1<-as.numeric(unlist(strsplit(temp[grep("N_proj_f",temp)+1],split=" "))[2:31]) # 1st year of projection
n.f.proj.2<-as.numeric(unlist(strsplit(temp[grep("N_proj_f",temp)+2],split=" "))[2:31]) # 2nd year of projection

n.m.proj.1<-as.numeric(unlist(strsplit(temp[grep("N_proj_m",temp)+1],split=" "))[2:31])
n.m.proj.2<-as.numeric(unlist(strsplit(temp[grep("N_proj_m",temp)+2],split=" "))[2:31])

### projected 2023 and 2024 age 4 bio
age4bio.yr1_proj<-n.f.proj.1*wt_f+n.m.proj.1*wt_m
age4bio.yr1_proj2<-age4bio.yr1_proj[-c(1,2)] #only age-4+
age4bio.yr1_proj.final<-sum(age4bio.yr1_proj2)

age4bio.yr2_proj<-n.f.proj.2*wt_f+n.m.proj.2*wt_m
age4bio.yr2_proj2<-age4bio.yr2_proj[-c(1,2)] #only age-4+
age4bio.yr2_proj.final<-sum(age4bio.yr2_proj2)


# now assign to region based on survey proportions

AI_prop <- term_apportionment[1]
BS_prop <- term_apportionment[2]
GOA_prop <- sum(term_apportionment[c(3:6)])

age4_bio_AI_proj_yr1 <-  age4bio.yr1_proj.final*AI_prop
age4_bio_AI_proj_yr2 <-  age4bio.yr2_proj.final*AI_prop

age4_bio_BS_proj_yr1 <-  age4bio.yr1_proj.final*BS_prop
age4_bio_BS_proj_yr2 <-  age4bio.yr2_proj.final*BS_prop

age4_bio_GOA_proj_yr1 <-  age4bio.yr1_proj.final*GOA_prop
age4_bio_GOA_proj_yr2 <-  age4bio.yr2_proj.final*GOA_prop



