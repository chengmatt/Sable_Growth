##################################################################################################################################################################
# Perform Francis Reweighting of Compositional Data Based on Francis (2011, 2017)
# Code adapted from Pete Hulson (pete.hulson@noaa.gov) by Dan Goethel (dallllniel.goethel@noaa.gov) for sablefish assessment model
# See associated README file for more details on methodology
##################################################################################################################################################################

#################################################################################
# File Setup
#################################################################################
### Requires a Model Folder, Input Files Folder, and Results Folder within the working directory
### Model folder should contain all files to run assessment (i.e., tpl, exe, pin, and any c++ files), except .ctl and .dat file
### The current year .ctl and .dat file should be in Input Files folder
#################################################################################

####################################################################################
# Inputs
###################################################################################
rm(list=(ls()))
library(here)
library(tidyverse)
YEAR <- 2023
ITERS<-10

ages<-seq(2,31)
nages<-length(ages)
lengths<-c(41, 43, 45, 47, 49, 51, 53, 55, 57, 59, 61, 63, 65, 67, 69, 71, 73, 75, 77, 79, 81, 83, 85, 87, 89, 91, 93, 95, 97, 99)
nleng<-length(lengths)


####################################################################################
# Setup directory structure
###################################################################################

path<-here("Prop_Within_WAA_Vary")
pathM<-paste(path,sep="")
pathD<-paste(path,"/Input Files",sep="")
pathR<-paste(path,"/Results",sep="")

###################################################################################

##################################################################################################################################################################
# Get data and report files and set up sections for ESS iteration
################################################################################################################################################################## 

DAT_init<-readLines(paste0(pathD,"/tem_",YEAR,"_na_wh.dat",sep=""),warn=FALSE)  # read in .dat file 
file.copy(from=paste0(pathD,"/tem_",YEAR,"_na_wh.dat",sep=""),to=paste0(pathM,"/tem_",YEAR,"_na_wh.dat",sep=""),overwrite=T)
file.copy(from=paste0(pathD,"/tem.pin",sep=""),to=paste0(pathM,"/tem.pin",sep=""),overwrite=T)


yrs_LLFA<-as.numeric(unlist(strsplit(DAT_init[grep("nyrs_fish_age: #-",DAT_init)+1],split=" "))) 
yrs_LLSA<-as.numeric(unlist(strsplit(DAT_init[grep("nyrs_srv1_age: #-",DAT_init)+1],split=" ")))
yrs_LLJSA<-as.numeric(unlist(strsplit(DAT_init[grep("nyrs_srv2_age: #-",DAT_init)+1],split=" ")))

yrs_LLFS<-as.numeric(unlist(strsplit(DAT_init[grep("nyrs_fish_size: #-",DAT_init)[1]+1],split=" ")))
yrs_JFS<-as.numeric(unlist(strsplit(DAT_init[grep("nyrs_fish4_size: #-",DAT_init)+1],split=" ")))
yrs_TFS<-as.numeric(unlist(strsplit(DAT_init[grep("nyrs_fish3_size: #-",DAT_init)+1],split=" ")))

yrs_JTFS<-as.numeric(unlist(strsplit(DAT_init[grep("nyrs_fish_size: #-",DAT_init)[2]+1],split=" ")))


yrs_LLSS<-as.numeric(unlist(strsplit(DAT_init[grep("nyrs_srv1_size: #-",DAT_init)+1],split=" ")))
yrs_LLJSS<-as.numeric(unlist(strsplit(DAT_init[grep("nyrs_srv2_size: #-",DAT_init)+1],split=" ")))
yrs_TSS<-as.numeric(unlist(strsplit(DAT_init[grep("nyrs_srv7_size: #-",DAT_init)+1],split=" ")))


ESS_init<-as.numeric(unlist(strsplit(DAT_init[grep("# Number of samples",DAT_init)+1],split=" ")))   # All comp sample sizes

ESS_LLFA_init<-ESS_init[c(1:yrs_LLFA)]
ESS_init<-ESS_init[-c(1:yrs_LLFA)]

ESS_LLSA_init<-ESS_init[c(1:yrs_LLSA)]
ESS_init<-ESS_init[-c(1:yrs_LLSA)]

ESS_LLJSA_init<-ESS_init[c(1:yrs_LLJSA)]
ESS_init<-ESS_init[-c(1:yrs_LLJSA)]

ESS_LLFS_init<-ESS_init[c(1:yrs_LLFS)]
ESS_init<-ESS_init[-c(1:yrs_LLFS)]

ESS_JFS_init<-ESS_init[c(1:yrs_JFS)]
ESS_init<-ESS_init[-c(1:yrs_JFS)]

ESS_TFS_init<-ESS_init[c(1:yrs_TFS)]
ESS_init<-ESS_init[-c(1:yrs_TFS)]

ESS_JTFS_init<-ESS_init[c(1:yrs_JTFS)]
ESS_init<-ESS_init[-c(1:yrs_JTFS)]

ESS_LLSS_init<-ESS_init[c(1:yrs_LLSS)]
ESS_init<-ESS_init[-c(1:yrs_LLSS)]

ESS_LLJSS_init<-ESS_init[c(1:yrs_LLJSS)]
ESS_init<-ESS_init[-c(1:yrs_LLJSS)]

ESS_TSS_init<-ESS_init[c(1:yrs_TSS)]
ESS_init<-ESS_init[-c(1:yrs_TSS)]


CTL_init<-readLines(paste(pathD,"/tem.ctl",sep=""),warn=FALSE)  # read in .ctl file 

wt_LLFA_line_f<-grep("#wt fish1 age comp female iter",CTL_init)
wt_LLFA_line_m<-grep("#wt fish1 age comp male iter",CTL_init)

wt_LLSA_line_f<-grep("#wt surv1 age comp female iter",CTL_init)
wt_LLSA_line_m<-grep("#wt surv1 age comp male iter",CTL_init)

wt_LLFS_line_m<-grep("#wt fish1 size comp male iter",CTL_init)
wt_LLFS_line_f<-grep("#wt fish1 size comp female iter",CTL_init)

wt_LLSS_line_m<-grep("#wt surv1 size comp male iter",CTL_init)
wt_LLSS_line_f<-grep("#wt surv1 size comp female iter",CTL_init)

wt_LLJSS_line_m<-grep("#wt surv2 size comp male iter",CTL_init)
wt_LLJSS_line_f<-grep("#wt surv2 size comp female iter",CTL_init)
wt_LLJSA_line_f <-grep("#wt surv2 age comp female iter",CTL_init)
wt_LLJSA_line_m <-grep("#wt surv2 age comp male iter",CTL_init)

wt_TFS_line_m<-grep("#wt fish 3 size comp male iter",CTL_init)
wt_TFS_line_f<-grep("#wt fish 3 size comp female iter",CTL_init)

wt_TSS_line_m<-grep("#wt srv7 size comp male iter",CTL_init)
wt_TSS_line_f<-grep("#wt srv7 size comp female iter",CTL_init)



# Set up results matrices
ESS_iter_LLFA_f<-matrix(nrow=ITERS+1,ncol=1)
rownames(ESS_iter_LLFA_f)<-c("Init",paste("Iter_",seq(1,ITERS),sep=""))
colnames(ESS_iter_LLFA_f)<-c("wt_LLFA_f")
ESS_iter_LLFA_f[1,1]<-1
SSQ_iter_LLFA_f<-matrix(nrow=ITERS,ncol=1)
rownames(SSQ_iter_LLFA_f)<-paste("Iter_",seq(1,ITERS),sep="")
colnames(SSQ_iter_LLFA_f)<-"SSQ_ESS_diff"

ESS_iter_LLFA_m<-matrix(nrow=ITERS+1,ncol=1)
rownames(ESS_iter_LLFA_m)<-c("Init",paste("Iter_",seq(1,ITERS),sep=""))
colnames(ESS_iter_LLFA_m)<-c("wt_LLFA_m")
ESS_iter_LLFA_m[1,1]<-1
SSQ_iter_LLFA_m<-matrix(nrow=ITERS,ncol=1)
rownames(SSQ_iter_LLFA_m)<-paste("Iter_",seq(1,ITERS),sep="")
colnames(SSQ_iter_LLFA_m)<-"SSQ_ESS_diff"

ESS_iter_LLSA_f<-matrix(nrow=ITERS+1,ncol=1)
rownames(ESS_iter_LLSA_f)<-c("Init",paste("Iter_",seq(1,ITERS),sep=""))
colnames(ESS_iter_LLSA_f)<-c("wt_LLSA_f")
ESS_iter_LLSA_f[1,1]<-1
SSQ_iter_LLSA_f<-matrix(nrow=ITERS,ncol=1)
rownames(SSQ_iter_LLSA_f)<-paste("Iter_",seq(1,ITERS),sep="")
colnames(SSQ_iter_LLSA_f)<-"SSQ_ESS_diff"

ESS_iter_LLSA_m<-matrix(nrow=ITERS+1,ncol=1)
rownames(ESS_iter_LLSA_m)<-c("Init",paste("Iter_",seq(1,ITERS),sep=""))
colnames(ESS_iter_LLSA_m)<-c("wt_LLSA_m")
ESS_iter_LLSA_m[1,1]<-1
SSQ_iter_LLSA_m<-matrix(nrow=ITERS,ncol=1)
rownames(SSQ_iter_LLSA_m)<-paste("Iter_",seq(1,ITERS),sep="")
colnames(SSQ_iter_LLSA_m)<-"SSQ_ESS_diff"

ESS_iter_LLFS_f<-matrix(nrow=ITERS+1,ncol=1)
rownames(ESS_iter_LLFS_f)<-c("Init",paste("Iter_",seq(1,ITERS),sep=""))
colnames(ESS_iter_LLFS_f)<-c("wt_LLFS_f")
ESS_iter_LLFS_f[1,1]<-1
SSQ_iter_LLFS_f<-matrix(nrow=ITERS,ncol=1)
rownames(SSQ_iter_LLFS_f)<-paste("Iter_",seq(1,ITERS),sep="")
colnames(SSQ_iter_LLFS_f)<-"SSQ_ESS_diff"

ESS_iter_LLFS_m<-matrix(nrow=ITERS+1,ncol=1)
rownames(ESS_iter_LLFS_m)<-c("Init",paste("Iter_",seq(1,ITERS),sep=""))
colnames(ESS_iter_LLFS_m)<-c("wt_LLFS_m")
ESS_iter_LLFS_m[1,1]<-1
SSQ_iter_LLFS_m<-matrix(nrow=ITERS,ncol=1)
rownames(SSQ_iter_LLFS_m)<-paste("Iter_",seq(1,ITERS),sep="")
colnames(SSQ_iter_LLFS_m)<-"SSQ_ESS_diff"

ESS_iter_LLSS_f<-matrix(nrow=ITERS+1,ncol=1)
rownames(ESS_iter_LLSS_f)<-c("Init",paste("Iter_",seq(1,ITERS),sep=""))
colnames(ESS_iter_LLSS_f)<-c("wt_LLSS_f")
ESS_iter_LLSS_f[1,1]<-1
SSQ_iter_LLSS_f<-matrix(nrow=ITERS,ncol=1)
rownames(SSQ_iter_LLSS_f)<-paste("Iter_",seq(1,ITERS),sep="")
colnames(SSQ_iter_LLSS_f)<-"SSQ_ESS_diff"

ESS_iter_LLSS_m<-matrix(nrow=ITERS+1,ncol=1)
rownames(ESS_iter_LLSS_m)<-c("Init",paste("Iter_",seq(1,ITERS),sep=""))
colnames(ESS_iter_LLSS_m)<-c("wt_LLSS_m")
ESS_iter_LLSS_m[1,1]<-1
SSQ_iter_LLSS_m<-matrix(nrow=ITERS,ncol=1)
rownames(SSQ_iter_LLSS_m)<-paste("Iter_",seq(1,ITERS),sep="")
colnames(SSQ_iter_LLSS_m)<-"SSQ_ESS_diff"

ESS_iter_LLJSS_f<-matrix(nrow=ITERS+1,ncol=1)
rownames(ESS_iter_LLJSS_f)<-c("Init",paste("Iter_",seq(1,ITERS),sep=""))
colnames(ESS_iter_LLJSS_f)<-c("wt_LLJSS_f")
ESS_iter_LLJSS_f[1,1]<-1
SSQ_iter_LLJSS_f<-matrix(nrow=ITERS,ncol=1)
rownames(SSQ_iter_LLJSS_f)<-paste("Iter_",seq(1,ITERS),sep="")
colnames(SSQ_iter_LLJSS_f)<-"SSQ_ESS_diff"

ESS_iter_LLJSS_m<-matrix(nrow=ITERS+1,ncol=1)
rownames(ESS_iter_LLJSS_m)<-c("Init",paste("Iter_",seq(1,ITERS),sep=""))
colnames(ESS_iter_LLJSS_m)<-c("wt_LLJSS_m")
ESS_iter_LLJSS_m[1,1]<-1
SSQ_iter_LLJSS_m<-matrix(nrow=ITERS,ncol=1)
rownames(SSQ_iter_LLJSS_m)<-paste("Iter_",seq(1,ITERS),sep="")
colnames(SSQ_iter_LLJSS_m)<-"SSQ_ESS_diff"

ESS_iter_LLJSA_f<-matrix(nrow=ITERS+1,ncol=1)
rownames(ESS_iter_LLJSA_f)<-c("Init",paste("Iter_",seq(1,ITERS),sep=""))
colnames(ESS_iter_LLJSA_f)<-c("wt_LLJSA_f")
ESS_iter_LLJSA_f[1,1]<-1
SSQ_iter_LLJSA_f<-matrix(nrow=ITERS,ncol=1)
rownames(SSQ_iter_LLJSA_f)<-paste("Iter_",seq(1,ITERS),sep="")
colnames(SSQ_iter_LLJSA_f)<-"SSQ_ESS_diff"

ESS_iter_LLJSA_m<-matrix(nrow=ITERS+1,ncol=1)
rownames(ESS_iter_LLJSA_m)<-c("Init",paste("Iter_",seq(1,ITERS),sep=""))
colnames(ESS_iter_LLJSA_m)<-c("wt_LLJSA_m")
ESS_iter_LLJSA_m[1,1]<-1
SSQ_iter_LLJSA_m<-matrix(nrow=ITERS,ncol=1)
rownames(SSQ_iter_LLJSA_m)<-paste("Iter_",seq(1,ITERS),sep="")
colnames(SSQ_iter_LLJSA_m)<-"SSQ_ESS_diff"

ESS_iter_TFS_f<-matrix(nrow=ITERS+1,ncol=1)
rownames(ESS_iter_TFS_f)<-c("Init",paste("Iter_",seq(1,ITERS),sep=""))
colnames(ESS_iter_TFS_f)<-c("wt_TFS_f")
ESS_iter_TFS_f[1,1]<-1
SSQ_iter_TFS_f<-matrix(nrow=ITERS,ncol=1)
rownames(SSQ_iter_TFS_f)<-paste("Iter_",seq(1,ITERS),sep="")
colnames(SSQ_iter_TFS_f)<-"SSQ_ESS_diff"

ESS_iter_TFS_m<-matrix(nrow=ITERS+1,ncol=1)
rownames(ESS_iter_TFS_m)<-c("Init",paste("Iter_",seq(1,ITERS),sep=""))
colnames(ESS_iter_TFS_m)<-c("wt_TFS_m")
ESS_iter_TFS_m[1,1]<-1
SSQ_iter_TFS_m<-matrix(nrow=ITERS,ncol=1)
rownames(SSQ_iter_TFS_m)<-paste("Iter_",seq(1,ITERS),sep="")
colnames(SSQ_iter_TFS_m)<-"SSQ_ESS_diff"


ESS_iter_TSS_f<-matrix(nrow=ITERS+1,ncol=1)
rownames(ESS_iter_TSS_f)<-c("Init",paste("Iter_",seq(1,ITERS),sep=""))
colnames(ESS_iter_TSS_f)<-c("wt_TSS_f")
ESS_iter_TSS_f[1,1]<-1
SSQ_iter_TSS_f<-matrix(nrow=ITERS,ncol=1)
rownames(SSQ_iter_TSS_f)<-paste("Iter_",seq(1,ITERS),sep="")
colnames(SSQ_iter_TSS_f)<-"SSQ_ESS_diff"

ESS_iter_TSS_m<-matrix(nrow=ITERS+1,ncol=1)
rownames(ESS_iter_TSS_m)<-c("Init",paste("Iter_",seq(1,ITERS),sep=""))
colnames(ESS_iter_TSS_m)<-c("wt_TSS_m")
ESS_iter_TSS_m[1,1]<-1
SSQ_iter_TSS_m<-matrix(nrow=ITERS,ncol=1)
rownames(SSQ_iter_TSS_m)<-paste("Iter_",seq(1,ITERS),sep="")
colnames(SSQ_iter_TSS_m)<-"SSQ_ESS_diff"


# Write new CTL file with iterated weights
CTL<-c(
  CTL_init[1:(wt_LLFA_line_f-1)],
  paste0(1," ","#wt fish1 age comp female iter",sep=""),
  paste0(1," ","#wt fish1 age comp male iter",sep=""),
  paste0(1," ","#wt surv1 age comp female iter",sep=""),
  paste0(1," ","#wt surv1 age comp male iter",sep=""),
  paste0(1," ","#wt surv2 age comp female iter",sep=""),
  paste0(1," ","#wt surv2 age comp male iter",sep=""),
  paste0(1," ","#wt fish1 size comp male iter",sep=""),
  paste0(1," ","#wt fish1 size comp female iter",sep=""),
  paste0(1," ","#wt surv1 size comp male iter",sep=""),
  paste0(1," ","#wt surv1 size comp female iter",sep=""),
  paste0(1," ","#wt surv2 size comp male iter",sep=""),
  paste0(1," ","#wt surv2 size comp female iter",sep=""),
  paste0(1," ","#wt fish 3 size comp male iter",sep=""),
  paste0(1," ","#wt fish 3 size comp female iter",sep=""),
  paste0(1," ","#wt srv7 size comp male iter",sep=""),
  paste0(1," ","#wt srv7 size comp female iter",sep=""),
  CTL_init[(wt_TSS_line_f+1):length(CTL_init)])


do_rwt<-grep("data_reweight_switch",CTL_init)              # make sure reweight strip is set to 1
CTL[do_rwt+1]<-1

write.table(CTL,file=paste(pathM,"/tem.ctl",sep=""),quote=F,row.names=F,col.names=F)


##################################################################################################################################################################
# Perform Analyses
################################################################################################################################################################## 

start_time<-Sys.time()

# Start iteration loop
for (i in 1:ITERS){

print(paste("i = ",i,sep=""))

# Run model
setwd(pathM)
# Run model
setwd(pathM)
# shell("tem.EXE -nohess")
# COMPILE ADMB HERE
# Set path to admb terminal tools!
Sys.setenv(PATH=paste(Sys.getenv("PATH"),
                      "/Users/matthewcheng/Desktop/ADMBTerminal.app/admb-12.3",
                      sep=":"))

# Compile admb template file
# system(paste("admb tem"))
system(paste("./tem"))
# Read in report file and get output ESSs
REP<-readLines(paste(pathM,"/tem.rep",sep=""),warn=FALSE)


line1<-grep("Obs_P_fish_age Female",REP)+1
line2<-grep("Pred_P_fish1_age Female",REP)-2
fy<-length(REP[line1:line2])
obs_fa<-REP[line1:line2]
obs_p_fa<-matrix(nrow=fy,ncol=nages+1)
line1<-grep("Pred_P_fish1_age Female",REP)+1
line2<-grep("Obs_P_fish_age Male",REP)-2
pred_fa<-REP[line1:line2]
pred_p_fa<-matrix(nrow=fy,ncol=nages+1)
v_y<-vector(length=fy)
w_denom<-vector(length=fy)
obs_bar<-vector(length=fy)
pred_bar<-vector(length=fy)

for(y in 1:fy){
pred_y<-as.numeric(strsplit(pred_fa[y]," ")[[1]])[3:(length(ages)+2)]
pred_bar[y]<-sum(ages*pred_y)
obs_y<-as.numeric(strsplit(obs_fa[y]," ")[[1]])[3:(length(ages)+2)]
obs_bar[y]<-sum(ages*obs_y)
v_y[y]<-sum(ages^2*pred_y)-pred_bar[y]^2
w_denom[y]<-(obs_bar[y]-pred_bar[y])/sqrt(v_y[y]/ESS_LLFA_init[y])
}

wt_LLFA_f<-1/var(w_denom)
ESS_iter_LLFA_f[i+1,1]<-wt_LLFA_f
SSQ_iter_LLFA_f[i]<-sum((ESS_iter_LLFA_f[i+1,1]-ESS_iter_LLFA_f[i,1])^2)


line1<-grep("Obs_P_fish_age Male",REP)+1
line2<-grep("Pred_P_fish1_age Male",REP)-2
fy<-length(REP[line1:line2])
obs_fa<-REP[line1:line2]
obs_p_fa<-matrix(nrow=fy,ncol=nages+1)
line1<-grep("Pred_P_fish1_age Male",REP)+1
line2<-grep("Obs_P_fish1_size Female",REP)-2
pred_fa<-REP[line1:line2]
pred_p_fa<-matrix(nrow=fy,ncol=nages+1)
v_y<-vector(length=fy)
w_denom<-vector(length=fy)
obs_bar<-vector(length=fy)
pred_bar<-vector(length=fy)

for(y in 1:fy){
  pred_y<-as.numeric(strsplit(pred_fa[y]," ")[[1]])[3:(length(ages)+2)]
  pred_bar[y]<-sum(ages*pred_y)
  obs_y<-as.numeric(strsplit(obs_fa[y]," ")[[1]])[3:(length(ages)+2)]
  obs_bar[y]<-sum(ages*obs_y)
  v_y[y]<-sum(ages^2*pred_y)-pred_bar[y]^2
  w_denom[y]<-(obs_bar[y]-pred_bar[y])/sqrt(v_y[y]/ESS_LLFA_init[y])
}

wt_LLFA_m<-1/var(w_denom)
ESS_iter_LLFA_m[i+1,1]<-wt_LLFA_m
SSQ_iter_LLFA_m[i]<-sum((ESS_iter_LLFA_m[i+1,1]-ESS_iter_LLFA_m[i,1])^2)

line1<-grep("Obs_P_fish1_size Female",REP)+1
line2<-grep("Pred_P_fish1_size",REP)-2
line2<-line2[1]
fy<-length(REP[line1:line2])
obs_LLfs_f<-REP[line1:line2]
obs_p_LLfs_f<-matrix(nrow=fy,ncol=nleng+1)
line1<-grep("Pred_P_fish1_size",REP)+1
line1<-line1[1]
line2<-grep("Obs_P_fish1_size Male",REP)-2
pred_LLfs_f<-REP[line1:line2]
pred_p_LLfs_f<-matrix(nrow=fy,ncol=nleng+1)
v_y<-vector(length=fy)
w_denom<-vector(length=fy)
obs_bar<-vector(length=fy)
pred_bar<-vector(length=fy)

for(y in 1:fy){
  pred_y<-as.numeric(strsplit(pred_LLfs_f[y]," ")[[1]])[3:(length(lengths)+2)]
  pred_bar[y]<-sum(lengths*pred_y)
  obs_y<-as.numeric(strsplit(obs_LLfs_f[y]," ")[[1]])[3:(length(lengths)+2)]
  obs_bar[y]<-sum(lengths*obs_y)
  v_y[y]<-sum(lengths^2*pred_y)-pred_bar[y]^2
  w_denom[y]<-(obs_bar[y]-pred_bar[y])/sqrt(v_y[y]/ESS_LLFS_init[y])
}

wt_LLFS_f<-1/var(w_denom)
ESS_iter_LLFS_f[i+1,1]<-wt_LLFS_f
SSQ_iter_LLFS_f[i]<-sum((ESS_iter_LLFS_f[i+1,1]-ESS_iter_LLFS_f[i,1])^2)


line1<-grep("Obs_P_fish1_size Male",REP)+1
line2<-grep("Pred_P_fish1_size",REP)-2
line2<-line2[2]
fy<-length(REP[line1:line2])
obs_LLfs_m<-REP[line1:line2]
obs_p_LLfs_m<-matrix(nrow=fy,ncol=nleng+1)
line1<-grep("Pred_P_fish1_size",REP)+1
line1<-line1[2]
line2<-grep("Obs_P_fish3_size Males",REP)-2
pred_LLfs_m<-REP[line1:line2]
pred_p_LLfs_m<-matrix(nrow=fy,ncol=nleng+1)
v_y<-vector(length=fy)
w_denom<-vector(length=fy)
obs_bar<-vector(length=fy)
pred_bar<-vector(length=fy)

for(y in 1:fy){
  pred_y<-as.numeric(strsplit(pred_LLfs_m[y]," ")[[1]])[3:(length(lengths)+2)]
  pred_bar[y]<-sum(lengths*pred_y)
  obs_y<-as.numeric(strsplit(obs_LLfs_m[y]," ")[[1]])[3:(length(lengths)+2)]
  obs_bar[y]<-sum(lengths*obs_y)
  v_y[y]<-sum(lengths^2*pred_y)-pred_bar[y]^2
  w_denom[y]<-(obs_bar[y]-pred_bar[y])/sqrt(v_y[y]/ESS_LLFS_init[y])
}

wt_LLFS_m<-1/var(w_denom)
ESS_iter_LLFS_m[i+1,1]<-wt_LLFS_m
SSQ_iter_LLFS_m[i]<-sum((ESS_iter_LLFS_m[i+1,1]-ESS_iter_LLFS_m[i,1])^2)


line1<-grep("Obs_P_fish3_size Female",REP)+1
line2<-grep("Pred_P_fish3_size",REP)-2
line2<-line2[2]
fy<-length(REP[line1:line2])
obs_Tfs_f<-REP[line1:line2]
obs_p_Tfs_f<-matrix(nrow=fy,ncol=nleng+1)
line1<-grep("Pred_P_fish3_size",REP)+1
line1<-line1[2]
line2<-grep("Obs_P_srv1_age",REP)-2
pred_Tfs_f<-REP[line1:line2]
pred_p_Tfs_f<-matrix(nrow=fy,ncol=nleng+1)
v_y<-vector(length=fy)
w_denom<-vector(length=fy)
obs_bar<-vector(length=fy)
pred_bar<-vector(length=fy)

for(y in 1:fy){
  pred_y<-as.numeric(strsplit(pred_Tfs_f[y]," ")[[1]])[3:(length(lengths)+2)]
  pred_bar[y]<-sum(lengths*pred_y)
  obs_y<-as.numeric(strsplit(obs_Tfs_f[y]," ")[[1]])[3:(length(lengths)+2)]
  obs_bar[y]<-sum(lengths*obs_y)
  v_y[y]<-sum(lengths^2*pred_y)-pred_bar[y]^2
  w_denom[y]<-(obs_bar[y]-pred_bar[y])/sqrt(v_y[y]/ESS_TFS_init[y])
}

wt_TFS_f<-1/var(w_denom)
ESS_iter_TFS_f[i+1,1]<-wt_TFS_f
SSQ_iter_TFS_f[i]<-sum((ESS_iter_TFS_f[i+1,1]-ESS_iter_TFS_f[i,1])^2)


line1<-grep("Obs_P_fish3_size Males",REP)+1
line2<-grep("Pred_P_fish3_size",REP)-2
line2<-line2[1]
fy<-length(REP[line1:line2])
obs_Tfs_m<-REP[line1:line2]
obs_p_Tfs_m<-matrix(nrow=fy,ncol=nleng+1)
line1<-grep("Pred_P_fish3_size",REP)+1
line1<-line1[1]
line2<-grep("Obs_P_fish3_size Female",REP)-2
pred_Tfs_m<-REP[line1:line2]
pred_p_Tfs_m<-matrix(nrow=fy,ncol=nleng+1)
v_y<-vector(length=fy)
w_denom<-vector(length=fy)
obs_bar<-vector(length=fy)
pred_bar<-vector(length=fy)

for(y in 1:fy){
  pred_y<-as.numeric(strsplit(pred_Tfs_m[y]," ")[[1]])[3:(length(lengths)+2)]
  pred_bar[y]<-sum(lengths*pred_y)
  obs_y<-as.numeric(strsplit(obs_Tfs_m[y]," ")[[1]])[3:(length(lengths)+2)]
  obs_bar[y]<-sum(lengths*obs_y)
  v_y[y]<-sum(lengths^2*pred_y)-pred_bar[y]^2
  w_denom[y]<-(obs_bar[y]-pred_bar[y])/sqrt(v_y[y]/ESS_TFS_init[y])
}

wt_TFS_m<-1/var(w_denom)
ESS_iter_TFS_m[i+1,1]<-wt_TFS_m
SSQ_iter_TFS_m[i]<-sum((ESS_iter_TFS_m[i+1,1]-ESS_iter_TFS_m[i,1])^2)


line1<-grep("Obs_P_srv1_age Female",REP)+1
line2<-grep("Pred_P_srv1_age Female",REP)-2
sy<-length(REP[line1:line2])
obs_LLsa<-REP[line1:line2]
obs_p_LLsa<-matrix(nrow=sy,ncol=nages+1)
line1<-grep("Pred_P_srv1_age Female",REP)+1
line2<-grep("Obs_P_srv1_age Male",REP)-2
pred_LLsa<-REP[line1:line2]
pred_p_LLsa<-matrix(nrow=sy,ncol=nages+1)
v_y<-vector(length=sy)
w_denom<-vector(length=sy)
obs_bar<-vector(length=sy)
pred_bar<-vector(length=sy)

for(y in 1:sy){
  pred_y<-as.numeric(strsplit(pred_LLsa[y]," ")[[1]])[3:(length(ages)+2)]
  pred_bar[y]<-sum(ages*pred_y)
  obs_y<-as.numeric(strsplit(obs_LLsa[y]," ")[[1]])[3:(length(ages)+2)]
  obs_bar[y]<-sum(ages*obs_y)
  v_y[y]<-sum(ages^2*pred_y)-pred_bar[y]^2
  w_denom[y]<-(obs_bar[y]-pred_bar[y])/sqrt(v_y[y]/ESS_LLSA_init[y])
}

wt_LLSA_f<-1/var(w_denom)
ESS_iter_LLSA_f[i+1,1]<-wt_LLSA_f
SSQ_iter_LLSA_f[i]<-sum((ESS_iter_LLSA_f[i+1,1]-ESS_iter_LLSA_f[i,1])^2)


line1<-grep("Obs_P_srv1_age Male",REP)+1
line2<-grep("Pred_P_srv1_age Male",REP)-2
sy<-length(REP[line1:line2])
obs_LLsa<-REP[line1:line2]
obs_p_LLsa<-matrix(nrow=sy,ncol=nages+1)
line1<-grep("Pred_P_srv1_age Male",REP)+1
line2<-grep("Obs_P_srv2_age Female",REP)-2
pred_LLsa<-REP[line1:line2]
pred_p_LLsa<-matrix(nrow=sy,ncol=nages+1)
v_y<-vector(length=sy)
w_denom<-vector(length=sy)
obs_bar<-vector(length=sy)
pred_bar<-vector(length=sy)

for(y in 1:sy){
  pred_y<-as.numeric(strsplit(pred_LLsa[y]," ")[[1]])[3:(length(ages)+2)]
  pred_bar[y]<-sum(ages*pred_y)
  obs_y<-as.numeric(strsplit(obs_LLsa[y]," ")[[1]])[3:(length(ages)+2)]
  obs_bar[y]<-sum(ages*obs_y)
  v_y[y]<-sum(ages^2*pred_y)-pred_bar[y]^2
  w_denom[y]<-(obs_bar[y]-pred_bar[y])/sqrt(v_y[y]/ESS_LLSA_init[y])
}

wt_LLSA_m<-1/var(w_denom)
ESS_iter_LLSA_m[i+1,1]<-wt_LLSA_m
SSQ_iter_LLSA_m[i]<-sum((ESS_iter_LLSA_m[i+1,1]-ESS_iter_LLSA_m[i,1])^2)

line1<-grep("Obs_P_srv2_age Female",REP)+1
line2<-grep("Pred_P_srv2_age Female",REP)-2
sy<-length(REP[line1:line2])
obs_LLJsa<-REP[line1:line2]
obs_p_LLJsa<-matrix(nrow=sy,ncol=nages+1)
line1<-grep("Pred_P_srv2_age Female",REP)+1
line2<-grep("Obs_P_srv2_age Male",REP)-2
line2<-line2[1]
pred_LLJsa<-REP[line1:line2]
pred_p_LLJsa<-matrix(nrow=sy,ncol=nages+1)
v_y<-vector(length=sy)
w_denom<-vector(length=sy)
obs_bar<-vector(length=sy)
pred_bar<-vector(length=sy)

for(y in 1:sy){
  pred_y<-as.numeric(strsplit(pred_LLJsa[y]," ")[[1]])[3:(length(ages)+2)]
  pred_bar[y]<-sum(ages*pred_y)
  obs_y<-as.numeric(strsplit(obs_LLJsa[y]," ")[[1]])[3:(length(ages)+2)]
  obs_bar[y]<-sum(ages*obs_y)
  v_y[y]<-sum(ages^2*pred_y)-pred_bar[y]^2
  w_denom[y]<-(obs_bar[y]-pred_bar[y])/sqrt(v_y[y]/ESS_LLJSA_init[y])
}

wt_LLJSA_f<-1/var(w_denom)
ESS_iter_LLJSA_f[i+1,1]<-wt_LLJSA_f
SSQ_iter_LLJSA_f[i]<-sum((ESS_iter_LLJSA_f[i+1,1]-ESS_iter_LLJSA_f[i,1])^2)

line1<-grep("Obs_P_srv2_age Male",REP)+1
line2<-grep("Pred_P_srv2_age Male",REP)-2
sy<-length(REP[line1:line2])
obs_LLJsa<-REP[line1:line2]
obs_p_LLJsa<-matrix(nrow=sy,ncol=nages+1)
line1<-grep("Pred_P_srv2_age Male",REP)+1
line2<-grep("Obs_P_srv1_size",REP)-2
line2<-line2[1]
pred_LLJsa<-REP[line1:line2]
pred_p_LLJsa<-matrix(nrow=sy,ncol=nages+1)
v_y<-vector(length=sy)
w_denom<-vector(length=sy)
obs_bar<-vector(length=sy)
pred_bar<-vector(length=sy)

for(y in 1:sy){
  pred_y<-as.numeric(strsplit(pred_LLJsa[y]," ")[[1]])[3:(length(ages)+2)]
  pred_bar[y]<-sum(ages*pred_y)
  obs_y<-as.numeric(strsplit(obs_LLJsa[y]," ")[[1]])[3:(length(ages)+2)]
  obs_bar[y]<-sum(ages*obs_y)
  v_y[y]<-sum(ages^2*pred_y)-pred_bar[y]^2
  w_denom[y]<-(obs_bar[y]-pred_bar[y])/sqrt(v_y[y]/ESS_LLJSA_init[y])
}

wt_LLJSA_m<-1/var(w_denom)
ESS_iter_LLJSA_m[i+1,1]<-wt_LLJSA_m
SSQ_iter_LLJSA_m[i]<-sum((ESS_iter_LLJSA_m[i+1,1]-ESS_iter_LLJSA_m[i,1])^2)


line1<-grep("Obs_P_srv1_size",REP)+1
line1<-line1[1]
line2<-grep("Pred_P_srv1_size",REP)-2
line2<-line2[1]
fy<-length(REP[line1:line2])
obs_LLss_f<-REP[line1:line2]
obs_p_LLss_f<-matrix(nrow=fy,ncol=nleng+1)
line1<-grep("Pred_P_srv1_size",REP)+1
line1<-line1[1]
line2<-grep("Obs_P_srv1_size Males",REP)-2
pred_LLss_f<-REP[line1:line2]
pred_p_LLss_f<-matrix(nrow=fy,ncol=nleng+1)
v_y<-vector(length=fy)
w_denom<-vector(length=fy)
obs_bar<-vector(length=fy)
pred_bar<-vector(length=fy)

for(y in 1:fy){
  pred_y<-as.numeric(strsplit(pred_LLss_f[y]," ")[[1]])[3:(length(lengths)+2)]
  pred_bar[y]<-sum(lengths*pred_y)
  obs_y<-as.numeric(strsplit(obs_LLss_f[y]," ")[[1]])[3:(length(lengths)+2)]
  obs_bar[y]<-sum(lengths*obs_y)
  v_y[y]<-sum(lengths^2*pred_y)-pred_bar[y]^2
  w_denom[y]<-(obs_bar[y]-pred_bar[y])/sqrt(v_y[y]/ESS_LLSS_init[y])
}

wt_LLSS_f<-1/var(w_denom)
ESS_iter_LLSS_f[i+1,1]<-wt_LLSS_f
SSQ_iter_LLSS_f[i]<-sum((ESS_iter_LLSS_f[i+1,1]-ESS_iter_LLSS_f[i,1])^2)


line1<-grep("Obs_P_srv1_size Males",REP)+1
line2<-grep("Pred_P_srv1_size",REP)-2
line2<-line2[2]
fy<-length(REP[line1:line2])
obs_LLss_m<-REP[line1:line2]
obs_p_LLss_m<-matrix(nrow=fy,ncol=nleng+1)
line1<-grep("Pred_P_srv1_size",REP)+1
line1<-line1[2]
line2<-grep("Obs_P_srv2_size",REP)-2
line2<-line2[1]
pred_LLss_m<-REP[line1:line2]
pred_p_LLss_m<-matrix(nrow=fy,ncol=nleng+1)
v_y<-vector(length=fy)
w_denom<-vector(length=fy)
obs_bar<-vector(length=fy)
pred_bar<-vector(length=fy)

for(y in 1:fy){
  pred_y<-as.numeric(strsplit(pred_LLss_m[y]," ")[[1]])[3:(length(lengths)+2)]
  pred_bar[y]<-sum(lengths*pred_y)
  obs_y<-as.numeric(strsplit(obs_LLss_m[y]," ")[[1]])[3:(length(lengths)+2)]
  obs_bar[y]<-sum(lengths*obs_y)
  v_y[y]<-sum(lengths^2*pred_y)-pred_bar[y]^2
  w_denom[y]<-(obs_bar[y]-pred_bar[y])/sqrt(v_y[y]/ESS_LLSS_init[y])
}

wt_LLSS_m<-1/var(w_denom)
ESS_iter_LLSS_m[i+1,1]<-wt_LLSS_m
SSQ_iter_LLSS_m[i]<-sum((ESS_iter_LLSS_m[i+1,1]-ESS_iter_LLSS_m[i,1])^2)


line1<-grep("Obs_P_srv2_size",REP)+1
line1<-line1[1]
line2<-grep("Pred_P_srv2_size",REP)-2
line2<-line2[1]
fy<-length(REP[line1:line2])
obs_LLJss_f<-REP[line1:line2]
obs_p_LLJss_f<-matrix(nrow=fy,ncol=nleng+1)
line1<-grep("Pred_P_srv2_size",REP)+1
line1<-line1[1]
line2<-grep("Obs_P_srv2_size Males",REP)-2
pred_LLJss_f<-REP[line1:line2]
pred_p_LLJss_f<-matrix(nrow=fy,ncol=nleng+1)
v_y<-vector(length=fy)
w_denom<-vector(length=fy)
obs_bar<-vector(length=fy)
pred_bar<-vector(length=fy)

for(y in 1:fy){
  pred_y<-as.numeric(strsplit(pred_LLJss_f[y]," ")[[1]])[3:(length(lengths)+2)]
  pred_bar[y]<-sum(lengths*pred_y)
  obs_y<-as.numeric(strsplit(obs_LLJss_f[y]," ")[[1]])[3:(length(lengths)+2)]
  obs_bar[y]<-sum(lengths*obs_y)
  v_y[y]<-sum(lengths^2*pred_y)-pred_bar[y]^2
  w_denom[y]<-(obs_bar[y]-pred_bar[y])/sqrt(v_y[y]/ESS_LLJSS_init[y])
}

wt_LLJSS_f<-1/var(w_denom)
ESS_iter_LLJSS_f[i+1,1]<-wt_LLJSS_f
SSQ_iter_LLJSS_f[i]<-sum((ESS_iter_LLJSS_f[i+1,1]-ESS_iter_LLJSS_f[i,1])^2)


line1<-grep("Obs_P_srv2_size Males",REP)+1
line2<-grep("Pred_P_srv2_size Males",REP)-2
fy<-length(REP[line1:line2])
obs_LLJss_m<-REP[line1:line2]
obs_p_LLJss_m<-matrix(nrow=fy,ncol=nleng+1)
line1<-grep("Pred_P_srv2_size Males",REP)+1
line2<-grep("Obs_P_srv7_size",REP)-2
line2<-line2[1]
pred_LLJss_m<-REP[line1:line2]
pred_p_LLJss_m<-matrix(nrow=fy,ncol=nleng+1)
v_y<-vector(length=fy)
w_denom<-vector(length=fy)
obs_bar<-vector(length=fy)
pred_bar<-vector(length=fy)

for(y in 1:fy){
  pred_y<-as.numeric(strsplit(pred_LLJss_m[y]," ")[[1]])[3:(length(lengths)+2)]
  pred_bar[y]<-sum(lengths*pred_y)
  obs_y<-as.numeric(strsplit(obs_LLJss_m[y]," ")[[1]])[3:(length(lengths)+2)]
  obs_bar[y]<-sum(lengths*obs_y)
  v_y[y]<-sum(lengths^2*pred_y)-pred_bar[y]^2
  w_denom[y]<-(obs_bar[y]-pred_bar[y])/sqrt(v_y[y]/ESS_LLJSS_init[y])
}

wt_LLJSS_m<-1/var(w_denom)
ESS_iter_LLJSS_m[i+1,1]<-wt_LLJSS_m
SSQ_iter_LLJSS_m[i]<-sum((ESS_iter_LLJSS_m[i+1,1]-ESS_iter_LLJSS_m[i,1])^2)


line1<-grep("Obs_P_srv7_size",REP)+1
line1<-line1[1]
line2<-grep("Pred_P_srv7_size",REP)-2
line2<-line2[1]
fy<-length(REP[line1:line2])
obs_Tss_f<-REP[line1:line2]
obs_p_Tss_f<-matrix(nrow=fy,ncol=nleng+1)
line1<-grep("Pred_P_srv7_size",REP)+1
line1<-line1[1]
line2<-grep("Obs_P_srv7_size Males",REP)-2
pred_Tss_f<-REP[line1:line2]
pred_p_Tss_f<-matrix(nrow=fy,ncol=nleng+1)
v_y<-vector(length=fy)
w_denom<-vector(length=fy)
obs_bar<-vector(length=fy)
pred_bar<-vector(length=fy)

for(y in 1:fy){
  pred_y<-as.numeric(strsplit(pred_Tss_f[y]," ")[[1]])[3:(length(lengths)+2)]
  pred_bar[y]<-sum(lengths*pred_y)
  obs_y<-as.numeric(strsplit(obs_Tss_f[y]," ")[[1]])[3:(length(lengths)+2)]
  obs_bar[y]<-sum(lengths*obs_y)
  v_y[y]<-sum(lengths^2*pred_y)-pred_bar[y]^2
  w_denom[y]<-(obs_bar[y]-pred_bar[y])/sqrt(v_y[y]/ESS_TSS_init[y])
}

wt_TSS_f<-1/var(w_denom)
ESS_iter_TSS_f[i+1,1]<-wt_TSS_f
SSQ_iter_TSS_f[i]<-sum((ESS_iter_TSS_f[i+1,1]-ESS_iter_TSS_f[i,1])^2)


line1<-grep("Obs_P_srv7_size Males",REP)+1
line2<-grep("Pred_P_srv7_size Males",REP)-2
fy<-length(REP[line1:line2])
obs_Tss_m<-REP[line1:line2]
obs_p_Tss_m<-matrix(nrow=fy,ncol=nleng+1)
line1<-grep("Pred_P_srv7_size Males",REP)+1
line2<-grep("Survey Biomass",REP)-2
line2<-line2[1]
pred_Tss_m<-REP[line1:line2]
pred_p_Tss_m<-matrix(nrow=fy,ncol=nleng+1)
v_y<-vector(length=fy)
w_denom<-vector(length=fy)
obs_bar<-vector(length=fy)
pred_bar<-vector(length=fy)

for(y in 1:fy){
  pred_y<-as.numeric(strsplit(pred_Tss_m[y]," ")[[1]])[3:(length(lengths)+2)]
  pred_bar[y]<-sum(lengths*pred_y)
  obs_y<-as.numeric(strsplit(obs_Tss_m[y]," ")[[1]])[3:(length(lengths)+2)]
  obs_bar[y]<-sum(lengths*obs_y)
  v_y[y]<-sum(lengths^2*pred_y)-pred_bar[y]^2
  w_denom[y]<-(obs_bar[y]-pred_bar[y])/sqrt(v_y[y]/ESS_TSS_init[y])
}

wt_TSS_m<-1/var(w_denom)
ESS_iter_TSS_m[i+1,1]<-wt_TSS_m
SSQ_iter_TSS_m[i]<-sum((ESS_iter_TSS_m[i+1,1]-ESS_iter_TSS_m[i,1])^2)

# Write new CTL file with iterated weights
CTL<-c(
CTL_init[1:(wt_LLFA_line_f-1)],
paste0(wt_LLFA_f," ","#wt fish1 age comp female iter",sep=""),
paste0(wt_LLFA_m," ","#wt fish1 age comp male iter",sep=""),
paste0(wt_LLSA_f," ","#wt surv1 age comp female iter",sep=""),
paste0(wt_LLSA_m," ","#wt surv1 age comp male iter",sep=""),
paste0(wt_LLJSA_f," ","#wt surv2 age comp female iter",sep=""),
paste0(wt_LLJSA_m," ","#wt surv2 age comp male iter",sep=""),
paste0(wt_LLFS_m," ","#wt fish1 size comp male iter",sep=""),
paste0(wt_LLFS_f," ","#wt fish1 size comp female iter",sep=""),
paste0(wt_LLSS_m," ","#wt surv1 size comp male iter",sep=""),
paste0(wt_LLSS_f," ","#wt surv1 size comp female iter",sep=""),
paste0(wt_LLJSS_m," ","#wt surv2 size comp male iter",sep=""),
paste0(wt_LLJSS_f," ","#wt surv2 size comp female iter",sep=""),
paste0(wt_TFS_m," ","#wt fish 3 size comp male iter",sep=""),
paste0(wt_TFS_f," ","#wt fish 3 size comp female iter",sep=""),
paste0(wt_TSS_m," ","#wt srv7 size comp male iter",sep=""),
paste0(wt_TSS_f," ","#wt srv7 size comp female iter",sep=""),
CTL_init[(wt_TSS_line_f+1):length(CTL_init)])

do_rwt<-grep("data_reweight_switch",CTL_init)              # make sure reweight strip is set to 1
CTL[do_rwt+1]<-1

write.table(CTL,file=paste(pathM,"/tem.ctl",sep=""),quote=F,row.names=F,col.names=F)

# file.copy(from=paste(pathM,"/tem.par",sep=""),to=paste(pathM,"/tem.pin",sep=""),overwrite=T)   # use .par as .pin to improve estimation (usually)

# print selectivity
sab_curr <- dget(paste0(path,"//tem.rdat",sep='')) 
selex_names = c("Derby fishery female","Derby fishery male","Trawl fishery female" ,"Trawl fishery male" ,"IFQ fishery female" ,"IFQ fishery male", "IFQ Recent fishery female" ,
                "IFQ Recent fishery male", "Domestic LL survey female","Domestic LL survey male", "Domestic LL Recent survey female","Domestic LL Recent survey male","Cooperative LL survey female" ,
                "Cooperative LL survey male" ,"GOA trawl survey female","GOA trawl survey male" ) # selectivity names

# pivot longer for plotting
selex = data.frame(sab_curr$agesel, ages = ages) 
names(selex) = c(selex_names, "ages") # rename columns for plotting

# pivot longer for plotting purposes
selex_df = selex %>%  
  pivot_longer(!ages, names_to = "type", values_to = "selex") %>% 
  mutate(type = factor(type, levels = selex_names),
         sex = case_when( # differentiate sexes
           str_detect(type, "female") ~ "Female",
           str_detect(type, "male") ~ "Male" ))

# plot selectivities!
pdf(paste0(pathR,"//Selectivity.pdf"), width = 10)
print(ggplot(selex_df, aes(x = ages, y = selex, color = sex)) +
        geom_line() +
        geom_point() +
        facet_wrap(~type) +
        labs(x = "Ages", y = "Selectivity") +
        theme(legend.position = "none"))
dev.off()


# end iteration loop
}

end_time<-Sys.time()

end_time-start_time


##################################################################################################################################################################
# Write results
################################################################################################################################################################## 

write.csv(ESS_iter_LLFA_f,paste(pathR,"/ESS_iter_LLFA_f.csv",sep=""))
write.csv(ESS_iter_LLSA_f,paste(pathR,"/ESS_iter_LLSA_f.csv",sep=""))
write.csv(ESS_iter_LLJSA_f,paste(pathR,"/ESS_iter_LLJSA_f.csv",sep=""))
write.csv(ESS_iter_LLFA_m,paste(pathR,"/ESS_iter_LLFA_m.csv",sep=""))
write.csv(ESS_iter_LLSA_m,paste(pathR,"/ESS_iter_LLSA_m.csv",sep=""))
write.csv(ESS_iter_LLJSA_m,paste(pathR,"/ESS_iter_LLJSA_m.csv",sep=""))
write.csv(ESS_iter_LLFS_m,paste(pathR,"/ESS_iter_LLFS_m.csv",sep=""))
write.csv(ESS_iter_LLFS_f,paste(pathR,"/ESS_iter_LLFS_f.csv",sep=""))
write.csv(ESS_iter_TFS_m,paste(pathR,"/ESS_iter_TFS_m.csv",sep=""))
write.csv(ESS_iter_TFS_f,paste(pathR,"/ESS_iter_TFS_f.csv",sep=""))
write.csv(ESS_iter_LLSS_m,paste(pathR,"/ESS_iter_LLSS_m.csv",sep=""))
write.csv(ESS_iter_LLSS_f,paste(pathR,"/ESS_iter_LLSS_f.csv",sep=""))
write.csv(ESS_iter_LLJSS_m,paste(pathR,"/ESS_iter_LLJSS_m.csv",sep=""))
write.csv(ESS_iter_LLJSS_f,paste(pathR,"/ESS_iter_LLJSS_f.csv",sep=""))
write.csv(ESS_iter_TSS_m,paste(pathR,"/ESS_iter_TSS_m.csv",sep=""))
write.csv(ESS_iter_TSS_f,paste(pathR,"/ESS_iter_TSS_f.csv",sep=""))

ESS_tot<-cbind(ESS_iter_LLFA_f,ESS_iter_LLSA_f,ESS_iter_LLJSA_f,
               ESS_iter_LLFA_m,ESS_iter_LLSA_m,ESS_iter_LLJSA_m, ESS_iter_LLFS_m,ESS_iter_LLFS_f,ESS_iter_TFS_m,ESS_iter_TFS_f,ESS_iter_LLSS_m,
               ESS_iter_LLSS_f,ESS_iter_LLJSS_m,ESS_iter_LLJSS_f,ESS_iter_TSS_m,ESS_iter_TSS_f)
write.csv(ESS_tot,paste(pathR,"/ESS_TOT.csv",sep=""))


write.csv(SSQ_iter_LLFA_f,paste(pathR,"/SSQ_iter_LLFA_f.csv",sep=""))
write.csv(SSQ_iter_LLSA_f,paste(pathR,"/SSQ_iter_LLSA_f.csv",sep=""))
write.csv(SSQ_iter_LLJSA_f,paste(pathR,"/SSQ_iter_LLJSA_f.csv",sep=""))
write.csv(SSQ_iter_LLFA_m,paste(pathR,"/SSQ_iter_LLFA_m.csv",sep=""))
write.csv(SSQ_iter_LLSA_m,paste(pathR,"/SSQ_iter_LLSA_m.csv",sep=""))
write.csv(SSQ_iter_LLJSA_m,paste(pathR,"/SSQ_iter_LLJSA_m.csv",sep=""))
write.csv(SSQ_iter_LLFS_m,paste(pathR,"/SSQ_iter_LLFS_m.csv",sep=""))
write.csv(SSQ_iter_LLFS_f,paste(pathR,"/SSQ_iter_LLFS_f.csv",sep=""))
write.csv(SSQ_iter_TFS_m,paste(pathR,"/SSQ_iter_TFS_m.csv",sep=""))
write.csv(SSQ_iter_TFS_f,paste(pathR,"/SSQ_iter_TFS_f.csv",sep=""))
write.csv(SSQ_iter_LLSS_m,paste(pathR,"/SSQ_iter_LLSS_m.csv",sep=""))
write.csv(SSQ_iter_LLSS_f,paste(pathR,"/SSQ_iter_LLSS_f.csv",sep=""))
write.csv(SSQ_iter_LLJSS_m,paste(pathR,"/SSQ_iter_LLJSS_m.csv",sep=""))
write.csv(SSQ_iter_LLJSS_f,paste(pathR,"/SSQ_iter_LLJSS_f.csv",sep=""))
write.csv(SSQ_iter_TSS_m,paste(pathR,"/SSQ_iter_TSS_m.csv",sep=""))
write.csv(SSQ_iter_TSS_f,paste(pathR,"/SSQ_iter_TSS_f.csv",sep=""))

SSQ_tot<-cbind(SSQ_iter_LLFA_f,SSQ_iter_LLSA_f,SSQ_iter_LLJSA_f,
               SSQ_iter_LLFA_m,SSQ_iter_LLSA_m,SSQ_iter_LLJSA_m, 
               SSQ_iter_LLFS_m,SSQ_iter_LLFS_f,SSQ_iter_TFS_m,SSQ_iter_TFS_f,SSQ_iter_LLSS_m,
               SSQ_iter_LLSS_f,SSQ_iter_LLJSS_m,SSQ_iter_LLJSS_f,SSQ_iter_TSS_m,SSQ_iter_TSS_f)
write.csv(SSQ_tot,paste(pathR,"/SSQ_TOT.csv",sep=""))
