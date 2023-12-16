# Purpose: To Run Francis Reweighting as proportions across
# Creator: Matthew LH. Cheng
# Date 11/20/23


# Set up ------------------------------------------------------------------

library(here)
library(tidyverse)

YEAR <- 2023
ITERS<-10
ages<-rep(seq(2,31),2) # repeat because 2 sex (need to be cautious when doing re-reweighting on age-aggregated comps...)
nages<-length(ages)
lengths<-rep(seq(41,99,2),2) # repeat because 2 sex
nleng<-length(lengths)

path<-here()
pathM<-paste(path,"/Growth_Time_Partial_Length",sep="")
pathD<-paste(pathM,"/Input Files",sep="")
pathR<-paste(pathM,"/Results",sep="")

# Read in data file
DAT_init<-readLines(paste0(pathD,"/tem_",YEAR,"_na_wh.dat",sep=""),warn=FALSE)  # read in .dat file 
file.copy(from=paste0(pathD,"/tem_",YEAR,"_na_wh.dat",sep=""),to=paste0(pathM,"/tem_",YEAR,"_na_wh.dat",sep=""),overwrite=T)
file.copy(from=paste0(pathD,"/tem.pin",sep=""),to=paste0(pathM,"/tem.pin",sep=""),overwrite=T)



# Get dataset years -------------------------------------------------------
yrs_LLFA<-as.numeric(unlist(strsplit(DAT_init[grep("# Fixed Gear Fishery Age Composition",DAT_init)+3],split=" "))) 
yrs_LLSA<-as.numeric(unlist(strsplit(DAT_init[grep("# Domestic LL survey Age Composition",DAT_init)+3],split=" ")))
yrs_LLJSA<-as.numeric(unlist(strsplit(DAT_init[grep("# Japanese LL survey Age Composition",DAT_init)+3],split=" ")))
yrs_LLFS<-as.numeric(unlist(strsplit(DAT_init[grep("# U.S. Fixed gear fishery length compositions",DAT_init)[1]+3],split=" ")))
yrs_JFS<-as.numeric(unlist(strsplit(DAT_init[grep("# Japanese longline fishery unsexed lengths",DAT_init)+3],split=" ")))
yrs_TFS<-as.numeric(unlist(strsplit(DAT_init[grep("# U.S. Trawl gear fishery length compositions",DAT_init)+3],split=" ")))
yrs_JTFS<-as.numeric(unlist(strsplit(DAT_init[grep("# Japanese Trawl survey unsexed size Composition",DAT_init)[1]+3],split=" ")))
yrs_LLSS<-as.numeric(unlist(strsplit(DAT_init[grep("# U.S. Domestic LL Survey Size Composition",DAT_init)+3],split=" ")))
yrs_LLJSS<-as.numeric(unlist(strsplit(DAT_init[grep("# Japanese LL survey size Composition",DAT_init)+3],split=" ")))
yrs_TSS<-as.numeric(unlist(strsplit(DAT_init[grep("# GOA trawl survey size composition",DAT_init)+3],split=" ")))
yrs_JPLLFS<-as.numeric(unlist(strsplit(DAT_init[grep("# Japanese longline fishery unsexed lengths",DAT_init)+3],split=" ")))
yrs_TSA<-as.numeric(unlist(strsplit(DAT_init[grep("# GOA trawl survey historical age composition",DAT_init)+3],split=" ")))
yrs_TSA<-as.numeric(unlist(strsplit(DAT_init[grep("# GOA trawl survey historical age composition",DAT_init)+3],split=" ")))
yrs_wo_LLFS_age = 9
yrs_wo_LLSS_age = 9

# Get Sample Sizes --------------------------------------------------------
# This code is gross....
# ESS_init<-as.numeric(unlist(strsplit(DAT_init[grep("# Number of samples",DAT_init)+1],split=" ")))   # All comp sample sizes
ESS_LLFA_init<-rep(20, yrs_LLFA)
ESS_LLSA_init<-rep(20, yrs_LLSA)
ESS_LLJSA_init<-rep(20, yrs_LLJSA)
ESS_LLFS_init<- c(rep(20, yrs_wo_LLFS_age), rep(0, yrs_LLFS - yrs_wo_LLFS_age))
ESS_JFS_init<-rep(20, yrs_JFS)
ESS_TFS_init<-rep(20, yrs_TFS)
ESS_JTFS_init<-rep(20, yrs_JTFS)
ESS_LLSS_init<-c(rep(20, yrs_wo_LLSS_age), rep(0, yrs_LLSS - yrs_wo_LLSS_age - 1), 20) # last year has survey size comps
ESS_LLJSS_init<- rep(20, yrs_LLJSS)
ESS_TSS_init<-rep(20, yrs_TSS)
ESS_JPLLFS_init = rep(20, yrs_JPLLFS)
ESS_TSA_init = rep(20, yrs_TSA)


# Get weights (CTL file) --------------------------------------------------
CTL_init<-readLines(paste(pathD,"/tem.ctl",sep=""),warn=FALSE)  # read in .ctl file 
wt_LLFA_line<-grep("#wt fish1 age comp iter",CTL_init)
wt_LLSA_line<-grep("#wt surv1 age comp iter",CTL_init)
wt_LLJSA_line<-grep("#wt surv2 age comp iter",CTL_init)
wt_LLFS_line<-grep("#wt fish1 size comp iter",CTL_init)
wt_LLSS_line<-grep("#wt surv1 size comp iter",CTL_init)
wt_LLJSS_line<-grep("#wt surv2 size comp iter",CTL_init)
wt_TFS_line<-grep("#wt fish 3 size comp iter",CTL_init)
wt_JPLLFS_line<-grep("#wt fish 2 size comp iter",CTL_init)
wt_TSS_line<-grep("#wt srv 7 size comp iter",CTL_init)

# Set up result matrices --------------------------------------------------
ESS_iter_LLFA<-matrix(nrow=ITERS+1,ncol=1)
rownames(ESS_iter_LLFA)<-c("Init",paste("Iter_",seq(1,ITERS),sep=""))
colnames(ESS_iter_LLFA)<-c("wt_LLFA")
ESS_iter_LLFA[1,1]<-1
SSQ_iter_LLFA<-matrix(nrow=ITERS,ncol=1)
rownames(SSQ_iter_LLFA)<-paste("Iter_",seq(1,ITERS),sep="")
colnames(SSQ_iter_LLFA)<-"SSQ_ESS_diff"

ESS_iter_LLSA<-matrix(nrow=ITERS+1,ncol=1)
rownames(ESS_iter_LLSA)<-c("Init",paste("Iter_",seq(1,ITERS),sep=""))
colnames(ESS_iter_LLSA)<-c("wt_LLSA")
ESS_iter_LLSA[1,1]<-1
SSQ_iter_LLSA<-matrix(nrow=ITERS,ncol=1)
rownames(SSQ_iter_LLSA)<-paste("Iter_",seq(1,ITERS),sep="")
colnames(SSQ_iter_LLSA)<-"SSQ_ESS_diff"

ESS_iter_LLJSA<-matrix(nrow=ITERS+1,ncol=1)
rownames(ESS_iter_LLJSA)<-c("Init",paste("Iter_",seq(1,ITERS),sep=""))
colnames(ESS_iter_LLJSA)<-c("wt_LLJSA")
ESS_iter_LLJSA[1,1]<-1
SSQ_iter_LLJSA<-matrix(nrow=ITERS,ncol=1)
rownames(SSQ_iter_LLJSA)<-paste("Iter_",seq(1,ITERS),sep="")
colnames(SSQ_iter_LLJSA)<-"SSQ_ESS_diff"

ESS_iter_LLFS<-matrix(nrow=ITERS+1,ncol=1)
rownames(ESS_iter_LLFS)<-c("Init",paste("Iter_",seq(1,ITERS),sep=""))
colnames(ESS_iter_LLFS)<-c("wt_LLFS")
ESS_iter_LLFS[1,1]<-1
SSQ_iter_LLFS<-matrix(nrow=ITERS,ncol=1)
rownames(SSQ_iter_LLFS)<-paste("Iter_",seq(1,ITERS),sep="")
colnames(SSQ_iter_LLFS)<-"SSQ_ESS_diff"

ESS_iter_LLSS<-matrix(nrow=ITERS+1,ncol=1)
rownames(ESS_iter_LLSS)<-c("Init",paste("Iter_",seq(1,ITERS),sep=""))
colnames(ESS_iter_LLSS)<-c("wt_LLSS")
ESS_iter_LLSS[1,1]<-1
SSQ_iter_LLSS<-matrix(nrow=ITERS,ncol=1)
rownames(SSQ_iter_LLSS)<-paste("Iter_",seq(1,ITERS),sep="")
colnames(SSQ_iter_LLSS)<-"SSQ_ESS_diff"

ESS_iter_LLJSS<-matrix(nrow=ITERS+1,ncol=1)
rownames(ESS_iter_LLJSS)<-c("Init",paste("Iter_",seq(1,ITERS),sep=""))
colnames(ESS_iter_LLJSS)<-c("wt_LLJSS")
ESS_iter_LLJSS[1,1]<-1
SSQ_iter_LLJSS<-matrix(nrow=ITERS,ncol=1)
rownames(SSQ_iter_LLJSS)<-paste("Iter_",seq(1,ITERS),sep="")
colnames(SSQ_iter_LLJSS)<-"SSQ_ESS_diff"

ESS_iter_TFS<-matrix(nrow=ITERS+1,ncol=1)
rownames(ESS_iter_TFS)<-c("Init",paste("Iter_",seq(1,ITERS),sep=""))
colnames(ESS_iter_TFS)<-c("wt_TFS")
ESS_iter_TFS[1,1]<-1
SSQ_iter_TFS<-matrix(nrow=ITERS,ncol=1)
rownames(SSQ_iter_TFS)<-paste("Iter_",seq(1,ITERS),sep="")
colnames(SSQ_iter_TFS)<-"SSQ_ESS_diff"

ESS_iter_TSS<-matrix(nrow=ITERS+1,ncol=1)
rownames(ESS_iter_TSS)<-c("Init",paste("Iter_",seq(1,ITERS),sep=""))
colnames(ESS_iter_TSS)<-c("wt_TSS")
ESS_iter_TSS[1,1]<-1
SSQ_iter_TSS<-matrix(nrow=ITERS,ncol=1)
rownames(SSQ_iter_TSS)<-paste("Iter_",seq(1,ITERS),sep="")
colnames(SSQ_iter_TSS)<-"SSQ_ESS_diff"

ESS_iter_JPLLFS<-matrix(nrow=ITERS+1,ncol=1)
rownames(ESS_iter_JPLLFS)<-c("Init",paste("Iter_",seq(1,ITERS),sep=""))
colnames(ESS_iter_JPLLFS)<-c("wt_JPLLFS")
ESS_iter_JPLLFS[1,1]<-1
SSQ_iter_JPLLFS<-matrix(nrow=ITERS,ncol=1)
rownames(SSQ_iter_JPLLFS)<-paste("Iter_",seq(1,ITERS),sep="")
colnames(SSQ_iter_JPLLFS)<-"SSQ_ESS_diff"

ESS_iter_TSA<-matrix(nrow=ITERS+1,ncol=1)
rownames(ESS_iter_TSA)<-c("Init",paste("Iter_",seq(1,ITERS),sep=""))
colnames(ESS_iter_TSA)<-c("wt_TSA")
ESS_iter_TSA[1,1]<-1
SSQ_iter_TSA<-matrix(nrow=ITERS,ncol=1)
rownames(SSQ_iter_TSA)<-paste("Iter_",seq(1,ITERS),sep="")
colnames(SSQ_iter_TSA)<-"SSQ_ESS_diff"

# Write new CTL file with iterated weights
CTL<-c(
  CTL_init[1:(wt_LLFA_line-1)],
  paste0(1," ","#wt fish1 age comp iter",sep=""),
  paste0(1," ","#wt surv1 age comp iter",sep=""),
  paste0(1," ","#wt surv2 age comp iter",sep=""),
  paste0(1," ","#wt fish1 size comp iter",sep=""),
  paste0(1," ","#wt surv1 size comp iter",sep=""),
  paste0(1," ","#wt surv2 size comp iter",sep=""),
  paste0(1," ","#wt fish 3 size comp iter",sep=""),
  paste0(0," ","#wt fish 2 size comp iter",sep=""), # 0 weight, not fitting
  paste0(0," ","#wt srv 7 age comp iter",sep=""), # 0 weight, not fitting
  paste0(1," ","#wt srv 7 size comp iter",sep=""),
  CTL_init[(wt_TSS_line+1):length(CTL_init)])

do_rwt<-grep("data_reweight_switch",CTL_init)              # make sure reweight strip is set to 1
CTL[do_rwt+1]<-1
write.table(CTL,file=paste(pathM,"/tem.ctl",sep=""),quote=F,row.names=F,col.names=F)


start_time<-Sys.time()

# Start iteration loop
for (i in 1:ITERS){
  
  print(paste("i = ",i,sep=""))
  
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
  
  
  
  # Fixed Gear Age Comps ----------------------------------------------------
  
  line1<-grep("Obs_P_fish_age",REP)+1
  line2<-grep("Pred_P_fish1_age",REP)-2
  fy<-length(REP[line1:line2])
  obs_fa<-REP[line1:line2]
  obs_p_fa<-matrix(nrow=fy,ncol=nages+1)
  line1<-grep("Pred_P_fish1_age",REP)+1
  line2<-grep("Obs_P_fish1_size",REP)-2
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
  
  wt_LLFA<-1/var(w_denom)
  ESS_iter_LLFA[i+1,1]<-wt_LLFA
  SSQ_iter_LLFA[i]<-sum((ESS_iter_LLFA[i+1,1]-ESS_iter_LLFA[i,1])^2)
  
  # Fixed Gear Size Comps ---------------------------------------------------
  
  line1<-grep("Obs_P_fish1_size ",REP)+1
  line2<-grep("Pred_P_fish1_size ",REP)-2
  fy<-length(REP[line1:line2])
  obs_LLfs<-REP[line1:line2]
  obs_p_LLfs<-matrix(nrow=fy,ncol=nleng+1)
  line1<-grep("Pred_P_fish1_size",REP)+1
  line2<-grep("Obs_P_fish3_size",REP)-2
  pred_LLfs<-REP[line1:line2]
  pred_p_LLfs<-matrix(nrow=fy,ncol=nleng+1)
  v_y<-vector(length=fy)
  w_denom<-vector(length=fy)
  obs_bar<-vector(length=fy)
  pred_bar<-vector(length=fy)
  
  for(y in 1:fy){
    pred_y<-as.numeric(strsplit(pred_LLfs[y]," ")[[1]])[3:(length(lengths)+2)]
    pred_bar[y]<-sum(lengths*pred_y)
    obs_y<-as.numeric(strsplit(obs_LLfs[y]," ")[[1]])[3:(length(lengths)+2)]
    obs_bar[y]<-sum(lengths*obs_y)
    v_y[y]<-sum(lengths^2*pred_y)-pred_bar[y]^2
    w_denom[y]<-(obs_bar[y]-pred_bar[y])/sqrt(v_y[y]/ESS_LLFS_init[y])
  }
  
  wt_LLFS<-1/var(w_denom[which(ESS_LLFS_init != 0)])
  ESS_iter_LLFS[i+1,1]<-wt_LLFS
  SSQ_iter_LLFS[i]<-sum((ESS_iter_LLFS[i+1,1]-ESS_iter_LLFS[i,1])^2)
  
  # Trawl Fleet Size Comps --------------------------------------------------
  
  line1<-grep("Obs_P_fish3_size",REP)+1
  line2<-grep("Pred_P_fish3_size",REP)-2
  fy<-length(REP[line1:line2])
  obs_Tfs<-REP[line1:line2]
  obs_p_Tfs<-matrix(nrow=fy,ncol=nleng+1)
  line1<-grep("Pred_P_fish3_size",REP)+1
  line2<-grep("Obs_P_srv1_age",REP)-2
  pred_Tfs<-REP[line1:line2]
  pred_p_Tfs<-matrix(nrow=fy,ncol=nleng+1)
  v_y<-vector(length=fy)
  w_denom<-vector(length=fy)
  obs_bar<-vector(length=fy)
  pred_bar<-vector(length=fy)
  
  for(y in 1:fy){
    pred_y<-as.numeric(strsplit(pred_Tfs[y]," ")[[1]])[3:(length(lengths)+2)]
    pred_bar[y]<-sum(lengths*pred_y)
    obs_y<-as.numeric(strsplit(obs_Tfs[y]," ")[[1]])[3:(length(lengths)+2)]
    obs_bar[y]<-sum(lengths*obs_y)
    v_y[y]<-sum(lengths^2*pred_y)-pred_bar[y]^2
    w_denom[y]<-(obs_bar[y]-pred_bar[y])/sqrt(v_y[y]/ESS_TFS_init[y])
  }
  
  wt_TFS<-1/var(w_denom)
  ESS_iter_TFS[i+1,1]<-wt_TFS
  SSQ_iter_TFS[i]<-sum((ESS_iter_TFS[i+1,1]-ESS_iter_TFS[i,1])^2)
  
  
  # LL Survey Ages ----------------------------------------------------------
  
  line1<-grep("Obs_P_srv1_age",REP)+1
  line2<-grep("Pred_P_srv1_age",REP)-2
  sy<-length(REP[line1:line2])
  obs_LLsa<-REP[line1:line2]
  obs_p_LLsa<-matrix(nrow=sy,ncol=nages+1)
  line1<-grep("Pred_P_srv1_age",REP)+1
  line2<-grep("Obs_P_srv2_age",REP)-2
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
  
  wt_LLSA<-1/var(w_denom)
  ESS_iter_LLSA[i+1,1]<-wt_LLSA
  SSQ_iter_LLSA[i]<-sum((ESS_iter_LLSA[i+1,1]-ESS_iter_LLSA[i,1])^2)
  
  
  # JP LL Survey Age -------------------------------------------------------
  
  line1<-grep("Obs_P_srv2_age",REP)+1
  line2<-grep("Pred_P_srv2_age",REP)-2
  sy<-length(REP[line1:line2])
  obs_LLJsa<-REP[line1:line2]
  obs_p_LLJsa<-matrix(nrow=sy,ncol=(length(ages))+1)
  line1<-grep("Pred_P_srv2_age",REP)+1
  line2<-grep("Obs_P_srv1_size",REP)-2
  pred_LLJsa<-REP[line1:line2]
  pred_p_LLJsa<-matrix(nrow=sy,ncol=(length(ages))+1)
  v_y<-vector(length=sy)
  w_denom<-vector(length=sy)
  obs_bar<-vector(length=sy)
  pred_bar<-vector(length=sy)             
  
  for(y in 1:sy){ # not sex-specific
    pred_y<-as.numeric(strsplit(pred_LLJsa[y]," ")[[1]])[3:((length(ages))+2)]
    pred_bar[y]<-sum(ages*pred_y)
    obs_y<-as.numeric(strsplit(obs_LLJsa[y]," ")[[1]])[3:((length(ages))+2)]
    obs_bar[y]<-sum(ages*obs_y)
    v_y[y]<-sum(ages^2*pred_y)-pred_bar[y]^2
    w_denom[y]<-(obs_bar[y]-pred_bar[y])/sqrt(v_y[y]/ESS_LLJSA_init[y])
  }
  
  wt_LLJSA<-(1/var(w_denom))
  ESS_iter_LLJSA[i+1,1]<-wt_LLJSA
  SSQ_iter_LLJSA[i]<-sum((ESS_iter_LLJSA[i+1,1]-ESS_iter_LLJSA[i,1])^2)
  
  
  # Survey Size Comps -------------------------------------------------------
  
  line1<-grep("Obs_P_srv1_size",REP)+1
  line2<-grep("Pred_P_srv1_size",REP)-2
  fy<-length(REP[line1:line2])
  obs_LLss<-REP[line1:line2]
  obs_p_LLss<-matrix(nrow=fy,ncol=nleng+1)
  line1<-grep("Pred_P_srv1_size",REP)+1
  line2<-grep("Obs_P_srv2_size",REP)-2
  pred_LLss<-REP[line1:line2]
  pred_p_LLss<-matrix(nrow=fy,ncol=nleng+1)
  v_y<-vector(length=fy)
  w_denom<-vector(length=fy)
  obs_bar<-vector(length=fy)
  pred_bar<-vector(length=fy)
  
  for(y in 1:fy){
    pred_y<-as.numeric(strsplit(pred_LLss[y]," ")[[1]])[3:(length(lengths)+2)]
    pred_bar[y]<-sum(lengths*pred_y)
    obs_y<-as.numeric(strsplit(obs_LLss[y]," ")[[1]])[3:(length(lengths)+2)]
    obs_bar[y]<-sum(lengths*obs_y)
    v_y[y]<-sum(lengths^2*pred_y)-pred_bar[y]^2
    w_denom[y]<-(obs_bar[y]-pred_bar[y])/sqrt(v_y[y]/ESS_LLSS_init[y])
  }
  
  wt_LLSS<-1/var(w_denom[which(ESS_LLSS_init != 0)])
  ESS_iter_LLSS[i+1,1]<-wt_LLSS
  SSQ_iter_LLSS[i]<-sum((ESS_iter_LLSS[i+1,1]-ESS_iter_LLSS[i,1])^2)
  
  
  # Japanese LL Survey Size Comps -------------------------------------------
  
  line1<-grep("Obs_P_srv2_size",REP)+1
  line2<-grep("Pred_P_srv2_size",REP)-2
  fy<-length(REP[line1:line2])
  obs_LLJss<-REP[line1:line2]
  obs_p_LLJss<-matrix(nrow=fy,ncol=nleng+1)
  line1<-grep("Pred_P_srv2_size",REP)+1
  line2<-grep("Obs_P_srv7_size",REP)-2
  pred_LLJss<-REP[line1:line2]
  pred_p_LLJss<-matrix(nrow=fy,ncol=nleng+1)
  v_y<-vector(length=fy)
  w_denom<-vector(length=fy)
  obs_bar<-vector(length=fy)
  pred_bar<-vector(length=fy)
  
  for(y in 1:fy){
    pred_y<-as.numeric(strsplit(pred_LLJss[y]," ")[[1]])[3:(length(lengths)+2)]
    pred_bar[y]<-sum(lengths*pred_y)
    obs_y<-as.numeric(strsplit(obs_LLJss[y]," ")[[1]])[3:(length(lengths)+2)]
    obs_bar[y]<-sum(lengths*obs_y)
    v_y[y]<-sum(lengths^2*pred_y)-pred_bar[y]^2
    w_denom[y]<-(obs_bar[y]-pred_bar[y])/sqrt(v_y[y]/ESS_LLJSS_init[y])
  }
  
  wt_LLJSS<-1/var(w_denom)
  ESS_iter_LLJSS[i+1,1]<-wt_LLJSS
  SSQ_iter_LLJSS[i]<-sum((ESS_iter_LLJSS[i+1,1]-ESS_iter_LLJSS[i,1])^2)
  
  
  # Trawl Survey Size -------------------------------------------------------
  
  line1<-grep("Obs_P_srv7_size",REP)+1
  line2<-grep("Pred_P_srv7_size",REP)-2
  fy<-length(REP[line1:line2])
  obs_Tss<-REP[line1:line2]
  obs_p_Tss<-matrix(nrow=fy,ncol=nleng+1)
  line1<-grep("Pred_P_srv7_size",REP)+1
  line2<-grep("Obs_P_fish2_size",REP)-2
  pred_Tss<-REP[line1:line2]
  pred_p_Tss<-matrix(nrow=fy,ncol=nleng+1)
  v_y<-vector(length=fy)
  w_denom<-vector(length=fy)
  obs_bar<-vector(length=fy)
  pred_bar<-vector(length=fy)
  
  for(y in 1:fy){
    pred_y<-as.numeric(strsplit(pred_Tss[y]," ")[[1]])[3:(length(lengths)+2)]
    pred_bar[y]<-sum(lengths*pred_y)
    obs_y<-as.numeric(strsplit(obs_Tss[y]," ")[[1]])[3:(length(lengths)+2)]
    obs_bar[y]<-sum(lengths*obs_y)
    v_y[y]<-sum(lengths^2*pred_y)-pred_bar[y]^2
    w_denom[y]<-(obs_bar[y]-pred_bar[y])/sqrt(v_y[y]/ESS_TSS_init[y])
  }
  
  wt_TSS<-1/var(w_denom)
  ESS_iter_TSS[i+1,1]<-wt_TSS
  SSQ_iter_TSS[i]<-sum((ESS_iter_TSS[i+1,1]-ESS_iter_TSS[i,1])^2)
  
  
  # JP LL Fishery Size Comps ------------------------------------------------
  
  line1<-grep("Obs_P_fish2_size",REP)+1
  line2<-grep("Pred_P_fish2_size",REP)-2
  fy<-length(REP[line1:line2])
  obs_JPLLFS<-REP[line1:line2]
  obs_p_JPLLFS<-matrix(nrow=fy,ncol=(nleng/2)+1)
  line1<-grep("Pred_P_fish2_size",REP)+1
  line2<-grep("Obs_P_srv7_age",REP)-2
  pred_JPLLFS<-REP[line1:line2]
  pred_p_JPLLFS<-matrix(nrow=fy,ncol=(nleng/2)+1)
  v_y<-vector(length=fy)
  w_denom<-vector(length=fy)
  obs_bar<-vector(length=fy)
  pred_bar<-vector(length=fy)
  
  for(y in 1:fy){
    pred_y<-as.numeric(strsplit(pred_JPLLFS[y]," ")[[1]])[3:((length(lengths)/2)+2)]
    pred_bar[y]<-sum((lengths[1:(length(lengths)/2)])*pred_y)
    obs_y<-as.numeric(strsplit(obs_JPLLFS[y]," ")[[1]])[3:((length(lengths)/2)+2)]
    obs_bar[y]<-sum((lengths[1:(length(lengths)/2)])*obs_y)
    v_y[y]<-sum((lengths[1:(length(lengths)/2)])^2*pred_y)-pred_bar[y]^2
    w_denom[y]<-(obs_bar[y]-pred_bar[y])/sqrt(v_y[y]/ESS_JPLLFS_init[y])
  }
  
  wt_JPLLFS<-1/var(w_denom)
  ESS_iter_JPLLFS[i+1,1]<-wt_JPLLFS
  SSQ_iter_JPLLFS[i]<-sum((ESS_iter_JPLLFS[i+1,1]-ESS_iter_JPLLFS[i,1])^2)
  
  # Trawl Survey Ages -------------------------------------------------------
  
  # Not sex-specific
  line1<-grep("Obs_P_srv7_age",REP)+1
  line2<-grep("Pred_P_srv7_age",REP)-2
  sy<-length(REP[line1:line2])
  obs_TSA<-REP[line1:line2]
  obs_p_TSA<-matrix(nrow=sy,ncol=(length(ages)/2)+1)
  line1<-grep("Pred_P_srv7_age",REP)+1
  line2<-grep("Obs_full_P_fish_age",REP)-2
  pred_TSA<-REP[line1:line2]
  pred_p_TSA<-matrix(nrow=sy,ncol=(length(ages)/2)+1)
  v_y<-vector(length=sy)
  w_denom<-vector(length=sy)
  obs_bar<-vector(length=sy)
  pred_bar<-vector(length=sy)             
  
  for(y in 1:sy){ # not sex-specific
    pred_y<-as.numeric(strsplit(pred_TSA[y]," ")[[1]])[3:((length(ages)/2)+2)]
    pred_bar[y]<-sum(ages[1:(length(ages)/2)]*pred_y)
    obs_y<-as.numeric(strsplit(obs_TSA[y]," ")[[1]])[3:((length(ages)/2)+2)]
    obs_bar[y]<-sum(ages[1:(length(ages)/2)]*obs_y)
    v_y[y]<-sum(ages[1:(length(ages)/2)]^2*pred_y)-pred_bar[y]^2
    w_denom[y]<-(obs_bar[y]-pred_bar[y])/sqrt(v_y[y]/ESS_TSA_init[y])
  }
  
  wt_TSA<-(1/var(w_denom))
  ESS_iter_TSA[i+1,1]<-wt_TSA
  SSQ_iter_TSA[i]<-sum((ESS_iter_TSA[i+1,1]-ESS_iter_TSA[i,1])^2)  
  # Write new CTL file with iterated weights
  CTL<-c(
    CTL_init[1:(wt_LLFA_line-1)],
    paste0(wt_LLFA," ","#wt fish1 age comp iter",sep=""),
    paste0(wt_LLSA," ","#wt surv1 age comp iter",sep=""),
    paste0(wt_LLJSA," ","#wt surv2 age comp iter",sep=""),
    paste0(wt_LLFS," ","#wt fish1 size comp iter",sep=""),
    paste0(wt_LLSS," ","#wt surv1 size comp iter",sep=""),
    paste0(wt_LLJSS," ","#wt surv2 size comp iter",sep=""),
    paste0(wt_TFS," ","#wt fish 3 size comp iter",sep=""),
    paste0(0," ","#wt fish 2 size comp iter",sep=""), # 0 weight, not fitting
    paste0(0," ","#wt srv 7  age comp iter",sep=""), # 0 weight, not fitting
    paste0(wt_TSS," ","#wt srv 7 size comp iter",sep=""),
    CTL_init[(wt_TSS_line+1):length(CTL_init)])
  
  
  do_rwt<-grep("data_reweight_switch",CTL_init)              # make sure reweight strip is set to 1
  CTL[do_rwt+1]<-1
  
  write.table(CTL,file=paste(pathM,"/tem.ctl",sep=""),quote=F,row.names=F,col.names=F)
  
  # file.copy(from=paste(pathM,"/tem.par",sep=""),to=paste(pathM,"/tem.pin",sep=""),overwrite=T)   # use .par as .pin to improve estimation (usually)
  
  
  # print selectivity runs
  sab_curr <- dget(paste0(pathM,"//tem.rdat",sep='')) 
  sab_rep <- readLines(paste0(pathM,"//sable.rep",sep=''))
  
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


write.csv(ESS_iter_LLFA,paste(pathR,"/ESS_iter_LLFA.csv",sep=""))
write.csv(ESS_iter_LLSA,paste(pathR,"/ESS_iter_LLSA.csv",sep=""))
write.csv(ESS_iter_LLJSA,paste(pathR,"/ESS_iter_LLJSA.csv",sep=""))
write.csv(ESS_iter_LLFS,paste(pathR,"/ESS_iter_LLFS.csv",sep=""))
write.csv(ESS_iter_TFS,paste(pathR,"/ESS_iter_TFS.csv",sep=""))
write.csv(ESS_iter_LLSS,paste(pathR,"/ESS_iter_LLSS.csv",sep=""))
write.csv(ESS_iter_LLJSS,paste(pathR,"/ESS_iter_LLJSS.csv",sep=""))
write.csv(ESS_iter_TSS,paste(pathR,"/ESS_iter_TSS.csv",sep=""))
write.csv(ESS_iter_JPLLFS,paste(pathR,"/ESS_iter_JPLLFS.csv",sep=""))
write.csv(ESS_iter_TSA,paste(pathR,"/ESS_iter_TSA.csv",sep=""))

ESS_tot<-cbind(ESS_iter_LLFA,ESS_iter_LLSA,ESS_iter_LLJSA,ESS_iter_LLFS,ESS_iter_TFS,ESS_iter_LLSS,
               ESS_iter_JPLLFS, ESS_iter_LLJSS,ESS_iter_TSS, ESS_iter_TSA)
write.csv(ESS_tot,paste(pathR,"/ESS_TOT.csv",sep=""))

write.csv(SSQ_iter_LLFA,paste(pathR,"/SSQ_iter_LLFA.csv",sep=""))
write.csv(SSQ_iter_LLSA,paste(pathR,"/SSQ_iter_LLSA.csv",sep=""))
write.csv(SSQ_iter_LLJSA,paste(pathR,"/SSQ_iter_LLJSA.csv",sep=""))
write.csv(SSQ_iter_LLFS,paste(pathR,"/SSQ_iter_LLFS.csv",sep=""))
write.csv(SSQ_iter_TFS,paste(pathR,"/SSQ_iter_TFS.csv",sep=""))
write.csv(SSQ_iter_LLSS,paste(pathR,"/SSQ_iter_LLSS.csv",sep=""))
write.csv(SSQ_iter_LLJSS,paste(pathR,"/SSQ_iter_LLJSS.csv",sep=""))
write.csv(SSQ_iter_TSS,paste(pathR,"/SSQ_iter_TSS.csv",sep=""))
write.csv(SSQ_iter_JPLLFS,paste(pathR,"/SSQ_iter_JPLLFS.csv",sep=""))
write.csv(SSQ_iter_TSA,paste(pathR,"/SSQ_iter_TSA.csv",sep=""))

SSQ_tot<-cbind(SSQ_iter_LLFA,SSQ_iter_LLSA,SSQ_iter_LLJSA,SSQ_iter_LLFS,SSQ_iter_TFS,SSQ_iter_LLSS,
               SSQ_iter_JPLLFS, SSQ_iter_LLJSS,SSQ_iter_TSS,SSQ_iter_TSA)
write.csv(SSQ_tot,paste(pathR,"/SSQ_TOT.csv",sep=""))
