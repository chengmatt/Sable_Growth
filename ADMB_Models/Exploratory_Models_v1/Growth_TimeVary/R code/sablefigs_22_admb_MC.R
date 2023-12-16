###########################################################################################################################
# Purpose: Graphics for NOAA AK Sablefish Assessment
# Creator: Matthew LH. Cheng (UAF-CFOS) (building on code developed by D. Hanselman and D. Goethel (NOAA-AFSC-ABL)); Updated by D. Goethel
# Date 9/27/23
###########################################################################################################################

rm(list=(ls()))

# Set up ------------------------------------------------------------------
library(here)
library(tidyverse)
library(reshape2)
library(data.table)
library(R2admb)
library(ggsci)


############### INPUTS TO BE ALTERED ######################################

SA_curr_YR<-2023 #### enter terminal year for previous stock assessment model
SA_prev_YR<-2022 #### enter terminal year for previous stock assessment model
#####################################################################################################


dir_R<-here()
dir_master<-here(dir_R, "Growth_TimeVary")
dir_results<-paste0(dir_master,"//Results",sep='')
dir.create(dir_results)


# Read in .rdat
sab_curr <- dget(paste0(dir_master,"//tem.rdat",sep='')) 
sab_rep <- readLines(paste0(dir_master,"//sable.rep",sep=''))
ages = 2:31

load(paste0(dir_R,"\\Sablefish_Data_",SA_curr_YR,".RData",sep=''))


# Set up theme for ggplot
theme_reg = function() {
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 14),
        legend.text = element_text(color = "black", size = 12),
        legend.title = element_text(color = "black", size = 14),
        legend.background = element_blank(),
        strip.text = element_text(size = 12))
}


#############################################################################################################################################

# Plots

#############################################################################################################################################



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

ggsave(paste0(dir_results , "//Fig. 3.2. Catch by gear.png"),plot=catch_gear_plot) 



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






# Likelihood components ---------------------------------------------------
# Likelihood component names
like_names<-c("LL Fish Age","LL Srv Age", "LL Fish Size_F","LL Fish Size_M","TRWL Fish Size_F","TRWL Fish Size_M",
              "LL Srv Size_F","LL Srv Size_M","Coop Srv Size_F","Coop Srv Size_M", "TRWL Srv Size_F","TRWL Srv Size_M",
              "LL Srv RPN","Coop Srv RPN","LL CPUE RPN", "JPN CPUE RPN", "Trawl Survey RPW","Catch",
              "Recruit_Pen","F_Pen","M_Prior")     

likecomps = sab_curr$likecomp # extract out likelihood components
likecomps_df = data.frame(nLL = likecomps, component = names(likecomps)) # dataframe for plotting

# Do some residual munging (getting rid of data not fit to and objfun, and assign groups to components)
likecomps_df = likecomps_df %>% 
  filter(nLL != 0, component != "obj.fun") %>% 
  mutate(component = like_names,
         # Setting up grouping structure for different likelihood types
         groups = case_when(
           str_detect(component, "Catch") ~ "Catch",
           str_detect(component, "Age") ~ "Age Comps",
           str_detect(component, "Size") ~ "Length Comps",
           str_detect(component, "RPW") ~ "Survey Indices",
           str_detect(component, "RPN") ~ "Survey Indices",
           str_detect(component, "CPUE") ~ "CPUE",
           str_detect(component, "Pen") ~ "Penalties",
           str_detect(component, "Prior") ~ "Penalties"),
         groups = factor(groups, levels = c("Age Comps", "Length Comps", "Survey Indices", "CPUE",
                                            "Catch", "Penalties"))) # set up order of plotting here

# Plot likelihood components
pdf(paste0(dir_results,"//Like_Comps.pdf",sep=''))
ggplot(likecomps_df, aes(x = component, y = nLL, fill = groups)) +
  geom_col() +
  labs(x = "Likelihood Component", y = "Likelihood", fill = "") +
  scale_x_discrete(guide = guide_axis(angle = 90), limits = unique(likecomps_df$component))  +
  theme_reg() +
  theme(legend.position = c(0.85, 0.85))
dev.off()
  
# Fits to indices ---------------------------------------------------------
idx_names = c("Domestic LL Survey Relative Population Weight", "Japanese LL Survey Relative Population Weight",
              "Domestic LL Survey Relative Population Numbers", "Japanese LL Survey Relative Population Numbers",
              "Domestic Fishery CPUE Index", "Japanese Fishery CPUE Index", "GOA Trawl Survey Biomass (kt)") # index names

# extract out index data
idx_df = data.frame()
idx_list = sab_curr[str_detect(names(sab_curr), "obssrv")]

# Loop through to extract stuff out
for(i in 1:length(idx_list)) {
  # extract out index dataframe (just extracting out components from a list)
  idx_tmp = data.frame(year = as.numeric(rownames(idx_list[[i]])), obs = idx_list[[i]][[1]],
                       lci = idx_list[[i]][[2]], uci = idx_list[[i]][[3]],
                       pred = idx_list[[i]][[4]], type = idx_names[i])
  # bind together
  idx_df = rbind(idx_df, idx_tmp)
} # end i loop

# relvel by index name
idx_df$type = factor(idx_df$type, levels = idx_names)

# Now plot index data
pdf(paste0(dir_results,"//Index_Fits.pdf"), width = 13, height = 8)
ggplot(idx_df) +
  geom_line(mapping = aes(x = year, y = pred), size = 1.5, col = "red") +
  geom_pointrange(mapping = aes(x = year, y = obs, ymin = lci, ymax = uci), 
                  alpha = 0.75, size = 0.75, col = "blue") +
  facet_wrap(~type, scales = "free") +
  scale_y_continuous(limits = c(0, NA)) +  # NA for no upper limit
  labs(x = "Year", y = "Index") +
  theme_reg() +
  theme(legend.position = "top")
dev.off()

# Fits to compositions ---------------------------------------------------------
pdf(paste0(dir_results,"//Composition_Fits.pdf"), width = 10, height = 20)
### Domestic LL Fishery Age Compositions -------------------------------
obs_ac_fish1 = data.frame(reshape2::melt(sab_curr$oac.fish1), type = "obs")
pred_ac_fish1 = data.frame(reshape2::melt(sab_curr$eac.fish1), type = "pred")

# Put these into a dataframe
ac_fish1 = rbind(obs_ac_fish1, pred_ac_fish1)
names(ac_fish1) = c("year", "age", "prop", "type")
ac_fish1$broodYear = ac_fish1$year - ac_fish1$age # create brood year to track cohort over time
ac_fish1 = ac_fish1 %>% mutate(broodYear = ifelse(age == 31, "31", broodYear)) # make plus group consistent color

ggplot() +
  geom_col(ac_fish1 %>% filter(type == "obs"), mapping = aes(x = age, y = prop, fill = factor(broodYear))) +
  geom_line(ac_fish1 %>% filter(type == "pred"), mapping = aes(x = age, y = prop), size = 1) +
  facet_wrap(~year, scales = "free_x", dir = "v") +
  labs(x = "Ages", y = "Proportion", title = "Domestic LL Fishery Age Compositions") +
  theme_reg() +
  theme(legend.position = "none")

### Domestic LL Survey Age Compositions -------------------------------
obs_ac_srv1 = data.frame(reshape2::melt(sab_curr$oac.srv1), type = "obs")
pred_ac_srv1 = data.frame(reshape2::melt(sab_curr$eac.srv1), type = "pred")

# Put these into a dataframe
ac_srv1 = rbind(obs_ac_srv1, pred_ac_srv1)
names(ac_srv1) = c("year", "age", "prop", "type")
ac_srv1$broodYear = ac_srv1$year - ac_srv1$age # create brood year to track cohort over time
ac_srv1 = ac_srv1 %>% mutate(broodYear = ifelse(age == 31, "31", broodYear)) # make plus group consistent color

ggplot() +
  geom_col(ac_srv1 %>% filter(type == "obs"), mapping = aes(x = age, y = prop, fill = factor(broodYear))) +
  geom_line(ac_srv1 %>% filter(type == "pred"), mapping = aes(x = age, y = prop), size = 1) +
  facet_wrap(~year, scales = "free_x",  dir = "v") +
  labs(x = "Ages", y = "Proportion", title = "Domestic LL Survey Age Compositions") +
  theme_reg() +
  theme(legend.position = "none")

### Japanese LL Survey Age Compositions -------------------------------
obs_ac_srv2 = data.frame(reshape2::melt(sab_curr$oac.srv2), type = "obs")
pred_ac_srv2 = data.frame(reshape2::melt(sab_curr$eac.srv2), type = "pred")

# Put these into a dataframe
ac_srv2 = rbind(obs_ac_srv2, pred_ac_srv2)
names(ac_srv2) = c("year", "age", "prop", "type")
ac_srv2$broodYear = ac_srv2$year - ac_srv2$age # create brood year to track cohort over time
ac_srv2 = ac_srv2 %>% mutate(broodYear = ifelse(age == 31, "31", broodYear)) # make plus group consistent color

ggplot() +
  geom_col(ac_srv2 %>% filter(type == "obs"), mapping = aes(x = age, y = prop, fill = factor(broodYear))) +
  geom_line(ac_srv2 %>% filter(type == "pred"), mapping = aes(x = age, y = prop), size = 1) +
  facet_wrap(~year, scales = "free_x", ncol = 2, dir = "v") +
  labs(x = "Ages", y = "Proportion", title = "Japanese LL Survey Age Compositions") +
  theme_reg() +
  theme(legend.position = "none")

### Domestic LL Fishery Length Compositions (Female) -------------------------------
obs_lc_fish1_female = data.frame(reshape2::melt(sab_curr$olc.fish1.f), type = "obs")
pred_lc_fish1_female = data.frame(reshape2::melt(sab_curr$elc.fish1.f), type = "pred")

# Put these into a dataframe
lc_fish1_female = rbind(obs_lc_fish1_female, pred_lc_fish1_female)
names(lc_fish1_female) = c("year", "age", "prop", "type")

ggplot() +
  geom_col(lc_fish1_female %>% filter(type == "obs"), mapping = aes(x = age, y = prop), fill = "darkgreen", alpha = 0.85) +
  geom_line(lc_fish1_female %>% filter(type == "pred"), mapping = aes(x = age, y = prop), size = 1) +
  facet_wrap(~year, scales = "free_x",  dir = "v") +
  labs(x = "Ages", y = "Proportion", title = "Domestic LL Fishery Length Compositions (Female)") +
  theme_reg() +
  theme(legend.position = "none")

### Domestic LL Fishery Length Compositions (Male) -------------------------------
obs_lc_fish1_male = data.frame(reshape2::melt(sab_curr$olc.fish1.m), type = "obs")
pred_lc_fish1_male = data.frame(reshape2::melt(sab_curr$elc.fish1.m), type = "pred")

# Put these into a dataframe
lc_fish1_male = rbind(obs_lc_fish1_male, pred_lc_fish1_male)
names(lc_fish1_male) = c("year", "age", "prop", "type")

ggplot() +
  geom_col(lc_fish1_male %>% filter(type == "obs"), mapping = aes(x = age, y = prop), fill = "darkgreen", alpha = 0.85) +
  geom_line(lc_fish1_male %>% filter(type == "pred"), mapping = aes(x = age, y = prop), size = 1) +
  facet_wrap(~year, scales = "free_x", dir = "v") +
  labs(x = "Ages", y = "Proportion", title = "Domestic LL Fishery Length Compositions (Male)") +
  theme_reg() +
  theme(legend.position = "none")

### Domestic Trawl Fishery Length Compositions (Female) -------------------------------
obs_lc_fish3_female = data.frame(reshape2::melt(sab_curr$olc.fish3.f), type = "obs")
pred_lc_fish3_female = data.frame(reshape2::melt(sab_curr$elc.fish3.f), type = "pred")

# Put these into a dataframe
lc_fish3_female = rbind(obs_lc_fish3_female, pred_lc_fish3_female)
names(lc_fish3_female) = c("year", "age", "prop", "type")

ggplot() +
  geom_col(lc_fish3_female %>% filter(type == "obs"), mapping = aes(x = age, y = prop), fill = "darkgreen", alpha = 0.85) +
  geom_line(lc_fish3_female %>% filter(type == "pred"), mapping = aes(x = age, y = prop), size = 1) +
  facet_wrap(~year, scales = "free_x", dir = "v") +
  labs(x = "Ages", y = "Proportion", title = "Domestic Trawl Fishery Length Compositions (Female)") +
  theme_reg() +
  theme(legend.position = "none")

### Domestic Trawl Fishery Length Compositions (Male) -------------------------------
obs_lc_fish3_male = data.frame(reshape2::melt(sab_curr$olc.fish3.m), type = "obs")
pred_lc_fish3_male = data.frame(reshape2::melt(sab_curr$elc.fish3.m), type = "pred")

# Put these into a dataframe
lc_fish3_male = rbind(obs_lc_fish3_male, pred_lc_fish3_male)
names(lc_fish3_male) = c("year", "age", "prop", "type")

ggplot() +
  geom_col(lc_fish3_male %>% filter(type == "obs"), mapping = aes(x = age, y = prop), fill = "darkgreen", alpha = 0.85) +
  geom_line(lc_fish3_male %>% filter(type == "pred"), mapping = aes(x = age, y = prop), size = 1) +
  facet_wrap(~year, scales = "free_x", ncol = 2, dir = "v") +
  labs(x = "Ages", y = "Proportion", title = "Domestic Trawl Fishery Length Compositions (Male)") +
  theme_reg() +
  theme(legend.position = "none")

### Domestic LL Survey Length Compositions (Female) -------------------------------
obs_lc_srv1_female = data.frame(reshape2::melt(sab_curr$olc.srv1.f), type = "obs")
pred_lc_srv1_female = data.frame(reshape2::melt(sab_curr$elc.srv1.f), type = "pred")

# Put these into a dataframe
lc_srv1_female = rbind(obs_lc_srv1_female, pred_lc_srv1_female)
names(lc_srv1_female) = c("year", "age", "prop", "type")

ggplot() +
  geom_col(lc_srv1_female %>% filter(type == "obs"), mapping = aes(x = age, y = prop), fill = "darkgreen", alpha = 0.85) +
  geom_line(lc_srv1_female %>% filter(type == "pred"), mapping = aes(x = age, y = prop), size = 1) +
  facet_wrap(~year, scales = "free_x", ncol = 2, dir = "v") +
  labs(x = "Ages", y = "Proportion", title = "Domestic LL Survey Length Compositions (Female)") +
  theme_reg() +
  theme(legend.position = "none")

### Domestic LL Survey Length Compositions (Male) -------------------------------
obs_lc_srv1_male = data.frame(reshape2::melt(sab_curr$olc.srv1.m), type = "obs")
pred_lc_srv1_male = data.frame(reshape2::melt(sab_curr$elc.srv1.m), type = "pred")

# Put these into a dataframe
lc_srv1_male = rbind(obs_lc_srv1_male, pred_lc_srv1_male)
names(lc_srv1_male) = c("year", "age", "prop", "type")

ggplot() +
  geom_col(lc_srv1_male %>% filter(type == "obs"), mapping = aes(x = age, y = prop), fill = "darkgreen", alpha = 0.85) +
  geom_line(lc_srv1_male %>% filter(type == "pred"), mapping = aes(x = age, y = prop), size = 1) +
  facet_wrap(~year, scales = "free_x", ncol = 2, dir = "v") +
  labs(x = "Ages", y = "Proportion", title = "Domestic LL Survey Length Compositions (Male)") +
  theme_reg() +
  theme(legend.position = "none")

### Japanese LL Survey Length Compositions (Female) -------------------------------
obs_lc_srv2_female = data.frame(reshape2::melt(sab_curr$olc.srv2.f), type = "obs")
pred_lc_srv2_female = data.frame(reshape2::melt(sab_curr$elc.srv2.f), type = "pred")

# Put these into a dataframe
lc_srv2_female = rbind(obs_lc_srv2_female, pred_lc_srv2_female)
names(lc_srv2_female) = c("year", "age", "prop", "type")

ggplot() +
  geom_col(lc_srv2_female %>% filter(type == "obs"), mapping = aes(x = age, y = prop), fill = "darkgreen", alpha = 0.85) +
  geom_line(lc_srv2_female %>% filter(type == "pred"), mapping = aes(x = age, y = prop), size = 1) +
  facet_wrap(~year, scales = "free_x", ncol = 2, dir = "v") +
  labs(x = "Ages", y = "Proportion", title = "Japanese LL Survey Length Compositions (Female)") +
  theme_reg() +
  theme(legend.position = "none")

### Japanese LL Survey Length Compositions (Male) -------------------------------
obs_lc_srv2_male = data.frame(reshape2::melt(sab_curr$olc.srv2.m), type = "obs")
pred_lc_srv2_male = data.frame(reshape2::melt(sab_curr$elc.srv2.m), type = "pred")

# Put these into a dataframe
lc_srv2_male = rbind(obs_lc_srv2_male, pred_lc_srv2_male)
names(lc_srv2_male) = c("year", "age", "prop", "type")

ggplot() +
  geom_col(lc_srv2_male %>% filter(type == "obs"), mapping = aes(x = age, y = prop), fill = "darkgreen", alpha = 0.85) +
  geom_line(lc_srv2_male %>% filter(type == "pred"), mapping = aes(x = age, y = prop), size = 1) +
  facet_wrap(~year, scales = "free_x", ncol = 2, dir = "v") +
  labs(x = "Ages", y = "Proportion", title = "Japanese LL Survey Length Compositions (Male)") +
  theme_reg() +
  theme(legend.position = "none")

### Domestic Trawl Survey Length Compositions (Female) -------------------------------
obs_lc_srv7_female = data.frame(reshape2::melt(sab_curr$olc.srv7.f), type = "obs")
pred_lc_srv7_female = data.frame(reshape2::melt(sab_curr$elc.srv7.f), type = "pred")

# Put these into a dataframe
lc_srv7_female = rbind(obs_lc_srv7_female, pred_lc_srv7_female)
names(lc_srv7_female) = c("year", "age", "prop", "type")

ggplot() +
  geom_col(lc_srv7_female %>% filter(type == "obs"), mapping = aes(x = age, y = prop), fill = "darkgreen", alpha = 0.85) +
  geom_line(lc_srv7_female %>% filter(type == "pred"), mapping = aes(x = age, y = prop), size = 1) +
  facet_wrap(~year, scales = "free_x", ncol = 2, dir = "v") +
  labs(x = "Ages", y = "Proportion", title = "Domestic Trawl Survey Length Compositions (Female)") +
  theme_reg() +
  theme(legend.position = "none")

### Domestic Trawl Survey Length Compositions (Male) -------------------------------
obs_lc_srv7_male = data.frame(reshape2::melt(sab_curr$olc.srv7.m), type = "obs")
pred_lc_srv7_male = data.frame(reshape2::melt(sab_curr$elc.srv7.m), type = "pred")

# Put these into a dataframe
lc_srv7_male = rbind(obs_lc_srv7_male, pred_lc_srv7_male)
names(lc_srv7_male) = c("year", "age", "prop", "type")

ggplot() +
  geom_col(lc_srv7_male %>% filter(type == "obs"), mapping = aes(x = age, y = prop), fill = "darkgreen", alpha = 0.85) +
  geom_line(lc_srv7_male %>% filter(type == "pred"), mapping = aes(x = age, y = prop), size = 1) +
  facet_wrap(~year, scales = "free_x", ncol = 2, dir = "v") +
  labs(x = "Ages", y = "Proportion", title = "Domestic Trawl Survey Length Compositions (Male)") +
  theme_reg() +
  theme(legend.position = "none")

dev.off()


# Selectivity -------------------------------------------------------------
sab_curr <- dget(paste0(dir_master,"//tem.rdat",sep='')) 
sab_rep <- readLines(paste0(dir_master,"//sable.rep",sep=''))

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
# pdf(paste0(dir_results,"//Selectivity.pdf"), width = 10)
ggplot(selex_df, aes(x = ages, y = selex, color = sex)) +
  geom_line() +
  geom_point() +
  facet_wrap(~type) +
  labs(x = "Ages", y = "Selectivity") +
  theme_reg() +
  theme(legend.position = "none")
# dev.off()

# Recruitment ~ SSB -------------------------------------------------------
rec_ssb_df = data.frame(year = sab_curr$t.series$year[-c(1:2)], rec = sab_curr$t.series$Recr[-c(1:2)],
                        ssb = sab_curr$t.series$spbiom[-c(1:2)]) # removing first 2 years from time series

pdf(paste0(dir_results,"//Rec_SSB.pdf"))
ggplot(rec_ssb_df, aes(x = ssb, y = rec, label = year - 2)) +
  geom_text(color = "blue") +
  scale_x_continuous(limits = c(0, NA)) +  # NA for no upper limit
  theme_reg() +
  labs(y = "Age 2 Recruits (millions)", x = "SSB (kt)")
dev.off()


# Recruitment, SSB, Catch -------------------------------------------------
rec_ssb_catch = data.frame(year = sab_curr$t.series$year, rec = sab_curr$t.series$Recr,
                           ssb = sab_curr$t.series$spbiom, catch = sab_curr$t.series$Catch_HAL + sab_curr$t.series$Catch_TWL)

# Plot!
pdf(paste0(dir_results,"//Rec_SSB_Catch.pdf"), height = 5)
ggplot(rec_ssb_catch) +
  geom_col(mapping = aes(x = year, y = rec, color = "Recruitment")) +
  geom_line(mapping = aes(x = year, y = ssb / 5, color = "SSB"), size = 1.75) +
  geom_line(mapping = aes(x = year, y = catch / 5, color = "Catch"), size = 1.75) +
  scale_y_continuous(sec.axis = sec_axis(~.*5, name = "SSB or Catch (kt)") ) +
  scale_color_manual(values = c("Recruitment" = "lightblue3", "SSB" = "orange", 
                                "Catch" = "yellow2"),  name = "") + 
  labs(x = "Year", y = "Recruitment (millions of fish)") +
  theme_reg() + theme(legend.position = c(0.15, 0.9)) 
dev.off()


# Numbers at age ----------------------------------------------------------
# Females
n_at_age_f = data.frame(year = rownames(sab_curr$natage.female), sab_curr$natage.female) %>% 
  pivot_longer(!year, names_to = "age", values_to = "N") %>% 
  mutate(age = as.numeric(str_remove(age, "X")), year = as.numeric(year),
         sex = "Female") 

# Males
n_at_age_m = data.frame(year = rownames(sab_curr$natage.male),  sab_curr$natage.male) %>% 
  pivot_longer(!year, names_to = "age", values_to = "N") %>% 
  mutate(age = as.numeric(str_remove(age, "X")), year = as.numeric(year),
         sex = "Male")

n_at_age = rbind(n_at_age_f, n_at_age_m) # bind together

# Create proportions to look at age structure
n_at_age = n_at_age %>% 
  group_by(year) %>% 
  mutate(prop = N / sum(N))

pdf(paste0(dir_results,"//NAA_plots.pdf"), width = 10)

# Plot numbers at age (females)
ggplot(n_at_age %>% filter(sex == "Female"), aes(x = year, y = N)) +
  geom_line(size = 1) +
  facet_wrap(~age, scales = "free") +
  theme_reg() +
  labs(x = "Year", y = "Numbers at age (Females)")

# Plot numbers at age (males)
ggplot(n_at_age %>% filter(sex == "Male"), aes(x = year, y = N)) +
  geom_line(size = 1) +
  facet_wrap(~age, scales = "free") +
  theme_reg() +
  labs(x = "Year", y = "Numbers at age (Males)")

# Plot age-structure of the population (filter to recent years)
ggplot(n_at_age %>% filter(year %in% c(2019)), aes(x = age, y = prop, fill = sex)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~year) +
  theme_reg() +
  theme(legend.position = "top") +
  labs(x = "Age", y = "Proportion", fill = "Sex")

dev.off()

