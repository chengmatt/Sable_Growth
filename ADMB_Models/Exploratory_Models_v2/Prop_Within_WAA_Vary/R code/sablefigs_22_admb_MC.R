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
# setwd('Prop_Within_Growth_SQ')
dir_master<-here("Prop_Within_Growth_Vary")
dir_results<-paste0(dir_master,"//Results",sep='')
dir.create(dir_results)


# Read in .rdat
sab_curr <- dget(paste0(dir_master,"//tem.rdat",sep='')) 
sab_rep <- readLines(paste0(dir_master,"//sable.rep",sep=''))
ages = 2:31

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
pdf(paste0(dir_results,"//Composition_Fits.pdf"), width = 15, height = 15)
### Domestic LL Fishery Age Compositions -------------------------------
obs_ac_fish1_f = data.frame(reshape2::melt(sab_curr$oac.fish1.f), type = "obs")
pred_ac_fish1_f = data.frame(reshape2::melt(sab_curr$eac.fish1.f), type = "pred")

# Put these into a dataframe
ac_fish1_f = rbind(obs_ac_fish1_f, pred_ac_fish1_f)
names(ac_fish1_f) = c("year", "age", "prop", "type")
ac_fish1_f$broodYear = ac_fish1_f$year - ac_fish1_f$age # create brood year to track cohort over time
ac_fish1_f = ac_fish1_f %>% mutate(broodYear = ifelse(age == 31, "31", broodYear)) # make plus group consistent color

ggplot() +
  geom_col(ac_fish1_f %>% filter(type == "obs"), mapping = aes(x = age, y = prop, fill = factor(broodYear))) +
  geom_line(ac_fish1_f %>% filter(type == "pred"), mapping = aes(x = age, y = prop), size = 1) +
  facet_wrap(~year, scales = "free_x", dir = "v") +
  labs(x = "Ages", y = "Proportion", title = "Domestic LL Fishery Age Compositions Females") +
  theme_reg() +
  theme(legend.position = "none")

obs_ac_fish1_m = data.frame(reshape2::melt(sab_curr$oac.fish1.m), type = "obs")
pred_ac_fish1_m = data.frame(reshape2::melt(sab_curr$eac.fish1.m), type = "pred")

# Put these into a dataframe
ac_fish1_m = rbind(obs_ac_fish1_m, pred_ac_fish1_m)
names(ac_fish1_m) = c("year", "age", "prop", "type")
ac_fish1_m$broodYear = ac_fish1_m$year - ac_fish1_m$age # create brood year to track cohort over time
ac_fish1_m = ac_fish1_m %>% mutate(broodYear = ifelse(age == 31, "31", broodYear)) # make plus group consistent color

ggplot() +
  geom_col(ac_fish1_m %>% filter(type == "obs"), mapping = aes(x = age, y = prop, fill = factor(broodYear))) +
  geom_line(ac_fish1_m %>% filter(type == "pred"), mapping = aes(x = age, y = prop), size = 1) +
  facet_wrap(~year, scales = "free_x", dir = "v") +
  labs(x = "Ages", y = "Proportion", title = "Domestic LL Fishery Age Compositions Males") +
  theme_reg() +
  theme(legend.position = "none")

# plot average fits as well
avg_ac_fish1_m = ac_fish1_m %>%
  group_by(age, type) %>% 
  summarize(prop = mean(prop)) %>% 
  mutate(sex = "Male")

# plot average fits as well
avg_ac_fish1_f = ac_fish1_f %>%
  group_by(age, type) %>% 
  summarize(prop = mean(prop)) %>% 
  mutate(sex = "Female")

avg_ac_fish1 = rbind(avg_ac_fish1_m, avg_ac_fish1_f)

ggplot() +
  geom_col(avg_ac_fish1 %>% filter(type == "obs"), mapping = aes(x = age, y = prop), 
           fill = "orange", alpha = 0.35, color = "black") +  
  geom_line(avg_ac_fish1 %>% filter(type == "pred"), mapping = aes(x = age, y = prop), size = 1) +
  labs(x = "Ages", y = "Proportion", title = "Average Domestic LL Fishery Age Compositions") +
  theme_reg() +
  facet_wrap(~sex) + 
  theme(legend.position = "none")


### Domestic LL Survey Age Compositions -------------------------------
obs_ac_srv1_f = data.frame(reshape2::melt(sab_curr$oac.srv1.f), type = "obs")
pred_ac_srv1_f = data.frame(reshape2::melt(sab_curr$eac.srv1.f), type = "pred")

# Put these into a dataframe
ac_srv1_f = rbind(obs_ac_srv1_f, pred_ac_srv1_f)
names(ac_srv1_f) = c("year", "age", "prop", "type")
ac_srv1_f$broodYear = ac_srv1_f$year - ac_srv1_f$age # create brood year to track cohort over time
ac_srv1_f = ac_srv1_f %>% mutate(broodYear = ifelse(age == 31, "31", broodYear)) # make plus group consistent color

ggplot() +
  geom_col(ac_srv1_f %>% filter(type == "obs"), mapping = aes(x = age, y = prop, fill = factor(broodYear))) +
  geom_line(ac_srv1_f %>% filter(type == "pred"), mapping = aes(x = age, y = prop), size = 1) +
  facet_wrap(~year, scales = "free_x", dir = "v") +
  labs(x = "Ages", y = "Proportion", title = "Domestic Survey Fishery Age Compositions Females") +
  theme_reg() +
  theme(legend.position = "none")

obs_ac_srv1_m = data.frame(reshape2::melt(sab_curr$oac.srv1.m), type = "obs")
pred_ac_srv1_m = data.frame(reshape2::melt(sab_curr$eac.srv1.m), type = "pred")

# Put these into a dataframe
ac_srv1_m = rbind(obs_ac_srv1_m, pred_ac_srv1_m)
names(ac_srv1_m) = c("year", "age", "prop", "type")
ac_srv1_m$broodYear = ac_srv1_m$year - ac_srv1_m$age # create brood year to track cohort over time
ac_srv1_m = ac_srv1_m %>% mutate(broodYear = ifelse(age == 31, "31", broodYear)) # make plus group consistent color

ggplot() +
  geom_col(ac_srv1_m %>% filter(type == "obs"), mapping = aes(x = age, y = prop, fill = factor(broodYear))) +
  geom_line(ac_srv1_m %>% filter(type == "pred"), mapping = aes(x = age, y = prop), size = 1) +
  facet_wrap(~year, scales = "free_x", dir = "v") +
  labs(x = "Ages", y = "Proportion", title = "Domestic Survey Fishery Age Compositions Males") +
  theme_reg() +
  theme(legend.position = "none")

# plot average fits as well
avg_ac_srv1_m = ac_srv1_m %>%
  group_by(age, type) %>% 
  summarize(prop = mean(prop)) %>% 
  mutate(sex = "Male")

# plot average fits as well
avg_ac_srv1_f = ac_srv1_f %>%
  group_by(age, type) %>% 
  summarize(prop = mean(prop)) %>% 
  mutate(sex = "Female")

avg_ac_fish1 = rbind(avg_ac_srv1_m, avg_ac_srv1_f)

ggplot() +
  geom_col(avg_ac_fish1 %>% filter(type == "obs"), mapping = aes(x = age, y = prop), 
           fill = "orange", alpha = 0.35, color = "black") +  
  geom_line(avg_ac_fish1 %>% filter(type == "pred"), mapping = aes(x = age, y = prop), size = 1) +
  labs(x = "Ages", y = "Proportion", title = "Average Domestic Survey Fishery Age Compositions") +
  theme_reg() +
  facet_wrap(~sex) + 
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
  labs(x = "Length", y = "Proportion", title = "Domestic LL Fishery Length Compositions (Male)") +
  theme_reg() +
  theme(legend.position = "none")

obs_lc_fish1_female = data.frame(reshape2::melt(sab_curr$olc.fish1.f), type = "obs")
pred_lc_fish1_female = data.frame(reshape2::melt(sab_curr$elc.fish1.f), type = "pred")

# Put these into a dataframe
lc_fish1_female = rbind(obs_lc_fish1_female, pred_lc_fish1_female)
names(lc_fish1_female) = c("year", "age", "prop", "type")

ggplot() +
  geom_col(lc_fish1_female %>% filter(type == "obs"), mapping = aes(x = age, y = prop), fill = "darkgreen", alpha = 0.85) +
  geom_line(lc_fish1_female %>% filter(type == "pred"), mapping = aes(x = age, y = prop), size = 1) +
  facet_wrap(~year, scales = "free_x", dir = "v") +
  labs(x = "Length", y = "Proportion", title = "Domestic LL Fishery Length Compositions (Female)") +
  theme_reg() +
  theme(legend.position = "none")

# Plot average fits
avg_lc_fish1_m = lc_fish1_male %>%
  group_by(age, type) %>% 
  summarize(prop = mean(prop)) %>% 
  mutate(sex = "Male")

# plot average fits as well
avg_lc_fish1_f = lc_fish1_female %>%
  group_by(age, type) %>% 
  summarize(prop = mean(prop)) %>% 
  mutate(sex = "Female")

avg_lc_fish1 = rbind(avg_lc_fish1_f, avg_lc_fish1_m)

ggplot() +
  geom_col(avg_lc_fish1 %>% filter(type == "obs"), mapping = aes(x = age, y = prop), 
           fill = "darkgreen", alpha = 0.2, color = "black") +  
  geom_line(avg_lc_fish1 %>% filter(type == "pred"), mapping = aes(x = age, y = prop), size = 1) +
  labs(x = "Length", y = "Proportion", title = "Average Domestic LL Fishery Length Compositions") +
  theme_reg() +
  facet_wrap(~sex) + 
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
  labs(x = "Length", y = "Proportion", title = "Domestic Trawl Fishery Length Compositions (Female)") +
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
  facet_wrap(~year, scales = "free_x", dir = "v") +
  labs(x = "Length", y = "Proportion", title = "Domestic Trawl Fishery Length Compositions (Male)") +
  theme_reg() +
  theme(legend.position = "none")

# Plot average fits
avg_lc_fish3_m = lc_fish3_male %>%
  group_by(age, type) %>% 
  summarize(prop = mean(prop)) %>% 
  mutate(sex = "Male")

# plot average fits as well
avg_lc_fish3_f = lc_fish3_female %>%
  group_by(age, type) %>% 
  summarize(prop = mean(prop)) %>% 
  mutate(sex = "Female")

avg_lc_fish3 = rbind(avg_lc_fish3_f, avg_lc_fish3_m)

ggplot() +
  geom_col(avg_lc_fish3 %>% filter(type == "obs"), mapping = aes(x = age, y = prop), 
           fill = "darkgreen", alpha = 0.2, color = "black") +  
  geom_line(avg_lc_fish3 %>% filter(type == "pred"), mapping = aes(x = age, y = prop), size = 1) +
  labs(x = "Ages", y = "Proportion", title = "Average Domestic LL Trawl Length Compositions") +
  theme_reg() +
  facet_wrap(~sex) + 
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
  facet_wrap(~year, scales = "free_x", dir = "v") +
  labs(x = "Length", y = "Proportion", title = "Domestic LL Survey Length Compositions (Female)") +
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
  facet_wrap(~year, scales = "free_x", dir = "v") +
  labs(x = "Ages", y = "Proportion", title = "Domestic LL Survey Length Compositions (Male)") +
  theme_reg() +
  theme(legend.position = "none")

# Plot average fits
avg_lc_srv1_m = lc_srv1_male %>%
  group_by(age, type) %>% 
  summarize(prop = mean(prop)) %>% 
  mutate(sex = "Male")

# plot average fits as well
avg_lc_srv1_f = lc_srv1_female %>%
  group_by(age, type) %>% 
  summarize(prop = mean(prop)) %>% 
  mutate(sex = "Female")

avg_lc_srv1 = rbind(avg_lc_srv1_f, avg_lc_srv1_m)

ggplot() +
  geom_col(avg_lc_srv1 %>% filter(type == "obs"), mapping = aes(x = age, y = prop), 
           fill = "darkgreen", alpha = 0.2, color = "black") +  
  geom_line(avg_lc_srv1 %>% filter(type == "pred"), mapping = aes(x = age, y = prop), size = 1) +
  labs(x = "Length", y = "Proportion", title = "Average Domestic LL Survey Length Compositions") +
  theme_reg() +
  facet_wrap(~sex) + 
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
  facet_wrap(~year, scales = "free_x",  dir = "v") +
  labs(x = "Length", y = "Proportion", title = "Japanese LL Survey Length Compositions (Female)") +
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
  facet_wrap(~year, scales = "free_x",  dir = "v") +
  labs(x = "Length", y = "Proportion", title = "Japanese LL Survey Length Compositions (Male)") +
  theme_reg() +
  theme(legend.position = "none")

# Plot average fits
avg_lc_srv2_m = lc_srv2_male %>%
  group_by(age, type) %>% 
  summarize(prop = mean(prop)) %>% 
  mutate(sex = "Male")

# plot average fits as well
avg_lc_srv2_f = lc_srv2_female %>%
  group_by(age, type) %>% 
  summarize(prop = mean(prop)) %>% 
  mutate(sex = "Female")

avg_lc_srv2 = rbind(avg_lc_srv2_f, avg_lc_srv2_m)

ggplot() +
  geom_col(avg_lc_srv2 %>% filter(type == "obs"), mapping = aes(x = age, y = prop), 
           fill = "darkgreen", alpha = 0.2, color = "black") +  
  geom_line(avg_lc_srv2 %>% filter(type == "pred"), mapping = aes(x = age, y = prop), size = 1) +
  labs(x = "Length", y = "Proportion", title = "Average Japanese LL Trawl Length Compositions") +
  theme_reg() +
  facet_wrap(~sex) + 
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
  labs(x = "Length", y = "Proportion", title = "Domestic Trawl Survey Length Compositions (Female)") +
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

# Plot average fits
avg_lc_srv7_m = lc_srv7_male %>%
  group_by(age, type) %>% 
  summarize(prop = mean(prop)) %>% 
  mutate(sex = "Male")

# plot average fits as well
avg_lc_srv7_f = lc_srv7_female %>%
  group_by(age, type) %>% 
  summarize(prop = mean(prop)) %>% 
  mutate(sex = "Female")

avg_lc_srv7 = rbind(avg_lc_srv7_f, avg_lc_srv7_m)

ggplot() +
  geom_col(avg_lc_srv7 %>% filter(type == "obs"), mapping = aes(x = age, y = prop), 
           fill = "darkgreen", alpha = 0.2, color = "black") +  
  geom_line(avg_lc_srv7 %>% filter(type == "pred"), mapping = aes(x = age, y = prop), size = 1) +
  labs(x = "Ages", y = "Proportion", title = "Average Domestic Trawl Survey Length Compositions") +
  theme_reg() +
  facet_wrap(~sex) + 
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
pdf(paste0(dir_results,"//Selectivity.pdf"), width = 10)
ggplot(selex_df, aes(x = ages, y = selex, color = sex)) +
  geom_line() +
  geom_point() +
  facet_wrap(~type) +
  labs(x = "Ages", y = "Selectivity") +
  theme_reg() +
  theme(legend.position = "none")
dev.off()

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

pdf(paste0(dir_results,"//NAA_plots.pdf"), width = 15)

# Plot numbers at age (females)
ggplot(n_at_age, aes(x = year, y = N, color = sex)) +
  geom_line(size = 1) +
  facet_wrap(~age, scales = "free") +
  theme_reg() +
  labs(x = "Year", y = "Numbers at age")

# Plot age-structure of the population (filter to recent years)
ggplot(n_at_age %>% filter(year %in% c(1965)), aes(x = age, y = prop, fill = sex)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~year) +
  theme_reg() +
  theme(legend.position = "top") +
  labs(x = "Age", y = "Proportion", fill = "Sex")

dev.off()

