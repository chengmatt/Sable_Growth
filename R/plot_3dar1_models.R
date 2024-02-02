# Purpose: To plot 3DAR1 models
# Creator: Matthew LH. Cheng (UAF-CFOS)
# Date: 11/16/23

theme_trevor_jordan <- function() {
  theme_bw() +
    theme(panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "white", color = NA),
          axis.text = element_text(size = 13, color = "black"),
          strip.text = element_text(size = 17),
          axis.title = element_text(size = 17),
          legend.position = "top",
          legend.title = element_text(size = 17),
          legend.text = element_text(size = 13))
}

# Set up ------------------------------------------------------------------

library(here)
library(tidyverse)
library(ggthemes)
source(here("R", "functions.R"))

dir.figs = here("figs", "Growth_Figs")
dir.create(dir.figs)

# growth estimates here
laa_ts = read.csv(here("output", "LAA_Models", "Growth_estimates.csv")) 
laa_sd = read.csv(here("output", "LAA_Models", "SDRep_estimates.csv")) 
waa_ts = read.csv(here("output", "WAA_Models", "Growth_estimates.csv"))
waa_sd = read.csv(here("output", "WAA_Models", "SDRep_estimates.csv"))
sd_all = rbind(waa_sd, laa_sd)

# Read in ageing error matrix
ageing_error_df = readLines(here("data", "SQ_Nominal.dat")) # read in dat file
ageage_mat = ageing_error_df[[which(str_detect(ageing_error_df, "#Age age transition matrix")) + 1]] # read in age-age error matrix
ageage_mat = matrix(as.numeric(str_split(ageage_mat, pattern = " ")[[1]] ), ncol = 30, nrow = 30, dimnames = list(c(2:31), c(2:31))) # coerce strings into matrix

# Plot Ageing Error -------------------------------------------------------
ageage_mat = reshape2::melt(ageage_mat) # melt into df
names(ageage_mat) = c("Assigned_Age", "True_Age", "Prob" )

pdf(here('figs', "ageing_error_matrix.pdf"), width = 10, height = 8)
# Now plot!
ggplot(ageage_mat, aes(x = factor(True_Age), y = factor(Assigned_Age), fill = Prob)) +
  geom_tile() +
  scale_fill_distiller(palette = "YlGnBu", direction = 1) +
  labs(x = "True Ages", y = "Assigned Ages", fill = "Probability of Assignment") +
  theme_bw() +
  theme(legend.text = element_text(size = 15),
        legend.title = element_text(size = 20),
        axis.text = element_text(size = 15, color = "black"),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 20),
        legend.key.width=unit(1,"cm"),
        legend.position = "top")
dev.off()

# Plot most recent year (time series) ---------------------------------------------------

# Length models
pdf(here(dir.figs, "LAA_TS.pdf"), width = 22, height = 13)
laa_ts %>% 
  filter(peel == 0, Sex == "female") %>% 
  ggplot(aes(x = Year, y = Mean, ymin = lwr_95, ymax = upr_95, color = re_model, fill = re_model)) +
  geom_line(size = 1) +
  geom_ribbon(alpha = 0.3, color = NA) +
  facet_wrap(~Age, scales = "free") +
  # theme_trevor_jordan() +
  scale_color_manual(values = c("orange", "black")) +
  scale_fill_manual(values = c("orange", "black")) +
  labs(x = "Year", y = "Length-at-age (Females)", fill = "Model", color = "Model")

laa_ts %>% 
  filter(peel == 0, Sex == "Male") %>% 
  ggplot(aes(x = Year, y = Mean, ymin = lwr_95, ymax = upr_95, color = re_model, fill = re_model)) +
  geom_line(size = 1) +
  geom_ribbon(alpha = 0.3, color = NA) +
  facet_wrap(~Age, scales = "free") +
  theme_trevor_jordan() +
  scale_color_manual(values = c("orange", "black")) +
  scale_fill_manual(values = c("orange", "black")) +
  labs(x = "Year", y = "Length-at-age (Males)", fill = "Model", color = "Model")

dev.off()

# Weight models
pdf(here(dir.figs, "WAA_TS.pdf"), width = 22, height = 13)
waa_ts %>% 
  filter(peel == 0, Sex == "female") %>% 
  ggplot(aes(x = Year, y = Mean, ymin = lwr_95, ymax = upr_95, color = re_model, fill = re_model)) +
  geom_line(size = 1) +
  geom_ribbon(alpha = 0.3, color = NA) +
  facet_wrap(~Age, scales = "free") +
  theme_trevor_jordan() +
  scale_color_manual(values = c("orange", "black")) +
  scale_fill_manual(values = c("orange", "black")) +
  labs(x = "Year", y = "Weight-at-age (Females)", fill = "Model", color = "Model")

waa_ts %>% 
  filter(peel == 0, Sex == "Male") %>% 
  ggplot(aes(x = Year, y = Mean, ymin = lwr_95, ymax = upr_95, color = re_model, fill = re_model)) +
  geom_line(size = 1) +
  geom_ribbon(alpha = 0.3, color = NA) +
  facet_wrap(~Age, scales = "free") +
  theme_trevor_jordan() +
  scale_color_manual(values = c("orange", "black")) +
  scale_fill_manual(values = c("orange", "black")) +
  labs(x = "Year", y = "Weight-at-age (Males)", fill = "Model", color = "Model")

dev.off()

pdf(here(dir.figs, "TS_OnePlot.pdf"), width = 10, height = 10)
laa_ts %>% 
  filter(peel == 0, re_model == "3dar1") %>% 
  ggplot(aes(x = Year, y = Mean, ymin = lwr_95, ymax = upr_95, color = factor(Age))) +
  geom_line(size = 1) +
  facet_wrap(~Sex) +
  theme_trevor_jordan() +
  labs(x = "Year", y = "Length-at-age", fill = "Age", color = "Age")

waa_ts %>% 
  filter(peel == 0, re_model == "3dar1") %>% 
  ggplot(aes(x = Year, y = Mean, ymin = lwr_95, ymax = upr_95, color = factor(Age))) +
  geom_line(size = 1) +
  theme_trevor_jordan() +
  facet_wrap(~Sex) +
  labs(x = "Year", y = "Weight-at-age", fill = "Model", color = "Model")
dev.off()

pdf(here(dir.figs, "Grwth_Curve.pdf"), width = 15, height = 13)
laa_ts %>% 
  filter(peel == 0, re_model == "3dar1") %>% 
  ggplot(aes(x = Age, y = Mean, ymin = lwr_95, ymax = upr_95, color = factor(Year))) +
  geom_line(size = 1) +
  facet_wrap(~Sex) +
  theme_trevor_jordan() +
  labs(x = "Age", y = "Length-at-age", fill = "Age", color = "Age")

waa_ts %>% 
  filter(peel == 0, re_model == "3dar1") %>% 
  ggplot(aes(x = Age, y = Mean, ymin = lwr_95, ymax = upr_95, color = factor(Year))) +
  geom_line(size = 1) +
  theme_trevor_jordan() +
  facet_wrap(~Sex) +
  labs(x = "Age", y = "Weight-at-age", fill = "Model", color = "Model")

dev.off()

# Cohort Plot (most recent year) -------------------------------------------------------------

# length models
pdf(here(dir.figs, "LAA_Cohort.pdf"), width = 18, height = 13)
laa_ts %>% filter(Sex == "Male",
                  broodYear %in% c(2002:2017), 
                  Year != 2023,
                  peel == 0) %>% 
ggplot(aes(x = Age, y = Mean, color = re_model, ymin = lwr_95, ymax = upr_95, fill = re_model), size = 1.3) +
  geom_line(size = 1.3) +
  geom_ribbon(alpha = 0.2, color = NA) +
  facet_wrap(~broodYear, scale = "free") +
  labs(x = "Age", y = "Male Length (cm)", fill = "Model", color = "Model") +
  theme(legend.position = "top") +
  scale_color_manual(values = c("orange", "black")) +
  scale_fill_manual(values = c("orange", "black")) +
  theme_trevor_jordan() 

laa_ts %>% filter(Sex == "female",
                  broodYear %in% c(2002:2017), 
                  Year != 2023,
                  peel == 0) %>% 
  ggplot(aes(x = Age, y = Mean, color = re_model, ymin = lwr_95, ymax = upr_95, fill = re_model), size = 1.3) +
  geom_line(size = 1.3) +
  geom_ribbon(alpha = 0.2, color = NA) +
  facet_wrap(~broodYear, scale = "free") +
  labs(x = "Age", y = "Female Length (cm)", fill = "Model", color = "Model") +
  scale_color_manual(values = c("orange", "black")) +
  scale_fill_manual(values = c("orange", "black")) +
  theme(legend.position = "top") +
  theme_trevor_jordan() 
dev.off()


# weight models
pdf(here(dir.figs, "WAA_Cohort.pdf"), width = 18, height = 13)
waa_ts %>% filter(Sex == "Male",
                  broodYear %in% c(2002:2017), 
                  Year != 2023,
                  peel == 0) %>% 
  ggplot(aes(x = Age, y = Mean, color = re_model, ymin = lwr_95, ymax = upr_95, fill = re_model), size = 1.3) +
  geom_line(size = 1.3) +
  geom_ribbon(alpha = 0.2, color = NA) +
  facet_wrap(~broodYear, scale = "free") +
  labs(x = "Age", y = "Male Weight (g)", fill = "Model", color = "Model") +
  theme(legend.position = "top") +
  scale_color_manual(values = c("orange", "black")) +
  scale_fill_manual(values = c("orange", "black")) +
  theme_trevor_jordan() 

waa_ts %>% filter(Sex == "female",
                  broodYear %in% c(2002:2017), 
                  Year != 2023,
                  peel == 0) %>% 
  ggplot(aes(x = Age, y = Mean, color = re_model, ymin = lwr_95, ymax = upr_95, fill = re_model), size = 1.3) +
  geom_line(size = 1.3) +
  geom_ribbon(alpha = 0.2, color = NA) +
  facet_wrap(~broodYear, scale = "free") +
  labs(x = "Age", y = "Female Weight (g)", fill = "Model", color = "Model") +
  scale_color_manual(values = c("orange", "black")) +
  scale_fill_manual(values = c("orange", "black")) +
  theme(legend.position = "top") +
  theme_trevor_jordan() 
dev.off()


# Cohort Plot Retrospective -----------------------------------------------

pdf(here(dir.figs, "LAA_Cohort_Retro.pdf"), width = 18, height = 13)

laa_ts %>% filter(Sex == "Male",
                  broodYear %in% c(2009:2017),  re_model == "3dar1") %>% 
  group_by(peel) %>% 
  filter(Year <= 2023 - peel - 1) %>% 
  ggplot(aes(x = Age, y = Mean, color = factor(2022 - peel), ymin = lwr_95, 
             ymax = upr_95, fill = factor(2022 - peel)), size = 1.3) +
  geom_line(size = 1.3, alpha = 1) +
  facet_wrap(~broodYear, scales = "free") +
  labs(x = "Age", y = "Male Length (cm)", fill = "Model Peel", color = "Model Peel") +
  theme(legend.position = "top") +
  theme_trevor_jordan() 

laa_ts %>% filter(Sex == "female",
                  broodYear %in% c(2009:2017),  re_model == "3dar1") %>% 
  group_by(peel) %>% 
  filter(Year <= 2023 - peel - 1) %>% 
  ggplot(aes(x = Age, y = Mean, color = factor(2022 - peel), 
             ymin = lwr_95, ymax = upr_95, fill = factor(2022 - peel)), size = 1.3) +
  geom_line(size = 1.3, alpha = 1) +
  facet_wrap(~broodYear, scales = "free") +
  labs(x = "Age", y = "Female Length (cm)", fill = "Model Peel", color = "Model Peel") +
  theme(legend.position = "top") +
  theme_trevor_jordan() 

dev.off()

# waa cohort plots
pdf(here(dir.figs, "WAA_Cohort_Retro.pdf"), width = 18, height = 13)
waa_ts %>% filter(Sex == "Male",
                  broodYear %in% c(2009:2017),  re_model == "3dar1") %>% 
  group_by(peel) %>% 
  filter(Year <= 2023 - peel - 1) %>% 
  ggplot(aes(x = Age, y = Mean, color = factor(2022 - peel), 
             ymin = lwr_95, ymax = upr_95, fill = factor(2022 - peel)), size = 1.3) +
  geom_line(size = 1.3, alpha = 1) +
  facet_wrap(~broodYear, scales = "free") +
  labs(x = "Age", y = "Male Weight (g)", fill = "Model Peel", color = "Model Peel") +
  theme(legend.position = "top") +
  theme_trevor_jordan() 

waa_ts %>% filter(Sex == "female",
                  broodYear %in% c(2009:2017), re_model == "3dar1") %>% 
  group_by(peel) %>% 
  filter(Year <= 2023 - peel - 1) %>% 
  ggplot(aes(x = Age, y = Mean, color = factor(2022 - peel), 
             ymin = lwr_95, ymax = upr_95, fill = factor(2022 - peel)), size = 1.3) +
  geom_line(size = 1.3, alpha = 1) +
  facet_wrap(~broodYear, scales = "free") +
  labs(x = "Age", y = "Female Weight (g)", fill = "Model Peel", color = "Model Peel") +
  theme(legend.position = "top") +
  theme_trevor_jordan() 
dev.off()


# Retrospective across time -----------------------------------------------

# Calculate LAA relative retrospective
ref_laa_ts = laa_ts %>% filter(re_model == "3dar1", peel == 0) %>% 
  select(Year, Age, Mean, Sex) %>% rename(Ref = Mean)

retro_laa = laa_ts %>% 
  left_join(ref_laa_ts, by = c("Age", "Sex", "Year")) %>%
  distinct() %>% 
  mutate(Anom = (Mean-Ref) / Ref)

pdf(here(dir.figs, "LAA_Rel_Retro.pdf"), width = 18, height = 13)
ggplot(retro_laa %>% filter(re_model == "3dar1",
                            Sex == "female"), 
       aes(x = Year, y = Anom, color = factor(2022 - peel))) +
  geom_line(size = 1) +
  facet_wrap(~Age) +
  labs(x = "Year", y = "Female Length Retrospective", color = "Model Peel") +
  theme(legend.position = "top") +
  theme_trevor_jordan()

ggplot(retro_laa %>% filter(re_model == "3dar1",
                            Sex == "Male"), 
       aes(x = Year, y = Anom, color = factor(2022 - peel))) +
  geom_line(size = 1) +
  facet_wrap(~Age) +
  labs(x = "Year", y = "Male Length Retrospective", color = "Model Peel") +
  theme(legend.position = "top") +
  theme_trevor_jordan()
dev.off()

# laa retrospectie across time
pdf(here(dir.figs, "LAA_Year_Retro.pdf"), width = 18, height = 13)
laa_ts %>% filter(Sex == "female", re_model == "3dar1") %>% 
  group_by(peel) %>% 
  filter(Year <= 2023 - peel - 1) %>% 
  ggplot(aes(x = Year, y = Mean, color = factor(2022 - peel), 
             fill = factor(2022 - peel)), size = 1.3) +
  geom_line(size = 1.3, alpha = 1) +
  facet_wrap(~Age, scale = "free") +
  labs(x = "Year", y = "Female Length (cm)", fill = "Model Peel", color = "Model Peel") +
  theme(legend.position = "top") +
  theme_trevor_jordan()

laa_ts %>% filter(Sex == "Male", re_model == "3dar1") %>% 
  group_by(peel) %>% 
  filter(Year <= 2023 - peel - 1) %>% 
  ggplot(aes(x = Year, y = Mean, color = factor(2022 - peel), 
             fill = factor(2022 - peel)), size = 1.3) +
  geom_line(size = 1.3, alpha = 1) +
  facet_wrap(~Age, scale = "free") +
  labs(x = "Year", y = "Male Length (cm)", fill = "Model Peel", color = "Model Peel") +
  theme(legend.position = "top") +
  theme_trevor_jordan()
dev.off()

# Calculate LAA relative retrospective
ref_waa_ts = waa_ts %>% filter(re_model == "3dar1", peel == 0) %>% 
  select(Year, Age, Mean, Sex) %>% rename(Ref = Mean)

retro_waa = waa_ts %>% 
  left_join(ref_waa_ts, by = c("Age", "Sex", "Year")) %>%
  distinct() %>% 
  mutate(Anom = (Mean-Ref) / Ref)

pdf(here(dir.figs, "WAA_Rel_Retro.pdf"), width = 18, height = 13)
ggplot(retro_waa %>% filter(re_model == "3dar1",
                            Sex == "female"), 
       aes(x = Year, y = Anom, color = factor(2022 - peel))) +
  geom_line(size = 1) +
  facet_wrap(~Age) +
  labs(x = "Year", y = "Female Weight Retrospective", color = "Model Peel") +
  theme(legend.position = "top") +
  theme_trevor_jordan()

ggplot(retro_waa %>% filter(re_model == "3dar1",
                            Sex == "Male"), 
       aes(x = Year, y = Anom, color = factor(2022 - peel))) +
  geom_line(size = 1) +
  facet_wrap(~Age) +
  labs(x = "Year", y = "Male Weight Retrospective", color = "Model Peel") +
  theme(legend.position = "top") +
  theme_trevor_jordan()
dev.off()

# waa retrospectie across time
pdf(here(dir.figs, "WAA_Year_Retro.pdf"), width = 18, height = 13)
waa_ts %>% filter(Sex == "female", re_model == "3dar1") %>% 
  group_by(peel) %>% 
  filter(Year <= 2023 - peel - 1) %>% 
  ggplot(aes(x = Year, y = Mean, color = factor(2022 - peel), 
             fill = factor(2022 - peel)), size = 1.3) +
  geom_line(size = 1.3, alpha = 1) +
  facet_wrap(~Age, scale = "free") +
  labs(x = "Year", y = "Female Weight (g)", fill = "Model Peel", color = "Model Peel") +
  theme(legend.position = "top") +
  theme_trevor_jordan()

waa_ts %>% filter(Sex == "Male", re_model == "3dar1") %>% 
  group_by(peel) %>% 
  filter(Year <= 2023 - peel - 1) %>% 
  ggplot(aes(x = Year, y = Mean, color = factor(2022 - peel), 
             fill = factor(2022 - peel)), size = 1.3) +
  geom_line(size = 1.3, alpha = 1) +
  facet_wrap(~Age, scale = "free") +
  labs(x = "Year", y = "Male Weight (g)", fill = "Model Peel", color = "Model Peel") +
  theme(legend.position = "top") +
  theme_trevor_jordan()
dev.off()


# Retrospective Correlations ----------------------------------------------
pdf(here(dir.figs, "Correlations_Retro.pdf"), width = 15)
sd_all %>% filter(str_detect(Parameter, "rho"),
                  Growth_model == "weight") %>% 
  ggplot(aes(x = 2022 - peel, y = Trans_Estimate, ymin = lwr_95, ymax = upr_95, color = Sex)) +
  geom_pointrange() +
  facet_grid(Sex~Parameter) +
  geom_hline(yintercept = 0, lty = 2, size = 1.3) +
  theme_trevor_jordan() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 90)) +
  labs(y = "Weight correlation values", x = "Model Peel") 
  
sd_all %>% filter(str_detect(Parameter, "rho"),
                  Growth_model == "length") %>% 
  ggplot(aes(x = 2022 - peel, y = Trans_Estimate, ymin = lwr_95, ymax = upr_95, color = Sex)) +
  geom_pointrange() +
  facet_grid(Sex~Parameter) +
  geom_hline(yintercept = 0, lty = 2, size = 1.3) +
  theme_trevor_jordan() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 90)) +
  labs(y = " Length correlation values", x = "Model Peel") 
dev.off()


# Anomaly plot ------------------------------------------------------------

# Calculate LAA anomaly
const_laa_ts = laa_ts %>% filter(re_model == "constant",
                                 Year == 1996, peel == 0) %>% 
  select(Age, Mean, Sex) %>% rename(Ref = Mean)
anom_laa = laa_ts %>% filter(re_model == "3dar1") %>% 
  left_join(const_laa_ts, by = c("Age", "Sex")) %>% 
  mutate(Anom = (Mean-Ref) / Ref)

# Calculate WAA anomaly
const_waa_ts = waa_ts %>% filter(re_model == "constant", Year == 1996, peel == 0) %>% 
  select(Age, Mean, Sex) %>% rename(Ref = Mean)
anom_waa = waa_ts %>% filter(re_model == "3dar1") %>% 
  left_join(const_waa_ts, by = c("Age", "Sex")) %>% 
  mutate(Anom = (Mean-Ref) / Ref)

pdf(here(dir.figs, "LAA_Anomaly.pdf"), width = 15)
ggplot(anom_laa %>% filter(peel == 0,
                           ), aes(x = factor(Year), y = factor(Age), fill = Anom)) +
  geom_tile(alpha = 1) +
  scale_fill_gradient2(midpoint = mean(anom_laa$Anom), 
                       high = scales::muted("red"), 
                       low = scales::muted("blue")) +
  facet_wrap(~Sex) +
  scale_y_discrete(breaks = seq(3, 31, 3)) +
  scale_x_discrete(breaks = seq(1996, 2023, 5)) +
  theme_trevor_jordan() +
  theme(legend.key.width = unit(3, "cm")) +
  labs(x = "Year", y = "Age", fill = "Length Anomalies")
dev.off()

pdf(here(dir.figs, "WAA_Anomaly.pdf"), width = 15)
ggplot(anom_waa %>% filter(peel == 0), 
       aes(x = factor(Year), y = factor(Age), fill = Anom)) +
  geom_tile(alpha = 1) +
  scale_fill_gradient2(midpoint = mean(anom_waa$Anom), 
                       high = scales::muted("red"), 
                       low = scales::muted("blue")) +
  facet_wrap(~Sex) +
  scale_y_discrete(breaks = seq(3, 31, 3)) +
  scale_x_discrete(breaks = seq(1996, 2023, 5)) +
  theme_trevor_jordan() +
  theme(legend.key.width = unit(3, "cm")) +
  labs(x = "Year", y = "Age", fill = "Weight Anomalies")
dev.off()

# Plot Size-Age Transition Matrix -----------------------------------------

age_bins = 2:31
len_bins = seq(41, 99, 2)
al_mat_all_df = data.frame()
years = unique(laa_ts$Year)
n_proj_years = 1
# extract out length models
length_models = list.files(here("output", "LAA_Models"))
length_models = length_models[str_detect(length_models, ".RData")]

for(i in 1:length(length_models)) {
  
    load(here("output", "LAA_Models", paste(length_models[i], sep = "")))

    # Extract out quantities
    laa_sigma = model$rep$obs_sigma_at
    laa_mat = matrix(model$sd_rep$value[names(model$sd_rep$value) == "mu_at"],
                     nrow = length(age_bins), ncol = length(years) )
    
    # Set up array
    al_array = array(data = 0, dim = c(length(age_bins), length(len_bins), ncol(laa_mat)),
                     dimnames = list(c(age_bins), c(len_bins), c(sort(years)) ))
    
    
    for(t in 1:ncol(laa_mat)) { # loop through time to get matrix
      al_array[,,t] = get_al_trans_matrix(age_bins = 2:31, len_bins = seq(41, 99, 2), 
                                          mean_length = laa_mat[,t], sd = laa_sigma[,t])
    } # end t loop
    
    plot(al_array[1,,1])
    
    # naming and mungingz
    al_df = reshape2::melt(al_array)
    names(al_df) = c("Age", "Length", "Year", "Value")
    al_df$sex = str_split(length_models[i], "_")[[1]][2]
    al_mat_all_df = rbind(al_mat_all_df, al_df)
    
} # end t loop

pdf(here(dir.figs, "AL_Matrix_3DAR1.pdf"), width = 15, height = 13)
val_range = range(al_mat_all_df$Value) # get range of values
print(
  ggplot(al_mat_all_df %>% filter(sex == "female"), aes(x = factor(Age), y = factor(Length), fill = Value)) +
    geom_tile() + 
    scale_fill_distiller(palette = "YlGnBu", direction = 1, limits = val_range) +
    guides(fill = guide_colourbar(barwidth = 12)) +
    facet_wrap(~Year) +
    scale_y_discrete(breaks = seq(41, 99, 5)) +
    scale_x_discrete(breaks = seq(3, 31, 3)) +
    labs(x = "Age", y = "Female Length (cm)", fill = "Probability") +
    theme_trevor_jordan() 
)

print(
  ggplot(al_mat_all_df %>% filter(sex == "Male"), aes(x = factor(Age), y = factor(Length), fill = Value)) +
    geom_tile() + 
    scale_fill_distiller(palette = "YlGnBu", direction = 1, limits = val_range) +
    guides(fill = guide_colourbar(barwidth = 12)) +
    facet_wrap(~Year) +
    scale_y_discrete(breaks = seq(41, 99, 5)) +
    scale_x_discrete(breaks = seq(3, 31, 3)) +
    labs(x = "Age", y = "Male Length (cm)", fill = "Probability") +
    theme_trevor_jordan() 
)
dev.off()

