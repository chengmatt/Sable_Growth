# Purpose: To run growth models
# Creator: Matthew LH. Cheng (UAF-CFOS)
# Date 11/5/23

theme_trevor_jordan <- function() {
  theme_bw() +
    theme(panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "white", color = NA),
          axis.text = element_text(size = 15, color = "black"),
          strip.text = element_text(size = 19),
          axis.title = element_text(size = 19),
          legend.position = "top",
          legend.title = element_text(size = 17),
          legend.text = element_text(size = 14))
}

# Set up ------------------------------------------------------------------

library(TMB)
library(TMBhelper)
library(here)
library(tidyverse)

setwd("src")
TMB::compile("Growth_Model.cpp")
dyn.unload(dynlib('Growth_Model'))
dyn.load(dynlib('Growth_Model'))
source(here("R", "functions.R"))
dir.out = here(dir.out, "exploratory_models")
dir.create(dir.out)

# Do some quick data munging and load in data
age_dat = read.csv(here("data", "age_view.csv"), check.names = FALSE) %>% 
  filter(!is.na(Age), !is.na(`Length (cm)`), !is.na(`Weight (g)`),
         `Sex Description` != "Unknown", `Error Flag` == 0) %>% 
  rename(Length = `Length (cm)`,
         Sex_name = `Sex Description`,
         Weight = `Weight (g)`) 

# set up stuff here
n_sexes = 2 # number of sexes
growth_model = c("length", "weight") # growth models
sexes = c("female", "Male") # sexes
re_model = c("iid", "2dar1", "3dar1", "constant") # random effects models
var_param = 0 # conditional variance for gmrf
age_offset = 2 # start age
mod_outputs = data.frame()
mod_list = list()
mod_counter = 1
n_proj_years = 10

# Loop through to run models
for(s in 1:n_sexes) { # sex loop
  
  # Filter to model data
  mod_df = age_dat %>% 
    filter(Sex_name == sexes[s]) %>% 
    filter(!Age %in% c(0, 1)) %>% 
    mutate(Age = ifelse(Age >= 31, 31, Age), # collapse ages to + group
           Age = Age - min(Age),
           Year = Year - min(Year)) %>% 
    select(Length, Weight, Age, Year) 
  
  for(g in 1:length(growth_model)) { # growth model loop
    for(r in 1:length(re_model)) { # random effects loop
      
      # set up model inputs
      model_inputs = prepare_data(data = mod_df, 
                                  growth_model_name = growth_model[g], 
                                  re_model = re_model[r], 
                                  var_param = var_param, 
                                  age_offset = age_offset,
                                  n_proj_years = n_proj_years)
      
      # set up AD object
      model = MakeADFun(data = model_inputs$data, 
                        parameters = model_inputs$parameters, 
                        DLL = "Growth_Model", random = "ln_eps_at",
                        map = model_inputs$map, silent = FALSE)
      
      # optimize
      optim = stats::nlminb(model$par, model$fn, model$gr,  
                            control = list(iter.max = 1e6, eval.max = 1e6))
      
      # take additional newton steps
      add_newton(n.newton = 3, ad_model = model, mle_optim = optim)
      model$optim = optim # optimizer information
      model$report = model$report() # get report
      model$sd_rep = sdreport(model) # get standard errors
      model$name = paste(sexes[s], growth_model[g], re_model[r]) # model name
      
      # extract out estimates
      growth = reshape2::melt(matrix(model$sd_rep$value[names(model$sd_rep$value) == "mu_at"], 
                                  nrow = length(unique(mod_df$Age)),
                                  ncol = length(unique(mod_df$Year)) + n_proj_years )) # save mean estimates
      
      growth_sd = reshape2::melt(matrix(model$sd_rep$sd[names(model$sd_rep$value) == "mu_at"], 
                                 nrow = length(unique(mod_df$Age)),
                                 ncol = length(unique(mod_df$Year)) + n_proj_years )) # standard errors from delta method
      
      # save estimates, etc.
      growth_df = growth %>% 
        left_join(growth_sd, by = c("Var1", "Var2")) %>% 
        rename(Age = Var1, Year = Var2,
               Mean = value.x, SD = value.y) %>% 
        mutate(Age = Age + 1, Year = Year + 1995,
               broodYear = Year-Age,
               AIC = TMBhelper::TMBAIC(model$optim)) 
      
      if(growth_model[g] == "length") {
        growth_df = growth_df %>% 
          mutate(lwr_95 = Mean - 1.96*SD,
                 upr_95 = Mean + 1.96*SD,
                 sex = sexes[s], re_model = re_model[r], 
                 type = growth_model[g] )
      } # length-at-age
      
      if(growth_model[g] == "weight") {
        growth_df = growth_df %>% 
          mutate(lwr_95 = exp(Mean - 1.96*SD),
                 upr_95 = exp(Mean + 1.96*SD),
                 Mean = exp(Mean),
                 sex = sexes[s], re_model = re_model[r], 
                 type = growth_model[g] )
      } # weight-at-age
      
      mod_outputs = rbind(mod_outputs, growth_df) # bind all outputs
      
      mod_list[[mod_counter]] = model # save model here
      mod_counter = mod_counter + 1 # counter for models
      
    } # end r loop
  } # end g loop
} # end s loop

save(mod_list, file =  here(dir.out, "models.RData"))
write.csv(mod_outputs, here(dir.out, "growth_estimates.csv"))


# Summary and model exploratory plots -----------------------------------------------------------

load(here(dir.out, "models.RData"))
mod_outputs = read.csv(here(dir.out, "growth_estimates.csv")) %>% 
  filter(re_model %in% c("constant", "3dar1"),
         Year <= 2022)

# Filter to model data
mod_df = age_dat %>% 
  filter(Sex_name == "female") %>% 
  filter(!Age %in% c(0, 1)) %>% 
  mutate(Age = ifelse(Age >= 31, 31, Age), # collapse ages to + group
         Age = Age - min(Age),
         Year = Year - min(Year)) %>% 
  select(Length, Weight, Age, Year) 


pdf(here(dir.out, "Growth_trends.pdf"), width = 11.5, height = 8)
# WAA over time (Females)
mod_outputs %>% 
  filter(sex == "female", type == "weight") %>% 
  ggplot(aes(x = Year, y = Mean, ymin = lwr_95, ymax = upr_95)) +
  geom_ribbon(alpha = 0.3, aes(fill = re_model)) +
  geom_line(aes(color = re_model)) +
  facet_wrap(~Age, scales = "free") +
  labs(x = "Year", y = "Weight (g)", color = "Model",
       fill = "Model", title = "Female WAA") +
  theme(legend.position = "top") 
  # geom_vline(xintercept = 2022)

# WAA over time (Males)
mod_outputs %>% 
  filter(sex == "Male", type == "weight") %>% 
  ggplot(aes(x = Year, y = Mean, ymin = lwr_95, ymax = upr_95)) +
  geom_ribbon(alpha = 0.3, aes(fill = re_model)) +
  geom_line(aes(color = re_model)) +
  facet_wrap(~Age, scales = "free") +
  labs(x = "Year", y = "Weight (g)", color = "Model",
       fill = "Model", title = "Male WAA") +
  theme(legend.position = "top") 
  # geom_vline(xintercept = 2022)

# LAA over time (Females)
mod_outputs %>% 
  filter(sex == "female", type == "length") %>% 
  ggplot(aes(x = Year, y = Mean, ymin = lwr_95, ymax = upr_95)) +
  geom_ribbon(alpha = 0.3, aes(fill = re_model)) +
  geom_line(aes(color = re_model)) +
  facet_wrap(~Age, scales = "free") +
  labs(x = "Year", y = "Length (cm)", color = "Model",
       fill = "Model", title = "Female LAA") +
  theme(legend.position = "top") 
  # geom_vline(xintercept = 2022)

# LAA over time (Males)
mod_outputs %>% 
  filter(sex == "Male", type == "length") %>% 
  ggplot(aes(x = Year, y = Mean, ymin = lwr_95, ymax = upr_95)) +
  geom_ribbon(alpha = 0.3, aes(fill = re_model)) +
  geom_line(aes(color = re_model)) +
  facet_wrap(~Age, scales = "free") +
  labs(x = "Year", y = "Length (cm)", color = "Model",
       fill = "Model", title = "Male LAA") +
  theme(legend.position = "top") +
  geom_vline(xintercept = 2022)

# LAA curve over time (Females)
mod_outputs %>% 
  filter(sex == "female", type == "length",
         Year > 2013) %>% 
  ggplot(aes(x = Age, y = Mean, ymin = lwr_95, ymax = upr_95)) +
  geom_ribbon(alpha = 0.3, aes(fill = re_model)) +
  geom_line(aes(color = re_model)) +
  facet_wrap(~Year, scales = "free") +
  labs(x = "Year", y = "Length (cm)", color = "Model",
       fill = "Model", title = "Female LAA") +
  theme(legend.position = "top")

# LAA curve over time (Males)
mod_outputs %>% 
  filter(sex == "Male", type == "length",
         re_model == "constant",
         Year > 2013) %>% 
  ggplot(aes(x = Age, y = Mean, ymin = lwr_95, ymax = upr_95)) +
  geom_ribbon(alpha = 0.3, aes(fill = re_model)) +
  geom_line(aes(color = re_model)) +
  facet_wrap(~Year, scales = "free") +
  labs(x = "Year", y = "Length (cm)", color = "Model",
       fill = "Model", title = "Male LAA") +
  theme(legend.position = "top")

# WAA curve over time (Females)
mod_outputs %>% 
  filter(sex == "female", type == "weight",
         Year > 2013) %>% 
  ggplot(aes(x = Age, y = Mean, ymin = lwr_95, ymax = upr_95)) +
  geom_ribbon(alpha = 0.3, aes(fill = re_model)) +
  geom_line(aes(color = re_model)) +
  facet_wrap(~Year, scales = "free") +
  labs(x = "Year", y = "Weight (g)", color = "Model",
       fill = "Model", title = "Female WAA") +
  theme(legend.position = "top")

# WAA curve over time (Males)
mod_outputs %>% 
  filter(sex == "Male", type == "weight",
         Year > 2013) %>% 
  ggplot(aes(x = Age, y = Mean, ymin = lwr_95, ymax = upr_95)) +
  geom_ribbon(alpha = 0.3, aes(fill = re_model)) +
  geom_line(aes(color = re_model)) +
  facet_wrap(~Year, scales = "free") +
  labs(x = "Year", y = "Weight (g)", color = "Model",
       fill = "Model", title = "Male WAA") +
  theme(legend.position = "top")

mod_outputs %>% 
  filter(sex == "female", type == "weight", broodYear %in% c(2010:2016)) %>% 
  ggplot(aes(x = Age, y = Mean, ymin = lwr_95, ymax = upr_95)) +
  geom_ribbon(alpha = 0.3, aes(fill = re_model)) +
  geom_line(aes(color = re_model)) +
  geom_vline(aes(xintercept = 2018 - broodYear)) +
  facet_wrap(~broodYear, scales = "free") +
  labs(x = "Age", y = "Weight (g)", color = "Model",
       fill = "Model", title = "Female WAA Cohorts") +
  theme(legend.position = "top") 

# WAA curve over time (Males)
mod_outputs %>% 
  filter(sex == "Male", type == "weight", broodYear %in% c(2010:2016)) %>% 
  ggplot(aes(x = Age, y = Mean, ymin = lwr_95, ymax = upr_95)) +
  geom_ribbon(alpha = 0.3, aes(fill = re_model)) +
  geom_line(aes(color = re_model)) +
  geom_vline(aes(xintercept = 2018 - broodYear)) +
  facet_wrap(~broodYear, scales = "free") +
  labs(x = "Age", y = "Weight (g)", color = "Model",
       fill = "Model", title = "Male WAA Cohorts") +
  theme(legend.position = "top") 

mod_outputs %>% 
  filter(sex == "female", type == "length", broodYear %in% c(2010:2016)) %>% 
  ggplot(aes(x = Age, y = Mean, ymin = lwr_95, ymax = upr_95)) +
  geom_ribbon(alpha = 0.3, aes(fill = re_model)) +
  geom_line(aes(color = re_model)) +
  geom_vline(aes(xintercept = 2018 - broodYear)) +
  facet_wrap(~broodYear, scales = "free") +
  labs(x = "Age", y = "Weight (g)", color = "Model",
       fill = "Model", title = "Female LAA Cohorts") +
  theme(legend.position = "top") 

# LAA curve over time (Males)
mod_outputs %>% 
  filter(sex == "Male", type == "length", broodYear %in% c(2010:2016)) %>% 
  ggplot(aes(x = Age, y = Mean, ymin = lwr_95, ymax = upr_95)) +
  geom_ribbon(alpha = 0.3, aes(fill = re_model)) +
  # geom_label(aes(label = Year, fill = re_model)) +
  geom_line(aes(color = re_model), size = 1.3) +
  facet_wrap(~broodYear, scales = "free") +
  geom_vline(aes(xintercept = 2018 - broodYear)) +
  labs(x = "Age", y = "Weight (g)", color = "Model",
       fill = "Model", title = "Male LAA Cohorts") +
  theme(legend.position = "top") 

dev.off()


# Multi panel plot --------------------------------------------------------

# LAA curve over time (Males)
laa_male_cohort <- mod_outputs %>% 
  filter(sex == "Male", type == "length", broodYear %in% c(2010:2016)) %>% 
  ggplot(aes(x = Age, y = Mean, ymin = lwr_95, ymax = upr_95)) +
  geom_ribbon(alpha = 0.3, aes(fill = re_model)) +
  geom_line(aes(color = re_model), size = 1.3) +
  facet_wrap(~broodYear, scales = "free") +
  geom_vline(aes(xintercept = 2018 - broodYear)) +
  labs(x = "Age", y = "Length (cm)", color = "Model",
       fill = "Model", title = "Male LAA Cohorts") +
  theme(legend.position = "top") +
  theme_trevor_jordan()

female_laa_cohort <- mod_outputs %>% 
  filter(sex == "female", type == "length", broodYear %in% c(2010:2016)) %>% 
  ggplot(aes(x = Age, y = Mean, ymin = lwr_95, ymax = upr_95)) +
  geom_ribbon(alpha = 0.3, aes(fill = re_model)) +
  geom_line(aes(color = re_model), size = 1.3) +
  geom_vline(aes(xintercept = 2018 - broodYear)) +
  facet_wrap(~broodYear, scales = "free") +
  labs(x = "Age", y = "Length (cm)", color = "Model",
       fill = "Model", title = "Female LAA Cohorts") +
  theme(legend.position = "top") +
  theme_trevor_jordan()

waa_male_cohort <- mod_outputs %>% 
  filter(sex == "Male", type == "weight", broodYear %in% c(2010:2016)) %>% 
  ggplot(aes(x = Age, y = Mean, ymin = lwr_95, ymax = upr_95)) +
  geom_ribbon(alpha = 0.3, aes(fill = re_model)) +
  geom_line(aes(color = re_model), size = 1.3) +
  facet_wrap(~broodYear, scales = "free") +
  geom_vline(aes(xintercept = 2018 - broodYear)) +
  labs(x = "Age", y = "Weight (g)", color = "Model",
       fill = "Model", title = "Male WAA Cohorts") +
  theme(legend.position = "top") +
  theme_trevor_jordan()

female_waa_cohort <- mod_outputs %>% 
  filter(sex == "female", type == "weight", broodYear %in% c(2010:2016)) %>% 
  ggplot(aes(x = Age, y = Mean, ymin = lwr_95, ymax = upr_95)) +
  geom_ribbon(alpha = 0.3, aes(fill = re_model)) +
  geom_line(aes(color = re_model), size = 1.3) +
  geom_vline(aes(xintercept = 2018 - broodYear)) +
  facet_wrap(~broodYear, scales = "free") +
  labs(x = "Age", y = "Weight (cm)", color = "Model",
       fill = "Model", title = "Female WAA Cohorts") +
  theme(legend.position = "top") +
  theme_trevor_jordan()

png(here(dir.out, "Growth_Cohort_Panel.png"), width = 1500, height = 1300)
ggpubr::ggarrange(laa_male_cohort, female_laa_cohort, 
                  waa_male_cohort, female_waa_cohort)
dev.off()

png(here(dir.out, "female_waa_trend.png"), width = 1500, height = 1300)
mod_outputs %>% 
  filter(sex == "female", type == "weight") %>% 
  ggplot(aes(x = Year, y = Mean, ymin = lwr_95, ymax = upr_95)) +
  geom_ribbon(alpha = 0.3, aes(fill = re_model)) +
  geom_line(aes(color = re_model)) +
  facet_wrap(~Age, scales = "free") +
  labs(x = "Year", y = "Weight (g)", color = "Model",
       fill = "Model", title = "Female WAA") +
  theme(legend.position = "top") +
  theme_trevor_jordan()
# geom_vline(xintercept = 2022)
dev.off()


# Get Size-Age Transition -------------------------------------------------

age_bins = 2:31
len_bins = seq(41, 99, 2)
al_mat_all_df = data.frame()

for(i in 1:length(mod_list)) {
  if(str_detect(mod_list[[i]]$name, "length")) { # only do this for length modelsmmmm
    # Extract out quantities
    laa_sigma = exp(mod_list[[i]]$sd_rep$par.fixed[names(mod_list[[i]]$sd_rep$par.fixed) == "ln_obs_sigma2"])
    laa_mat = matrix(mod_list[[i]]$sd_rep$value[names(mod_list[[i]]$sd_rep$value) == "mu_at"],
                     nrow = length(unique(mod_df$Age)), ncol = length(unique(mod_df$Year)) + n_proj_years )
    
    # Set up array
    al_array = array(data = 0, dim = c(length(age_bins), length(len_bins), ncol(laa_mat)),
                     dimnames = list(c(age_bins), c(len_bins), c(sort(unique(age_dat$Year)), 2023:2032) ))

    
    for(t in 1:ncol(laa_mat)) { # loop through time to get matrix
      al_array[,,t] = get_al_trans_matrix(age_bins = 2:31, len_bins = seq(41, 99, 2), 
                                          mean_length = laa_mat[,t], sd = laa_sigma)
    } # end t loop
     
    # naming and mungingz
    al_df = reshape2::melt(al_array)
    names(al_df) = c("Age", "Length", "Year", "Value")
    al_df$sex = str_split(mod_list[[i]]$name, " ")[[1]][1] # differentiate sex
    al_df$re_model = str_split(mod_list[[i]]$name, " ")[[1]][3] # differentiate model
    
    al_mat_all_df = rbind(al_mat_all_df, al_df)
  } # end if
} # end t loop


pdf(here(dir.out, "size-age-transition.pdf"), width = 15, height = 13)
val_range = range(al_mat_all_df$Value) # get range of values
for(s in 1:length(sexes)) {
  for(r in 1:length(re_model)) {
    
    # filter for given elements
    plot_al_df = al_mat_all_df[al_mat_all_df$sex == sexes[s] &
                               al_mat_all_df$re_model == re_model[r] ,]

    print(
      ggplot(plot_al_df, aes(x = Age, y = Length, color = Value, size = Value, alpha = Value)) +
        geom_point() + 
        scale_color_viridis_c(limits = val_range) +
        facet_wrap(~Year) +
        labs(x = "Age", y = "Length (cm)", title = paste(sexes[s], re_model[r])) +
        theme_classic() +
        coord_cartesian(ylim = c(36, 99))
    )
    
    print(
      ggplot(plot_al_df, aes(x = Length, y = Value, color = Year, group = Year))+
        geom_line(size = 1.5, alpha = 0.75) + 
        scale_color_viridis_c() +
        facet_wrap(~Age) +
        labs(y = "P(Length|Age)", x = "Length (cm)", title = paste(sexes[s], re_model[r])) +
        theme_classic() +
        coord_cartesian(ylim = c(0, 0.3))
    )

  } # end r loop
} # end s loop
dev.off()


# Extract out estimated parameters
parameter_df = data.frame()
for(i in 1:length(mod_list)) {
  
  # Extract out names
  names = str_split(mod_list[[i]]$name, " ")[[1]]
  # extract out parameters here
  par_df = reshape2::melt(summary(mod_list[[i]]$sd_rep)) %>% 
    filter(!Var1 %in% c('mu_at', "ln_eps_at") ) %>% 
    pivot_wider(values_from = value, names_from = Var2) %>% 
    rename(Parameters = Var1, Std_Error = `Std. Error`) %>% 
    mutate(model = names[3], sex = names[1], type = names[2],
           lwr_95 = Estimate - 1.96*Std_Error, upr_95 = Estimate + 1.96*Std_Error)
  
  parameter_df = rbind(parameter_df, par_df)
    
} # end i loop

ggplot(parameter_df %>% 
         filter(Parameters %in% c("rho_a", "rho_y","rho_c")), 
       aes(x = sex, y = Estimate, ymin = lwr_95, ymax = upr_95, color = Parameters)) +
  facet_grid(model~type) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_pointrange(position = position_dodge(width = 0.2), size = 1) +
  theme_bw() +
  labs(x = "Model", y = "Estimate", title = 'Correlation')

ggplot(parameter_df %>% 
         filter(Parameters %in% c("ln_obs_sigma2")), 
       aes(x = model, y = Estimate, ymin = lwr_95, ymax = upr_95, color = sex)) +
  facet_wrap(~type, scales = "free") +
  geom_pointrange(position = position_dodge(width = 0.2), size = 1) +
  theme_bw() +
  labs(x = "Model", y = "Estimate", title = 'Observation Error')

