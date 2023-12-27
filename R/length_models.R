# Purpose: To run 3dar1 length-at-age models (including retrospective)
# Creator: Matthew LH. Cheng (UAF CFOS)
# Date: 12/15/23


# Set up ------------------------------------------------------------------

library(TMB)
library(TMBhelper)
library(here)
library(tidyverse)

# setwd("src")
TMB::compile("Growth_Model.cpp")
dyn.unload(dynlib("Growth_Model"))
dyn.load(dynlib('Growth_Model'))

source(here("R", "functions.R"))

dir.out = here("output", "LAA_Models")
dir.create(dir.out)

# Read in length-at-age data here
age_dat = read.csv(here("data", "age_view.csv"), check.names = FALSE) %>% 
  filter(!is.na(Age), !is.na(`Length (cm)`), !is.na(`Weight (g)`),
         `Sex Description` != "Unknown", `Error Flag` == 0,
         !Age %in% c(0, 1)) %>% 
  rename(Length = `Length (cm)`,
         Sex_name = `Sex Description`,
         Weight = `Weight (g)`) %>% 
  mutate(Age = ifelse(Age >= 31, 31, Age), # collapse ages to + group
         Year = Year - min(Year)) %>% 
  select(Length, Weight, Age, Year, Sex_name) %>% 
  mutate(Read_Age = Age) %>% 
  rename(Tester_Age = Age) # rename

age_bins = sort(unique(age_dat$Tester_Age)) # set up number of age bins
re_model = c("constant", "3dar1")
sex_names = unique(age_dat$Sex_name)
years = sort(unique(age_dat$Year))
n_retro = 15
n_proj = 1 # projection years

# Run Models --------------------------------------------------------------
sd_rep_all = data.frame() # store sdrep
growth_all = data.frame() # store growth estimatesl

# Order ot operations
# RE model loop
# Sex loop
# Year loop (retrospectives)

for(r in 1:length(re_model)) {
  for(s in 1:length(sex_names)) {
    for(y in 0:n_retro) {
      
      # Filter to model data
      mod_df = age_dat %>% filter(Sex_name == sex_names[s], # filter to sex
                         Year %in% c(0:(length(years)-y-1))) %>%  # peel back years
        mutate(Tester_Age = Tester_Age - min(Tester_Age)) %>% # subtract age for tmb (prob not necessary)
        select(-Sex_name) %>% 
        filter(Length >= 41)
      
      # set up model inputs
      model_inputs = prepare_data(data = mod_df, # model dataframe
                                  obs_sd_mat = matrix(0), # no observed sd for LAA models
                                  age_bins = age_bins, # age bins
                                  years = sort(unique(mod_df$Year)), # years to model
                                  growth_model_name = "length", # growth model
                                  re_model = re_model[r], # random effects model
                                  var_param = 0, # variance parameterization for 3dar1
                                  n_proj_years = n_proj) # number of projection years
      
      if(sex_names[s] == "Male") model_inputs$parameters$ln_X_inf = log(65)

      # set up AD object
      model = MakeADFun(data = model_inputs$data, 
                        parameters = model_inputs$parameters, 
                        DLL = "Growth_Model", random = c("ln_eps_at"), # specify what to integrate out
                        map = model_inputs$map, silent = TRUE)
      
      # optimize
      optim = stats::nlminb(model$par, model$fn, model$gr,  
                            control = list(iter.max = 5e5, eval.max = 5e5))
      
      add_newton(n.newton = 3, ad_model = model, mle_optim = optim) # add newton steps
      model$sd_rep = sdreport(model) # get standard errors
      model$optim = optim # optimizer information
      model$report = model$report() # get report
      
      # Assess model convergence here
      if(model$sd_rep$pdHess & max(abs(model$sd_rep$gradient.fixed)) <= 1e-3 &
         sum(is.nan( model$sd_rep$sd)) == 0) {
        model$Convergence = "Converged"
      } else {
        model$Convergence = "Not Converged"
      } # if else for model convergence
      
# Summarize and collate results -------------------------------------------
      
      # Extract out parameter estimates
      sd_rep_summary = reshape2::melt(summary(model$sd_rep))
      
      # some residual data munging
      sd_rep_summary = sd_rep_summary %>% filter(!Var1 %in% c("ln_eps_at", "mu_at")) %>% 
        rename(Parameter = Var1, Type = Var2) %>% 
        mutate(peel = y, Convergence = model$Convergence,
               re_Model = re_model[r], Sex = sex_names[s],
               Growth_model = "length") %>% 
        mutate(Parameter = ifelse(Parameter == "ln_obs_sigma2",
                                         paste("ln_obs_sigma2", c("min", "max")),
                                         as.character(Parameter))) %>% 
        pivot_wider(names_from = "Type", values_from = "value") %>% 
        # Calculating uncertainty here
        mutate(lwr_95 = case_when(
          str_detect(Parameter, "ln") ~ exp(Estimate - (1.96*`Std. Error`)),
          !str_detect(Parameter, "ln") ~ Estimate - (1.96*`Std. Error`)
        ),
        upr_95 = case_when(
          str_detect(Parameter, "ln") ~ exp(Estimate + (1.96*`Std. Error`)),
          !str_detect(Parameter, "ln") ~ Estimate + (1.96*`Std. Error`)
        ),
        Trans_Estimate = case_when(
          str_detect(Parameter, "ln") ~ exp(Estimate),
          !str_detect(Parameter, "ln") ~ Estimate # whether to log transform
        ))
      
      # extract out estimates
      growth = reshape2::melt(matrix(model$sd_rep$value[names(model$sd_rep$value) == "mu_at"], 
                                     nrow = length(unique(age_bins)),
                                     ncol = length(unique(mod_df$Year)) + n_proj )) # save mean estimates
      
      growth_sd = reshape2::melt(matrix(model$sd_rep$sd[names(model$sd_rep$value) == "mu_at"], 
                                        nrow = length(unique(age_bins)),
                                        ncol = length(unique(mod_df$Year)) + n_proj  )) # standard errors from delta method
      
      # save estimates, etc.
      growth_df = growth %>% 
        left_join(growth_sd, by = c("Var1", "Var2")) %>% 
        rename(Age = Var1, Year = Var2,
               Mean = value.x, SD = value.y) %>% 
        mutate(Age = Age + 1, Year = Year + 1995,
               broodYear = Year-Age, peel = y,
               Convergence = model$Convergence,
               Growth_Model = "length", Sex = sex_names[s],
               re_model = re_model[r], # name model elements
               lwr_95 = Mean - 1.96*SD,
               upr_95 = Mean + 1.96*SD) 
      
      growth_all = rbind(growth_all, growth_df)
      sd_rep_all = rbind(sd_rep_summary, sd_rep_all)
      if(y == 0 & re_model[r] == "3dar1") {
        save(model, file =  here(dir.out, paste("length", sex_names[s], "3DAR1_Model.RData", sep = "_"))) # save model
      } # save model for terminal
      print(paste("Peel", y))
    } # end y loop
  } # end s loop
} # end r loop

write.csv(growth_all, here(dir.out, "Growth_estimates.csv"))
write.csv(sd_rep_all, here(dir.out, "SDRep_estimates.csv"))

