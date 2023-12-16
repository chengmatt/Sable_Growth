# Purpose: To run 3dar1 growth models (exploring how to incorporate age-ing error)
# Creator: Matthew LH. Cheng (UAF-CFOS)
# Date 11/5/23

library(TMB)
library(TMBhelper)
library(here)
library(tidyverse)
library(truncnorm)
library(doSNOW)
library(parallel)

setwd("src")
TMB::compile("Growth_Model.cpp")
dyn.unload(dynlib("Growth_Model"))
dyn.load(dynlib('Growth_Model'))
source(here("R", "functions.R"))

# set up cores
ncores <- detectCores() 
# Register cluster here
cl <- makeCluster(ncores - 2)
registerDoSNOW(cl)

# Read in ageing error matrix
ageing_error_df = readLines(here("data", "SQ_Nominal.dat")) # read in dat file
ageage_mat = ageing_error_df[[which(str_detect(ageing_error_df, "#Age age transition matrix")) + 1]] # read in age-age error matrix
ageage_mat = matrix(as.numeric(str_split(ageage_mat, pattern = " ")[[1]] ), ncol = 30, nrow = 30, dimnames = list(c(2:31), c(2:31))) # coerce strings into matrix
ewaa = read.csv(here("output", "ewaa.csv"))

# Do some quick data munging and load in data
age_dat = read.csv(here("data", "age_view.csv"), check.names = FALSE) %>% 
  filter(!is.na(Age), !is.na(`Length (cm)`), !is.na(`Weight (g)`),
         `Sex Description` != "Unknown", `Error Flag` == 0) %>% 
  rename(Length = `Length (cm)`,
         Sex_name = `Sex Description`,
         Weight = `Weight (g)`) 

n_iter = 100 # number of iterations to run
age_bins = 2:31 # set up number of age bins
sexes = unique(age_dat$Sex_name) # sex names
grwth_model = c("length", "weight") # length and weight growth models
n_proj_years = 1 # number of projection years
growth_all = data.frame() # empty dataframe to store everything
sd_rep_all = data.frame() # empty dataframe to store all sdreps
sexes = c("female", "Male")
# Run Models (incorporating aging variability) --------------------------------------------------------------

for(s in 1:length(sexes)) { # start sex loop
  
  # Filter to model data
  mod_df = age_dat %>% 
    filter(Sex_name == sexes[s]) %>%  # filter to a given sex
    filter(!Age %in% c(0, 1)) %>% 
    mutate(Age = ifelse(Age >= 31, 31, Age), # collapse ages to + group
           Year = Year - min(Year)) %>% 
    select(Length, Weight, Age, Year) %>% 
    mutate(Read_Age = Age) %>% 
    rename(Impute_Age = Age) # rename
  
  for(g in 1:length(grwth_model)) { # start growth loop
    
    dir.out = here("output", "Ageing_Error", paste(grwth_model[g], sexes[s], "3DAR1_Model", sep = "_")) # directory
    dir.create(dir.out) # create directory
    
    run_models <- foreach(iter = 1:n_iter, .packages = c("TMB", "here", "tidyverse")) %dopar% {
      # Compile TMB
      TMB::compile("Growth_Model.cpp")
      dyn.load(dynlib('Growth_Model'))
      
      # Impute ageing error/within age reader variability
      for(a in 1:nrow(mod_df)) {
        Read_Age = mod_df$Read_Age[a] # extract out read age
        Impute_Age = reshape2::melt(rmultinom(1,1,ageage_mat[,Read_Age-1])) # multinomial draw from ageing error to imput
        mod_df$Impute_Age[a] <- Impute_Age$Var1[Impute_Age$value == 1] # impute the age in
      } # end a loop
      
      # Subtract age for indexing in TMB
      mod_df <- mod_df %>% mutate(Impute_Age = Impute_Age - min(Impute_Age))
      
      # set up model inputs
      model_inputs = prepare_data(data = mod_df, # model dataframe
                                  age_bins = age_bins,
                                  growth_model_name = grwth_model[g], # growth model
                                  re_model = "3dar1", # random effects model
                                  var_param = 0, # variance parameterization for 3dar1
                                  age_offset = 2, # age offset for calcualting growth
                                  n_proj_years = n_proj_years) # number of projection years
      
      # set up AD object
      model = MakeADFun(data = model_inputs$data, 
                        parameters = model_inputs$parameters, 
                        DLL = "Growth_Model", random = c("ln_eps_at"),
                        map = model_inputs$map, silent = FALSE)
      
      # optimize
      optim = stats::nlminb(model$par, model$fn, model$gr,  
                            control = list(iter.max = 1e5, eval.max = 1e5))
      
      add_newton(n.newton = 4, ad_model = model, mle_optim = optim) # add newton steps
      model$optim = optim # optimizer information
      model$report = model$report() # get report
      model$sd_rep = sdreport(model) # get standard errors
      # Figure out model convergence
      if(model$sd_rep$pdHess & max(abs(model$sd_rep$gradient.fixed)) <= 1e-3 &
         sum(is.nan( model$sd_rep$sd)) == 0) model$Convergence = "Converged"
      else model$Convergence = "Not Converged"
      model # save this as our output for the parallelized loop
    } # end foreach loop
    
    for(iter2 in 1:n_iter) {
      
      # Extract out parameter estimates
      sd_rep_summary = reshape2::melt(summary(run_models[[iter2]]$sd_rep))
      
      # some residual data munging
      sd_rep_summary = sd_rep_summary %>% filter(!Var1 %in% c("ln_eps_at", "mu_at")) %>% 
        pivot_wider(names_from = "Var2") %>% 
        rename(Parameter = Var1, SE = `Std. Error`) %>% 
        mutate(iter = iter2, Convergence = run_models[[iter2]]$Convergence,
               Growth_Model = grwth_model[g], Sex = sexes[s],
               Transformed_Estimate = case_when(
                 str_detect(Parameter, "ln_") ~ exp(Estimate), 
                 !str_detect(Parameter, "ln_") ~ Estimate),
               lwr_95 = case_when(
                 str_detect(Parameter, "ln_") ~ exp(Estimate - 1.96 * SE), 
                 !str_detect(Parameter, "ln_") ~ Estimate - 1.96 * SE),
               upr_95 = case_when(
                 str_detect(Parameter, "ln_") ~ exp(Estimate + 1.96 * SE), 
                 !str_detect(Parameter, "ln_") ~ Estimate + 1.96 * SE))
      
      sd_rep_all = rbind(sd_rep_all, sd_rep_summary)
      
      # extract out estimates
      growth = reshape2::melt(matrix(run_models[[iter2]]$sd_rep$value[names(run_models[[iter2]]$sd_rep$value) == "mu_at"], 
                                     nrow = length(unique(mod_df$Read_Age)),
                                     ncol = length(unique(mod_df$Year)) + n_proj_years )) # save mean estimates
      
      growth_sd = reshape2::melt(matrix(run_models[[iter2]]$sd_rep$sd[names(run_models[[iter2]]$sd_rep$value) == "mu_at"], 
                                        nrow = length(unique(mod_df$Read_Age)),
                                        ncol = length(unique(mod_df$Year)) + n_proj_years )) # standard errors from delta method
      
      # save estimates, etc.
      growth_df = growth %>% 
        left_join(growth_sd, by = c("Var1", "Var2")) %>% 
        rename(Age = Var1, Year = Var2,
               Mean = value.x, SD = value.y) %>% 
        mutate(Age = Age + 1, Year = Year + 1995,
               broodYear = Year-Age, iter = iter2,
               Convergence = run_models[[iter2]]$Convergence,
               Growth_Model = grwth_model[g], Sex = sexes[s]) # name model elements
      
      if(grwth_model[g] == "length") {
        growth_df = growth_df %>% 
          mutate(lwr_95 = Mean - 1.96*SD,
                 upr_95 = Mean + 1.96*SD)
      } # length-at-age
      
      if(grwth_model[g] == "weight") {
        growth_df = growth_df %>% 
          mutate(lwr_95 = exp(Mean - 1.96*SD),
                 upr_95 = exp(Mean + 1.96*SD),
                 Mean = exp(Mean))
      } # weight-at-age
      growth_all = rbind(growth_all, growth_df)
    } # end second iteration loop to output stuff
    
    save(run_models, file =  here(dir.out, paste(grwth_model[g], sexes[s], "3DAR1_Model.RData", sep = "_"))) # save model
    print(paste("Done with growth model", grwth_model[g]))
    
  } # end g loop
  print(paste("Done with sex", sexes[s]))
} # end s loop

write.csv(growth_all, here("output", "Growth_All_3DAR1_Models.csv"))
write.csv(sd_rep_all, here("output", "SDREP_3DAR1_Models.csv"))


# Run Models (without aging variability) ----------------------------------
ewaa = read.csv(here("output", "ewaa.csv"))

growth_noage = data.frame() # empty dataframe to store everything (no aging error)
sd_rep_no_age = data.frame() # empty dataframe to store all sdreps (no aging error)

for(s in 1:length(sexes)) { # start sex loop
  
  # Filter to model data
  mod_df = age_dat %>% 
    filter(Sex_name == sexes[s]) %>%  # filter to a given sex
    filter(!Age %in% c(0, 1)) %>% 
    mutate(Age = ifelse(Age >= 31, 31, Age), # collapse ages to + group
           Year = Year - min(Year)) %>% 
    select(Length, Weight, Age, Year) %>% 
    mutate(Read_Age = Age) %>% 
    rename(Impute_Age = Age) %>%  # rename
    mutate(Impute_Age = Impute_Age - min(Impute_Age))
  
  for(g in 1:length(grwth_model)) { # start growth loop
    
    dir.out = here("output", "No_Ageing", paste(grwth_model[g], sexes[s], "3DAR1_Model", sep = "_")) # directory
    dir.create(dir.out) # create directory
    
    TMB::compile("Growth_Model.cpp")
    dyn.unload(dynlib("Growth_Model"))
    dyn.load(dynlib('Growth_Model'))

    # set up model inputs
    model_inputs = prepare_data(data = mod_df, # model dataframe
                                age_bins = age_bins,
                                growth_model_name = grwth_model[g], # growth model
                                re_model = "3dar1", # random effects model
                                var_param = 0, # variance parameterization for 3dar1
                                n_proj_years = n_proj_years) # number of projection years
    
    if(s == 2) model_inputs$parameters$ln_X_inf = log(60)

    model_inputs$parameters$ln_obs_sigma2 = rep(0, 30)
    
    log_sd = ewaa %>% filter(Sex_name == sexes[s]) %>%
      select(Year, Tester_Age, log_sd) %>%
      pivot_wider(names_from = "Year", values_from = "log_sd")
    log_sd = cbind(log_sd[,as.character(0:26)], `27` = NA)
    log_sd[log_sd == 0] = NA
    model_inputs$data$obs_sd_mat = as.matrix(log_sd)
    image(t(as.matrix(log_sd)))
    model_inputs$map$ln_obs_sigma2 = factor(rep(NA, 30))

    mean_waa = ewaa %>% filter(Sex_name == sexes[s]) %>%
      select(Year, Tester_Age, mean) %>%
      pivot_wider(names_from = "Year", values_from = "mean")
    mean_waa = cbind(mean_waa[,as.character(0:26)], `27` = NA)
    model_inputs$data$obs_mat = as.matrix(mean_waa)
    image(t(as.matrix(mean_waa)))

    # set up AD object
    model = MakeADFun(data = model_inputs$data, 
                      parameters = model_inputs$parameters, 
                      DLL = "Growth_Model", random = c("ln_eps_at"),
                      map = model_inputs$map, silent = FALSE)
    
    # optimize
    optim = stats::nlminb(model$par, model$fn, model$gr,  
                          control = list(iter.max = 5e5, eval.max = 5e5))
    
    add_newton(n.newton = 4, ad_model = model, mle_optim = optim) # add newton steps
    model$sd_rep = sdreport(model) # get standard errors
    model$optim = optim # optimizer information
    model$report = model$report() # get report
    
    # Figure out model convergence
    if(model$sd_rep$pdHess & max(abs(model$sd_rep$gradient.fixed)) <= 1e-3 &
       sum(is.nan( model$sd_rep$sd)) == 0) {model$Convergence = "Converged"} else model$Convergence = "Not Converged"
      
    # Extract out parameter estimates
    sd_rep_summary = reshape2::melt(summary(model$sd_rep))

    # some residual data munging
    sd_rep_summary = sd_rep_summary %>% filter(!Var1 %in% c("ln_eps_at", "mu_at")) %>% 
      rename(Parameter = Var1, Type = Var2) %>% 
      mutate(iter = 1, Convergence = model$Convergence,
             Growth_Model = grwth_model[g], Sex = sexes[s])
    
    sd_rep_no_age = rbind(sd_rep_no_age, sd_rep_summary)
      
    # Filter to model data
    mod_df = age_dat %>% 
      filter(Sex_name == sexes[s]) %>%  # filter to a given sex
      filter(!Age %in% c(0, 1)) %>% 
      mutate(Age = ifelse(Age >= 31, 31, Age), # collapse ages to + group
             Year = Year - min(Year)) %>% 
      select(Length, Weight, Age, Year) %>% 
      mutate(Read_Age = Age) %>% 
      rename(Impute_Age = Age) %>%  # rename
      mutate(Impute_Age = Impute_Age - min(Impute_Age))
    
    # extract out estimates
    growth = reshape2::melt(matrix(model$sd_rep$value[names(model$sd_rep$value) == "mu_at"], 
                                   nrow = length(unique(mod_df$Read_Age)),
                                   ncol = length(unique(mod_df$Year)) + n_proj_years )) # save mean estimates
    
    growth_sd = reshape2::melt(matrix(model$sd_rep$sd[names(model$sd_rep$value) == "mu_at"], 
                                      nrow = length(unique(mod_df$Read_Age)),
                                      ncol = length(unique(mod_df$Year)) + n_proj_years  )) # standard errors from delta method
    
    # save estimates, etc.
    growth_df = growth %>% 
      left_join(growth_sd, by = c("Var1", "Var2")) %>% 
      rename(Age = Var1, Year = Var2,
             Mean = value.x, SD = value.y) %>% 
      mutate(Age = Age + 1, Year = Year + 1995,
             broodYear = Year-Age, iter = 1,
             Convergence = model$Convergence,
             Growth_Model = grwth_model[g], Sex = sexes[s]) # name model elements
    
      growth_df = growth_df %>% 
        mutate(lwr_95 = Mean - 1.96*SD,
               upr_95 = Mean + 1.96*SD)
      
      growth_df %>% 
        filter(Year != 2023) %>% 
        ggplot(aes(x = Year, y = Mean, color = factor(Age))) +
        geom_line() 
      
      growth_df %>% 
        filter(broodYear %in% c(2000:2019), Year != 2023) %>% 
        ggplot(aes(x = Age, y = Mean, color = factor(broodYear))) +
        geom_line()  +
        geom_line(const %>%  filter(broodYear %in% c(2000:2019), Year != 2023), 
                  mapping = aes(x = Age, y = Mean), color = "black") +
        facet_wrap(~broodYear, scales = "free")
      
      ggplot() +
        geom_line(const, mapping = aes(x = Age, y = Mean), color = "black")
      
      growth_df %>%
        ggplot(aes(x = Age, y = Mean, color = factor(Year))) +
        geom_line()

  
  growth_noage = rbind(growth_noage, growth_df)
  save(model, file =  here(dir.out, paste(grwth_model[g], sexes[s], "3DAR1_Model.RData", sep = "_"))) # save model
  print(paste("Done with growth model", grwth_model[g]))
    
  } # end g loop
  print(paste("Done with sex", sexes[s]))
} # end s loop

write.csv(growth_noage, here("output", "Growth_NoAge_3DAR1_Models.csv"))
write.csv(sd_rep_no_age, here("output", "SDREP_NoAge_3DAR1_Models.csv"))
