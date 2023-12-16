library(TMB)
library(TMBhelper)
library(here)
library(tidyverse)

setwd("src")
TMB::compile("Growth_Model.cpp")
dyn.unload(dynlib('Growth_Model'))
dyn.load(dynlib('Growth_Model'))
source(here("R", "functions.R"))

# Do some quick data munging and load in data
age_dat = read.csv(here("data", "age_view.csv"), check.names = FALSE) %>% 
  filter(!is.na(Age), !is.na(`Length (cm)`), !is.na(`Weight (g)`),
         `Sex Description` != "Unknown", `Error Flag` == 0)  # only keep age samples

n_sexes = 2
sexes = c("female", "Male")
mod_name = c("2dar1", "3dar1", "iid")
plot_df = data.frame()
mod_list = list()

# filter to model data
mod_age_dat = age_dat %>% 
  rename(Length = `Length (cm)`,
         Sex_name = `Sex Description`,
         Weight = `Weight (g)`) %>% 
  filter(Sex_name == sexes[s]) %>% 
  select(Length, Age, Year, Weight) %>% 
  filter(!Age %in% c(0, 1)) %>% 
  mutate(Age = ifelse(Age >= 31, 31, Age),
         Age = Age - min(Age),
         Year = Year - min(Year))

  
  # get age year index
  ay_Index <- as.matrix(expand.grid("age" = seq_len(length(unique(mod_age_dat$Age))), 
                                    "year" = seq_len(length(unique(mod_age_dat$Year)) ) ))
  
  parameters = list(ln_X_inf = log(5),
                    ln_k = log(0.25),
                    t0 = -0.2,
                    ln_alpha = log(1e-05),
                    ln_beta = log(3.02),
                    ln_obs_sigma2 = log(0.1),
                    ln_eps_sigma2 = log(0.05),
                    rho_a = 0.3,
                    rho_y = 0.3,
                    rho_c = 0.3,
                    ln_eps_at = array(data = log(0.01), dim = c(length(unique(mod_age_dat$Age)),
                                                                length(unique(mod_age_dat$Year)) ) ) )
  
  # set up data
  data = list(obs_mat = as.matrix(mod_age_dat),
              ages = as.vector(sort(unique(mod_age_dat$Age)) + 2),
              ay_Index = ay_Index, var_param = 0, dev_param = 3, growth_model = 2)
  map = list(ln_beta = factor(NA), ln_alpha = factor(NA))
  
  
  model = MakeADFun(data = data, parameters = parameters, 
                    DLL = "Growth_Model", random = "ln_eps_at",
                    map = map, silent = FALSE)
  
  optim = stats::nlminb(model$par, model$fn, model$gr,  
                        control = list(iter.max = 1e6, eval.max = 1e6))
  
  add_newton(n.newton = 3, ad_model = model, mle_optim = optim)
  report = model$report()
  model_sd_rep <- sdreport(model)
  
  # save data
  waa = reshape2::melt(matrix(model_sd_rep$value[names(model_sd_rep$value) == "mu_at"], 
                              nrow = length(unique(mod_age_dat$Age)),
                              ncol = length(unique(mod_age_dat$Year)) ))
  
  sd = reshape2::melt(matrix(model_sd_rep$sd[names(model_sd_rep$value) == "mu_at"], 
                             nrow = length(unique(mod_age_dat$Age)),
                             ncol = length(unique(mod_age_dat$Year)) ))
  
  waa_df = waa %>% 
    left_join(sd, by = c("Var1", "Var2")) %>% 
    mutate(lwr_95 = value.x - 1.96*value.y,
           upr_95 = value.x + 1.96*value.y,
           sex = sexes[s], model = mod_name[mod] ) %>% 
    rename(Age = Var1, Year = Var2) %>% 
    mutate(Age = Age + 1, Year = Year + 1995, broodYear = Year-Age,
           AIC = TMBhelper::TMBAIC(optim)) 
  
  ggplot(waa_df, aes(x = Year, y = exp(value.x))) +
    geom_point(mod_age_dat_all %>% filter(sex == "Male"),
               mapping = aes(x = Year, y = Weight), alpha = 0.15, color = "black",
               size = 1) +
    geom_line(aes(color = factor(model)), size = 1.5) +
    geom_ribbon(alpha = 0.25, aes(ymin = lwr_95, ymax = upr_95, fill = factor(model))) +
    geom_hline(aes(yintercept = mu), lty = 2) +
    facet_wrap(~Age, scales = 'free') +
    labs(x = "Year", y = 'Weight at age (g)', color = "Model", fill = "Model",
         title = "Males") +
    theme(legend.position = "top")
  
  

for(s in 1:n_sexes) {
  
  # filter to model data
  mod_age_dat = age_dat %>% 
    rename(Length = `Length (cm)`,
           Sex_name = `Sex Description`,
           Weight = `Weight (g)`) %>% 
    filter(Sex_name == sexes[s]) %>% 
    select(Length, Age, Year, Weight) %>% 
    filter(!Age %in% c(0, 1)) %>% 
    mutate(Age = ifelse(Age >= 31, 31, Age),
           Age = Age - min(Age),
           Year = Year - min(Year))
  
  
  # get age year index
  ay_Index <- as.matrix(expand.grid("age" = seq_len(length(unique(mod_age_dat$Age))), 
                                    "year" = seq_len(length(unique(mod_age_dat$Year)) ) ))
  
  
  parameters = list(ln_Linf = log(80),
                    ln_k = log(0.25),
                    ln_L0 = log(35),
                    ln_alpha = log(1e-05),
                    ln_beta = log(3),
                    ln_waa_sigma2 = log(0.1),
                    ln_laa_sigma2 = log(5),
                    ln_eps_sigma2 = log(0.05),
                    rho_a = 0.3,
                    rho_y = 0.3,
                    rho_c = 0.3,
                    ln_eps_at = array(data = log(0.01), dim = c(length(unique(mod_age_dat$Age)),
                                                                length(unique(mod_age_dat$Year)) ) ) )
  
  for(mod in 1:length(mod_name)) {
    
    if(mod_name[mod] == "3dar1") {
      # set up data
      data = list(obs_mat = as.matrix(mod_age_dat),
                  ages = as.vector(sort(unique(mod_age_dat$Age)) + 2),
                  ay_Index = ay_Index, var_param = 0, dev_param = 3)
      map = list(ln_beta = factor(NA))
    }
    if(mod_name[mod] == "2dar1") {
      # set up data
      data = list(obs_mat = as.matrix(mod_age_dat),
                  ages = as.vector(sort(unique(mod_age_dat$Age)) + 2),
                  ay_Index = ay_Index, var_param = 0, dev_param = 2)
      map = list(ln_beta = factor(NA),
                 rho_c = factor(NA))
    }
    
    if(mod_name[mod] == "iid") {
      # set up data
      data = list(obs_mat = as.matrix(mod_age_dat),
                  ages = as.vector(sort(unique(mod_age_dat$Age)) + 2),
                  ay_Index = ay_Index, var_param = 0, dev_param = 0)
      map = list(ln_beta = factor(NA),
                 rho_c = factor(NA),
                 rho_a = factor(NA),
                 rho_y = factor(NA))
    }
    
    model = MakeADFun(data = data, parameters = parameters, 
                       DLL = "Growth_Model", random = "ln_eps_at",
                       map = map, silent = FALSE)
    
    optim = stats::nlminb(model$par, model$fn, model$gr,  
                          control = list(iter.max = 1e6, eval.max = 1e6))
    
    add_newton(n.newton = 3, ad_model = model, mle_optim = optim)
    report = model$report()
    model_sd_rep <- sdreport(model)

    # save data
    waa = reshape2::melt(matrix(model_sd_rep$value[names(model_sd_rep$value) == "WAA_at"], 
                                     nrow = length(unique(mod_age_dat$Age)),
                                     ncol = length(unique(mod_age_dat$Year)) ))
    
    sd = reshape2::melt(matrix(model_sd_rep$sd[names(model_sd_rep$value) == "WAA_at"], 
                                    nrow = length(unique(mod_age_dat$Age)),
                                    ncol = length(unique(mod_age_dat$Year)) ))
    
    waa_df = waa %>% 
      left_join(sd, by = c("Var1", "Var2")) %>% 
      mutate(lwr_95 = value.x - 1.96*value.y,
             upr_95 = value.x + 1.96*value.y,
             sex = sexes[s], model = mod_name[mod] ) %>% 
      rename(Age = Var1, Year = Var2) %>% 
      mutate(Age = Age + 1, Year = Year + 1995, broodYear = Year-Age,
             AIC = TMBhelper::TMBAIC(optim)) 
    
    plot_df = rbind(waa_df, plot_df)
    
  } # end mod loop
  
  
} # end sex loop

# Get assessment weight at age relationships
assess = data.frame(Age = 2:31,
                    female = exp(log(5.87) + 3.02 * log(1 - exp(-0.17 * (2:31 + 2.98)))) ,
                    Male = exp(log(3.22) + 3.02 * log(1 - exp(-0.27 * (2:31 + 2.41)))) ) %>% 
  pivot_longer(!Age, names_to = "sex", values_to = "mu")

# all raw data
mod_age_dat_all = age_dat %>% 
  rename(Length = `Length (cm)`,
         sex = `Sex Description`,
         Weight = `Weight (g)`) %>% 
  select(Length, Age, Year, Weight, sex) %>% 
  filter(!Age %in% c(0, 1)) %>% 
  mutate(Age = ifelse(Age >= 31, 31, Age))

  ggplot(plot_df %>% filter(sex == "Male") %>% 
           left_join(assess, by = c("sex", "Age")), 
         aes(x = Year, y = value.x)) +
    geom_point(mod_age_dat_all %>% filter(sex == "Male"),
               mapping = aes(x = Year, y = Weight), alpha = 0.15, color = "black",
               size = 1) +
  geom_line(aes(color = factor(model)), size = 1.5) +
  geom_ribbon(alpha = 0.25, aes(ymin = lwr_95, ymax = upr_95, fill = factor(model))) +
  geom_hline(aes(yintercept = mu), lty = 2) +
  facet_wrap(~Age, scales = 'free') +
  labs(x = "Year", y = 'Weight at age (g)', color = "Model", fill = "Model",
       title = "Males") +
    theme(legend.position = "top")
  
  ggplot(plot_df %>% filter(sex == "Male") %>% 
           left_join(assess, by = c("sex", "Age")), aes(x = Year, y = value.x)) +
    geom_line(aes(color = factor(model)), size = 1.5) +
    geom_ribbon(alpha = 0.25, aes(ymin = lwr_95, ymax = upr_95, fill = factor(model))) +
    geom_hline(aes(yintercept = mu), lty = 2) +
    facet_wrap(~Age, scales = 'free') +
    labs(x = "Year", y = 'Weight at age (g)', color = "Model", fill = "Model",
         title = "Males") +
    theme(legend.position = "top")
  
  ggplot(plot_df %>% filter(sex == "female") %>% 
           left_join(assess, by = c("sex", "Age")), 
         aes(x = Year, y = value.x)) +
    geom_line(aes(color = factor(model)), size = 1.5) +
    geom_ribbon(alpha = 0.25, aes(ymin = lwr_95, ymax = upr_95, fill = factor(model))) +
    geom_hline(aes(yintercept = mu), lty = 2) +
    facet_wrap(~Age, scales = 'free') +
    labs(x = "Year", y = 'Weight at age (g)', color = "Model", fill = "Model",
         title = "Males") +
    theme_bw() +
    theme(legend.position = "top") 
  
  ggplot(plot_df %>% filter(sex == "female") %>% 
           left_join(assess, by = c("sex", "Age")), 
         aes(x = Year, y = value.x)) +
    geom_point(mod_age_dat_all %>% filter(sex == "female"),
               mapping = aes(x = Year, y = Weight), alpha = 0.15, color = "black",
               size = 1) +
    geom_line(aes(color = factor(model)), size = 1.5) +
    geom_ribbon(alpha = 0.25, aes(ymin = lwr_95, ymax = upr_95, fill = factor(model))) +
    geom_hline(aes(yintercept = mu), lty = 2) +
    facet_wrap(~Age, scales = 'free') +
    labs(x = "Year", y = 'Weight at age (g)', color = "Model", fill = "Model",
         title = "Males") +
    theme_bw() +
    theme(legend.position = "top") 
  
  
  ggplot(plot_df %>% filter(sex == "Male") %>% 
           left_join(assess, by = c("sex", "Age")), 
         aes(x = Age, y = value.x, ymin = lwr_95, 
             ymax = upr_95, fill = factor(model))) +
    geom_line(aes(color = factor(model)), size = 1.3) +
    geom_line(aes(y = mu), lty = 2) +
    geom_ribbon(alpha = 0.25) +
    facet_wrap(~Year, scales = 'free') +
    labs(x = "Age", y = 'Weight at age (g)', color = "Model", fill = "Model",
         title = "Males") +
    theme(legend.position = "top")
  
  ggplot(plot_df %>% filter(sex == "female") %>% 
           left_join(assess, by = c("sex", "Age")), 
         aes(x = Age, y = value.x, ymin = lwr_95, 
             ymax = upr_95, fill = factor(model))) +
    geom_line(aes(color = factor(model)), size = 1.3) +
    geom_line(aes(y = mu), lty = 2) +
    geom_ribbon(alpha = 0.25) +
    facet_wrap(~Year, scales = 'free') +
    labs(x = "Age", y = 'Weight at age (g)', color = "Model", fill = "Model",
         title = "Females") +
    theme(legend.position = "top")
  
  aic = plot_df %>% 
    group_by(model, sex) %>% 
    select(AIC, model, sex) %>% 
    unique()

  

# Get Size age transition -------------------------------------------------

  
  laa_sigma = exp(model_sd_rep$par.fixed[names(model_sd_rep$par.fixed) == "ln_laa_sigma2"])
  laa_mat = matrix(model_sd_rep$value[names(model_sd_rep$value) == "LAA_at"], 
                   nrow = length(unique(mod_age_dat$Age)),
                   ncol = length(unique(mod_age_dat$Year)) )
  
  al_array = array(data = 0, dim = c(30, 30, ncol(laa_mat)))
  
for(i in 1:ncol(laa_mat)) {
  al_array[,,i] = get_al_trans_matrix(age_bins = 2:31, len_bins = seq(41, 101, 2), 
                      mean_length = laa_mat[,i], sd = laa_sigma)
}

al_df_male = reshape2::melt(al_array)
names(al_df_male) = c("Age", "Length", "Year", "Value")
al_df_male$sex = "male"
  
al_df_female = reshape2::melt(al_array)
names(al_df_female) = c("Age", "Length", "Year", "Value")
al_df_female$sex = "female"

al_df = rbind(al_df_female, al_df_male)

ggplot(al_df %>% filter(Year %in% c(23:27)), 
       aes(x = Age, y = Length, color = Value, size = Value, alpha = Value)) +
  geom_point() + 
  scale_color_viridis_c() +
  facet_grid(Year~sex)
