# Plot updated admb assessment outputs and associated diagnostics
# Creator: Matthew LH. Cheng (UAF-CFOS)
# Date: 11/21/23

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

par_coeff <- function(df, coeff_ex, type) {
  
  # Extract them out!
  coeff <- df$coefficients[names(df$coefficients) %in% c(coeff_ex)]
  se <- df$se[names(df$coefficients) %in% c(coeff_ex)]
  
  par_df <- data.frame(coeff = coeff, se = se,
                       upr = (coeff + (1.96*se)),
                       downr = (coeff - (1.96*se)),
                       type = type,
                       coeff_type = names(coeff))
  
  return(par_df)
}

#' Title Get OSA residuals
#'
#' @param obs observed matrix (n_years, n_bins)
#' @param pred predicted matrix (n_years, n_bins)
#' @param iss Input sample size
#' @param iter_wt Iterative weights (if any)
#' @param index Index (age bins or length bins)
#' @param drop_bin Which age or length bin to drop
#'
#' @return
#' @export
#'
#' @examples
get_osa_res = function(obs, pred, iss, iter_wt, index, drop_bin) {
  require(compResidual)
  # Do some munging
  years = rownames(pred) # get years
  ess_wt = round(iss * iter_wt, 0) # get ess
  obs1 = round((obs + 0.001) * (ess_wt)) # calculate observations in numbers
  pred = pred + 0.001 # add constant so we can divide finite
  # determine which bin to drop (software drops last bin automatically)
  obs1 <- cbind(obs1[,-drop_bin], obs1[,drop_bin])
  pred1 <- cbind(pred[,-drop_bin], pred[,drop_bin])
  
  # Get OSA here
  osa_res <- NULL
  while (is.null(osa_res) || !is.finite(sum(osa_res))) {
    osa_res <- resMulti(t(obs1), t(pred1)) 
  } # caluclate residual until it does not return Nans
  
  mat <- t(matrix(osa_res, nrow=nrow(osa_res), ncol=ncol(osa_res))) # munge into matrix
  dimnames(mat) <- list(year=years, index=index) # input dimensions 
  reslong <- reshape2::melt(mat, value.name='resid') # reshape
  return(list(reslong))
}

#' Title Plot residuals
#'
#' @param data Dataframe with columns year, index, resid
#' @param comp_type Composition type 
#' @param res_type Residual type
#' @param ymin lower ylimit on axis
#'
#' @return
#' @export
#'
#' @examples
res_plot = function(data, comp_type, res_type, ymin) {
  require(ggplot2)
  ggplot(data, aes(year, index, size=abs(resid), color=resid>0)) + geom_point(alpha = 0.75) +
    labs(x = "Year", y = comp_type, title = res_type, 
         size = "Absolute Residual", color = "Residual > 0") +
    theme_reg() +
    scale_size_area(limits = c(0, 3)) +
    scale_y_continuous(limits = c(ymin, NA)) +
    guides(color = guide_legend(order = 0),
           size = guide_legend(order = 1))
}

# Set up ------------------------------------------------------------------

library(here)
library(tidyverse)

ages = 2:31
lengths = seq(41,99,2)
iss = 20

# Read in report file -----------------------------------------------------
model_path = here("ADMB_Models", "Final_Models")

grwth_vary_rep <- readLines(here(model_path, "Prop_Within_Growth_Vary", "tem.rep"))
grwth_vary_sab_curr <- dget(here(model_path, "Prop_Within_Growth_Vary", "tem.rdat")) 
grw_vary_ctl_wts <- readLines(here(model_path, "Prop_Within_Growth_Vary", "tem.ctl"))
grwth_sq_rep <- readLines(here(model_path, "Prop_Within_Growth_SQ", "tem.rep"))
grwth_sq_sab_curr <- dget(here(model_path, "Prop_Within_Growth_SQ", "tem.rdat")) 
grw_sq_ctl_wts <- readLines(here(model_path, "Prop_Within_Growth_SQ", "tem.ctl"))

# Extract out comps -------------------------------------------------------
## Fixed Gear Ages ---------------------------------------------------------

# Growth Vary
# Get fixed gear ages
gv_obs_ac_fish1_f = data.frame(reshape2::melt(grwth_vary_sab_curr$oac.fish1.f), type = "obs") %>% mutate(sex = "Female")
gv_pred_ac_fish1_f = data.frame(reshape2::melt(grwth_vary_sab_curr$eac.fish1.f), type = "pred") %>% mutate(sex = "Female")
gv_obs_ac_fish1_m = data.frame(reshape2::melt(grwth_vary_sab_curr$oac.fish1.m), type = "obs") %>% mutate(sex = "Male")
gv_pred_ac_fish1_m = data.frame(reshape2::melt(grwth_vary_sab_curr$eac.fish1.m), type = "pred") %>% mutate(sex = "Male")

# Put these into a dataframe
gv_ac_fish1 = rbind(gv_obs_ac_fish1_f, gv_pred_ac_fish1_f, gv_obs_ac_fish1_m, gv_pred_ac_fish1_m)
names(gv_ac_fish1) = c("year", "age", "prop", "type", "sex")
gv_ac_fish1$broodYear = gv_ac_fish1$year - gv_ac_fish1$age # create brood year to track cohort over time
gv_ac_fish1 = gv_ac_fish1 %>% mutate(broodYear = ifelse(age == 31, "31", broodYear), Model = "Growth Vary") # make plus group consistent color

# Status-quo model
# Get fixed gear ages
sq_obs_ac_fish1_f = data.frame(reshape2::melt(grwth_sq_sab_curr$oac.fish1.f), type = "obs") %>% mutate(sex = "Female")
sq_pred_ac_fish1_f = data.frame(reshape2::melt(grwth_sq_sab_curr$eac.fish1.f), type = "pred") %>% mutate(sex = "Female")
sq_obs_ac_fish1_m = data.frame(reshape2::melt(grwth_sq_sab_curr$oac.fish1.m), type = "obs") %>% mutate(sex = "Male")
sq_pred_ac_fish1_m = data.frame(reshape2::melt(grwth_sq_sab_curr$eac.fish1.m), type = "pred") %>% mutate(sex = "Male")

# Put these into a dataframe
sq_ac_fish1 = rbind(sq_obs_ac_fish1_f, sq_pred_ac_fish1_f, sq_obs_ac_fish1_m, sq_pred_ac_fish1_m)
names(sq_ac_fish1) = c("year", "age", "prop", "type", "sex")
sq_ac_fish1$broodYear = sq_ac_fish1$year - sq_ac_fish1$age # create brood year to track cohort over time
sq_ac_fish1 = sq_ac_fish1 %>% mutate(broodYear = ifelse(age == 31, "31", broodYear), Model = "Growth SQ") # make plus group consistent color

# plot average fits as well
avg_gv_ac_fish1 = gv_ac_fish1 %>%
  group_by(age, type, sex, Model) %>% 
  summarize(prop = mean(prop))

avg_sq_ac_fish1 = sq_ac_fish1 %>%
  group_by(age, type, sex, Model) %>% 
  summarize(prop = mean(prop)) 

# Bind these together
all_fixed_gear_plot_df = rbind(gv_ac_fish1, sq_ac_fish1)
all_avg_fixed_gear_plot_df = rbind(avg_gv_ac_fish1, avg_sq_ac_fish1)


# Residuals ---------------------------------------------------------------
sq_fish1f = get_osa_res(obs = grwth_sq_sab_curr$oac.fish1.f, 
                        pred = grwth_sq_sab_curr$eac.fish1.f, 
                        iss = rep(20, nrow(grwth_sq_sab_curr$eac.fish1.f)), 
                        iter_wt = as.numeric(strsplit(grw_sq_ctl_wts[str_detect(grw_sq_ctl_wts, "#wt fish1 age comp female iter")], " ")[[1]][1]), 
                        index = seq(3,31,1),
                        drop_bin = 1)

(sq_fish1f_plot = res_plot(data = sq_fish1f[[1]], 
                           comp_type = "Age", res_type = "(Growth SQ Fishery LL Age Comps Females)",
                           ymin = 2))

gv_fish1f = get_osa_res(obs = grwth_vary_sab_curr$oac.fish1.f, 
                        pred = grwth_vary_sab_curr$eac.fish1.f, 
                        iss = rep(20, nrow(grwth_vary_sab_curr$eac.fish1.f)), 
                        iter_wt = as.numeric(strsplit(grw_vary_ctl_wts[str_detect(grw_vary_ctl_wts, "#wt fish1 age comp female iter")], " ")[[1]][1]), 
                        index = seq(3,31,1),
                        drop_bin = 1)

(gv_fish1f_plot = res_plot(data = gv_fish1f[[1]], 
                           comp_type = "Age", res_type = "(Growth Vary Fishery LL Age Comps Females)",
                           ymin = 2))

sq_fish1m = get_osa_res(obs = grwth_sq_sab_curr$oac.fish1.m, 
                        pred = grwth_sq_sab_curr$eac.fish1.m, 
                        iss = rep(20, nrow(grwth_sq_sab_curr$eac.fish1.m)), 
                        iter_wt = as.numeric(strsplit(grw_sq_ctl_wts[str_detect(grw_sq_ctl_wts, "#wt fish1 age comp male iter")], " ")[[1]][1]), 
                        index = seq(3,31,1),
                        drop_bin = 1)

sq_fish1m_plot = res_plot(data = sq_fish1m[[1]], 
                          comp_type = "Age", res_type = "(Growth SQ Fishery LL Age Comps Males)",
                          ymin = 2)

gv_fish1m = get_osa_res(obs = grwth_vary_sab_curr$oac.fish1.m, 
                        pred = grwth_vary_sab_curr$eac.fish1.m, 
                        iss = rep(20, nrow(grwth_vary_sab_curr$eac.fish1.m)), 
                        iter_wt = as.numeric(strsplit(grw_vary_ctl_wts[str_detect(grw_vary_ctl_wts, "#wt fish1 age comp male iter")], " ")[[1]][1]), 
                        index = seq(3,31,1),
                        drop_bin = 1)

gv_fish1m_plot = res_plot(data = gv_fish1m[[1]], 
                          comp_type = "Age", res_type = "(Growth Vary Fishery LL Age Comps Males)",
                          ymin = 2)

pdf(here("figs", "Model Comparison", "FixedGear_Age_Resid.pdf"), width = 15)
ggpubr::ggarrange(sq_fish1f_plot, gv_fish1f_plot)
ggpubr::ggarrange(sq_fish1m_plot, gv_fish1m_plot)
dev.off()

### Fits --------------------------------------------------------------------

# Plot this out!
# fits across years
pdf(here("figs", "Model Comparison", "FixedGear_Age_Years.pdf"), width = 13, height = 10)
ggplot() +
  geom_col(all_fixed_gear_plot_df %>% filter(type == "obs", Model == "Growth SQ", sex == "Female"),
           mapping = aes(x = age, y = prop), alpha = 0.25, color = "black", fill = "orange") +
  geom_line(all_fixed_gear_plot_df %>% filter(type == "pred", sex == "Female"), 
            mapping = aes(x = age, y = prop, color = Model, group = Model), size = 0.85) +
  facet_wrap(~year) +
  theme_reg() +
  labs(x = "Age", y = "Female Proportions", fill = "Sex", color = "Model")

ggplot() +
  geom_col(all_fixed_gear_plot_df %>% filter(type == "obs", Model == "Growth SQ", sex == "Male"),
           mapping = aes(x = age, y = prop), alpha = 0.25, color = "black", fill = "orange") +
  geom_line(all_fixed_gear_plot_df %>% filter(type == "pred", sex == "Male"), 
            mapping = aes(x = age, y = prop, color = Model, group = Model), size = 0.85) +
  facet_wrap(~year) +
  theme_reg() +
  labs(x = "Age", y = "Male Proportions", fill = "Sex", color = "Model")

ggplot() +
  geom_col(all_fixed_gear_plot_df %>% filter(type == "obs", Model == "Growth SQ", sex == "Female", 
                                             year %in% 2014:2019),
           mapping = aes(x = age, y = prop), alpha = 0.25, color = "black", fill = "orange") +
  geom_line(all_fixed_gear_plot_df %>% filter(type == "pred", sex == "Female", 
                                              year %in% 2014:2019), 
            mapping = aes(x = age, y = prop, color = Model, group = Model), size = 0.85) +
  facet_wrap(~year) +
  theme_reg() +
  labs(x = "Age", y = "Female Proportions", fill = "Sex", color = "Model")

ggplot() +
  geom_col(all_fixed_gear_plot_df %>% filter(type == "obs", Model == "Growth SQ", sex == "Male",
                                             year %in% 2014:2019),
           mapping = aes(x = age, y = prop), alpha = 0.25, color = "black", fill = "orange") +
  geom_line(all_fixed_gear_plot_df %>% filter(type == "pred", sex == "Male",
                                              year %in% 2014:2019), 
            mapping = aes(x = age, y = prop, color = Model, group = Model), size = 0.85) +
  facet_wrap(~year) +
  theme_reg() +
  labs(x = "Age", y = "Male Proportions", fill = "Sex", color = "Model")

dev.off()

# Average fits
pdf(here("figs", "Model Comparison", "FixedGear_Age_Avg.pdf"), height = 5)
ggplot() +
  geom_col(all_avg_fixed_gear_plot_df %>% filter(type == "obs", Model  == "Growth SQ"),
           mapping = aes(x = age, y = prop), alpha = 0.25, color = "black", fill = "orange") +
  geom_line(all_avg_fixed_gear_plot_df %>% filter(type == "pred"), 
            mapping = aes(x = age, y = prop, color = Model, group = Model), size = 1) +
  theme_reg() +
  facet_wrap(~sex) +
  labs(x = "Age", y = "Proportions", fill = "Sex", Model = "Sex")
dev.off()

## Fixed Gear lens ---------------------------------------------------------

# Growth Vary
# Get fixed gear lens
gv_obs_lc_fish1_f = data.frame(reshape2::melt(grwth_vary_sab_curr$olc.fish1.f), type = "obs") %>% mutate(sex = "Female")
gv_pred_lc_fish1_f = data.frame(reshape2::melt(grwth_vary_sab_curr$elc.fish1.f), type = "pred") %>% mutate(sex = "Female")
gv_obs_lc_fish1_m = data.frame(reshape2::melt(grwth_vary_sab_curr$olc.fish1.m), type = "obs") %>% mutate(sex = "Male")
gv_pred_lc_fish1_m = data.frame(reshape2::melt(grwth_vary_sab_curr$elc.fish1.m), type = "pred") %>% mutate(sex = "Male")

# Put these into a dataframe
gv_lc_fish1 = rbind(gv_obs_lc_fish1_f, gv_pred_lc_fish1_f, gv_obs_lc_fish1_m, gv_pred_lc_fish1_m)
names(gv_lc_fish1) = c("year", "len", "prop", "type", "sex")
gv_lc_fish1 = gv_lc_fish1 %>% mutate(Model = "Growth Vary") 

# Status-quo model
# Get fixed gear lens
sq_obs_lc_fish1_f = data.frame(reshape2::melt(grwth_sq_sab_curr$olc.fish1.f), type = "obs") %>% mutate(sex = "Female")
sq_pred_lc_fish1_f = data.frame(reshape2::melt(grwth_sq_sab_curr$elc.fish1.f), type = "pred") %>% mutate(sex = "Female")
sq_obs_lc_fish1_m = data.frame(reshape2::melt(grwth_sq_sab_curr$olc.fish1.m), type = "obs") %>% mutate(sex = "Male")
sq_pred_lc_fish1_m = data.frame(reshape2::melt(grwth_sq_sab_curr$elc.fish1.m), type = "pred") %>% mutate(sex = "Male")

# Put these into a dataframe
sq_lc_fish1 = rbind(sq_obs_lc_fish1_f, sq_pred_lc_fish1_f, sq_obs_lc_fish1_m, sq_pred_lc_fish1_m)
names(sq_lc_fish1) = c("year", "len", "prop", "type", "sex")
sq_lc_fish1 = sq_lc_fish1 %>% mutate(Model = "Growth SQ") # make plus group consistent color

# plot average fits as well
avg_gv_lc_fish1 = gv_lc_fish1 %>%
  group_by(len, type, sex, Model) %>% 
  summarize(prop = mean(prop))

avg_sq_lc_fish1 = sq_lc_fish1 %>%
  group_by(len, type, sex, Model) %>% 
  summarize(prop = mean(prop)) 

# Bind these together
all_fixed_gear_plot_df = rbind(gv_lc_fish1, sq_lc_fish1)
all_avg_fixed_gear_plot_df = rbind(avg_gv_lc_fish1, avg_sq_lc_fish1)


# Residuals ---------------------------------------------------------------

sq_lenfish1f = get_osa_res(obs = grwth_sq_sab_curr$olc.fish1.f, 
                           pred = grwth_sq_sab_curr$elc.fish1.f, 
                           iss = rep(20, nrow(grwth_sq_sab_curr$elc.fish1.f)), 
                           iter_wt = as.numeric(strsplit(grw_sq_ctl_wts[str_detect(grw_sq_ctl_wts, "#wt fish1 size comp female iter")], " ")[[1]][1]), 
                           index = seq(43,99,2),
                           drop_bin = 1)

sq_lenfish1f_plot = res_plot(data = sq_lenfish1f[[1]], 
                             comp_type = "Length", res_type = "(Growth SQ Fishery LL size Comps Females)",
                             ymin = 41)

gv_lenfish1f = get_osa_res(obs = grwth_vary_sab_curr$olc.fish1.f, 
                           pred = grwth_vary_sab_curr$elc.fish1.f, 
                           iss = rep(20, nrow(grwth_vary_sab_curr$elc.fish1.f)), 
                           iter_wt = as.numeric(strsplit(grw_vary_ctl_wts[str_detect(grw_vary_ctl_wts, "#wt fish1 size comp female iter")], " ")[[1]][1]), 
                           index = seq(43,99,2),
                           drop_bin = 1)

gv_lenfish1f_plot = res_plot(data = gv_lenfish1f[[1]], 
                             comp_type = "Length", res_type = "(Growth Vary Fishery LL size Comps Females)",
                             ymin = 43)

sq_fish1m = get_osa_res(obs = grwth_sq_sab_curr$olc.fish1.m, 
                        pred = grwth_sq_sab_curr$elc.fish1.m, 
                        iss = rep(20, nrow(grwth_sq_sab_curr$elc.fish1.m)), 
                        iter_wt = as.numeric(strsplit(grw_sq_ctl_wts[str_detect(grw_sq_ctl_wts, "#wt fish1 size comp male iter")], " ")[[1]][1]), 
                        index = seq(43,99,2),
                        drop_bin = 1)

sq_fish1m_plot = res_plot(data = sq_fish1m[[1]], 
                          comp_type = "size", res_type = "(Growth SQ Fishery LL size Comps Males)",
                          ymin = 43)

gv_fish1m = get_osa_res(obs = grwth_vary_sab_curr$olc.fish1.m, 
                        pred = grwth_vary_sab_curr$elc.fish1.m, 
                        iss = rep(20, nrow(grwth_vary_sab_curr$elc.fish1.m)), 
                        iter_wt = as.numeric(strsplit(grw_vary_ctl_wts[str_detect(grw_vary_ctl_wts, "#wt fish1 size comp male iter")], " ")[[1]][1]), 
                        index = seq(43,99,2),
                        drop_bin = 1)

(gv_fish1m_plot = res_plot(data = gv_fish1m[[1]], 
                           comp_type = "size", res_type = "(Growth Vary Fishery LL size Comps Males)",
                           ymin = 43))

pdf(here("figs", "Model Comparison", "FisheryLL_size_Resid.pdf"), width = 15)
ggpubr::ggarrange(sq_lenfish1f_plot, gv_lenfish1f_plot, sq_fish1m_plot, gv_fish1m_plot)
dev.off()


### Fits --------------------------------------------------------------------

# Plot this out!
# fits across years
pdf(here("figs", "Model Comparison", "FixedGear_Len_Years.pdf"), width = 13, height = 10)
ggplot() +
  geom_col(all_fixed_gear_plot_df %>% filter(type == "obs", Model == "Growth SQ", sex == "Female"),
           mapping = aes(x = len, y = prop), alpha = 0.35, color = "black", fill = "darkgreen") +
  geom_line(all_fixed_gear_plot_df %>% filter(type == "pred", sex == "Female"), 
            mapping = aes(x = len, y = prop, color = Model, group = Model), size = 0.85) +
  facet_wrap(~year) +
  theme_reg() +
  labs(x = "Length", y = "Female Proportions", fill = "Sex", color = "Model")

ggplot() +
  geom_col(all_fixed_gear_plot_df %>% filter(type == "obs", Model == "Growth SQ", sex == "Male"),
           mapping = aes(x = len, y = prop), alpha = 0.35, color = "black", fill = "darkgreen") +
  geom_line(all_fixed_gear_plot_df %>% filter(type == "pred", sex == "Male"), 
            mapping = aes(x = len, y = prop, color = Model, group = Model), size = 0.85) +
  facet_wrap(~year) +
  theme_reg() +
  labs(x = "Length", y = "Male Proportions", fill = "Sex", color = "Model")
ggplot() +
  geom_col(all_fixed_gear_plot_df %>% filter(type == "obs", Model == "Growth SQ", sex == "Female",
                                             year %in% 2014:2019),
           mapping = aes(x = len, y = prop), alpha = 0.35, color = "black", fill = "darkgreen") +
  geom_line(all_fixed_gear_plot_df %>% filter(type == "pred", sex == "Female",
                                              year %in% 2014:2019), 
            mapping = aes(x = len, y = prop, color = Model, group = Model), size = 0.85) +
  facet_wrap(~year) +
  theme_reg() +
  labs(x = "Length", y = "Female Proportions", fill = "Sex", color = "Model")

ggplot() +
  geom_col(all_fixed_gear_plot_df %>% filter(type == "obs", Model == "Growth SQ", sex == "Male",
                                             year %in% 2014:2019),
           mapping = aes(x = len, y = prop), alpha = 0.35, color = "black", fill = "darkgreen") +
  geom_line(all_fixed_gear_plot_df %>% filter(type == "pred", sex == "Male",
                                              year %in% 2014:2019), 
            mapping = aes(x = len, y = prop, color = Model, group = Model), size = 0.85) +
  facet_wrap(~year) +
  theme_reg() +
  labs(x = "Length", y = "Male Proportions", fill = "Sex", color = "Model")

dev.off()

# Average fits
pdf(here("figs", "Model Comparison", "FixedGear_Len_Avg.pdf"), height = 5)
ggplot() +
  geom_col(all_avg_fixed_gear_plot_df %>% filter(type == "obs", Model  == "Growth SQ"),
           mapping = aes(x = len, y = prop), alpha = 0.25, color = "black", fill = "orange") +
  geom_line(all_avg_fixed_gear_plot_df %>% filter(type == "pred"), 
            mapping = aes(x = len, y = prop, color = Model, group = Model), size = 0.85) +
  theme_reg() +
  facet_wrap(~sex) +
  labs(x = "Length", y = "Proportions", fill = "Sex", Model = "Sex")
dev.off()

# Trawl Gear --------------------------------------------------------------
## Trawl Gear lens ---------------------------------------------------------

# Growth Vary
# Get Trawl Gear lens
gv_obs_lc_fish3_f = data.frame(reshape2::melt(grwth_vary_sab_curr$olc.fish3.f), type = "obs") %>% mutate(sex = "Female")
gv_pred_lc_fish3_f = data.frame(reshape2::melt(grwth_vary_sab_curr$elc.fish3.f), type = "pred") %>% mutate(sex = "Female")
gv_obs_lc_fish3_m = data.frame(reshape2::melt(grwth_vary_sab_curr$olc.fish3.m), type = "obs") %>% mutate(sex = "Male")
gv_pred_lc_fish3_m = data.frame(reshape2::melt(grwth_vary_sab_curr$elc.fish3.m), type = "pred") %>% mutate(sex = "Male")

# Put these into a dataframe
gv_lc_fish3 = rbind(gv_obs_lc_fish3_f, gv_pred_lc_fish3_f, gv_obs_lc_fish3_m, gv_pred_lc_fish3_m)
names(gv_lc_fish3) = c("year", "len", "prop", "type", "sex")
gv_lc_fish3 = gv_lc_fish3 %>% mutate(Model = "Growth Vary") 

# Status-quo model
# Get Trawl Gear lens
sq_obs_lc_fish3_f = data.frame(reshape2::melt(grwth_sq_sab_curr$olc.fish3.f), type = "obs") %>% mutate(sex = "Female")
sq_pred_lc_fish3_f = data.frame(reshape2::melt(grwth_sq_sab_curr$elc.fish3.f), type = "pred") %>% mutate(sex = "Female")
sq_obs_lc_fish3_m = data.frame(reshape2::melt(grwth_sq_sab_curr$olc.fish3.m), type = "obs") %>% mutate(sex = "Male")
sq_pred_lc_fish3_m = data.frame(reshape2::melt(grwth_sq_sab_curr$elc.fish3.m), type = "pred") %>% mutate(sex = "Male")

# Put these into a dataframe
sq_lc_fish3 = rbind(sq_obs_lc_fish3_f, sq_pred_lc_fish3_f, sq_obs_lc_fish3_m, sq_pred_lc_fish3_m)
names(sq_lc_fish3) = c("year", "len", "prop", "type", "sex")
sq_lc_fish3 = sq_lc_fish3 %>% mutate(Model = "Growth SQ") # make plus group consistent color


# plot average fits as well
avg_gv_lc_fish3 = gv_lc_fish3 %>%
  group_by(len, type, sex, Model) %>% 
  summarize(prop = mean(prop))

avg_sq_lc_fish3 = sq_lc_fish3 %>%
  group_by(len, type, sex, Model) %>% 
  summarize(prop = mean(prop)) 

# Bind these together
all_trawl_gear_plot_df = rbind(gv_lc_fish3, sq_lc_fish3)
all_avg_trawl_gear_plot_df = rbind(avg_gv_lc_fish3, avg_sq_lc_fish3)


# Residuals ---------------------------------------------------------------

sq_lenfish3f = get_osa_res(obs = grwth_sq_sab_curr$olc.fish3.f, 
                           pred = grwth_sq_sab_curr$elc.fish3.f, 
                           iss = rep(20, nrow(grwth_sq_sab_curr$elc.fish3.f)), 
                           iter_wt = as.numeric(strsplit(grw_sq_ctl_wts[str_detect(grw_sq_ctl_wts, "#wt fish 3 size comp female iter")], " ")[[1]][1]), 
                           index = seq(43,99,2),
                           drop_bin = 1)

sq_lenfish3f_plot = res_plot(data = sq_lenfish3f[[1]], 
                             comp_type = "Length", res_type = "(Growth SQ Fishery Trawl size Comps Females)",
                             ymin = 41)

gv_lenfish3f = get_osa_res(obs = grwth_vary_sab_curr$olc.fish3.f, 
                           pred = grwth_vary_sab_curr$elc.fish3.f, 
                           iss = rep(20, nrow(grwth_vary_sab_curr$elc.fish3.f)), 
                           iter_wt = as.numeric(strsplit(grw_vary_ctl_wts[str_detect(grw_vary_ctl_wts, "#wt fish 3 size comp female iter")], " ")[[1]][1]), 
                           index = seq(43,99,2),
                           drop_bin = 1)

gv_lenfish3f_plot = res_plot(data = gv_lenfish3f[[1]], 
                             comp_type = "Length", res_type = "(Growth Vary Fishery Trawl size Comps Females)",
                             ymin = 43)

sq_fish3m = get_osa_res(obs = grwth_sq_sab_curr$olc.fish3.m, 
                        pred = grwth_sq_sab_curr$elc.fish3.m, 
                        iss = rep(20, nrow(grwth_sq_sab_curr$elc.fish3.m)), 
                        iter_wt = as.numeric(strsplit(grw_sq_ctl_wts[str_detect(grw_sq_ctl_wts, "#wt fish 3 size comp male iter")], " ")[[1]][1]), 
                        index = seq(43,99,2),
                        drop_bin = 1)

sq_fish3m_plot = res_plot(data = sq_fish3m[[1]], 
                          comp_type = "size", res_type = "(Growth SQ Fishery Trawl size Comps Males)",
                          ymin = 43)

gv_fish3m = get_osa_res(obs = grwth_vary_sab_curr$olc.fish3.m, 
                        pred = grwth_vary_sab_curr$elc.fish3.m, 
                        iss = rep(20, nrow(grwth_vary_sab_curr$elc.fish3.m)), 
                        iter_wt = as.numeric(strsplit(grw_vary_ctl_wts[str_detect(grw_vary_ctl_wts, "#wt fish 3 size comp male iter")], " ")[[1]][1]), 
                        index = seq(43,99,2),
                        drop_bin = 1)

gv_fish3m_plot = res_plot(data = gv_fish3m[[1]], 
                          comp_type = "size", res_type = "(Growth Vary Fishery Trawl size Comps Males)",
                          ymin = 43)

pdf(here("figs", "Model Comparison", "FisheryTrawl_size_Resid.pdf"), width = 15)
ggpubr::ggarrange(sq_lenfish3f_plot, gv_lenfish3f_plot, sq_fish3m_plot, gv_fish3m_plot)
dev.off()


### Fits --------------------------------------------------------------------

# Plot this out!
# fits across years
pdf(here("figs", "Model Comparison", "TrawlGear_Len_Years.pdf"), width = 13, height = 10)
ggplot() +
  geom_col(all_trawl_gear_plot_df %>% filter(type == "obs", Model == "Growth SQ", sex == "Female"),
           mapping = aes(x = len, y = prop), alpha = 0.35, color = "black", fill = "darkgreen") +
  geom_line(all_trawl_gear_plot_df %>% filter(type == "pred", sex == "Female"), 
            mapping = aes(x = len, y = prop, color = Model, group = Model), size = 0.85) +
  facet_wrap(~year) +
  theme_reg() +
  labs(x = "Length", y = "Female Proportions", fill = "Sex", color = "Model")

ggplot() +
  geom_col(all_trawl_gear_plot_df %>% filter(type == "obs", Model == "Growth SQ", sex == "Male"),
           mapping = aes(x = len, y = prop), alpha = 0.35, color = "black", fill = "darkgreen") +
  geom_line(all_trawl_gear_plot_df %>% filter(type == "pred", sex == "Male"), 
            mapping = aes(x = len, y = prop, color = Model, group = Model), size = 0.85) +
  facet_wrap(~year) +
  theme_reg() +
  labs(x = "Length", y = "Male Proportions", fill = "Sex", color = "Model")
ggplot() +
  geom_col(all_trawl_gear_plot_df %>% filter(type == "obs", Model == "Growth SQ", sex == "Female",
                                             year %in% 2014:2019),
           mapping = aes(x = len, y = prop), alpha = 0.35, color = "black", fill = "darkgreen") +
  geom_line(all_trawl_gear_plot_df %>% filter(type == "pred", sex == "Female",
                                              year %in% 2014:2019), 
            mapping = aes(x = len, y = prop, color = Model, group = Model), size = 0.85) +
  facet_wrap(~year) +
  theme_reg() +
  labs(x = "Length", y = "Female Proportions", fill = "Sex", color = "Model")

ggplot() +
  geom_col(all_trawl_gear_plot_df %>% filter(type == "obs", Model == "Growth SQ", sex == "Male",
                                             year %in% 2014:2019),
           mapping = aes(x = len, y = prop), alpha = 0.35, color = "black", fill = "darkgreen") +
  geom_line(all_trawl_gear_plot_df %>% filter(type == "pred", sex == "Male",
                                              year %in% 2014:2019), 
            mapping = aes(x = len, y = prop, color = Model, group = Model), size = 0.85) +
  facet_wrap(~year) +
  theme_reg() +
  labs(x = "Length", y = "Male Proportions", fill = "Sex", color = "Model")

dev.off()

# Average fits
pdf(here("figs", "Model Comparison", "TrawlGear_Len_Avg.pdf"), height = 5)
ggplot() +
  geom_col(all_avg_trawl_gear_plot_df %>% filter(type == "obs", Model  == "Growth SQ"),
           mapping = aes(x = len, y = prop), alpha = 0.25, color = "black", fill = "orange") +
  geom_line(all_avg_trawl_gear_plot_df %>% filter(type == "pred"), 
            mapping = aes(x = len, y = prop, color = Model, group = Model), size = 0.85) +
  theme_reg() +
  facet_wrap(~sex) +
  labs(x = "Length", y = "Proportions", fill = "Sex", Model = "Sex")
dev.off()


# Survey LL Ages -----------------------------------------------------------
# Status-quo model 
# Get srv ll ages
sq_obs_ac_srv1_f = data.frame(reshape2::melt(grwth_sq_sab_curr$oac.srv1.f), type = "obs") %>% mutate(sex = "Female")
sq_pred_ac_srv1_f = data.frame(reshape2::melt(grwth_sq_sab_curr$eac.srv1.f), type = "pred") %>% mutate(sex = "Female")
sq_obs_ac_srv1_m = data.frame(reshape2::melt(grwth_sq_sab_curr$oac.srv1.m), type = "obs") %>% mutate(sex = "Male")
sq_pred_ac_srv1_m = data.frame(reshape2::melt(grwth_sq_sab_curr$eac.srv1.m), type = "pred") %>% mutate(sex = "Male")

# Put these into a dataframe
sq_ac_srv1 = rbind(sq_obs_ac_srv1_f, sq_pred_ac_srv1_f, sq_obs_ac_srv1_m, sq_pred_ac_srv1_m)
names(sq_ac_srv1) = c("year", "age", "prop", "type", "sex")
sq_ac_srv1$broodYear = sq_ac_srv1$year - sq_ac_srv1$age # create brood year to track cohort over time
sq_ac_srv1 = sq_ac_srv1 %>% mutate(broodYear = ifelse(age == 31, "31", broodYear), Model = "Growth SQ") # make plus group consistent color

# Growth model vary 
# Get srv ll ages
gv_obs_ac_srv1_f = data.frame(reshape2::melt(grwth_vary_sab_curr$oac.srv1.f), type = "obs") %>% mutate(sex = "Female")
gv_pred_ac_srv1_f = data.frame(reshape2::melt(grwth_vary_sab_curr$eac.srv1.f), type = "pred") %>% mutate(sex = "Female")
gv_obs_ac_srv1_m = data.frame(reshape2::melt(grwth_vary_sab_curr$oac.srv1.m), type = "obs") %>% mutate(sex = "Male")
gv_pred_ac_srv1_m = data.frame(reshape2::melt(grwth_vary_sab_curr$eac.srv1.m), type = "pred") %>% mutate(sex = "Male")

# Put these into a dataframe
gv_ac_srv1 = rbind(gv_pred_ac_srv1_m, gv_obs_ac_srv1_m, gv_pred_ac_srv1_f, gv_obs_ac_srv1_f)
names(gv_ac_srv1) = c("year", "age", "prop", "type", "sex")
gv_ac_srv1$broodYear = sq_ac_srv1$year - sq_ac_srv1$age # create brood year to track cohort over time
gv_ac_srv1 = gv_ac_srv1 %>% mutate(broodYear = ifelse(age == 31, "31", broodYear), Model = "Growth Vary") # make plus group consistent color

# plot average fits as well
avg_gv_ac_srv1 = gv_ac_srv1 %>%
  group_by(age, type, sex, Model) %>% 
  summarize(prop = mean(prop))

avg_sq_ac_srv1 = sq_ac_srv1 %>%
  group_by(age, type, sex, Model) %>% 
  summarize(prop = mean(prop)) 

# Bind these together
all_srv_ll_plot_df = rbind(gv_ac_srv1, sq_ac_srv1)
all_avg_srv_ll_plot_df = rbind(avg_gv_ac_srv1, avg_sq_ac_srv1)


# Residuals ---------------------------------------------------------------

sq_srv1f = get_osa_res(obs = grwth_sq_sab_curr$oac.srv1.f, 
                       pred = grwth_sq_sab_curr$eac.srv1.f, 
                       iss = rep(20, nrow(grwth_sq_sab_curr$eac.srv1.f)), 
                       iter_wt = as.numeric(strsplit(grw_sq_ctl_wts[str_detect(grw_sq_ctl_wts, "#wt surv1 age comp female iter")], " ")[[1]][1]), 
                       index = seq(3,31,1),
                       drop_bin = 1)

sq_srv1f_plot = res_plot(data = sq_srv1f[[1]], 
                         comp_type = "Age", res_type = "(Growth SQ Survey LL Age Comps Females)",
                         ymin = 2)

gv_srv1f = get_osa_res(obs = grwth_vary_sab_curr$oac.srv1.f, 
                       pred = grwth_vary_sab_curr$eac.srv1.f, 
                       iss = rep(20, nrow(grwth_vary_sab_curr$eac.srv1.f)), 
                       iter_wt = as.numeric(strsplit(grw_vary_ctl_wts[str_detect(grw_vary_ctl_wts, "#wt surv1 age comp female iter")], " ")[[1]][1]), 
                       index = seq(3,31,1),
                       drop_bin = 1)

gv_srv1f_plot = res_plot(data = gv_srv1f[[1]], 
                         comp_type = "Age", res_type = "(Growth Vary Survey LL Age Comps Females)",
                         ymin = 2)

sq_srv1m = get_osa_res(obs = grwth_sq_sab_curr$oac.srv1.m, 
                       pred = grwth_sq_sab_curr$eac.srv1.m, 
                       iss = rep(20, nrow(grwth_sq_sab_curr$eac.srv1.m)), 
                       iter_wt = as.numeric(strsplit(grw_sq_ctl_wts[str_detect(grw_sq_ctl_wts, "#wt surv1 age comp male iter")], " ")[[1]][1]), 
                       index = seq(3,31,1),
                       drop_bin = 1)

sq_srv1m_plot = res_plot(data = sq_srv1m[[1]], 
                         comp_type = "Age", res_type = "(Growth SQ Survey LL Age Comps Males)",
                         ymin = 2)

gv_srv1m = get_osa_res(obs = grwth_vary_sab_curr$oac.srv1.m, 
                       pred = grwth_vary_sab_curr$eac.srv1.m, 
                       iss = rep(20, nrow(grwth_vary_sab_curr$eac.srv1.m)), 
                       iter_wt = as.numeric(strsplit(grw_vary_ctl_wts[str_detect(grw_vary_ctl_wts, "#wt surv1 age comp male iter")], " ")[[1]][1]), 
                       index = seq(3,31,1),
                       drop_bin = 1)

gv_srv1m_plot = res_plot(data = gv_srv1m[[1]], 
                         comp_type = "Age", res_type = "(Growth Vary Survey LL Age Comps Males)",
                         ymin = 2)

pdf(here("figs", "Model Comparison", "SurveyLL_Age_Resid.pdf"), width = 15)
ggpubr::ggarrange(sq_srv1f_plot, gv_srv1f_plot, sq_srv1m_plot, gv_srv1m_plot)
dev.off()


### Fits --------------------------------------------------------------------

# Plot this out!
# fits across years
pdf(here("figs", "Model Comparison", "SurveyLL_Age_Years.pdf"), width = 13, height = 10)
ggplot() +
  geom_col(all_srv_ll_plot_df %>% filter(type == "obs", Model == "Growth SQ", sex == "Female"),
           mapping = aes(x = age, y = prop), alpha = 0.25, color = "black", fill = "orange") +
  geom_line(all_srv_ll_plot_df %>% filter(type == "pred", sex == "Female"), 
            mapping = aes(x = age, y = prop, color = Model, group = Model), size = 0.85) +
  facet_wrap(~year) +
  theme_reg() +
  labs(x = "Age", y = "Female Proportions", fill = "Sex", color = "Model")

ggplot() +
  geom_col(all_srv_ll_plot_df %>% filter(type == "obs", Model == "Growth SQ", sex == "Male"),
           mapping = aes(x = age, y = prop), alpha = 0.25, color = "black", fill = "orange") +
  geom_line(all_srv_ll_plot_df %>% filter(type == "pred", sex == "Male"), 
            mapping = aes(x = age, y = prop, color = Model, group = Model), size = 0.85) +
  facet_wrap(~year) +
  theme_reg() +
  labs(x = "Age", y = "Male Proportions", fill = "Sex", color = "Model")

ggplot() +
  geom_col(all_srv_ll_plot_df %>% filter(type == "obs", Model == "Growth SQ", sex == "Female", 
                                         year %in% 2014:2019),
           mapping = aes(x = age, y = prop), alpha = 0.25, color = "black", fill = "orange") +
  geom_line(all_srv_ll_plot_df %>% filter(type == "pred", sex == "Female", 
                                          year %in% 2014:2019), 
            mapping = aes(x = age, y = prop, color = Model, group = Model), size = 0.85) +
  facet_wrap(~year) +
  theme_reg() +
  labs(x = "Age", y = "Female Proportions", fill = "Sex", color = "Model")

ggplot() +
  geom_col(all_srv_ll_plot_df %>% filter(type == "obs", Model == "Growth SQ", sex == "Male",
                                         year %in% 2014:2019),
           mapping = aes(x = age, y = prop), alpha = 0.25, color = "black", fill = "orange") +
  geom_line(all_srv_ll_plot_df %>% filter(type == "pred", sex == "Male",
                                          year %in% 2014:2019), 
            mapping = aes(x = age, y = prop, color = Model, group = Model), size = 0.85) +
  facet_wrap(~year) +
  theme_reg() +
  labs(x = "Age", y = "Male Proportions", fill = "Sex", color = "Model")
dev.off()

# Average fits
pdf(here("figs", "Model Comparison", "SurveyLL_Age_Avg.pdf"), height = 5)
ggplot() +
  geom_col(all_avg_srv_ll_plot_df %>% filter(type == "obs", Model  == "Growth SQ"),
           mapping = aes(x = age, y = prop), alpha = 0.25, color = "black", fill = "orange") +
  geom_line(all_avg_srv_ll_plot_df %>% filter(type == "pred"), 
            mapping = aes(x = age, y = prop, color = Model, group = Model), size = 0.85) +
  theme_reg() +
  facet_wrap(~sex) +
  labs(x = "Age", y = "Proportions", fill = "Sex", Model = "Sex")
dev.off()

# Survey LL Sizes ---------------------------------------------------------

# Growth Vary
# Get srv ll lens
gv_obs_lc_srv1_f = data.frame(reshape2::melt(grwth_vary_sab_curr$olc.srv1.f), type = "obs") %>% mutate(sex = "Female")
gv_pred_lc_srv1_f = data.frame(reshape2::melt(grwth_vary_sab_curr$elc.srv1.f), type = "pred") %>% mutate(sex = "Female")
gv_obs_lc_srv1_m = data.frame(reshape2::melt(grwth_vary_sab_curr$olc.srv1.m), type = "obs") %>% mutate(sex = "Male")
gv_pred_lc_srv1_m = data.frame(reshape2::melt(grwth_vary_sab_curr$elc.srv1.m), type = "pred") %>% mutate(sex = "Male")

# Put these into a dataframe
gv_lc_srv1 = rbind(gv_obs_lc_srv1_f, gv_pred_lc_srv1_f, gv_obs_lc_srv1_m, gv_pred_lc_srv1_m)
names(gv_lc_srv1) = c("year", "len", "prop", "type", "sex")
gv_lc_srv1 = gv_lc_srv1 %>% mutate(Model = "Growth Vary") 

# Status-quo model
# Get srv ll lens
sq_obs_lc_srv1_f = data.frame(reshape2::melt(grwth_sq_sab_curr$olc.srv1.f), type = "obs") %>% mutate(sex = "Female")
sq_pred_lc_srv1_f = data.frame(reshape2::melt(grwth_sq_sab_curr$elc.srv1.f), type = "pred") %>% mutate(sex = "Female")
sq_obs_lc_srv1_m = data.frame(reshape2::melt(grwth_sq_sab_curr$olc.srv1.m), type = "obs") %>% mutate(sex = "Male")
sq_pred_lc_srv1_m = data.frame(reshape2::melt(grwth_sq_sab_curr$elc.srv1.m), type = "pred") %>% mutate(sex = "Male")

# Put these into a dataframe
sq_lc_srv1 = rbind(sq_obs_lc_srv1_f, sq_pred_lc_srv1_f, sq_obs_lc_srv1_m, sq_pred_lc_srv1_m)
names(sq_lc_srv1) = c("year", "len", "prop", "type", "sex")
sq_lc_srv1 = sq_lc_srv1 %>% mutate(Model = "Growth SQ") # make plus group consistent color

# plot average fits as well
avg_gv_lc_srv1 = gv_lc_srv1 %>%
  group_by(len, type, sex, Model) %>% 
  summarize(prop = mean(prop))

avg_sq_lc_srv1 = sq_lc_srv1 %>%
  group_by(len, type, sex, Model) %>% 
  summarize(prop = mean(prop)) 

# Bind these together
all_srv_ll_plot_df = rbind(gv_lc_srv1, sq_lc_srv1)
all_avg_srv_ll_plot_df = rbind(avg_gv_lc_srv1, avg_sq_lc_srv1)


# Residuals ---------------------------------------------------------------

sq_lensrv1f = get_osa_res(obs = grwth_sq_sab_curr$olc.srv1.f, 
                          pred = grwth_sq_sab_curr$elc.srv1.f, 
                          iss = rep(20, nrow(grwth_sq_sab_curr$elc.srv1.f)), 
                          iter_wt = as.numeric(strsplit(grw_sq_ctl_wts[str_detect(grw_sq_ctl_wts, "#wt surv1 size comp female iter")], " ")[[1]][1]), 
                          index = seq(43,99,2),
                          drop_bin = 1)

sq_lensrv1f_plot = res_plot(data = sq_lensrv1f[[1]], 
                            comp_type = "Length", res_type = "(Growth SQ Survey LL size Comps Females)",
                            ymin = 41)

gv_lensrv1f = get_osa_res(obs = grwth_vary_sab_curr$olc.srv1.f, 
                          pred = grwth_vary_sab_curr$elc.srv1.f, 
                          iss = rep(20, nrow(grwth_vary_sab_curr$elc.srv1.f)), 
                          iter_wt = as.numeric(strsplit(grw_vary_ctl_wts[str_detect(grw_vary_ctl_wts, "#wt surv1 size comp female iter")], " ")[[1]][1]), 
                          index = seq(43,99,2),
                          drop_bin = 1)

gv_lensrv1f_plot = res_plot(data = gv_lensrv1f[[1]], 
                            comp_type = "Length", res_type = "(Growth Vary Survey LL size Comps Females)",
                            ymin = 43)

sq_srv1m = get_osa_res(obs = grwth_sq_sab_curr$olc.srv1.m, 
                       pred = grwth_sq_sab_curr$elc.srv1.m, 
                       iss = rep(20, nrow(grwth_sq_sab_curr$elc.srv1.m)), 
                       iter_wt = as.numeric(strsplit(grw_sq_ctl_wts[str_detect(grw_sq_ctl_wts, "#wt surv1 size comp male iter")], " ")[[1]][1]), 
                       index = seq(43,99,2),
                       drop_bin = 1)

sq_srv1m_plot = res_plot(data = sq_srv1m[[1]], 
                         comp_type = "size", res_type = "(Growth SQ Survey LL size Comps Males)",
                         ymin = 43)

gv_srv1m = get_osa_res(obs = grwth_vary_sab_curr$olc.srv1.m, 
                       pred = grwth_vary_sab_curr$elc.srv1.m, 
                       iss = rep(20, nrow(grwth_vary_sab_curr$elc.srv1.m)), 
                       iter_wt = as.numeric(strsplit(grw_vary_ctl_wts[str_detect(grw_vary_ctl_wts, "#wt surv1 size comp male iter")], " ")[[1]][1]), 
                       index = seq(43,99,2),
                       drop_bin = 1)

gv_srv1m_plot = res_plot(data = gv_srv1m[[1]], 
                         comp_type = "size", res_type = "(Growth Vary Survey LL size Comps Males)",
                         ymin = 43)

pdf(here("figs", "Model Comparison", "SurveyLL_size_Resid.pdf"), width = 15)
ggpubr::ggarrange(sq_lensrv1f_plot, gv_lensrv1f_plot, sq_srv1m_plot, gv_srv1m_plot)
dev.off()



### Fits --------------------------------------------------------------------

# Plot this out!
# fits across years
pdf(here("figs", "Model Comparison", "SurveyLL_len_Years.pdf"), width = 13, height = 10)
ggplot() +
  geom_col(all_srv_ll_plot_df %>% filter(type == "obs", Model == "Growth SQ", sex == "Female"),
           mapping = aes(x = len, y = prop), alpha = 0.25, color = "black", fill = "orange") +
  geom_line(all_srv_ll_plot_df %>% filter(type == "pred", sex == "Female"), 
            mapping = aes(x = len, y = prop, color = Model, group = Model), size = 0.85) +
  facet_wrap(~year) +
  theme_reg() +
  labs(x = "Length", y = "Female Proportions", fill = "Sex", color = "Model")

ggplot() +
  geom_col(all_srv_ll_plot_df %>% filter(type == "obs", Model == "Growth SQ", sex == "Male"),
           mapping = aes(x = len, y = prop), alpha = 0.25, color = "black", fill = "orange") +
  geom_line(all_srv_ll_plot_df %>% filter(type == "pred", sex == "Male"), 
            mapping = aes(x = len, y = prop, color = Model, group = Model), size = 0.85) +
  facet_wrap(~year) +
  theme_reg() +
  labs(x = "Length", y = "Male Proportions", fill = "Sex", color = "Model")

ggplot() +
  geom_col(all_srv_ll_plot_df %>% filter(type == "obs", Model == "Growth SQ", sex == "Female",
                                         year %in% 2014:2019),
           mapping = aes(x = len, y = prop), alpha = 0.35, color = "black", fill = "darkgreen") +
  geom_line(all_srv_ll_plot_df %>% filter(type == "pred", sex == "Female",
                                          year %in% 2014:2019), 
            mapping = aes(x = len, y = prop, color = Model, group = Model), size = 0.85) +
  facet_wrap(~year) +
  theme_reg() +
  labs(x = "Length", y = "Female Proportions", fill = "Sex", color = "Model")

ggplot() +
  geom_col(all_srv_ll_plot_df %>% filter(type == "obs", Model == "Growth SQ", sex == "Male",
                                         year %in% 2014:2019),
           mapping = aes(x = len, y = prop), alpha = 0.35, color = "black", fill = "darkgreen") +
  geom_line(all_srv_ll_plot_df %>% filter(type == "pred", sex == "Male",
                                          year %in% 2014:2019), 
            mapping = aes(x = len, y = prop, color = Model, group = Model), size = 0.85) +
  facet_wrap(~year) +
  theme_reg() +
  labs(x = "Length", y = "Male Proportions", fill = "Sex", color = "Model")

dev.off()

# Averlen fits
pdf(here("figs", "Model Comparison", "SurveyLL_len_Avg.pdf"), height = 5)
ggplot() +
  geom_col(all_avg_srv_ll_plot_df %>% filter(type == "obs", Model  == "Growth SQ"),
           mapping = aes(x = len, y = prop), alpha = 0.25, color = "black", fill = "orange") +
  geom_line(all_avg_srv_ll_plot_df %>% filter(type == "pred"), 
            mapping = aes(x = len, y = prop, color = Model, group = Model), size = 0.85) +
  theme_reg() +
  facet_wrap(~sex) +
  labs(x = "Length", y = "Proportions", fill = "Sex", Model = "Sex")
dev.off()

# Trawl Size --------------------------------------------------------------

# Growth Vary
# Get srv trawl lens
gv_obs_lc_srv7_f = data.frame(reshape2::melt(grwth_vary_sab_curr$olc.srv7.f), type = "obs") %>% mutate(sex = "Female")
gv_pred_lc_srv7_f = data.frame(reshape2::melt(grwth_vary_sab_curr$elc.srv7.f), type = "pred") %>% mutate(sex = "Female")
gv_obs_lc_srv7_m = data.frame(reshape2::melt(grwth_vary_sab_curr$olc.srv7.m), type = "obs") %>% mutate(sex = "Male")
gv_pred_lc_srv7_m = data.frame(reshape2::melt(grwth_vary_sab_curr$elc.srv7.m), type = "pred") %>% mutate(sex = "Male")

# Put these into a dataframe
gv_lc_srv7 = rbind(gv_obs_lc_srv7_f, gv_pred_lc_srv7_f, gv_obs_lc_srv7_m, gv_pred_lc_srv7_m)
names(gv_lc_srv7) = c("year", "len", "prop", "type", "sex")
gv_lc_srv7 = gv_lc_srv7 %>% mutate(Model = "Growth Vary") 

# Status-quo model
# Get srv trawl lens
sq_obs_lc_srv7_f = data.frame(reshape2::melt(grwth_sq_sab_curr$olc.srv7.f), type = "obs") %>% mutate(sex = "Female")
sq_pred_lc_srv7_f = data.frame(reshape2::melt(grwth_sq_sab_curr$elc.srv7.f), type = "pred") %>% mutate(sex = "Female")
sq_obs_lc_srv7_m = data.frame(reshape2::melt(grwth_sq_sab_curr$olc.srv7.m), type = "obs") %>% mutate(sex = "Male")
sq_pred_lc_srv7_m = data.frame(reshape2::melt(grwth_sq_sab_curr$elc.srv7.m), type = "pred") %>% mutate(sex = "Male")

# Put these into a dataframe
sq_lc_srv7 = rbind(sq_obs_lc_srv7_f, sq_pred_lc_srv7_f, sq_obs_lc_srv7_m, sq_pred_lc_srv7_m)
names(sq_lc_srv7) = c("year", "len", "prop", "type", "sex")
sq_lc_srv7 = sq_lc_srv7 %>% mutate(Model = "Growth SQ") # make plus group consistent color

# plot average fits as well
avg_gv_lc_srv7 = gv_lc_srv7 %>%
  group_by(len, type, sex, Model) %>% 
  summarize(prop = mean(prop))

avg_sq_lc_srv7 = sq_lc_srv7 %>%
  group_by(len, type, sex, Model) %>% 
  summarize(prop = mean(prop)) 

# Bind these together
all_srv_trawl_plot_df = rbind(gv_lc_srv7, sq_lc_srv7)
all_avg_srv_trawl_plot_df = rbind(avg_gv_lc_srv7, avg_sq_lc_srv7)


# Residuals ---------------------------------------------------------------

sq_lensrv7f = get_osa_res(obs = grwth_sq_sab_curr$olc.srv7.f, 
                          pred = grwth_sq_sab_curr$elc.srv7.f, 
                          iss = rep(20, nrow(grwth_sq_sab_curr$elc.srv7.f)), 
                          iter_wt = as.numeric(strsplit(grw_sq_ctl_wts[str_detect(grw_sq_ctl_wts, "#wt srv7 size comp female iter")], " ")[[1]][1]), 
                          index = seq(43,99,2),
                          drop_bin = 1)

sq_lensrv7f_plot = res_plot(data = sq_lensrv7f[[1]], 
                            comp_type = "Length", res_type = "(Growth SQ Survey Trawl size Comps Females)",
                            ymin = 41)

gv_lensrv7f = get_osa_res(obs = grwth_vary_sab_curr$olc.srv7.f, 
                          pred = grwth_vary_sab_curr$elc.srv7.f, 
                          iss = rep(20, nrow(grwth_vary_sab_curr$elc.srv7.f)), 
                          iter_wt = as.numeric(strsplit(grw_vary_ctl_wts[str_detect(grw_vary_ctl_wts, "#wt srv7 size comp female iter")], " ")[[1]][1]), 
                          index = seq(43,99,2),
                          drop_bin = 1)

gv_lensrv7f_plot = res_plot(data = gv_lensrv7f[[1]], 
                            comp_type = "Length", res_type = "(Growth Vary Survey Trawl size Comps Females)",
                            ymin = 43)

sq_srv7m = get_osa_res(obs = grwth_sq_sab_curr$olc.srv7.m, 
                       pred = grwth_sq_sab_curr$elc.srv7.m, 
                       iss = rep(20, nrow(grwth_sq_sab_curr$elc.srv7.m)), 
                       iter_wt = as.numeric(strsplit(grw_sq_ctl_wts[str_detect(grw_sq_ctl_wts, "#wt srv7 size comp male iter")], " ")[[1]][1]), 
                       index = seq(43,99,2),
                       drop_bin = 1)

sq_srv7m_plot = res_plot(data = sq_srv7m[[1]], 
                         comp_type = "size", res_type = "(Growth SQ Survey Trawl size Comps Males)",
                         ymin = 43)

gv_srv7m = get_osa_res(obs = grwth_vary_sab_curr$olc.srv7.m, 
                       pred = grwth_vary_sab_curr$elc.srv7.m, 
                       iss = rep(20, nrow(grwth_vary_sab_curr$elc.srv7.m)), 
                       iter_wt = as.numeric(strsplit(grw_vary_ctl_wts[str_detect(grw_vary_ctl_wts, "#wt srv7 size comp male iter")], " ")[[1]][1]), 
                       index = seq(43,99,2),
                       drop_bin = 1)

gv_srv7m_plot = res_plot(data = gv_srv7m[[1]], 
                         comp_type = "size", res_type = "(Growth Vary Survey Trawl size Comps Males)",
                         ymin = 43)

pdf(here("figs", "Model Comparison", "SurveyTrawl_size_Resid.pdf"), width = 15)
ggpubr::ggarrange(sq_lensrv7f_plot, gv_lensrv7f_plot, sq_srv7m_plot, gv_srv7m_plot)
dev.off()

### Fits --------------------------------------------------------------------

# Plot this out!
# fits across years
pdf(here("figs", "Model Comparison", "SurveyTrawl_len_Years.pdf"), width = 13, height = 10)
ggplot() +
  geom_col(all_srv_trawl_plot_df %>% filter(type == "obs", Model == "Growth SQ", sex == "Female"),
           mapping = aes(x = len, y = prop), alpha = 0.25, color = "black", fill = "orange") +
  geom_line(all_srv_trawl_plot_df %>% filter(type == "pred", sex == "Female"), 
            mapping = aes(x = len, y = prop, color = Model, group = Model), size = 0.85) +
  facet_wrap(~year) +
  theme_reg() +
  labs(x = "Length", y = "Female Proportions", fill = "Sex", color = "Model")

ggplot() +
  geom_col(all_srv_trawl_plot_df %>% filter(type == "obs", Model == "Growth SQ", sex == "Male"),
           mapping = aes(x = len, y = prop), alpha = 0.25, color = "black", fill = "orange") +
  geom_line(all_srv_trawl_plot_df %>% filter(type == "pred", sex == "Male"), 
            mapping = aes(x = len, y = prop, color = Model, group = Model), size = 0.85) +
  facet_wrap(~year) +
  theme_reg() +
  labs(x = "Length", y = "Male Proportions", fill = "Sex", color = "Model")
ggplot() +
  geom_col(all_srv_trawl_plot_df %>% filter(type == "obs", Model == "Growth SQ", sex == "Female",
                                            year %in% 2014:2019),
           mapping = aes(x = len, y = prop), alpha = 0.35, color = "black", fill = "darkgreen") +
  geom_line(all_srv_trawl_plot_df %>% filter(type == "pred", sex == "Female",
                                             year %in% 2014:2019), 
            mapping = aes(x = len, y = prop, color = Model, group = Model), size = 0.85) +
  facet_wrap(~year) +
  theme_reg() +
  labs(x = "Length", y = "Female Proportions", fill = "Sex", color = "Model")

ggplot() +
  geom_col(all_srv_trawl_plot_df %>% filter(type == "obs", Model == "Growth SQ", sex == "Male",
                                            year %in% 2014:2019),
           mapping = aes(x = len, y = prop), alpha = 0.35, color = "black", fill = "darkgreen") +
  geom_line(all_srv_trawl_plot_df %>% filter(type == "pred", sex == "Male",
                                             year %in% 2014:2019), 
            mapping = aes(x = len, y = prop, color = Model, group = Model), size = 0.85) +
  facet_wrap(~year) +
  theme_reg() +
  labs(x = "Length", y = "Male Proportions", fill = "Sex", color = "Model")

dev.off()

# Averlen fits
pdf(here("figs", "Model Comparison", "SurveyTrawl_len_Avg.pdf"), height = 5)
ggplot() +
  geom_col(all_avg_srv_trawl_plot_df %>% filter(type == "obs", Model  == "Growth SQ"),
           mapping = aes(x = len, y = prop), alpha = 0.25, color = "black", fill = "orange") +
  geom_line(all_avg_srv_trawl_plot_df %>% filter(type == "pred"), 
            mapping = aes(x = len, y = prop, color = Model, group = Model), size = 0.85) +
  theme_reg() +
  facet_wrap(~sex) +
  labs(x = "Length", y = "Proportions", fill = "Sex", Model = "Sex")
dev.off()



# Compare SSB -------------------------------------------------------------
gv_par <- R2admb::read_pars(fn = here(model_path, "Prop_Within_Growth_Vary","tem"))
sq_par <- R2admb::read_pars(fn = here(model_path, "Prop_Within_Growth_SQ","tem"))

# Get par estiamtes for terminaly year ssb
gv_ssb <- par_coeff(gv_par, "ssbsd", type = "Growth Vary") %>% mutate(Year = 1960:2023)
sq_ssb <- par_coeff(sq_par, "ssbsd", type = "Growth SQ") %>% mutate(Year = 1960:2023)

# Get par estiamtes for projections too
gv_ssb_proj <- par_coeff(gv_par, "spawn_biom_proj", type = "Growth Vary") %>% mutate(Year = 2024:2038)
sq_ssb_proj <- par_coeff(sq_par, "spawn_biom_proj", type = "Growth SQ") %>% mutate(Year = 2024:2038)

# Get b40 as well
gv_b40 <- par_coeff(gv_par, "B40", type = "Growth Vary") 
sq_b40 <- par_coeff(sq_par, "B40", type = "Growth SQ") 
b40_df = rbind(gv_b40, sq_b40)

# residual munging
ssb_df = rbind(gv_ssb, sq_ssb, gv_ssb_proj, sq_ssb_proj) %>% 
  rename(ssb = coeff) # rename
ssb_df = ssb_df %>% left_join(b40_df %>% rename(b40 = coeff) %>% 
                                select(type, b40), by = "type") %>% 
  mutate(Bratio = ssb/b40) %>% 
  mutate(Model_Comps = ifelse(str_detect(type, "Partial"), "Partial Length", "Full"))

pdf(here("figs", "Model Comparison", "SSB_Comparison.pdf"), width = 13)
ggplot(ssb_df %>% filter(Year <= 2023), 
       aes(x = Year, y = ssb, ymin = downr, ymax = upr, fill = type)) +
  geom_line(aes(color = type), size = 1.5) +
  geom_ribbon(alpha = 0.3) +
  geom_hline(aes(yintercept = b40, color = type), lty = 2, size = 1.5) +
  geom_text(label = "B40%", aes(x = 2013, y = 135), size = 13, check_overlap = TRUE) +
  theme_reg() +
  theme(legend.position = "top") +
  labs(x = "Year", y = "SSB (kt)", color = "Model", fill = "Model")
dev.off()

# Compare SSB vs B40 ------------------------------------------------------
pdf(here("figs", "Model Comparison", "Bratio_Comparison.pdf"), width = 13)
ggplot(ssb_df %>% filter(Year >= 1996, Year<= 2023), 
       aes(x = Year, y = Bratio, color = type)) +
  geom_line(size = 1.5, alpha = 0.55) +
  geom_point(size = 4) +
  geom_hline(yintercept = 1, lty = 2, size = 1.3) +
  theme_reg() +
  theme(legend.position = "top") +
  labs(x = "Year", y = "SSB/B40%", color = "Model", fill = "Model")
dev.off()

# Compare Recruitment -----------------------------------------------------
gv_rec <- par_coeff(gv_par, "pred_rec", type = "Growth Vary") %>% mutate(Year = 1960:2023)
sq_rec <- par_coeff(sq_par, "pred_rec", type = "Growth SQ") %>%  mutate(Year = 1960:2023)

# Get mean recruitment as well
gv_meanrec <- par_coeff(gv_par, "log_mean_rec", type = "Growth Vary") 
sq_meanrec <- par_coeff(sq_par, "log_mean_rec", type = "Growth SQ") 

meanrec_df = rbind(gv_meanrec, sq_meanrec)
rec_df = rbind(gv_rec, sq_rec) %>%  mutate(downr = ifelse(downr <=0, 0, downr)) %>% rename(rec = coeff)
rec_df = rec_df %>% left_join(meanrec_df %>% rename(meanrec = coeff) %>% mutate(meanrec = exp(meanrec)) %>% 
                                select(type, meanrec), by = "type") %>% 
  mutate(Model_Comps = ifelse(str_detect(type, "Partial"), "Partial Length", "Full"))

pdf(here("figs", "Model Comparison", "Recruitment_Comparison.pdf"), width = 13)
ggplot(rec_df %>% filter(Year != 2023), aes(x = Year - 2, y = rec, ymin = downr, ymax = upr, fill = type)) +
  geom_line(aes(color = type), size = 1.2) +
  geom_ribbon(alpha = 0.4) +
  geom_hline(aes(yintercept = meanrec, color = type), lty = 2, size = 1.3) +
  theme_reg() +
  theme(legend.position = "top") +
  ylim(0, NA) +
  labs(x = "Year", y = "Age-2 Recruitment (millions)", color = "Model", fill = "Model")
dev.off()


# Selectivity Comparisons -------------------------------------------------
# names for ordering
sel_rename = c("Derby fishery F","Derby fishery M","Trawl fishery F" ,"Trawl fishery M" ,"IFQ fishery F" ,"IFQ fishery M", "IFQ Recent fishery F" ,
               "IFQ Recent fishery M", "US LL survey F","US LL survey M", "US LL Recent survey F","US LL Recent survey M",
               "Cooperative LL survey F" ,"Cooperative LL survey M" ,"GOA trawl survey F","GOA trawl survey M" ) # selectivity names

# extract unique selectivity names
selex_names = unique(names(sq_par$coefficients)[str_detect(names(sq_par$coefficients), "sel")])
all_sel_df = data.frame()
# Loop through to extract these quantities
for(i in 1:length(selex_names)) {
  # status-quo model
  sq_sel <- par_coeff(sq_par, selex_names[i], type = "Growth SQ") %>% 
    mutate(Age = 2:31, sel = selex_names[i]) %>% 
    mutate(downr = ifelse(downr <=0, 0, downr), upr = ifelse(upr >= 1, 1, upr)) 
  # time-varying model
  gv_sel <- par_coeff(gv_par, selex_names[i], type = "Growth Vary") %>% 
    mutate(Age = 2:31, sel = selex_names[i]) %>% 
    mutate(downr = ifelse(downr <=0, 0, downr), upr = ifelse(upr >= 1, 1, upr)) 
  all_sel_df = rbind(all_sel_df, sq_sel, gv_sel)
} # end i

# Rename selectivity stuff
all_sel_df = all_sel_df %>% 
  mutate(sel = factor(sel, labels = sel_rename)) 

pdf(here("figs", "Model Comparison", "Selectivity_Comparison.pdf"),  height = 10, width = 10)
ggplot(all_sel_df %>% filter(
  str_detect(sel, "trawl") | 
    str_detect(sel, "Trawl")),
  aes(x = Age, y = coeff, ymin = downr, ymax = upr, fill = type)) +
  geom_line(aes(color = type), size = 1.2, alpha = 0.5) +
  geom_point(aes(color = type), size = 1.5) +
  geom_ribbon(alpha = 0.3) +
  theme_reg() +
  theme(legend.position = "top") +
  facet_wrap(~sel) +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Year", y = "Proportion Selected", color = "Model", fill = "Model")

ggplot(all_sel_df %>% filter(
  str_detect(sel, "Derby") |
    str_detect(sel, "Recent")),
  aes(x = Age, y = coeff, ymin = downr, ymax = upr, fill = type)) +
  geom_line(aes(color = type), size = 1.2, alpha = 0.5) +
  geom_point(aes(color = type), size = 1.5) +
  geom_ribbon(alpha = 0.3) +
  theme_reg() +
  theme(legend.position = "top") +
  facet_wrap(~sel) +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Year", y = "Proportion Selected", color = "Model", fill = "Model")
dev.off()


# Reference Points --------------------------------------------------------

gv_abc<- par_coeff(gv_par, "ABC", type = "Growth Vary")
sq_abc <- par_coeff(sq_par, "ABC", type = "Growth SQ") 
abc_df = rbind(gv_abc, sq_abc)
ref_df = rbind(b40_df, abc_df)

pdf(here("figs", "Model Comparison", "ReferencePoint_Comparison.pdf"), width = 8)
ggplot(ref_df, aes(x = type, y = coeff, ymin = downr, ymax = upr, color = type)) +
  geom_pointrange() +
  facet_wrap(~coeff_type, scales = "free") +
  labs(x = "Model", y = "Value") +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_reg() +
  theme(legend.position = "none")

dev.off()
# Index Fits --------------------------------------------------------------

idx_names = c("Domestic LL Survey Relative Population Weight", "Japanese LL Survey Relative Population Weight",
              "Domestic LL Survey Relative Population Numbers", "Japanese LL Survey Relative Population Numbers",
              "Domestic Fishery CPUE Index", "Japanese Fishery CPUE Index", "GOA Trawl Survey Biomass (kt)") # index names

# extract out index data
gv_idx_df = data.frame()
gv_idx_list = grwth_vary_sab_curr[str_detect(names(grwth_vary_sab_curr), "obssrv")]
sq_idx_df = data.frame()
sq_idx_list = grwth_sq_sab_curr[str_detect(names(grwth_sq_sab_curr), "obssrv")]

# Loop through to extract stuff out
for(i in 1:length(gv_idx_list)) {
  # extract out index dataframe (just extracting out components from a list) - growth vary model
  gv_idx_tmp = data.frame(year = as.numeric(rownames(gv_idx_list[[i]])), obs = gv_idx_list[[i]][[1]],
                          lci = gv_idx_list[[i]][[2]], uci = gv_idx_list[[i]][[3]],
                          pred = gv_idx_list[[i]][[4]], type = idx_names[i], Model = "Growth Vary")
  # extract out index dataframe (just extracting out components from a list) - sq model
  sq_idx_tmp = data.frame(year = as.numeric(rownames(sq_idx_list[[i]])), obs = sq_idx_list[[i]][[1]],
                          lci = sq_idx_list[[i]][[2]], uci = sq_idx_list[[i]][[3]],
                          pred = sq_idx_list[[i]][[4]], type = idx_names[i], Model = "Growth SQ")
  # bind together
  gv_idx_df = rbind(gv_idx_df, gv_idx_tmp)
  sq_idx_df = rbind(sq_idx_df, sq_idx_tmp)
} # end i loop

idx_df = rbind(gv_idx_df, sq_idx_df)
# relvel by index name
idx_df$type = factor(idx_df$type, levels = idx_names)

pdf(here("figs", "Model Comparison", "Index_Fits.pdf"), height = 11, width = 15)
ggplot(idx_df) +
  geom_line(mapping = aes(x = year, y = pred, col = Model), size = 1.5) +
  geom_pointrange(mapping = aes(x = year, y = obs, ymin = lci, ymax = uci), 
                  alpha = 0.35, size = 0.75, col = "black") +
  facet_wrap(~type, scales = "free") +
  scale_y_continuous(limits = c(0, NA)) +  # NA for no upper limit
  labs(x = "Year", y = "Index") +
  theme_reg() +
  theme(legend.position = "top")
dev.off()


# Fishing Mortality Rates -------------------------------------------------
gv_f1<- par_coeff(gv_par, "Fmort_fish1", type = "Growth Vary") %>% 
  mutate(Year = 1960:2023)
sq_f1 <- par_coeff(sq_par, "Fmort_fish1", type = "Growth SQ") %>% 
  mutate(Year = 1960:2023)
gv_f3<- par_coeff(gv_par, "Fmort_fish3", type = "Growth Vary") %>% 
  mutate(Year = 1960:2023)
sq_f3 <- par_coeff(sq_par, "Fmort_fish3", type = "Growth SQ") %>% 
  mutate(Year = 1960:2023)
fmort_all = rbind(gv_f1, sq_f1, sq_f3, gv_f3) %>% 
  mutate(coeff_type = case_when(
    coeff_type == "Fmort_fish1" ~ "Fixed-Gear Fishery",
    coeff_type == "Fmort_fish3" ~ "Trawl Gear Fishery",
  ))

pdf(here("figs", "Model Comparison", "FishingMortality.pdf"), width = 10)
ggplot(fmort_all, aes(x = Year, y = coeff, ymin = downr, ymax = upr, fill = type)) +
  geom_line(aes(color = type), size = 1) +
  geom_ribbon(alpha = 0.3) +
  facet_wrap(~coeff_type) +
  theme_reg() +
  theme(legend.position = "top") +
  labs(x = "Year", y = "Annual Fishing Mortality Multiplier", color = "Model", fill = "Model")
dev.off()


# Fits to catch -----------------------------------------------------------
# Fixed Gear Catch
gv_pred_fixed_catch = as.numeric(strsplit(grwth_vary_rep[str_detect(grwth_vary_rep, "Pred_Catch Fixed Gear")], " ")[[1]][-c(1:4)])
gv_obs_fixed_catch = as.numeric(strsplit(grwth_vary_rep[str_detect(grwth_vary_rep, "Obs_Catch Fixed Gear")], " ")[[1]][-c(1:4)]) 
gv_fixed_catch = data.frame(Obs = gv_obs_fixed_catch, Pred = gv_pred_fixed_catch, Model = "Growth Vary", Year = 1960:2023)
sq_pred_fixed_catch = as.numeric(strsplit(grwth_sq_rep[str_detect(grwth_sq_rep, "Pred_Catch Fixed Gear")], " ")[[1]][-c(1:4)])
sq_obs_fixed_catch  = as.numeric(strsplit(grwth_sq_rep[str_detect(grwth_sq_rep, "Obs_Catch Fixed Gear")], " ")[[1]][-c(1:4)])
sq_fixed_catch = data.frame(Obs = sq_obs_fixed_catch, Pred = sq_pred_fixed_catch, Model = "Growth SQ", Year = 1960:2023)

# Trawl Gear catch
gv_pred_trawl_catch = as.numeric(strsplit(grwth_vary_rep[str_detect(grwth_vary_rep, "Pred catch trawl")], " ")[[1]][-c(1:4)])
gv_obs_trawl_catch = as.numeric(strsplit(grwth_vary_rep[str_detect(grwth_vary_rep, "Obs_Catch Trawl")], " ")[[1]][-c(1:4)]) 
gv_trawl_catch = data.frame(Obs = gv_obs_trawl_catch, Pred = gv_pred_trawl_catch, Model = "Growth Vary", Year = 1960:2023)
sq_pred_trawl_catch = as.numeric(strsplit(grwth_sq_rep[str_detect(grwth_sq_rep, "Pred catch trawl")], " ")[[1]][-c(1:4)])
sq_obs_trawl_catch  = as.numeric(strsplit(grwth_sq_rep[str_detect(grwth_sq_rep, "Obs_Catch Trawl")], " ")[[1]][-c(1:4)])
sq_trawl_catch = data.frame(Obs = sq_obs_trawl_catch, Pred = sq_pred_trawl_catch, Model = "Growth SQ", Year = 1960:2023)

# Cleaning up
fixed_catch = rbind(gv_fixed_catch, sq_fixed_catch) %>% mutate(Gear = 'Fixed Gear')
trawl_catch = rbind(gv_trawl_catch, sq_trawl_catch) %>% mutate(Gear = 'Trawl Gear')
all_catch = rbind(fixed_catch, trawl_catch)

pdf(here("figs", "Model Comparison", "CatchFits.pdf"), width = 10)
ggplot(all_catch) +
  geom_line(mapping = aes(x = Year, y = Pred, color = Model), size = 1.3, alpha = 0.85) +
  geom_point(mapping = aes(x = Year, y = Obs), size = 2, alpha = 0.5) +
  theme_reg() +
  facet_wrap(~Gear) +
  theme(legend.position = "top") +
  labs(x = "Year", y = "Catch (kt)", color = "Model", fill = "Model")

dev.off()

