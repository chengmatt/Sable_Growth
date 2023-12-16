# Purpose: To employ OSA residuals, and compare to Pearson residuals
# Creator: Matthew LH. Cheng
# Date 9/27/23

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

# Installation
# https://github.com/fishfollower/compResidual
# TMB:::install.contrib("https://github.com/vtrijoulet/OSA_multivariate_dists/archive/main.zip")
# remotes::install_github("fishfollower/compResidual/compResidual", force=TRUE)
# Set up ------------------------------------------------------------------
library(here)
library(tidyverse)
library(compResidual) 
library(cowplot)

# Read in ctl file to get weights from multinomial
source(here("2. Base (23.5)_final model", "R code", "residual_functions_MC.r"))
ctl_wts <- readLines(here("2. Base (23.5)_final model", "tem.ctl"))
sab_curr <- dget(here("2. Base (23.5)_final model", "tem.rdat")) # read in observed and predicted comps

iss = 20 # input sample size
plot_lens = seq(41, 99, 2)  # lens for plotting
plot_ages = seq(2, 31, 1) # ages for plotting
osa_lens = seq(43, 99, 2)  # start at 43 bc dropping first bin
osa_ages = seq(3, 31, 1) # start at 3 bc dropping first bin

pdf(here("2. Base (23.5)_final model", "Results", "Residual_Comparison.pdf"))
# Fishery LL Age Comps ----------------------------------------------------
# get effective sample size weights
iter_wt = parse_number(ctl_wts[str_detect(ctl_wts, "#wt fish1 age comp iter")]) 

# Get fishery osa_ages OSA
fish1_age_osa = get_osa_res(obs = sab_curr$oac.fish1, pred = sab_curr$eac.fish1, 
                               iss = iss, iter_wt = iter_wt, index = osa_ages, drop_bin = 1)

# Get fishery osa_ages Pearson
fish1_age_pearson = get_pearson_res(iter_wt = iter_wt, iss = iss,
                                       obs = sab_curr$oac.fish1, pred = sab_curr$eac.fish1)

# Plot OSA and Pearson Comparison
osa_plot = res_plot(data = fish1_age_osa[[1]], 
                    comp_type = "Age", res_type = "OSA w/o first bin (Fishery LL Age Comps)",
                    ymin = min(plot_ages))
pearson_plot = res_plot(data = fish1_age_pearson, comp_type = "Age", 
                        res_type = "Pearson (Fishery LL Age Comps)",
                        ymin = min(plot_ages))
plot_grid(osa_plot , pearson_plot, ncol = 1) 

# Domestic Survey LL Age Comps ----------------------------------------------------
# get effective sample size weights
iter_wt = parse_number(ctl_wts[str_detect(ctl_wts, "#wt surv1 age comp iter")]) 

# Get domestic survey osa_ages OSA
srv1_age_osa = get_osa_res(obs = sab_curr$oac.srv1, pred = sab_curr$eac.srv1, 
                            iss = iss, iter_wt = iter_wt, index = osa_ages, drop_bin = 1)

# Get domestic survey osa_ages Pearson
srv1_age_pearson = get_pearson_res(iter_wt = iter_wt, iss = iss,
                                    obs = sab_curr$oac.srv1, pred = sab_curr$eac.srv1)

# Plot OSA and Pearson Comparison
osa_plot = res_plot(data = srv1_age_osa[[1]], comp_type = "Age", 
                    res_type = "OSA w/o first bin (Domestic LL Survey Age Comps)",
                    ymin = min(plot_ages))
pearson_plot = res_plot(data = srv1_age_pearson, comp_type = "Age",
                        res_type = "Pearson (Domestic LL Survey Age Comps)",
                        ymin = min(plot_ages))
plot_grid(osa_plot , pearson_plot, ncol = 1) 

# Japanese Survey LL Age Comps ----------------------------------------------------
# get effective sample size weights
iter_wt = parse_number(ctl_wts[str_detect(ctl_wts, "#wt surv2 age comp iter")]) 

# Get Japanese Survey osa_ages OSA
srv2_age_osa = get_osa_res(obs = sab_curr$oac.srv2, pred = sab_curr$eac.srv2, 
                           iss = iss, iter_wt = iter_wt, index = osa_ages, drop_bin = 1)

# Get Japanese Survey osa_ages Pearson
srv2_age_pearson = get_pearson_res(iter_wt = iter_wt, iss = iss,
                                   obs = sab_curr$oac.srv2, pred = sab_curr$eac.srv2)

# Plot OSA and Pearson Comparison
osa_plot = res_plot(data = srv2_age_osa[[1]], comp_type = "Age", 
                    res_type = "OSA w/o first bin (Japanese LL Survey Age Comps)",
                    ymin = min(plot_ages))
pearson_plot = res_plot(data = srv2_age_pearson, comp_type = "Age", 
                        res_type = "Pearson (Japanese LL Survey Age Comps)",
                        ymin = min(plot_ages))
plot_grid(osa_plot , pearson_plot, ncol = 1) 

# Domestic LL Survey (Length Females) -----------------------------------------------------
# get effective sample size weights
iter_wt = parse_number(ctl_wts[str_detect(ctl_wts, "#wt surv1 size comp female iter")]) 

# Get Survey length females residuals (OSA)
srv1_fem_len_osa = get_osa_res(obs = sab_curr$olc.srv1.f, pred = sab_curr$elc.srv1.f, 
                           iss = iss, iter_wt = iter_wt, index = osa_lens, drop_bin = 1)

# Get Survey length females residuals (Pearson)
srv1_fem_len_pearson = get_pearson_res(iter_wt = iter_wt, iss = iss,
                                       obs = sab_curr$olc.srv1.f, pred = sab_curr$elc.srv1.f)

# Plot OSA and Pearson Comparison
osa_plot = res_plot(data = srv1_fem_len_osa[[1]], 
                    comp_type = "Length", res_type = "OSA w/o first bin (Domestic LL Survey Female Length Comps)",
                    ymin = min(plot_lens))
pearson_plot = res_plot(data = srv1_fem_len_pearson,
                        comp_type = "Length", res_type = "Pearson (Domestic LL Survey Female Length Comps)",
                        ymin = min(plot_lens))
plot_grid(osa_plot, pearson_plot, ncol = 1)

# Domestic LL Survey (Length Males) -----------------------------------------------------
# get effective sample size weights
iter_wt = parse_number(ctl_wts[str_detect(ctl_wts, "#wt surv1 size comp male iter")]) 

# Get Survey length males residuals (OSA)
srv1_male_len_osa = get_osa_res(obs = sab_curr$olc.srv1.m, pred = sab_curr$elc.srv1.m, 
                               iss = iss, iter_wt = iter_wt, index = osa_lens, drop_bin = 1)

# Get Survey length males residuals (Pearson)
srv1_male_len_pearson = get_pearson_res(iter_wt = iter_wt, iss = iss,
                                       obs = sab_curr$olc.srv1.m, pred = sab_curr$elc.srv1.m)

# Plot OSA and Pearson Comparison
osa_plot = res_plot(data = srv1_male_len_osa[[1]], 
                    comp_type = "Length", res_type = "OSA w/o first bin (Domestic LL Survey Male Length Comps)",
                    ymin = min(plot_lens))
pearson_plot = res_plot(data = srv1_male_len_pearson, 
                        comp_type = "Length", res_type = "Pearson (Domestic LL Survey Male Length Comps)",
                        ymin = min(plot_lens))
plot_grid(osa_plot, pearson_plot, ncol = 1)

# Japanese LL Survey (Length Females) -----------------------------------------------------
# get effective sample size weights
iter_wt = parse_number(ctl_wts[str_detect(ctl_wts, "#wt surv2 size comp female iter")]) 

# Get Survey length females residuals (OSA)
srv2_fem_len_osa = get_osa_res(obs = sab_curr$olc.srv2.f, pred = sab_curr$elc.srv2.f, 
                               iss = iss, iter_wt = iter_wt, index = osa_lens, drop_bin = 1)

# Get Survey length females residuals (Pearson)
srv2_fem_len_pearson = get_pearson_res(iter_wt = iter_wt, iss = iss,
                                       obs = sab_curr$olc.srv2.f, pred = sab_curr$elc.srv2.f)

# Plot OSA and Pearson Comparison
osa_plot = res_plot(data = srv2_fem_len_osa[[1]], 
                    comp_type = "Length", res_type = "OSA w/o first bin (Japanese LL Survey Female Length Comps)",
                    ymin = min(plot_lens))
pearson_plot = res_plot(data = srv2_fem_len_pearson, 
                        comp_type = "Length", res_type = "Pearson (Japanese LL Survey Female Length Comps)",
                        ymin = min(plot_lens))
plot_grid(osa_plot, pearson_plot, ncol = 1)

# Japanese LL Survey (Length Males) -----------------------------------------------------
# get effective sample size weights
iter_wt = parse_number(ctl_wts[str_detect(ctl_wts, "#wt surv2 size comp male iter")]) 

# Get Survey length males residuals (OSA)
srv2_male_len_osa = get_osa_res(obs = sab_curr$olc.srv2.m, pred = sab_curr$elc.srv2.m, 
                                iss = iss, iter_wt = iter_wt, index = osa_lens, drop_bin = 1)

# Get Survey length males residuals (Pearson)
srv2_male_len_pearson = get_pearson_res(iter_wt = iter_wt, iss = iss,
                                        obs = sab_curr$olc.srv2.m, pred = sab_curr$elc.srv2.m)

# Plot OSA and Pearson Comparison
osa_plot = res_plot(data = srv2_male_len_osa[[1]], 
                    comp_type = "Length", res_type = "OSA w/o first bin (Japanese LL Survey Male Length Comps)",
                    ymin = min(plot_lens))
pearson_plot = res_plot(data = srv2_male_len_pearson, 
                        comp_type = "Length", res_type = "Pearson (Japanese LL Survey Male Length Comps)",
                        ymin = min(plot_lens))
plot_grid(osa_plot, pearson_plot, ncol = 1)

# Domestic LL Fishery (Length Females) -----------------------------------------------------
# get effective sample size weights
iter_wt = parse_number(ctl_wts[str_detect(ctl_wts, "#wt fish1 size comp female iter")]) 

# Get fishery length females residuals (OSA)
fish1_fem_len_osa = get_osa_res(obs = sab_curr$olc.fish1.f, pred = sab_curr$elc.fish1.f, 
                               iss = iss, iter_wt = iter_wt, index = osa_lens, drop_bin = 1)

# Get fishery length females residuals (Pearson)
fish1_fem_len_pearson = get_pearson_res(iter_wt = iter_wt, iss = iss,
                                       obs = sab_curr$olc.fish1.f, pred = sab_curr$elc.fish1.f)

# Plot OSA and Pearson Comparison
osa_plot = res_plot(data = fish1_fem_len_osa[[1]], 
                    comp_type = "Length", res_type = "OSA w/o first bin (Fishery LL Female Length Comps)",
                    ymin = min(plot_lens))
pearson_plot = res_plot(data = fish1_fem_len_pearson, 
                        comp_type = "Length", res_type = "Pearson (Fishery LL Female Length Comps)",
                        ymin = min(plot_lens))
plot_grid(osa_plot, pearson_plot, ncol = 1)

# Domestic LL Fishery (Length Males) -----------------------------------------------------
# get effective sample size weights
iter_wt = parse_number(ctl_wts[str_detect(ctl_wts, "#wt fish1 size comp male iter")]) 

# Get fishery length males residuals (OSA)
fish1_male_len_osa = get_osa_res(obs = sab_curr$olc.fish1.m, pred = sab_curr$elc.fish1.m, 
                                iss = iss, iter_wt = iter_wt, index = osa_lens, drop_bin = 1)

# Get fishery length males residuals (Pearson)
fish1_male_len_pearson = get_pearson_res(iter_wt = iter_wt, iss = iss,
                                        obs = sab_curr$olc.fish1.m, pred = sab_curr$elc.fish1.m)

# Plot OSA and Pearson Comparison
osa_plot = res_plot(data = fish1_male_len_osa[[1]], 
                    comp_type = "Length", res_type = "OSA w/o first bin (Fishery LL Male Length Comps)",
                    ymin = min(plot_lens))
pearson_plot = res_plot(data = fish1_male_len_pearson, 
                        comp_type = "Length", res_type = "Pearson (Fishery LL Male Length Comps)",
                        ymin = min(plot_lens))
plot_grid(osa_plot, pearson_plot, ncol = 1)

# Domestic Trawl Fishery (Length Females) -----------------------------------------------------
# get effective sample size weights
iter_wt = parse_number(ctl_wts[str_detect(ctl_wts, "#wt fish 3 size comp female iter")]) 

# Get fishery length females residuals (OSA)
fish3_fem_len_osa = get_osa_res(obs = sab_curr$olc.fish3.f, pred = sab_curr$elc.fish3.f, 
                                iss = iss, iter_wt = iter_wt, index = osa_lens, drop_bin = 1)

# Get fishery length females residuals (Pearson)
fish3_fem_len_pearson = get_pearson_res(iter_wt = iter_wt, iss = iss,
                                       obs = sab_curr$olc.fish3.f, pred = sab_curr$elc.fish3.f)

# Plot OSA and Pearson Comparison
osa_plot = res_plot(data = fish3_fem_len_osa[[1]], 
                    comp_type = "Length", res_type = "OSA w/o first bin (Fishery Trawl Female Length Comps)",
                    ymin = min(plot_lens))
pearson_plot = res_plot(data = fish3_fem_len_pearson, 
                        comp_type = "Length", res_type = "Pearson (Fishery Trawl Female Length Comps)",
                        ymin = min(plot_lens))
plot_grid(osa_plot, pearson_plot, ncol = 1)

# Domestic Trawl Fishery (Length Males) -----------------------------------------------------
# get effective sample size weights
iter_wt = parse_number(ctl_wts[str_detect(ctl_wts, "#wt fish 3 size comp male iter")]) 

# Get fishery length males residuals (OSA)
fish3_male_len_osa = get_osa_res(obs = sab_curr$olc.fish3.m, pred = sab_curr$elc.fish3.m, 
                                 iss = iss, iter_wt = iter_wt, index = osa_lens, drop_bin = 1)

# Get fishery length males residuals (Pearson)
fish3_male_len_pearson = get_pearson_res(iter_wt = iter_wt, iss = iss,
                                        obs = sab_curr$olc.fish3.m, pred = sab_curr$elc.fish3.m)

# Plot OSA and Pearson Comparison
osa_plot = res_plot(data = fish3_male_len_osa[[1]], 
                    comp_type = "Length", res_type = "OSA w/o first bin (Fishery Trawl Male Length Comps)",
                    ymin = min(plot_lens))
pearson_plot = res_plot(data = fish3_male_len_pearson, 
                        comp_type = "Length", res_type = "Pearson (Fishery Trawl Male Length Comps)",
                        ymin = min(plot_lens))
plot_grid(osa_plot, pearson_plot, ncol = 1)

# Domestic Trawl Survey (Length Females) -----------------------------------------------------
# get effective sample size weights
iter_wt = parse_number(ctl_wts[str_detect(ctl_wts, "#wt srv7 size comp female iter")]) 

# Get Survey length females residuals (OSA)
srv7_fem_len_osa = get_osa_res(obs = sab_curr$olc.srv7.f, pred = sab_curr$elc.srv7.f, 
                               iss = iss, iter_wt = iter_wt, index = osa_lens, drop_bin = 1)

# Get Survey length females residuals (Pearson)
srv7_fem_len_pearson = get_pearson_res(iter_wt = iter_wt, iss = iss,
                                       obs = sab_curr$olc.srv7.f, pred = sab_curr$elc.srv7.f)

# Plot OSA and Pearson Comparison
osa_plot = res_plot(data = srv7_fem_len_osa[[1]], 
                    comp_type = "Length", res_type = "OSA w/o first bin (Domestic Trawl Survey Female Length Comps)",
                    ymin = min(plot_lens))
pearson_plot = res_plot(data = srv7_fem_len_pearson, 
                        comp_type = "Length", res_type = "Pearson (Domestic Trawl Survey Female Length Comps)",
                        ymin = min(plot_lens))
plot_grid(osa_plot, pearson_plot, ncol = 1)

# Domestic Trawl Survey (Length Males) -----------------------------------------------------
# get effective sample size weights
iter_wt = parse_number(ctl_wts[str_detect(ctl_wts, "#wt srv7 size comp male iter")]) 

# Get Survey length males residuals (OSA)
srv7_male_len_osa = get_osa_res(obs = sab_curr$olc.srv7.m, pred = sab_curr$elc.srv7.m, 
                                iss = iss, iter_wt = iter_wt, index = osa_lens, drop_bin = 1)

# Get Survey length males residuals (Pearson)
srv7_male_len_pearson = get_pearson_res(iter_wt = iter_wt, iss = iss,
                                        obs = sab_curr$olc.srv7.m, pred = sab_curr$elc.srv7.m)

# Plot OSA and Pearson Comparison
osa_plot = res_plot(data = srv7_male_len_osa[[1]], 
                    comp_type = "Length", res_type = "OSA w/o first bin (Domestic Trawl Survey Male Length Comps)",
                    ymin = min(plot_lens))
pearson_plot = res_plot(data = srv7_male_len_pearson, 
                        comp_type = "Length", 
                        res_type = "Pearson (Domestic Trawl Survey Male Length Comps)",
                        ymin = min(plot_lens))
plot_grid(osa_plot, pearson_plot, ncol = 1)

dev.off()



















































