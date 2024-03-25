# Purpose: Plots for manuscript
# Creator: Matthew LH. Cheng (UAF-CFOS)
# 1/30/24

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

# Set up ------------------------------------------------------------------

library(here)
library(tidyverse)
library(ggpubr)

ages = 2:31
lengths = seq(41,99,2)

# Read in results
source(here("R", "functions.R"))
model_path <- here("ADMB_Models", "Final_Models")
sq_df <- dget(here(model_path, "Prop_Within_Growth_SQ", "tem.rdat"))
gv_df <- dget(here(model_path, "Prop_Within_Growth_Vary", "tem.rdat"))
gv_par <- R2admb::read_pars(fn = here(model_path, "Prop_Within_Growth_Vary","tem"))
sq_par <- R2admb::read_pars(fn = here(model_path, "Prop_Within_Growth_SQ","tem"))

# SSB and B40 -------------------------------------------------------------

# ABC stuff
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

# B40 Stuff
gv_b40_rat <- par_coeff(gv_par, "b40rat", type = "Growth Vary") %>% mutate(Year = 1960:2023)
sq_b40_rat <- par_coeff(sq_par, "b40rat", type = "Growth SQ") %>% mutate(Year = 1960:2023)
b40rat_df = rbind(gv_b40_rat, sq_b40_rat)

# ABC stuff
gv_abc<- par_coeff(gv_par, "ABC", type = "Growth Vary")
sq_abc <- par_coeff(sq_par, "ABC", type = "Growth SQ") 
abc_df = rbind(gv_abc, sq_abc)

# Recruitment Stuff
gv_rec <- par_coeff(gv_par, "pred_rec", type = "Growth Vary") %>% mutate(Year = 1960:2023)
sq_rec <- par_coeff(sq_par, "pred_rec", type = "Growth SQ") %>%  mutate(Year = 1960:2023)

# Get mean recruitment as well
gv_meanrec <- par_coeff(gv_par, "log_mean_rec", type = "Growth Vary") 
sq_meanrec <- par_coeff(sq_par, "log_mean_rec", type = "Growth SQ") 

meanrec_df = rbind(gv_meanrec, sq_meanrec)
rec_df = rbind(gv_rec, sq_rec) %>%  mutate(downr = ifelse(downr <=0, 0, downr)) %>% rename(rec = coeff)
rec_df = rec_df %>% left_join(meanrec_df %>% rename(meanrec = coeff) %>% mutate(meanrec = exp(meanrec)) %>% 
                                select(type, meanrec), by = "type")

## SSB, Rec, ABC plot ------------------------------------------------------

ssb_plot = ggplot(ssb_df %>% filter(Year <= 2023), aes(x = Year, y = ssb, color = type,
                                                       ymin = downr, ymax = upr, fill = type)) +
  geom_ribbon(alpha = 0.35, color = NA) +
  geom_line(aes(color = type), size = 1.3) +
  geom_hline(aes(yintercept = b40, color = type), lty = 2, size = 1.3) +
  theme_reg() +
  theme(legend.position = c(0.88, 0.85)) +
  scale_color_manual(values = c("#0072B2", "#D55E00")) +
  scale_fill_manual(values = c("#0072B2", "#D55E00")) +
  labs(x = "Year", y = "Spawning Biomass (kt)", color = "Model", fill = "Model")

rec_plot = ggplot(rec_df %>% filter(Year %in% c(1960:2023)), 
       aes(x = Year, y = rec, ymin = downr, ymax = upr, fill = type)) +
  geom_line(aes(color = type), size = 1.3) +
  geom_ribbon(alpha = 0.35) +
  geom_hline(aes(yintercept = meanrec, color = type), lty = 2, size = 1.3) +
  theme_reg() +
  scale_color_manual(values = c("#0072B2", "#D55E00")) +
  scale_fill_manual(values = c("#0072B2", "#D55E00")) +
  theme(legend.position = "none") +
  ylim(0, NA) +
  labs(x = "Year", y = "Age-2 Recruitment (millions)", color = "Model", fill = "Model")

abc_plot = ggplot(abc_df, aes(x = type, y = coeff, ymin = downr, ymax = upr, color = type)) +
  geom_pointrange(position = position_dodge(width = 0.2), lwd = 1.3, fatten = 15) +
  labs(x = "Model", y = "Acceptable Biological Catch (ABC)", color = "Model") +
  theme_reg() +
  scale_color_manual(values = c("#0072B2", "#D55E00")) +
  scale_fill_manual(values = c("#0072B2", "#D55E00")) +
  theme(legend.position = "none")

pdf(here("figs", "ms_figs", "stock_status.pdf"), width = 9, height = 8)
ab = ggarrange(ssb_plot, rec_plot, ncol = 1, labels = c('', "B"), 
               font.label = list(size = 25), hjust = -0.3, vjust = 1.3)
ggarrange(ab, abc_plot, widths = c(0.7, 0.3), labels = c('A', "C"),
          font.label = list(size = 25), hjust = -0.3, vjust = 1.3)
dev.off()

# Look at mean differences in recruitment
rec_df %>% 
  # mutate(Year = Year - 2) %>% 
  filter(Year %in% c(2016:2022)) %>% 
  select(Year, type, rec) %>% 
  pivot_wider(names_from = "type", values_from = "rec") %>% 
  mutate(Diff = (`Growth Vary` - `Growth SQ`)  ) %>% 
  summarize(mean = mean(Diff))

# Anomaly Plot ------------------------------------------------------------

# growth estimates here
waa_ts = read.csv(here("output", "WAA_Models", "Growth_estimates.csv"))
waa_sd = read.csv(here("output", "WAA_Models", "SDRep_estimates.csv"))
sd_all = rbind(waa_sd)

# Calculate WAA anomaly
const_waa_ts = waa_ts %>% filter(re_model == "constant", Year == 1996, peel == 0) %>% 
  select(Age, Mean, Sex) %>% rename(Ref = Mean)
anom_waa = waa_ts %>% filter(re_model == "3dar1") %>% 
  left_join(const_waa_ts, by = c("Age", "Sex")) %>% 
  mutate(Anom = (Mean-Ref) / Ref,
         Sex = ifelse(Sex == "female", "Female", "Male"))


pdf(here("figs", "ms_figs", "waa_anom.pdf"), width = 15, height = 7)
ggplot(anom_waa %>% filter(peel == 0), 
       aes(x = factor(Year), y = factor(Age), fill = Anom)) +
  geom_tile(alpha = 1) +
  scale_fill_gradient2(midpoint = mean(anom_waa$Anom), 
                       high = scales::muted("red"), 
                       low = scales::muted("blue")) +
  facet_wrap(~Sex) +
  scale_y_discrete(breaks = seq(3, 31, 3)) +
  scale_x_discrete(breaks = seq(1996, 2023, 5)) +
  theme_reg() +
  theme(legend.key.width = unit(3, "cm")) +
  labs(x = "Year", y = "Age", fill = "Weight Anomalies") +
  theme(legend.position = "top")
dev.off()


# Cohort Plot -------------------------------------------------------------

# create ages at which large recruitment events were observed
waa_ts = waa_ts %>% mutate(large_rec_age1 = (2014 - broodYear) + 2,
                           large_rec_age2 = (2016 - broodYear) + 2,
                           large_rec_age1 = ifelse(large_rec_age1 <= 1, NA, large_rec_age1),
                           large_rec_age2 = ifelse(large_rec_age2 <= 1, NA, large_rec_age2))

male_plot = waa_ts %>% filter(Sex == "Male",
                  broodYear %in% c(2010:2018), 
                  Year != 2023,
                  peel == 0) %>% 
  ggplot(aes(x = Age, y = Mean, color = re_model, 
             ymin = lwr_95, ymax = upr_95, fill = re_model), size = 1.3) +
  geom_line(size = 1.3) +
  geom_ribbon(alpha = 0.3, color = NA) +
  facet_wrap(~broodYear, scale = "free") +
  geom_vline(aes(xintercept = large_rec_age1), size = 0.9) +
  geom_vline(aes(xintercept = large_rec_age2), lty = 2, size = 0.9) +
  labs(x = "Age", y = "Male Weight (g)", fill = "Model", color = "Model") +
  theme_reg() +
  theme(legend.position = c(0.925, 0.085))  +
  scale_color_manual(values = c("#005AB5", "black")) +
  scale_fill_manual(values = c("#005AB5", "black")) 

female_plot = waa_ts %>% filter(Sex == "female",
                                broodYear %in% c(2010:2018), 
                                Year != 2023,
                  peel == 0) %>% 
  ggplot(aes(x = Age, y = Mean, color = re_model, 
             ymin = lwr_95, ymax = upr_95, fill = re_model), size = 1.3) +
  geom_line(size = 1.3) +
  geom_ribbon(alpha = 0.35, color = NA) +
  facet_wrap(~broodYear, scale = "free") +
  geom_vline(aes(xintercept = large_rec_age1), size = 0.9) +
  geom_vline(aes(xintercept = large_rec_age2), lty = 2, size = 0.9) +
  labs(x = "Age", y = "Female Weight (g)", fill = "Model", color = "Model") +
  scale_color_manual(values = c("#DC3220", "black")) +
  scale_fill_manual(values = c("#DC3220", "black")) +
  theme_reg() +
  theme(legend.position = c(0.925, 0.085)) 
  
pdf(here("figs", "ms_figs", "cohort_plot.pdf"), width = 17, height = 8)
ggarrange(female_plot, male_plot, ncol = 2, labels = c('A', "B"), 
          font.label = list(size = 25), hjust = -0.6)
dev.off()


# Compare input differences -----------------------------------------------------
# Size Transition Differences ---------------------------------------------
age_bins = 2:31
len_bins = seq(41, 99, 2)
al_mat_all_df = data.frame()
years = unique(waa_ts$Year)
n_proj_years = 1
# extract out length models
length_models = list.files(here("output", "LAA_Models"))
length_models = length_models[str_detect(length_models, ".RData")]

# Get size transition matrix
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

  # naming and mungingz
  al_df = reshape2::melt(al_array)
  names(al_df) = c("Age", "Length", "Year", "Value")
  al_df$sex = str_split(length_models[i], "_")[[1]][2]
  al_mat_all_df = rbind(al_mat_all_df, al_df)
  
} # end t loop

# Get status-quo values
# females
al_sq_f_mat <- sq_df$sizeage.f.block2
colnames(al_sq_f_mat) <- seq(41,99,2)
rownames(al_sq_f_mat) <- 2:31
al_sq_f <- reshape2::melt(al_sq_f_mat)
colnames(al_sq_f) = c("Age", "Length", "Ref")
al_sq_f$sex = "female" 

# males
al_sq_m_mat <- sq_df$sizeage.m.block2
colnames(al_sq_m_mat) <- seq(41,99,2)
rownames(al_sq_m_mat) <- 2:31
al_sq_m <- reshape2::melt(al_sq_m_mat)
colnames(al_sq_m) = c("Age", "Length", "Ref")
al_sq_m$sex = "Male" 

# combine
al_sq <- rbind(al_sq_f, al_sq_m)

# compute relative differences
al_rel_diff <- al_mat_all_df %>% 
  left_join(al_sq, by = c("sex", "Length", "Age")) %>% 
  mutate(Anom = (Value - Ref))
val_range <- range(al_rel_diff$Anom)


### Plot AL Matrix Comparison -----------------------------------------------

pdf(here("figs", "ms_figs", "AL_comp_plot.pdf"), width = 14)
print(
  ggplot(al_rel_diff %>% filter(Age %in% c(2:6)) %>% 
           mutate(sex = ifelse(sex == "female", "Female", "Male"))) +
    geom_line(aes(x = Length,  y = Value, color = Year, group = Year), size = 1, alpha = 0.5) +
    geom_line(aes(x = Length,  y = Ref), lty = 2, size = 1.3) +
    scale_color_viridis_c() +
    facet_grid(sex~Age) +
    labs(x = "Length (cm)", y = "Probability", color = "Year") +
    theme_reg() 
)
dev.off()

