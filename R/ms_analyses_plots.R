# Purpose: Plots for manuscript
# Creator: Matthew LH. Cheng (UAF-CFOS)
# 1/30/24

# remove the dang yellow color from colorblind pallete
scale_fill_colorblind7 = function(.ColorList = c(1:4,6:8), ...){
  scale_fill_discrete(..., type = ggthemes::colorblind_pal()(8)[.ColorList])
}
# Color
scale_color_colorblind7 = function(.ColorList = c(1:4,6:8), ...){
  scale_color_discrete(..., type = ggthemes::colorblind_pal()(8)[.ColorList])
}

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

integer_breaks <- function(n = 5, ...) {
  fxn <- function(x) {
    breaks <- floor(pretty(x, n, ...))
    names(breaks) <- attr(breaks, "labels")
    breaks
  }
  return(fxn)
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
waa_dat <- read.csv(here("output", "ewaa.csv"))
laa_dat <- read.csv(here("data", "age_view.csv"), check.names = FALSE) %>% 
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

# Munging WAA data for plotting
waa_dat <- waa_dat %>% 
  rename(Age = Tester_Age, 
         Sex = Sex_name, 
         Mean = mean) %>% 
  mutate(lwr_95 = exp(log(Mean) - 1.96 * log_sd),
         upr_95 = exp(log(Mean) + 1.96 * log_sd),
         Year = Year + 1996,
         broodYear = Year-Age, 
         Sex = ifelse(Sex == "female", "Female", "Male"),
         large_rec_age1 = (2014 - broodYear) + 2,
         large_rec_age2 = (2016 - broodYear) + 2,
         large_rec_age1 = ifelse(large_rec_age1 <= 1, NA, large_rec_age1),
         large_rec_age2 = ifelse(large_rec_age2 <= 1, NA, large_rec_age2))

# Munging LAA data for plotting
laa_dat <- laa_dat %>% 
  rename(Age = Tester_Age, 
         Sex = Sex_name) %>% 
  # group_by(Age, Sex, Year) %>% 
  # summarize(Mean = mean(Length), sd = sd(Length)) %>% 
  mutate(
         # lwr_95 = Mean - 1.96 * sd,
         # upr_95 = Mean + 1.96 * sd,
         Year = Year + 1996,
         broodYear = Year-Age,
         Sex = ifelse(Sex == "female", "Female", "Male"),
         large_rec_age1 = (2014 - broodYear) + 2,
         large_rec_age2 = (2016 - broodYear) + 2,
         large_rec_age1 = ifelse(large_rec_age1 <= 1, NA, large_rec_age1),
         large_rec_age2 = ifelse(large_rec_age2 <= 1, NA, large_rec_age2))


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

rec_plot = ggplot(rec_df, aes(x = Year, y = rec, ymin = downr, ymax = upr, fill = type)) +
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


# Recruitment Comparisons -------------------------------------------------

# Look at differences in ssb and recruitment (probably just distributing
# age classes differently because of the age-length transition)
plot((gv_ssb$coeff - sq_ssb$coeff)/sq_ssb$coeff)
cor(gv_ssb$coeff, sq_ssb$coeff, method = "pearson")

plot((gv_rec$coeff - sq_rec$coeff)/sq_rec$coeff, type = "l")
cor(gv_rec$coeff, sq_rec$coeff, method = "pearson")

gv_abc$coeff
sq_abc$coeff

rec_df %>% 
  group_by(type) %>% 
  filter(Year %in% c(2016:2021)) %>%
  summarize(mean = sum(rec))

# Anomaly Plot ------------------------------------------------------------

# growth estimates here
waa_ts <- read.csv(here("output", "WAA_Models", "Growth_estimates.csv")) %>% 
  mutate(Sex = ifelse(Sex == "female", "Female", "Male"))
laa_ts <- read.csv(here("output", "LAA_Models", "Growth_estimates.csv"))  %>% 
  mutate(Sex = ifelse(Sex == "female", "Female", "Male"))

# Calculate WAA anomaly
const_waa_ts = waa_ts %>% filter(re_model == "constant", Year == 1996, peel == 0) %>% 
  select(Age, Mean, Sex) %>% rename(Ref = Mean)
anom_waa = waa_ts %>% filter(re_model == "3dar1") %>% 
  left_join(const_waa_ts, by = c("Age", "Sex")) %>% 
  mutate(Anom = (Mean-Ref) / Ref,
         Sex = ifelse(Sex == "female", "Female", "Male"))

# Calculate LAA anomaly
const_laa_ts = laa_ts %>% filter(re_model == "constant", Year == 1996, peel == 0) %>% 
  select(Age, Mean, Sex) %>% rename(Ref = Mean)
anom_laa = laa_ts %>% filter(re_model == "3dar1") %>% 
  left_join(const_laa_ts, by = c("Age", "Sex")) %>% 
  mutate(Anom = (Mean-Ref) / Ref,
         Sex = ifelse(Sex == "female", "Female", "Male"))

# anom_waa_plot <- ggplot(anom_waa %>% filter(peel == 0, Year != 2023), 
#                         aes(x = factor(Year), y = factor(Age), fill = Anom)) +
#   geom_tile(alpha = 1) +
#   scale_fill_gradient2(midpoint = mean(anom_waa$Anom), 
#                        high = scales::muted("red"), 
#                        low = scales::muted("blue")) +
#   facet_wrap(~Sex, ncol = 1) +
#   scale_y_discrete(breaks = seq(3, 31, 3)) +
#   scale_x_discrete(breaks = seq(1996, 2023, 5)) +
#   theme_reg() +
#   theme(legend.key.width = unit(1, "cm")) +
#   labs(x = "Year", y = "Age", fill = "Weight-at-age Anomalies") +
#   theme(legend.position = "top")
# 
# anom_laa_plot <- ggplot(anom_laa %>% filter(peel == 0, Year != 2023), 
#                         aes(x = factor(Year), y = factor(Age), fill = Anom)) +
#   geom_tile(alpha = 1) +
#   scale_fill_gradient2(midpoint = mean(anom_laa$Anom), 
#                        high = scales::muted("red"), 
#                        low = scales::muted("blue")) +
#   facet_wrap(~Sex, ncol = 1) +
#   scale_y_discrete(breaks = seq(3, 31, 3)) +
#   scale_x_discrete(breaks = seq(1996, 2023, 5)) +
#   theme_reg() +
#   theme(legend.key.width = unit(1, "cm")) +
#   labs(x = "Year", y = "Age", fill = "Length-at-age Anomalies") +
#   theme(legend.position = "top")
# 
# pdf(here("figs", "ms_figs", "waa_anom.pdf"), width = 12, height = 5)
# ggarrange(anom_waa_plot, anom_laa_plot, labels = c("A", "B"), 
#           font.label = list(size = 25), vjust = 3.5)
# dev.off()


# Compare Residuals -------------------------------------------------------
### Munging residual dataframe -----------------------------------------------------

# Weight at age residuals
waa_resids <- waa_ts %>% 
  filter(peel == 0) %>% 
  rename(Pred = Mean) %>% 
  left_join(waa_dat %>% select(Year, Age, Sex, Mean) %>% rename(Obs = Mean), by = c("Sex", "Age", "Year")) %>% 
  mutate(re_model = ifelse(re_model == "3dar1", "Growth Varying", "Constant Growth"))

# Length-at-age residuals
laa_dat_sub <- laa_dat %>% 
  select(Year, Age, Sex, Length) %>% 
  rename(Obs = Length)
  
# Subset to laa constant models
laa_constant_ts_sb <- laa_ts %>% 
  filter(peel == 0, re_model == "constant") %>%
  rename(Pred = Mean) %>% 
  select(Sex, Age, Year, Pred)

# subset to time varying 3dar1 model
laa_3dar1_ts_sb <- laa_ts %>% 
  filter(peel == 0, re_model == "3dar1") %>%
  rename(Pred = Mean) %>% 
  select(Sex, Age, Year, Pred)
  
# left join here constant model
laa_const_resids <- laa_constant_ts_sb %>% 
    left_join(laa_dat_sub, by = c("Sex", "Age", "Year")) %>% 
    mutate(Model = "Constant Growth")

# left join 3dar1 model
laa_3dar1_resids <- laa_3dar1_ts_sb %>% 
  left_join(laa_dat_sub, by = c("Sex", "Age", "Year")) %>% 
  mutate(Model = "Growth Varying")

# now rbind everything
laa_resids <- rbind(laa_const_resids, laa_3dar1_resids)

# Get standardized residuals WAA
waa_resids <- waa_resids %>% 
  mutate(raw = Obs - Pred) %>% 
  group_by(Sex, re_model) %>% 
  mutate(sd = sd(raw, na.rm = TRUE),
         std_resid = raw / sd)

# Get standardized residuals LAA
laa_resids <- laa_resids %>% 
  mutate(raw = Obs - Pred) %>% 
  group_by(Sex, Model) %>% 
  mutate(sd = sd(raw, na.rm = TRUE),
         std_resid = raw / sd)

### Plot all! ---------------------------------------------------------------

# annual residuals
waa_yr <- ggplot(waa_resids, aes(x = Year, y = std_resid)) +
  geom_point(alpha = 0.075) +
  geom_smooth(aes(fill = re_model, color = re_model), size = 2, method = "gam") +
  geom_hline(yintercept = 0, lty = 2, size = 1) +
  facet_grid(re_model~Sex) +
  scale_color_manual(values = c("#0072B2", "#D55E00")) +
  scale_fill_manual(values = c("#0072B2", "#D55E00")) +
  scale_x_continuous(breaks = seq(1995, 2020, 5)) +
  coord_cartesian(ylim = c(-4, 4)) +
  theme_reg() +
  theme(legend.position = "none") +
  labs(x = "Year", y = "WAA Standardized Residuals", fill = "Model")

# cohort residuals
waa_cohort <- ggplot(waa_resids, aes(x = broodYear, y = std_resid)) +
  geom_point(alpha = 0.075) +
  geom_smooth(aes(fill = re_model, color = re_model), size = 2, method = "gam") +
  geom_hline(yintercept = 0, lty = 2, size = 1) +
  facet_grid(re_model~Sex) +
  scale_x_continuous(breaks = seq(1965, 2020, 10)) +
  scale_color_manual(values = c("#0072B2", "#D55E00")) +
  scale_fill_manual(values = c("#0072B2", "#D55E00")) +
  coord_cartesian(ylim = c(-4, 4)) +
  theme_reg() +
  theme(legend.position = "none") +
  labs(x = "Cohort", y = "WAA Standardized Residuals", fill = "Model")

# length year resids
laa_yr <- ggplot(laa_resids, aes(x = Year, y = std_resid)) +
  geom_point(alpha = 0.008) +
  geom_smooth(aes(fill = Model, color = Model), size = 2, method = "gam") +
  geom_hline(yintercept = 0, lty = 2, size = 1) +
  facet_grid(Model~Sex) +
  coord_cartesian(ylim = c(-4, 4)) +
  scale_x_continuous(breaks = seq(1995, 2020, 5)) +
  scale_color_manual(values = c("#0072B2", "#D55E00")) +
  scale_fill_manual(values = c("#0072B2", "#D55E00")) +
  theme_reg() +
  theme(legend.position = "none") +
  labs(x = "Year", y = "LAA Standardized Residuals", fill = "Model")

# length cohort resids
laa_cohort <- ggplot(laa_resids %>% mutate(broodYear = Year - Age), 
       aes(x = broodYear, y = std_resid)) +
  geom_point(alpha = 0.008) +
  geom_smooth(aes(fill = Model, color = Model), size = 2, method = "gam") +
  geom_hline(yintercept = 0, lty = 2, size = 1) +
  facet_grid(Model~Sex) +
  coord_cartesian(ylim = c(-4, 4)) +
  scale_x_continuous(breaks = seq(1965, 2020, 10)) +
  scale_color_manual(values = c("#0072B2", "#D55E00")) +
  scale_fill_manual(values = c("#0072B2", "#D55E00")) +
  theme_reg() +
  theme(legend.position = "none") +
  labs(x = "Cohort", y = "LAA Standardized Residuals", fill = "Model")

pdf(here("figs", "ms_figs", "growth_residuals.pdf"), width = 13, height = 10)
ggarrange(waa_yr, laa_yr, 
          waa_cohort, laa_cohort, labels = c("A", "B", "C", "D"), font.label = list(size = 20))
dev.off()

# Cohort Plot -------------------------------------------------------------

# create ages at which large recruitment events were observed
waa_ts <- waa_ts %>% mutate(large_rec_age1 = (2014 - broodYear) + 2,
                           large_rec_age2 = (2016 - broodYear) + 2,
                           large_rec_age1 = ifelse(large_rec_age1 <= 1, NA, large_rec_age1),
                           large_rec_age2 = ifelse(large_rec_age2 <= 1, NA, large_rec_age2))

# create ages at which large recruitment events were observed
laa_ts <- laa_ts %>% mutate(large_rec_age1 = (2014 - broodYear) + 2,
                           large_rec_age2 = (2016 - broodYear) + 2,
                           large_rec_age1 = ifelse(large_rec_age1 <= 1, NA, large_rec_age1),
                           large_rec_age2 = ifelse(large_rec_age2 <= 1, NA, large_rec_age2))

waa_plot <- ggplot() +
  geom_line(waa_ts %>% filter(broodYear %in% seq(2010,2016,1), Year != 2023, peel == 0, re_model == "3dar1"),
            mapping = aes(x = Age, y = Mean, color = Sex), size = 1.3) +
  geom_ribbon(waa_ts %>% filter(broodYear %in% seq(2010,2016,1), Year != 2023, peel == 0, re_model == "3dar1"),
              mapping = aes(x = Age, y = Mean, ymin = lwr_95, ymax = upr_95, fill = Sex), alpha = 0.3) +
  geom_line(waa_ts %>% filter(broodYear %in% seq(2010,2016,1), Year != 2023, peel == 0, re_model == "constant"),
            mapping = aes(x = Age, y = Mean), size = 1.3, lty = 3) +
  geom_ribbon(waa_ts %>% filter(broodYear %in% seq(2010,2016,1), Year != 2023, peel == 0, re_model == "constant"),
              mapping = aes(x = Age, y = Mean, ymin = lwr_95, ymax = upr_95), alpha = 0.3) +
  geom_vline(waa_ts %>% filter(broodYear %in% seq(2010,2016,1), Year != 2023, peel == 0, re_model == "constant"), 
             mapping = aes(xintercept = large_rec_age1), size = 0.75) +
  geom_vline(waa_ts %>% filter(broodYear %in% seq(2010,2016,1), Year != 2023, peel == 0, re_model == "constant"), 
             mapping = aes(xintercept = large_rec_age2), lty = 2, size = 0.75) +
  scale_x_continuous(breaks = integer_breaks()) +
  facet_grid(Sex~broodYear, scale = "free") +
  scale_color_manual(values = c("#DC3220", "#005AB5")) +
  scale_fill_manual(values = c("#DC3220", "#005AB5")) +
  labs(x = "Age", y = "Weight (g)", fill = "Sex", color = "Sex") +
  theme_reg() +
  theme(legend.position = "none")

laa_plot <- ggplot() +
  geom_line(laa_ts %>% filter(broodYear %in% seq(2010,2016,1), Year != 2023, peel == 0, re_model == "3dar1"),
            mapping = aes(x = Age, y = Mean, color = Sex), size = 1.3) +
  geom_ribbon(laa_ts %>% filter(broodYear %in% seq(2010,2016,1), Year != 2023, peel == 0, re_model == "3dar1"),
              mapping = aes(x = Age, y = Mean, ymin = lwr_95, ymax = upr_95, fill = Sex), alpha = 0.3) +
  geom_line(laa_ts %>% filter(broodYear %in% seq(2010,2016,1), Year != 2023, peel == 0, re_model == "constant"),
            mapping = aes(x = Age, y = Mean), size = 1.3, lty = 3) +
  geom_vline(laa_ts %>% filter(broodYear %in% seq(2010,2016,1), Year != 2023, peel == 0, re_model == "constant"), 
             mapping = aes(xintercept = large_rec_age1), size = 0.75) +
  geom_ribbon(laa_ts %>% filter(broodYear %in% seq(2010,2016,1), Year != 2023, peel == 0, re_model == "constant"),
              mapping = aes(x = Age, y = Mean, ymin = lwr_95, ymax = upr_95), alpha = 0.3) +
  geom_vline(laa_ts %>% filter(broodYear %in% seq(2010,2016,1), Year != 2023, peel == 0, re_model == "constant"), 
             mapping = aes(xintercept = large_rec_age2), lty = 2, size = 0.75) +
  scale_x_continuous(breaks = integer_breaks()) +
  facet_grid(Sex~broodYear, scale = "free") +
  scale_color_manual(values = c("#DC3220", "#005AB5")) +
  scale_fill_manual(values = c("#DC3220", "#005AB5")) +
  labs(x = "Age", y = "Length (cm)", fill = "Sex", color = "Sex") +
  theme_reg() +
  theme(legend.position = "none")
  
pdf(here("figs", "ms_figs", "cohort_plot.pdf"), width = 13, height = 8)
ggarrange(waa_plot, laa_plot, ncol = 1,
          labels = c('A', "B"),
          font.label = list(size = 25), hjust = -0.6)
dev.off()

# Time Series Plot --------------------------------------------------------

waa_ts_plot <-
ggplot(waa_ts) +
  annotate("rect", xmin = 2014, xmax = Inf, ymin = -Inf, ymax = Inf, color = "grey", fill = "grey", alpha = 0.75) +
  geom_line(waa_ts %>% filter(Age %in% c(3,6,9,12,18,21), re_model == "3dar1", peel == 0, Year != 2023),
            mapping = aes(x = Year, y = Mean, color = factor(Age)), size = 1.3) +
  geom_ribbon(waa_ts %>% filter(Age %in% c(3,6,9,12,18,21), re_model == "3dar1", peel == 0, Year != 2023),
              mapping = aes(x = Year, y = Mean, ymin = lwr_95, ymax = upr_95, fill = factor(Age)), alpha = 0.2) +
  geom_line(waa_ts %>% filter(Age %in% c(3,6,9,12,18,21), re_model == "constant", peel == 0, Year != 2023),
            mapping = aes(x = Year, y = Mean, color = factor(Age)), size = 1.3, lty = 2) +
  geom_ribbon(waa_ts %>% filter(Age %in% c(3,6,9,12,18,21), re_model == "constant", peel == 0, Year != 2023),
            mapping = aes(x = Year, y = Mean, ymin = lwr_95, ymax = upr_95, fill = factor(Age)), alpha = 0.2) +
  facet_wrap(~Sex, ncol = 1, scales = "free") +
  scale_color_colorblind7() +
  scale_fill_colorblind7() +
  theme_reg() +
  labs(y = "Weight (g)", fill = "Ages", color = "Ages")

laa_ts_plot <- 
  ggplot(laa_ts) +
  annotate("rect", xmin = 2014, xmax = Inf, ymin = -Inf, ymax = Inf, color = "grey", fill = "grey", alpha = 0.75) +
  geom_line(laa_ts %>% filter(Age %in% c(3,6,9,12,18,21), re_model == "3dar1", peel == 0, Year != 2023),
            mapping = aes(x = Year, y = Mean, color = factor(Age)), size = 1.3) +
  geom_ribbon(laa_ts %>% filter(Age %in% c(3,6,9,12,18,21), re_model == "3dar1", peel == 0, Year != 2023),
              mapping = aes(x = Year, y = Mean, ymin = lwr_95, ymax = upr_95, fill = factor(Age)), alpha = 0.2) +
  geom_line(laa_ts %>% filter(Age %in% c(3,6,9,12,18,21), re_model == "constant", peel == 0, Year != 2023),
            mapping = aes(x = Year, y = Mean, color = factor(Age)), size = 1.3, lty = 2) +
  geom_ribbon(laa_ts %>% filter(Age %in% c(3,6,9,12,18,21), re_model == "constant", peel == 0, Year != 2023),
              mapping = aes(x = Year, y = Mean, ymin = lwr_95, ymax = upr_95, fill = factor(Age)), alpha = 0.2) +
  scale_color_colorblind7() +
  scale_fill_colorblind7() +
  facet_wrap(~Sex, ncol = 1, scales = "free") +
  theme_reg() +
  labs(y = "Length (cm)", fill = "Ages", color = "Ages")

pdf(here("figs", "ms_figs", "waa_ts.pdf"), width = 10, height = 10)
ggarrange(waa_ts_plot, laa_ts_plot, labels = c("A", "B"),  legend = "right",
          font.label = list(size = 20), common.legend = TRUE)
dev.off()


# Supplementary (All Ages) ------------------------------------------------

waa_ts_plot <-
  ggplot(waa_ts) +
  annotate("rect", xmin = 2014, xmax = Inf, ymin = -Inf, ymax = Inf, color = "grey", fill = "grey", alpha = 0.75) +
  geom_line(waa_ts %>% filter(re_model == "3dar1", peel == 0, Year != 2023),
            mapping = aes(x = Year, y = Mean, color = factor(Sex)), size = 1.3) +
  geom_ribbon(waa_ts %>% filter(re_model == "3dar1", peel == 0, Year != 2023),
              mapping = aes(x = Year, y = Mean, ymin = lwr_95, ymax = upr_95, fill = factor(Sex)), alpha = 0.2) +
  geom_line(waa_ts %>% filter(re_model == "constant", peel == 0, Year != 2023),
            mapping = aes(x = Year, y = Mean, color = factor(Sex)), size = 1.3, lty = 2) +
  geom_ribbon(waa_ts %>% filter(re_model == "constant", peel == 0, Year != 2023),
              mapping = aes(x = Year, y = Mean, ymin = lwr_95, ymax = upr_95, fill = factor(Sex)), alpha = 0.2) +
  facet_wrap(~Age, scales = "free") +
  scale_color_manual(values = c("#DC3220", "#005AB5")) +
  scale_fill_manual(values = c("#DC3220", "#005AB5")) +
  theme_reg() +
  labs(y = "Weight (g)", fill = "Sex", color = "Sex") +
  theme(legend.title = element_text(size = 23),
        legend.text = element_text(size = 20),
        axis.title = element_text(size = 23),
        axis.text = element_text(size = 15),
        strip.text = element_text(size = 23)) 

laa_ts_plot <-
  ggplot(laa_ts) +
  annotate("rect", xmin = 2014, xmax = Inf, ymin = -Inf, ymax = Inf, color = "grey", fill = "grey", alpha = 0.75) +
  geom_line(laa_ts %>% filter(re_model == "3dar1", peel == 0, Year != 2023),
            mapping = aes(x = Year, y = Mean, color = factor(Sex)), size = 1.3) +
  geom_ribbon(laa_ts %>% filter(re_model == "3dar1", peel == 0, Year != 2023),
              mapping = aes(x = Year, y = Mean, ymin = lwr_95, ymax = upr_95, fill = factor(Sex)), alpha = 0.2) +
  geom_line(laa_ts %>% filter(re_model == "constant", peel == 0, Year != 2023),
            mapping = aes(x = Year, y = Mean, color = factor(Sex)), size = 1.3, lty = 2) +
  geom_ribbon(laa_ts %>% filter(re_model == "constant", peel == 0, Year != 2023),
              mapping = aes(x = Year, y = Mean, ymin = lwr_95, ymax = upr_95, fill = factor(Sex)), alpha = 0.2) +
  facet_wrap(~Age, scales = "free") +
  scale_color_manual(values = c("#DC3220", "#005AB5")) +
  scale_fill_manual(values = c("#DC3220", "#005AB5")) +
  theme_reg() +
  labs(y = "Length (cm)", fill = "Sex", color = "Sex") +
    theme(legend.title = element_text(size = 23),
        legend.text = element_text(size = 20),
        axis.title = element_text(size = 23),
        axis.text = element_text(size = 15),
        strip.text = element_text(size = 23)) 

pdf(here("figs", "ms_figs", "supp_waa_ts.pdf"), width = 23, height = 15)
waa_ts_plot
dev.off()

pdf(here("figs", "ms_figs", "supp_laa_ts.pdf"), width = 23, height = 15)
laa_ts_plot
dev.off()

# Compare RPN and Growth --------------------------------------------------

# Weight-at-age
waa_sub <- waa_ts %>% filter(re_model == "3dar1", peel == 0, Year != 2023) # subset waa estimates
rpns <- data.frame(RPN = sq_df$obssrv3$obssrv3, Year = as.numeric(rownames(sq_df$obssrv3))) # extract survey RPNs
rpns_sub <- rpns %>% filter(Year %in% c(1996:2023))
waa_sub <- waa_sub %>% left_join(rpns_sub, by = c("Year")) # leftjoin rpns to WAA data

# Length-at-age
laa_sub <- laa_ts %>% filter(re_model == "3dar1", peel == 0, Year != 2023) # subset laa estimates
laa_sub <- laa_sub %>% left_join(rpns_sub, by = c("Year")) # leftjoin rpns to WAA data

# Plot Weights!
# Females
waa_rpn_plot <-
  ggplot() +
  geom_point(waa_sub %>% filter(Age %in% c(3,6,9,12,18,21)), 
             mapping = aes(x = log(RPN), y = log(Mean)), size = 2, alpha = 0.3) +
  geom_path(waa_sub %>% filter(Year %in% c(2014:2023), Age %in% c(3,6,9,12,18,21)),
            mapping = aes(x = log(RPN), y = log(Mean)), alpha = 0.55, size = 0.85, arrow = arrow()) +
  geom_smooth(waa_sub %>% filter(Age %in% c(3,6,9,12,18,21)), 
                mapping = aes(x = log(RPN), y = log(Mean), fill = Sex, color = Sex), method="lm", se=TRUE, 
                formula=y ~ poly(x, 2, raw=TRUE), alpha = 0.2, lty = 2) + 
  ggpmisc::stat_poly_eq(waa_sub %>% filter( Age %in% c(3,6,9,12,18,21)),  p.digits = 3,
                        mapping = aes(x = log(RPN), y = log(Mean)),
                        label.x = "right",
                        label.y = "top",
                        method="lm", formula=y ~ poly(x, 2, raw=TRUE)) +
  ggpmisc::stat_fit_glance(waa_sub %>% filter(Age %in% c(3,6,9,12,18,21)), method = "lm", 
                    method.args = list(formula = y ~ poly(x, 2, raw=TRUE)),
                    label.x = "right",
                    label.y = 0.8,
                    mapping = aes(x = log(RPN), y = log(Mean), 
                    label = sprintf("italic(p)*\"= \"*%.3g", 
                                    round(after_stat(p.value), 4))), parse = TRUE) +
  facet_grid(Age~Sex, scales = "free") +
  scale_color_manual(values = c("#DC3220", "#005AB5")) +
  scale_fill_manual(values = c("#DC3220", "#005AB5")) +
  theme_reg() +
  labs(x = "log(Survey RPN)", y = "log(Weight (g))")  +
  theme(legend.position = "none")

# Lengths
laa_rpn_plot <- ggplot() +
  geom_point(laa_sub %>% filter(Age %in% c(3,6,9,12,18,21)), 
             mapping = aes(x = log(RPN), y = log(Mean)), size = 2, alpha = 0.3) +
  geom_path(laa_sub %>% filter(Year %in% c(2014:2023), Age %in% c(3,6,9,12,18,21)),
            mapping = aes(x = log(RPN), y = log(Mean)), alpha = 0.55, size = 0.85, arrow = arrow()) +
  geom_smooth(laa_sub %>% filter(Age %in% c(3,6,9,12,18,21)), 
              mapping = aes(x = log(RPN), y = log(Mean), fill = Sex, color = Sex), method="lm", se=TRUE, 
              formula=y ~ poly(x, 2, raw=TRUE), alpha = 0.2, lty = 2) + 
  ggpmisc::stat_poly_eq(laa_sub %>% filter( Age %in% c(3,6,9,12,18,21)),  p.digits = 3,
                        mapping = aes(x = log(RPN), y = log(Mean)),
                        label.x = "right",
                        label.y = "top",
                        method="lm", formula=y ~ poly(x, 2, raw=TRUE)) +
  ggpmisc::stat_fit_glance(laa_sub %>% filter(Age %in% c(3,6,9,12,18,21)), method = "lm", 
                           method.args = list(formula = y ~ poly(x, 2, raw=TRUE)),
                           label.x = "right",
                           label.y = 0.8,
                           mapping = aes(x = log(RPN), y = log(Mean), 
                                         label = sprintf("italic(p)*\"= \"*%.3g", 
                                                         round(after_stat(p.value), 4))), parse = TRUE) +
  facet_grid(Age~Sex, scales = "free") +
  scale_color_manual(values = c("#DC3220", "#005AB5")) +
  scale_fill_manual(values = c("#DC3220", "#005AB5")) +
  theme_reg() +
  labs(x = "log(Survey RPN)", y = "log(Length (cm))") +
  theme(legend.position = "none")

pdf(here("figs", "ms_figs", "RPN_Growth_plot.pdf"), width = 10, height = 10)
ggarrange(waa_rpn_plot, laa_rpn_plot, ncol = 2,
          labels = c('A', "B"),
          font.label = list(size = 25), hjust = -0.6)
dev.off()


# Supplementary (All Ages) ------------------------------------------------

# all ages WAA
f_waa_rpn_plot <-
  ggplot() +
  geom_point(waa_sub %>% filter(Sex == "Female"), 
             mapping = aes(x = log(RPN), y = log(Mean)), size = 2, alpha = 0.3) +
  geom_path(waa_sub %>% filter(Year %in% c(2014:2023), Sex == "Female"),
            mapping = aes(x = log(RPN), y = log(Mean)), alpha = 0.55, size = 0.85, arrow = arrow()) +
  geom_smooth(waa_sub %>% filter(Sex == "Female"), 
              mapping = aes(x = log(RPN), y = log(Mean), fill = Sex, color = Sex), method="lm", se=TRUE, 
              formula=y ~ poly(x, 2, raw=TRUE), alpha = 0.2, lty = 2) + 
  ggpmisc::stat_poly_eq(waa_sub %>% filter(Sex == "Female"),  p.digits = 3,
                        mapping = aes(x = log(RPN), y = log(Mean)),
                        label.x = "right",
                        label.y = "top",
                        method="lm", formula=y ~ poly(x, 2, raw=TRUE)) +
  ggpmisc::stat_fit_glance(waa_sub %>% filter(Sex == "Female"), method = "lm", 
                           method.args = list(formula = y ~ poly(x, 2, raw=TRUE)),
                           label.x = "right",
                           label.y = 0.8,
                           mapping = aes(x = log(RPN), y = log(Mean), 
                                         label = sprintf("italic(p)*\"= \"*%.3g", 
                                                         round(after_stat(p.value), 4))), parse = TRUE) +
  facet_wrap(~Age, scales = "free") +
  scale_color_manual(values = c("#DC3220")) +
  scale_fill_manual(values = c("#DC3220")) +
  theme_reg() +
  labs(x = "log(Survey RPN)", y = "log(Weight (g))")  +
  theme(legend.position = "none")

m_waa_rpn_plot <-
  ggplot() +
  geom_point(waa_sub %>% filter(Sex == "Male"), 
             mapping = aes(x = log(RPN), y = log(Mean)), size = 2, alpha = 0.3) +
  geom_path(waa_sub %>% filter(Year %in% c(2014:2023), Sex == "Male"),
            mapping = aes(x = log(RPN), y = log(Mean)), alpha = 0.55, size = 0.85, arrow = arrow()) +
  geom_smooth(waa_sub %>% filter(Sex == "Male"), 
              mapping = aes(x = log(RPN), y = log(Mean), fill = Sex, color = Sex), method="lm", se=TRUE, 
              formula=y ~ poly(x, 2, raw=TRUE), alpha = 0.2, lty = 2) + 
  ggpmisc::stat_poly_eq(waa_sub %>% filter(Sex == "Male"),  p.digits = 3,
                        mapping = aes(x = log(RPN), y = log(Mean)),
                        label.x = "right",
                        label.y = "top",
                        method="lm", formula=y ~ poly(x, 2, raw=TRUE)) +
  ggpmisc::stat_fit_glance(waa_sub %>% filter(Sex == "Male"), method = "lm", 
                           method.args = list(formula = y ~ poly(x, 2, raw=TRUE)),
                           label.x = "right",
                           label.y = 0.8,
                           mapping = aes(x = log(RPN), y = log(Mean), 
                                         label = sprintf("italic(p)*\"= \"*%.3g", 
                                                         round(after_stat(p.value), 4))), parse = TRUE) +
  facet_wrap(~Age, scales = "free") +
  scale_color_manual(values = c("#005AB5")) +
  scale_fill_manual(values = c("#005AB5")) +
  theme_reg() +
  labs(x = "log(Survey RPN)", y = "log(Weight (g))")  +
  theme(legend.position = "none")

pdf(here("figs", "ms_figs", "supp_rpn_waaF.pdf"), width = 12, height = 10)
f_waa_rpn_plot
dev.off()

pdf(here("figs", "ms_figs", "supp_rpn_waaM.pdf"), width = 12, height = 10)
m_waa_rpn_plot
dev.off()

f_laa_rpn_plot <-
  ggplot() +
  geom_point(laa_sub %>% filter(Sex == "Female"), 
             mapping = aes(x = log(RPN), y = log(Mean)), size = 2, alpha = 0.3) +
  geom_path(laa_sub %>% filter(Year %in% c(2014:2023), Sex == "Female"),
            mapping = aes(x = log(RPN), y = log(Mean)), alpha = 0.55, size = 0.85, arrow = arrow()) +
  geom_smooth(laa_sub %>% filter(Sex == "Female"), 
              mapping = aes(x = log(RPN), y = log(Mean), fill = Sex, color = Sex), method="lm", se=TRUE, 
              formula=y ~ poly(x, 2, raw=TRUE), alpha = 0.2, lty = 2) + 
  ggpmisc::stat_poly_eq(laa_sub %>% filter(Sex == "Female"),  p.digits = 3,
                        mapping = aes(x = log(RPN), y = log(Mean)),
                        label.x = "right",
                        label.y = "top",
                        method="lm", formula=y ~ poly(x, 2, raw=TRUE)) +
  ggpmisc::stat_fit_glance(laa_sub %>% filter(Sex == "Female"), method = "lm", 
                           method.args = list(formula = y ~ poly(x, 2, raw=TRUE)),
                           label.x = "right",
                           label.y = 0.8,
                           mapping = aes(x = log(RPN), y = log(Mean), 
                                         label = sprintf("italic(p)*\"= \"*%.3g", 
                                                         round(after_stat(p.value), 4))), parse = TRUE) +
  facet_wrap(~Age, scales = "free") +
  scale_color_manual(values = c("#DC3220")) +
  scale_fill_manual(values = c("#DC3220")) +
  theme_reg() +
  labs(x = "log(Survey RPN)", y = "log(Length (cm))")  +
  theme(legend.position = "none")

m_laa_rpn_plot <-
  ggplot() +
  geom_point(laa_sub %>% filter(Sex == "Male"), 
             mapping = aes(x = log(RPN), y = log(Mean)), size = 2, alpha = 0.3) +
  geom_path(laa_sub %>% filter(Year %in% c(2014:2023), Sex == "Male"),
            mapping = aes(x = log(RPN), y = log(Mean)), alpha = 0.55, size = 0.85, arrow = arrow()) +
  geom_smooth(laa_sub %>% filter(Sex == "Male"), 
              mapping = aes(x = log(RPN), y = log(Mean), fill = Sex, color = Sex), method="lm", se=TRUE, 
              formula=y ~ poly(x, 2, raw=TRUE), alpha = 0.2, lty = 2) + 
  ggpmisc::stat_poly_eq(laa_sub %>% filter(Sex == "Male"),  p.digits = 3,
                        mapping = aes(x = log(RPN), y = log(Mean)),
                        label.x = "right",
                        label.y = "top",
                        method="lm", formula=y ~ poly(x, 2, raw=TRUE)) +
  ggpmisc::stat_fit_glance(laa_sub %>% filter(Sex == "Male"), method = "lm", 
                           method.args = list(formula = y ~ poly(x, 2, raw=TRUE)),
                           label.x = "right",
                           label.y = 0.8,
                           mapping = aes(x = log(RPN), y = log(Mean), 
                                         label = sprintf("italic(p)*\"= \"*%.3g", 
                                                         round(after_stat(p.value), 4))), parse = TRUE) +
  facet_wrap(~Age, scales = "free") +
  scale_color_manual(values = c("#005AB5")) +
  scale_fill_manual(values = c("#005AB5")) +
  theme_reg() +
  labs(x = "log(Survey RPN)", y = "log(Length (cm))")  +
  theme(legend.position = "none")

pdf(here("figs", "ms_figs", "supp_rpn_laaF.pdf"), width = 12, height = 10)
f_laa_rpn_plot
dev.off()

pdf(here("figs", "ms_figs", "supp_rpn_laaM.pdf"), width = 12, height = 10)
m_laa_rpn_plot
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
length_models = length_models[str_detect(length_models, "3DAR1_Model.RData")]

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

pdf(here("figs", "ms_figs", "AL_comp_plot.pdf"), width = 8, height = 8)
print(
  ggplot(al_rel_diff %>% filter(Age %in% c(3,6,9,12,18,21)) %>% 
           mutate(sex = ifelse(sex == "female", "Female", "Male"))) +
    geom_line(aes(x = Length,  y = Value, color = Year, group = Year), size = 0.85, alpha = 0.7) +
    geom_line(aes(x = Length,  y = Ref), lty = 2, size = 1.3) +
    scale_color_distiller(palette = "YlGnBu", direction = -1) +
    facet_grid(Age~sex) +
    labs(x = "Length (cm)", y = "P(Length|Age)", color = "Year") +
    theme_reg() 
  )
dev.off()

# Compare AIC -------------------------------------------------------------

# Length
# Females
load(here("output", "LAA_Models", "length_female_3DAR1_Model.RData"))
len_f_3dar1 <- model # rename models
load(here("output", "LAA_Models", "length_female_Constant_Model.RData"))
len_f_constant <- model # rename models


# Males
load(here("output", "LAA_Models", "length_Male_3DAR1_Model.RData"))
len_m_3dar1 <- model # rename models
load(here("output", "LAA_Models", "length_Male_Constant_Model.RData"))
len_m_constant <- model # rename models

# Weight
# Females
load(here("output", "WAA_Models", "weight_female_3DAR1_Model.RData"))
wt_f_3dar1 <- model # rename models
load(here("output", "WAA_Models", "weight_female_Constant_Model.RData"))
wt_f_constant <- model # rename models

# Males
load(here("output", "WAA_Models", "weight_Male_3DAR1_Model.RData"))
wt_m_3dar1 <- model # rename models
load(here("output", "WAA_Models", "weight_Male_Constant_Model.RData"))
wt_m_constant <- model # rename models

# put these models into a list
model_list <- list(
  len_f_3dar1 = len_f_3dar1,
  len_f_constant = len_f_constant,
  len_m_3dar1 = len_m_3dar1,
  len_m_constant = len_m_constant,
  wt_f_3dar1 = wt_f_3dar1,
  wt_f_constant = wt_f_constant,
  wt_m_3dar1 = wt_m_3dar1,
  wt_m_constant = wt_m_constant
)

# bind together to get AIC values across models
AIC_df <- data.frame()
for(i in 1:length(model_list)) {
  AIC_value <- TMBhelper::TMBAIC(model_list[[i]]$optim)
  tmp <- data.frame(AIC = AIC_value, Model = names(model_list)[i])
  AIC_df <- rbind(AIC_df, tmp)
} # end i

# Create grouping variable
AIC_df <- AIC_df %>% mutate(Type = ifelse(str_detect(Model, "constant"), "Constant", "3DAR1"),
                            Model = str_remove(Model, "_3dar1"), Model = str_remove(Model, "_constant"))

AIC_df %>% pivot_wider(names_from = "Type", values_from = "AIC") %>% 
  mutate(difference = `3DAR1` - Constant) %>% 
  summarize(mean(difference))
