# Purpose: To Compare outputs from WAA only runs
# Creator: Matthew LH. Cheng
# Date 12/1/23

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

ages = 2:31
lengths = seq(41,99,2)
iss = 20

# Read in report files
waa_par <- R2admb::read_pars(fn = here("Prop_Within_WAA_Vary","tem"))
sq_par <- R2admb::read_pars(fn = here("Prop_Within_WAA_SQ","tem"))
waa_vary_sab_curr <- dget(here("Prop_Within_WAA_Vary", "tem.rdat")) 
waa_sq_sab_curr <- dget(here("Prop_Within_WAA_SQ", "tem.rdat")) 

# SSB Stuff ---------------------------------------------------------------

# Get par estiamtes for terminal year ssb
waa_ssb <- par_coeff(waa_par, "ssbsd", type = "WAA Vary") %>% mutate(Year = 1960:2023)
sq_ssb <- par_coeff(sq_par, "ssbsd", type = "WAA SQ") %>% mutate(Year = 1960:2023)

# Get par estiamtes for projections too
waa_ssb_proj <- par_coeff(waa_par, "spawn_biom_proj", type = "WAA Vary") %>% mutate(Year = 2024:2038)
sq_ssb_proj <- par_coeff(sq_par, "spawn_biom_proj", type = "WAA SQ") %>% mutate(Year = 2024:2038)

# Get b40 as well
waa_b40 <- par_coeff(waa_par, "B40", type = "WAA Vary") 
sq_b40 <- par_coeff(sq_par, "B40", type = "WAA SQ") 
b40_df = rbind(waa_b40, sq_b40)

# residual munging
ssb_df = rbind(waa_ssb, sq_ssb, waa_ssb_proj, sq_ssb_proj) %>% rename(ssb = coeff) # rename
ssb_df = ssb_df %>% left_join(b40_df %>% rename(b40 = coeff) %>% 
                                select(type, b40), by = "type") %>% 
  mutate(Bratio = ssb/b40) 

ggplot(ssb_df %>% filter(Year <= 2023), 
       aes(x = Year, y = ssb, ymin = downr, ymax = upr, fill = type)) +
  geom_line(aes(color = type), size = 1.5) +
  geom_ribbon(alpha = 0.3) +
  geom_hline(aes(yintercept = b40, color = type), lty = 2, size = 1.5) +
  geom_text(label = "B40%", aes(x = 2013, y = 135), size = 13, check_overlap = TRUE) +
  theme_reg() +
  theme(legend.position = "top") +
  labs(x = "Year", y = "SSB (kt)", color = "Model", fill = "Model")

ggplot(ssb_df, 
       aes(x = Year, y = ssb, ymin = downr, ymax = upr, fill = type)) +
  geom_line(aes(color = type), size = 1.5) +
  geom_ribbon(alpha = 0.3) +
  geom_hline(aes(yintercept = b40, color = type), lty = 2, size = 1.5) +
  geom_text(label = "B40%", aes(x = 2013, y = 135), size = 13, check_overlap = TRUE) +
  theme_reg() +
  theme(legend.position = "top") +
  labs(x = "Year", y = "SSB (kt)", color = "Model", fill = "Model")

ggplot(ssb_df %>% filter(Year >= 1996, Year <= 2023), 
       aes(x = Year, y = ssb, ymin = downr, ymax = upr, fill = type)) +
  geom_line(aes(color = type), size = 1.5) +
  geom_ribbon(alpha = 0.3) +
  geom_hline(aes(yintercept = b40, color = type), lty = 2, size = 1.5) +
  geom_text(label = "B40%", aes(x = 2013, y = 135), size = 13, check_overlap = TRUE) +
  theme_reg() +
  theme(legend.position = "top") +
  labs(x = "Year", y = "SSB (kt)", color = "Model", fill = "Model")

ggplot(ssb_df %>% filter(Year >= 1996, Year<= 2023), 
       aes(x = Year, y = Bratio, color = type)) +
  geom_line(size = 1.5, alpha = 0.55) +
  geom_point(size = 4) +
  geom_hline(yintercept = 1, lty = 2, size = 1.3) +
  theme_reg() +
  theme(legend.position = "top") +
  labs(x = "Year", y = "SSB/B40%", color = "Model", fill = "Model")


# Total Biomass -----------------------------------------------------------
# Get par estiamtes for biomass
waa_biom <- par_coeff(waa_par, "tot_biom", type = "WAA Vary") %>% mutate(Year = 1960:2023)
sq_biom <- par_coeff(sq_par, "tot_biom", type = "WAA SQ") %>% mutate(Year = 1960:2023)
biom_df = rbind(waa_biom, sq_biom) %>% rename(biomass = coeff) # rename

ggplot(biom_df, 
       aes(x = Year, y = biomass, ymin = downr, ymax = upr, fill = type)) +
  geom_line(aes(color = type), size = 1.5) +
  geom_ribbon(alpha = 0.3) +
  theme_reg() +
  theme(legend.position = "top") +
  labs(x = "Year", y = "Total Biomass (kt)", color = "Model", fill = "Model")

# Recruitment Stuff -------------------------------------------------------
waa_rec <- par_coeff(waa_par, "pred_rec", type = "WAA Vary") %>% mutate(Year = 1960:2023)
sq_rec <- par_coeff(sq_par, "pred_rec", type = "WAA SQ") %>%  mutate(Year = 1960:2023)

# Get mean recruitment as well
waa_meanrec <- par_coeff(waa_par, "log_mean_rec", type = "WAA Vary") 
sq_meanrec <- par_coeff(sq_par, "log_mean_rec", type = "WAA SQ") 

meanrec_df = rbind(waa_meanrec, sq_meanrec)
rec_df = rbind(waa_rec, sq_rec) %>%  mutate(downr = ifelse(downr <=0, 0, downr)) %>% rename(rec = coeff)
rec_df = rec_df %>% left_join(meanrec_df %>% rename(meanrec = coeff) %>% 
                                mutate(meanrec = exp(meanrec)) %>%  select(type, meanrec), by = "type") 

ggplot(rec_df %>% filter(Year != 2023), aes(x = Year - 2, y = rec, ymin = downr, ymax = upr, fill = type)) +
  geom_line(aes(color = type), size = 1.2) +
  geom_ribbon(alpha = 0.4) +
  geom_hline(aes(yintercept = meanrec, color = type), lty = 2, size = 1.3) +
  theme_reg() +
  theme(legend.position = "top") +
  ylim(0, NA) +
  labs(x = "Year", y = "Age-2 Recruitment (millions)", color = "Model", fill = "Model")


# Reference Point Stuff ---------------------------------------------------

waa_abc<- par_coeff(waa_par, "ABC", type = "WAA Vary")
sq_abc <- par_coeff(sq_par, "ABC", type = "WAA SQ") 
abc_df = rbind(waa_abc, sq_abc)
ref_df = rbind(b40_df, abc_df)

ggplot(ref_df, aes(x = type, y = coeff, ymin = downr, ymax = upr, color = type)) +
  geom_pointrange() +
  facet_wrap(~coeff_type, scales = "free") +
  labs(x = "Model", y = "Value") +
  theme_reg() +
  theme(legend.position = "none")


# Index Fits --------------------------------------------------------------
idx_df = data.frame()
idx_names = c("Domestic LL Survey Relative Population Weight", "Japanese LL Survey Relative Population Weight",
              "Domestic LL Survey Relative Population Numbers", "Japanese LL Survey Relative Population Numbers",
              "Domestic Fishery CPUE Index", "Japanese Fishery CPUE Index", "GOA Trawl Survey Biomass (kt)") # index names

# extract out index data
waa_idx_df = data.frame()
waa_idx_list = grwth_vary_sab_curr[str_detect(names(waa_vary_sab_curr), "obssrv")]
sq_idx_df = data.frame()
sq_idx_list = grwth_sq_sab_curr[str_detect(names(waa_sq_sab_curr), "obssrv")]

# Loop through to extract stuff out
for(i in 1:length(waa_idx_list)) {
  # extract out index dataframe (just extracting out components from a list) - growth vary model
  waa_idx_tmp = data.frame(year = as.numeric(rownames(waa_idx_list[[i]])), obs = waa_idx_list[[i]][[1]],
                          lci = waa_idx_list[[i]][[2]], uci = waa_idx_list[[i]][[3]],
                          pred = waa_idx_list[[i]][[4]], type = idx_names[i], Model = "WAA Vary")
  # extract out index dataframe (just extracting out components from a list) - sq model
  sq_idx_tmp = data.frame(year = as.numeric(rownames(sq_idx_list[[i]])), obs = sq_idx_list[[i]][[1]],
                          lci = sq_idx_list[[i]][[2]], uci = sq_idx_list[[i]][[3]],
                          pred = sq_idx_list[[i]][[4]], type = idx_names[i], Model = "WAA SQ")
  # bind together
  waa_idx_df = rbind(waa_idx_df, waa_idx_tmp)
  sq_idx_df = rbind(sq_idx_df, sq_idx_tmp)
} # end i loop

idx_df = rbind(waa_idx_df, sq_idx_df)
# relvel by index name
idx_df$type = factor(idx_df$type, levels = idx_names)

ggplot(idx_df) +
  geom_line(mapping = aes(x = year, y = pred, col = Model), size = 1.5) +
  geom_pointrange(mapping = aes(x = year, y = obs, ymin = lci, ymax = uci), 
                  alpha = 0.35, size = 0.75, col = "black") +
  facet_wrap(~type, scales = "free") +
  scale_y_continuous(limits = c(0, NA)) +  # NA for no upper limit
  labs(x = "Year", y = "Index") +
  theme_reg() +
  theme(legend.position = "top")


# Fishing Stuff -----------------------------------------------------------
waa_f1<- par_coeff(waa_par, "Fmort_fish1", type = "WAA Vary") %>% mutate(Year = 1960:2023)
sq_f1 <- par_coeff(sq_par, "Fmort_fish1", type = "WAA SQ") %>% mutate(Year = 1960:2023)
waa_f3<- par_coeff(waa_par, "Fmort_fish3", type = "WAA Vary") %>%  mutate(Year = 1960:2023)
sq_f3 <- par_coeff(sq_par, "Fmort_fish3", type = "WAA SQ") %>% mutate(Year = 1960:2023)
fmort_all = rbind(waa_f1, sq_f1, sq_f3, waa_f3) %>% 
  mutate(coeff_type = case_when(
    coeff_type == "Fmort_fish1" ~ "Fixed-Gear Fishery",
    coeff_type == "Fmort_fish3" ~ "Trawl Gear Fishery",
  ))

ggplot(fmort_all, aes(x = Year, y = coeff, ymin = downr, ymax = upr, fill = type)) +
  geom_line(aes(color = type), size = 1) +
  geom_ribbon(alpha = 0.3) +
  facet_wrap(~coeff_type) +
  theme_reg() +
  theme(legend.position = "top") +
  labs(x = "Year", y = "Annual Fishing Mortality Multiplier", color = "Model", fill = "Model")
