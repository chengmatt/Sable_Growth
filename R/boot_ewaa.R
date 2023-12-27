# Purpose: To bootstrap hauls and specimen data to create EWAA matrix
# Creator: Matthew LH. Cheng
# Date: 12/12/ 23


# Set up ------------------------------------------------------------------

library(TMB)
library(TMBhelper)
library(here)
library(tidyverse)
library(truncnorm)
library(doSNOW)
library(parallel)
source(here("R", "functions.R"))
dir.figs = here("figs", "Bootstrap")
dir.create(dir.figs)

# set up cores
ncores <- detectCores() 
# Register cluster here
cl <- makeCluster(ncores - 2)
registerDoSNOW(cl)

# Read in ageing error matrix from ADMB assessment
ageing_error_df = readLines(here("data", "SQ_Nominal.dat")) # read in dat file
ageage_mat = ageing_error_df[[which(str_detect(ageing_error_df, "#Age age transition matrix")) + 1]] # read in age-age error matrix
ageage_mat = matrix(as.numeric(str_split(ageage_mat, pattern = " ")[[1]] ), ncol = 30, nrow = 30, dimnames = list(c(2:31), c(2:31))) # coerce strings into matrix

# Read in Survey Age Comp data
age_dat = read.csv(here("data", "age_view.csv"), check.names = FALSE) %>% 
  filter(!is.na(Age), !is.na(`Length (cm)`), !is.na(`Weight (g)`),
         `Sex Description` != "Unknown", `Error Flag` == 0) %>% 
  rename(Length = `Length (cm)`,
         Sex_name = `Sex Description`,
         Weight = `Weight (g)`,
         Haul = `Station Number`) 

sexes = unique(age_dat$Sex_name) # sex names

# Filter to model data
mod_df = age_dat %>% 
  filter(!Age %in% c(0, 1)) %>% 
  mutate(Age = ifelse(Age >= 31, 31, Age), # collapse ages to + group
         Year = Year - min(Year)) %>% 
  select(Sex_name, Length, Weight, Age, Year, Haul) %>% 
  mutate(Read_Age = Age) %>% 
  rename(Tester_Age = Age) 

# Bootstrap procedure -----------------------------------------------------

iters = 2000
age_bins = 2:31
years = sort(unique(mod_df$Year))
boot_all = data.frame()

# progress bar
pb <- txtProgressBar(max = iters, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

boot_all = foreach(iter = 1:iters, .packages = c("TMB", "here", "tidyverse"),
                   .options.snow = opts) %dopar% {
# for(iter in 1:iters) {
  resamp_all = data.frame()
  for(y in 1:length(years)) { # resample specimens
    yrhaul_df = mod_df %>% filter(Year == years[y]) %>% select(Haul) %>% distinct()
    # boot_hauls = sample(yrhaul_df$Haul, length(yrhaul_df$Haul), replace = TRUE)
    boot_hauls = yrhaul_df$Haul # specify hauls to saple from
    for(b in 1:length(boot_hauls)) { # second stage resampling specimens
      yrspec_haul_df =  mod_df %>% filter(Year == years[y], Haul == boot_hauls[b])
      resamp = sample_n(yrspec_haul_df, nrow(yrspec_haul_df), replace = TRUE)
      for(a in 1:nrow(resamp)) {
        read_age = resamp[a,]$Read_Age # get reader age
        prob_tester_ages = ageage_mat[,read_age-1] # get probabilities from tester ages
        resamp[a,]$Tester_Age = sample(age_bins, 1, prob = prob_tester_ages)
      } # end ageing error loop
      resamp_all = rbind(resamp_all, resamp)
    } # end b
  } # end y
  resamp_all %>% mutate(boot_iter = iter)
} # end foreach loop

boot_df = data.table::rbindlist(boot_all) # rbindlist everything


# Summarize bootstrap -----------------------------------------------------

ewaa = boot_df %>% 
  group_by(Tester_Age, Year, Sex_name, boot_iter) %>% 
  summarize(mean = mean(Weight)) 

ewaa_summary = ewaa %>% 
  group_by(Tester_Age, Year, Sex_name) %>%
  summarize(log_sd = sd(log(mean), na.rm = T),
            sd = sd(mean, na.rm = T),
            mean = mean(mean, na.rm = T))


# Plot Bootstrap ----------------------------------------------------------

pdf(here(dir.figs, "Mean_WAA.pdf"))
ggplot(ewaa_summary, aes(x = Year, y = mean, color = factor(Tester_Age))) +
  geom_line() +
  facet_wrap(~Sex_name, scales = "free") +
  theme_bw() +
  theme(legend.position = "top")+
  labs(x = "Year", y = "Mean WAA", color = "Age")
dev.off()

pdf(here(dir.figs, "logsd_meanWAA.pdf"))
ggplot(ewaa_summary, aes(x = Year, y = log_sd, color = Sex_name)) +
  geom_line() +
  facet_wrap(~Tester_Age) +
  theme_bw() +
  theme(legend.position = "top")+
  labs(x = "Year", y = "logSD on WAA", color = "Sex")
dev.off()

write.csv(ewaa_summary, here("output", "ewaa.csv"))
# write.csv(boot_df, here("output", "ewaa_boot.csv"))
