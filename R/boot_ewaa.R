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
  mutate(Age = ifelse(Age >= 31, 31, Age)) %>% 
  mutate(Read_Age = Age) %>% 
  select(Sex_name, Length, Weight, Age, Year, Haul, `Geographic Area Name`, Read_Age) %>% 
  rename(Tester_Age = Age, Area_Name = `Geographic Area Name`)

# Read in catch data
cat_dat = readxl::read_xlsx(here("data", "Area Stratum RPN.xlsx"), skip = 12)

# Munge to get RPN proportions to weight age-length key
cat_prop = cat_dat %>% 
  select(Year, RPN, `Geographic Area Name`) %>% 
  group_by(Year, `Geographic Area Name`) %>% 
  summarize(RPN = sum(RPN, na.rm = T)) %>% # summing RPNs
  rename(Area_Name = `Geographic Area Name`) %>% 
  mutate(RPN = RPN + 0.01) # adding 1 so some samples have some weights

# Left join weight proportion dataframe
cat_age_df = mod_df %>% 
  left_join(cat_prop, by = c("Year", "Area_Name")) %>% 
  mutate(Year = Year - min(Year))

# Bootstrap ---------------------------------------------------------------

iters = 2e3
age_bins = 2:31
years = sort(unique(cat_age_df$Year))

# progress bar
pb <- txtProgressBar(max = iters, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

boot_all = foreach(iter = 1:iters, .packages = c("TMB", "here", "tidyverse"),
                   .options.snow = opts) %dopar% {
                     
 boot_waa_all = data.frame() # reset after every bootstrap
 
 for(y in 1:length(years)) { # resample specimens
   
   resamp_df = data.frame() # reset resampling in each year
   # Get hauls for a given year
   yrhaul_df = cat_age_df %>% filter(Year == years[y]) %>% select(Haul) %>% distinct()
   boot_hauls = sample(yrhaul_df$Haul, length(yrhaul_df$Haul), replace = TRUE) # resample hauls
   boot_hauls = yrhaul_df$Haul # specify hauls to saple from
   
   # Bootstrap hauls for specimens
   for(b in 1:length(boot_hauls)) { # second stage resampling specimens
     
     yrspec_haul_df =  cat_age_df %>% filter(Year == years[y], Haul == boot_hauls[b]) # get specimens to resample from
     resamp = sample_n(yrspec_haul_df, nrow(yrspec_haul_df), replace = TRUE) # resample specimens
     
     for(a in 1:nrow(resamp)) { # resample using ageing error key
       
       read_age = resamp[a,]$Read_Age # get reader age
       prob_tester_ages = ageage_mat[,read_age-1] # get probabilities from tester ages
       resamp[a,]$Tester_Age = sample(age_bins, 1, prob = prob_tester_ages)
       
     } # end ageing error loop
     
     resamp_df = rbind(resamp_df, resamp)
     
   } # end b
   
   # Do density weighting process here
   waa_df = resamp_df %>% 
     select(Tester_Age, Weight, RPN, Sex_name) %>%
     group_by(Tester_Age, Sex_name) %>%
     mutate(prop = RPN/sum(RPN))  %>% # get weights for each length and age
     group_by(Tester_Age, Sex_name) %>% 
     summarize(mean_Weight = sum(Weight * prop)) %>%  # do weighted average to get the mean
     mutate(Year = years[y])

   boot_waa_all = rbind(boot_waa_all, waa_df)
 } # end y
 
 boot_waa_all %>% mutate(boot_iter = iter) # output this
                     
 } # end foreach loop

boot_df = data.table::rbindlist(boot_all) # rbindlist everything

ewaa_summary = boot_df %>% 
  group_by(Tester_Age, Year, Sex_name) %>%
  summarize(log_sd = sd(log(mean_Weight), na.rm = T),
            sd = sd(mean_Weight, na.rm = T),
            mean = mean(mean_Weight, na.rm = T))

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
