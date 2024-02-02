# Purpose: To create inputs for ADMB dat file
# Creator: Matthew LH. Cheng
# Date: 11/18/23


# Set up ------------------------------------------------------------------
library(here)
library(tidyverse)
source(here("R", "functions.R"))

dir.figs = here("figs", "Data Exploration")
dir.out = here("output", "ADMB Inputs")
dir.create(dir.out)
dir.create(dir.figs)

# Read in data
laa_est = read.csv(here("output", "LAA_Models", "Growth_estimates.csv"))
waa_est = read.csv(here("output", "WAA_Models", "Growth_estimates.csv")) 

# define years to put int constant
impute_years = 1960:1995 # (fill in with the mean of the time series)

# Create WAA inputs -------------------------------------------------------
### Females -----------------------------------------------------------------
# Get mean waa 
mean_growth_female = waa_est %>% filter(Sex == "female", re_model == "3dar1", peel == 0) %>% 
  filter(Year != 2023) %>% # don't fill in projection year for the mean
  select(Age, Year, Mean) %>% 
  group_by(Age) %>% 
  summarize(Mean = mean(Mean))

# Repeat dataframe 
rep_growth_female = mean_growth_female[rep(seq_len(nrow(mean_growth_female)), 
                                            each = length(impute_years)),]

# create column with imputed years
rep_growth_female$Year <- rep(impute_years, times = nrow(mean_growth_female))

# Rbind this dataframe
growth_female_tv = waa_est %>% filter(Sex == "female", re_model == "3dar1", peel == 0) %>% 
  select(Age, Year, Mean) 

# Output this
growth_female_all = rbind(rep_growth_female, growth_female_tv)
growth_female_all = growth_female_all %>% pivot_wider(names_from = "Year", values_from = "Mean")
write.table(t(as.matrix(growth_female_all)[, -1]), here(dir.out, "female_waa_tv.txt"), sep = " ", row.names = FALSE, col.names = FALSE)

### Males -----------------------------------------------------------------
# Get mean waa 
mean_growth_male = waa_est %>% filter(Sex == "Male", re_model == "3dar1", peel == 0) %>% 
  filter(Year != 2023) %>% # don't fill in projection year for the mean
  select(Age, Year, Mean) %>% 
  group_by(Age) %>% 
  summarize(Mean = mean(Mean))

# Repeat dataframe 
rep_growth_male = mean_growth_male[rep(seq_len(nrow(mean_growth_male)), 
                                           each = length(impute_years)),]

# create column with imputed years
rep_growth_male$Year <- rep(impute_years, times = nrow(mean_growth_male))

# Rbind this dataframe
growth_male_tv = waa_est %>% 
  filter(Sex == "Male", peel == 0, re_model == "3dar1") %>% 
  select(Age, Year, Mean)
  
# Output this
growth_male_all = rbind(rep_growth_male, growth_male_tv)
growth_male_all = growth_male_all %>% pivot_wider(names_from = "Year", values_from = "Mean")
write.table(t(as.matrix(growth_male_all)[, -1]), here(dir.out, "male_waa_tv.txt"), sep = " ", row.names = FALSE, col.names = FALSE)

# Length Processes --------------------------------------------------------
### Females -----------------------------------------------------------------
# Get mean growth 
mean_growth_female = laa_est %>% filter(re_model == "3dar1", Sex == "female", peel == 0) %>% 
  filter(Year != 2023) %>% 
  select(Age, Year, Mean) %>% 
  group_by(Age) %>% 
  summarize(Mean = mean(Mean))

# Repeat dataframe 
rep_growth_female = mean_growth_female[rep(seq_len(nrow(mean_growth_female)), 
                                           each = length(impute_years)),]

# create column with imputed years
rep_growth_female$Year <- rep(impute_years, times = nrow(mean_growth_female))

# Rbind this dataframe
growth_female_tv = laa_est %>% 
  filter(Sex == "female",re_model == "3dar1", peel == 0) %>% 
  select(Age, Year, Mean) 

# Output this
growth_female_all = rbind(rep_growth_female, growth_female_tv)
growth_female_all = growth_female_all %>% pivot_wider(names_from = "Year", values_from = "Mean")
write.table(t(as.matrix(growth_female_all)[, -1]), here(dir.out, "female_laa_tv.txt"), sep = " ", row.names = FALSE, col.names = FALSE)

### Males -----------------------------------------------------------------
# Get mean growth 
mean_growth_male = laa_est %>% filter(re_model == "3dar1", Sex == "Male", peel == 0) %>% 
  filter(Year != 2023) %>% 
  select(Age, Year, Mean) %>% 
  group_by(Age) %>% 
  summarize(Mean = mean(Mean))

# Repeat dataframe 
rep_growth_male = mean_growth_male[rep(seq_len(nrow(mean_growth_male)), 
                                       each = length(impute_years)),]

# create column with imputed years
rep_growth_male$Year <- rep(impute_years, times = nrow(mean_growth_male))

# Rbind this dataframe
growth_male_tv = laa_est %>% 
  filter(re_model == "3dar1", Sex == "Male", peel == 0) %>% 
  select(Age, Year, Mean) 

# Output this
growth_male_all = rbind(rep_growth_male, growth_male_tv)
growth_male_all = growth_male_all %>% pivot_wider(names_from = "Year", values_from = "Mean")
write.table(t(as.matrix(growth_male_all)[, -1]), here(dir.out, "male_laa_tv.txt"), sep = " ", row.names = FALSE, col.names = FALSE)

# Create AL matrices ------------------------------------------------------

### Females -------------------------------------------------------------------
len_bins = seq(41, 99, 2)
age_bins = 2:31
# Extract out quantities
load(here("output", "LAA_Models", paste("length_female_3DAR1_Model.RData", sep = "")))
laa_sigma = cbind(matrix(rowMeans(model$rep$obs_sigma_at[,-ncol(model$rep$obs_sigma_at)]), ncol = 36, nrow = 30), model$rep$obs_sigma_at)
# image(laa_sigma)
laa_mat = t(as.matrix(growth_female_all)[, -1]) # female LAA

# Set up array
female_al_array = array(data = 0, dim = c(nrow(laa_mat), length(len_bins), length(age_bins)),
                        dimnames = list(c(sort(names(growth_female_all[,-1]))), c(len_bins), c(age_bins) ))

for(t in 1:nrow(laa_mat)) { # loop through time to get matrix
  female_al_array[t,,] = t(
    get_al_trans_matrix(age_bins = 2:31, len_bins = seq(41, 99, 2), 
                        mean_length = laa_mat[t,], sd = laa_sigma[,t])
  )
} # end t loop

colSums(female_al_array[1,,], na.rm = T)
plot(female_al_array[60,,2], type = "l")
image(t(female_al_array[64,,]))

write.table(female_al_array, here(dir.out, "female_al_tv.txt"), sep = " ", row.names = FALSE, col.names = FALSE)

### Males -------------------------------------------------------------------
# Extract out quantities
load(here("output", "LAA_Models", paste("length_Male_3DAR1_Model.RData", sep = "")))
laa_sigma = cbind(matrix(rowMeans(model$rep$obs_sigma_at[,-ncol(model$rep$obs_sigma_at)]), ncol = 36, nrow = 30), model$rep$obs_sigma_at)
laa_mat = t(as.matrix(growth_male_all)[, -1]) # female LAA
image(laa_sigma)

# Set up array
male_al_array = array(data = 0, dim = c(nrow(laa_mat), length(len_bins), length(age_bins)),
                        dimnames = list(c(sort(names(growth_male_all[,-1]))), c(len_bins), c(age_bins) ))

for(t in 1:nrow(laa_mat)) { # loop through time to get matrix
  male_al_array[t,,] = t(
    get_al_trans_matrix(age_bins = 2:31, len_bins = seq(41, 99, 2), 
                        mean_length = laa_mat[t,], sd = laa_sigma[,t])
  )
} # end t loop

colSums(male_al_array[1,,])
plot(male_al_array[60,,1])
image(t(male_al_array[61,,]))

write.table(male_al_array, here(dir.out, "male_al_tv.txt"), sep = " ", row.names = FALSE, col.names = FALSE)

# Composition Inputs ----------------------------------------------------

# Load in data
ll_srv_ages = data.table::fread(here("data", "ll_survey_ages.csv"), skip = 5)
ll_srv_lengths = data.table::fread(here("data", "ll_survey_lengths.csv"), skip = 5)
jp_ll_srv_ages = data.table::fread(here('data', "jp_ll_survey_ages.csv"), skip = 5)
jp_ll_srv_lengths = data.table::fread(here('data', "jp_ll_survey_lengths.csv"), skip = 5)
goa_trawl_srv_lengths = data.table::fread(here("data", "goa_trawl_survey_lengths.csv"), skip = 7)
fish_ages = data.table::fread(here("data", "fishery_ages.csv"), skip = 6)
fish_lens = data.table::fread(here("data", "fishery_lengths.csv"), skip = 6)

# Fixed-Gear Fishery Age Compositions ------------------------------------------------

# Munging for fishery ages
fixed_gear_fishages = fish_ages %>% 
  rename_all(tolower) %>%  # Change to lower case
  # Filtering certai things out
  filter(
    `gear description` %in% c("LONGLINER", "POT OR TRAP"), # fixed-gear fleet
    year > 1998,
    !is.na(`fmp subarea`),
    !is.na(age),
    age > 1,
    sex != "U",
    `fmp area`!= 'INSID',
  ) %>% 
  # Create plus group
  mutate(age = ifelse(age > 31, 31, age)) 

# get input sample sizes for fixed-gear fleet
iss_fixedgear_ages = fixed_gear_fishages %>% 
  group_by(year) %>%
  summarize(SquareRoot_Hauls = sqrt(n_distinct(`haul join`))) %>% 
  mutate(Data = "Fixed-Gear Fishery Ages")

# munge to get proportions across sexes
fixed_gear_fishages_prop = fixed_gear_fishages %>% 
  select(age, year, sex) %>% 
  group_by(year) %>% 
  mutate(total = n()) %>% # get total sampled by year
  group_by(year, age, sex, total) %>% 
  summarize(count_ages = n()) %>% # get total ages sampeld for a given age, year, sex
  mutate(prop = count_ages/total) %>% # get proportionsa cross sex
  select(year, age, sex, prop)

# Get proportions across
fem_fixed_gear_fishages_prop = fixed_gear_fishages_prop %>% 
  filter(sex == "F") %>% 
  pivot_wider(names_from = "age", values_from = "prop") %>% 
  dplyr::select(year, as.character(2:31)) # pivot wider for proporitions across

male_fixed_gear_fishages_prop = fixed_gear_fishages_prop %>% 
  filter(sex == "M") %>% 
  pivot_wider(names_from = "age", values_from = "prop") %>% 
  dplyr::select(year, as.character(2:31)) # pivot wider for proporitions across

# Construct this into a matrix (rows = years, columns = ages, females then males)
fixed_gear_fishages_mat = as.matrix(cbind(fem_fixed_gear_fishages_prop[,-c(1:2)], male_fixed_gear_fishages_prop[,-c(1:2)]))
fixed_gear_fishages_mat[is.na(fixed_gear_fishages_mat)] = 0 # replace NAs with 0s

# Now, create our text file
n_years_fixed_gear_ages = nrow(fixed_gear_fishages_mat)
years_fixed_gear_ages = unique(fixed_gear_fishages_prop$year)
iss_fixed_gear_ages = round(iss_fixedgear_ages$SquareRoot_Hauls)

# Write this out into a text file matching the admb format
all_info_mat <- c(
  paste("# Number of years age comps: nyrs_fishage"),
  paste(n_years_fixed_gear_ages, collapse = " "),
  paste("# Fishery age comps years"),
  paste(years_fixed_gear_ages, collapse = " "),
  paste("# Number of samples: nsamples_fish_age(1,nyrs_fish_age)"),
  paste(iss_fixed_gear_ages, collapse = " "),
  paste("# Observed fishery age compositions (proportions at age): oac_fish(1,nyrs_fish_age,1,nages*2)"),
  apply(fixed_gear_fishages_mat, 1, function(row) paste(row, collapse = " ")),
  paste("#: #!")
)

write.table(all_info_mat, here(dir.out, "fixed_gear_fishages.txt"), sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE)
# plot comps
pdf(here(dir.figs, "Fixed_Gear_Comps.pdf"))
ggplot(fixed_gear_fishages_prop, aes(x = age, y = prop, color = sex)) +
  geom_line() +
  facet_wrap(~year)
dev.off()  


### Proportions within ------------------------------------------------------

# munge to get proportions across sexes
fixed_gear_fishages_prop = fixed_gear_fishages %>% 
  select(age, year, sex) %>% 
  group_by(year, sex) %>% 
  mutate(total = n()) %>% # get total sampled by year
  group_by(year, age, sex, total) %>% 
  summarize(count_ages = n()) %>% # get total ages sampeld for a given age, year, sex
  mutate(prop = count_ages/total) %>% # get proportionsa within sex
  select(year, age, sex, prop)

# Get proportions across
fem_fixed_gear_fishages_prop = fixed_gear_fishages_prop %>% 
  filter(sex == "F") %>% 
  pivot_wider(names_from = "age", values_from = "prop") %>% 
  dplyr::select(year, as.character(2:31)) # pivot wider for proporitions across

male_fixed_gear_fishages_prop = fixed_gear_fishages_prop %>% 
  filter(sex == "M") %>% 
  pivot_wider(names_from = "age", values_from = "prop") %>% 
  dplyr::select(year, as.character(2:31)) # pivot wider for proporitions across


# Construct this into a matrix (rows = years, columns = ages, females then males)
f_fixed_gear_fishages_mat = as.matrix(fem_fixed_gear_fishages_prop[,-c(1:2)])
m_fixed_gear_fishages_mat = as.matrix(male_fixed_gear_fishages_prop[,-c(1:2)])
f_fixed_gear_fishages_mat[is.na(f_fixed_gear_fishages_mat)] = 0 # replace NAs with 0s
m_fixed_gear_fishages_mat[is.na(m_fixed_gear_fishages_mat)] = 0 # replace NAs with 0s

# Now, create our text file
n_years_fixed_gear_ages = nrow(fixed_gear_fishages_mat)
years_fixed_gear_ages = unique(fixed_gear_fishages_prop$year)
iss_fixed_gear_ages = round(iss_fixedgear_ages$SquareRoot_Hauls)

# Write this out into a text file matching the admb format
all_info_mat <- c(
  paste("# Number of years age comps: nyrs_fishage"),
  paste(n_years_fixed_gear_ages, collapse = " "),
  paste("# Fishery age comps years"),
  paste(years_fixed_gear_ages, collapse = " "),
  paste("# Number of samples: nsamples_fish_age(1,nyrs_fish_age)"),
  paste(iss_fixed_gear_ages, collapse = " "),
  paste("# Observed male prop:"),
  apply(m_fixed_gear_fishages_mat, 1, function(row) paste(row, collapse = " ")),
  paste("# Observed female prop:"),
  apply(f_fixed_gear_fishages_mat, 1, function(row) paste(row, collapse = " ")),
  paste("#: #!")
)

write.table(all_info_mat, here(dir.out, "propwithin_fixed_gear_fishages.txt"), sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE)

# Fixed-Gear Fishery Length Compositions ---------------------------------------------
# Munging for fishery lengths
fixedgear_lengths = fish_lens %>% 
  rename_all(tolower) %>%  # Change to lower case
  # Filtering certai things out
  filter(
    `gear description` %in% c("LONGLINER", "POT OR TRAP"), # fixed-gear fleet
    !is.na(`fmp subarea`),
    !is.na(`length (cm)`),
    !`length (cm)` < 40,
    sex != "U",
    `fmp area`!= 'INSID',
  ) %>% 
  # Create plus group
  mutate(length_2 = ifelse(`length (cm)` < 41, 41,
                           ifelse(
                             `length (cm)` >= max(len_bins), max(len_bins), `length (cm)`
                           )),
         
         length_bins = cut(length_2, breaks = seq(39.99, 99.99, 2),
                           labels = paste(seq(41, 99, 2))))  

# get input sample sizes for ll survey
iss_fixedgear_lengths = fixedgear_lengths %>% 
  group_by(year) %>%
  summarize(SquareRoot_Hauls = sqrt(n_distinct(`haul join`))) %>% 
  mutate(Data = "Fixed-Gear Lengths")

# munge to get proportions across sexes
fixedgear_lengths_prop = fixedgear_lengths %>% 
  select(length_bins, year, sex) %>% 
  group_by(year) %>% 
  mutate(total = n()) %>% # get total sampled by year
  group_by(year, length_bins, sex, total) %>% 
  summarize(count_lengths = n()) %>% # get total ages sampeld for a given age, year, sex
  mutate(prop = count_lengths/total) %>% # get proportionsa cross sex
  select(year, length_bins, sex, prop)

# Get proportions across
fem_fixedgear_lengths_prop = fixedgear_lengths_prop %>% 
  filter(sex == "F") %>% 
  pivot_wider(names_from = "length_bins", values_from = "prop") %>% 
  dplyr::select(year, as.character(seq(41,99,2))) # pivot wider for proporitions across

male_fixedgear_lengths_prop = fixedgear_lengths_prop %>% 
  filter(sex == "M") %>% 
  pivot_wider(names_from = "length_bins", values_from = "prop") %>% 
  dplyr::select(year, as.character(seq(41,99,2))) # pivot wider for proporitions across

# Construct this into a matrix (rows = years, columns = ages, females then males)
fixedgear_lengths_mat = as.matrix(cbind(fem_fixedgear_lengths_prop[,-c(1:2)], male_fixedgear_lengths_prop[,-c(1:2)]))
fixedgear_lengths_mat[is.na(fixedgear_lengths_mat)] = 0 # replace NAs with 0s

# Now, create our text file
n_years_fixedgear_lengths = nrow(fixedgear_lengths_mat)
years_fixedgear_lengths = unique(fixedgear_lengths_prop$year)
iss_fixedgear_lengths_input = round(iss_fixedgear_lengths$SquareRoot_Hauls)

# Write this out into a text file matching the admb format
all_info_mat <- c(
  paste("# Number of years length comps: nyrs_fish1lens"),
  paste(n_years_fixedgear_lengths, collapse = " "),
  paste("# Survey length comps years"),
  paste(years_fixedgear_lengths, collapse = " "),
  paste("# Number of samples: nsamples_fish1_age(1,nyrs_fish1lens)"),
  paste(iss_fixedgear_lengths_input, collapse = " "),
  paste("# Observed survey length compositions (proportions at age): osc_srv(1,nyrs_fish1lens,1,nlens*2)"),
  apply(fixedgear_lengths_mat, 1, function(row) paste(row, collapse = " ")),
  paste("#: #!")
)

write.table(all_info_mat, here(dir.out, "fixed_gear_fishlens.txt"), sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE)

# plot comps
pdf(here(dir.figs, "Fixed_Gear_LenComps.pdf"))
ggplot(fixedgear_lengths_prop, aes(x = length_bins, y = prop, color = factor(sex), group = sex)) +
  geom_line() +
  facet_wrap(~year)
dev.off()  

# Trawl Gear Fishery Length Compositions ----------------------------------
# Munging for fishery lengths
trawlgear_lengths = fish_lens %>% 
  rename_all(tolower) %>%  # Change to lower case
  # Filtering certai things out
  filter(
    year %in% c(1990:1991, 1999, 2005:2022),
    `gear description` %in% c("PELAGIC", "NON PELAGIC"), # Trawl Fleet
    !is.na(`fmp subarea`),
    !is.na(`length (cm)`),
    !`length (cm)` < 40,
    sex != "U",
    `fmp area`!= 'INSID',
  ) %>% 
  # Create plus group
  mutate(length_2 = ifelse(`length (cm)` < 41, 41,
                           ifelse(
                             `length (cm)` >= max(len_bins), max(len_bins), `length (cm)`
                           )),
         
         length_bins = cut(length_2, breaks = seq(39.99, 99.99, 2),
                           labels = paste(seq(41, 99, 2))))  

# get input sample sizes for ll survey
iss_trawlgear_lengths = trawlgear_lengths %>% 
  group_by(year) %>%
  summarize(SquareRoot_Hauls = sqrt(n_distinct(`haul join`))) %>% 
  mutate(Data = "Trawl Gear Lengths")

# munge to get proportions across sexes
trawlgear_lengths_prop = trawlgear_lengths %>% 
  select(length_bins, year, sex) %>% 
  group_by(year) %>% 
  mutate(total = n()) %>% # get total sampled by year
  group_by(year, length_bins, sex, total) %>% 
  summarize(count_lengths = n()) %>% # get total ages sampeld for a given age, year, sex
  mutate(prop = count_lengths/total) %>% # get proportionsa cross sex
  select(year, length_bins, sex, prop)

# Get proportions across
fem_trawlgear_lengths_prop = trawlgear_lengths_prop %>% 
  filter(sex == "F") %>% 
  pivot_wider(names_from = "length_bins", values_from = "prop") %>% 
  dplyr::select(year, as.character(seq(41,99,2))) # pivot wider for proporitions across

male_trawlgear_lengths_prop = trawlgear_lengths_prop %>% 
  filter(sex == "M") %>% 
  pivot_wider(names_from = "length_bins", values_from = "prop") %>% 
  mutate(`95` = NA, `97` = NA, `99` = NA) %>% # appending lengths without any observaitons
  dplyr::select(year, as.character(seq(41,99,2))) # pivot wider for proporitions across

# Construct this into a matrix (rows = years, columns = ages, females then males)
trawlgear_lengths_mat = as.matrix(cbind(fem_trawlgear_lengths_prop[,-c(1:2)], male_trawlgear_lengths_prop[,-c(1:2)]))
trawlgear_lengths_mat[is.na(trawlgear_lengths_mat)] = 0 # replace NAs with 0s

# Now, create our text file
n_years_trawlgear_lengths = nrow(trawlgear_lengths_mat)
years_trawlgear_lengths = unique(trawlgear_lengths_prop$year)
iss_trawlgear_length_input = round(iss_trawlgear_lengths$SquareRoot_Hauls)

# Write this out into a text file matching the admb format
all_info_mat <- c(
  paste("# Number of years length comps: nyrs_fish3lens"),
  paste(n_years_trawlgear_lengths, collapse = " "),
  paste("# fish length comps years"),
  paste(years_trawlgear_lengths, collapse = " "),
  paste("# Number of samples: nsamples_fish3_age(1,nyrs_fish3lens)"),
  paste(iss_trawlgear_length_input, collapse = " "),
  paste("# Observed fish length compositions (proportions at age): osc_fish(1,nyrs_fish3lens,1,nlens*2)"),
  apply(trawlgear_lengths_mat, 1, function(row) paste(row, collapse = " ")),
  paste("#: #!")
)

write.table(all_info_mat, here(dir.out, "traw_gear_fishlens.txt"), sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE)

# plot comps
pdf(here(dir.figs, "Trawl_Gear_LenComps.pdf"))
ggplot(trawlgear_lengths_prop, aes(x = length_bins, y = prop, color = factor(sex), group = sex)) +
  geom_line() +
  facet_wrap(~year)
dev.off()  

# LL Survey Age Compositions -------------------------------------------------
# Munging for ll survey ages
llsurvey_ages = ll_srv_ages %>% 
  rename_all(tolower) %>%  # Change to lower case
  # Filtering certai things out
  filter(
    year >= 1996, 
    !str_detect(`geographic area name`, "Gully"),
    !is.na(age),
    age > 1,
    `sex description` != "Unknown"
  ) %>% 
  # Create plus group
  mutate(age = ifelse(age > 31, 31, age))

# get input sample sizes for ll survey
iss_llsrv_ages = llsurvey_ages %>% 
  group_by(year) %>%
  summarize(SquareRoot_Hauls = sqrt(n_distinct(`station number`))) %>% 
  mutate(Data = "LL Survey Ages")

# munge to get proportions across sexes
llsurvey_ages_prop = llsurvey_ages %>% 
  select(age, year, `sex description`) %>% 
  group_by(year) %>% 
  mutate(total = n()) %>% # get total sampled by year
  group_by(year, age, `sex description`, total) %>% 
  summarize(count_ages = n()) %>% # get total ages sampeld for a given age, year, sex
  mutate(prop = count_ages/total) %>% # get proportionsa cross sex
  select(year, age, `sex description`, prop)

# Get proportions across
fem_llsurvey_ages_prop = llsurvey_ages_prop %>% 
  filter(`sex description` == "female") %>% 
  pivot_wider(names_from = "age", values_from = "prop") %>% 
  dplyr::select(year, as.character(2:31)) # pivot wider for proporitions across

male_llsurvey_ages_prop = llsurvey_ages_prop %>% 
  filter(`sex description` == "Male") %>% 
  pivot_wider(names_from = "age", values_from = "prop") %>% 
  dplyr::select(year, as.character(2:31)) # pivot wider for proporitions across

# Construct this into a matrix (rows = years, columns = ages, females then males)
llsurvey_ages_mat = as.matrix(cbind(fem_llsurvey_ages_prop[,-c(1:2)], male_llsurvey_ages_prop[,-c(1:2)]))
llsurvey_ages_mat[is.na(llsurvey_ages_mat)] = 0 # replace NAs with 0s

# Now, create our text file
n_years_llsurvey_ages = nrow(llsurvey_ages_mat)
years_llsurvey_ages = unique(llsurvey_ages_prop$year)
iss_llsurvey_ages_input = round(iss_llsrv_ages$SquareRoot_Hauls)

# Write this out into a text file matching the admb format
all_info_mat <- c(
  paste("# Number of years age comps: nyrs_srv1age"),
  paste(n_years_llsurvey_ages, collapse = " "),
  paste("# Fishery age comps years"),
  paste(years_llsurvey_ages, collapse = " "),
  paste("# Number of samples: nsamples_srv1_age(1,nyrs_srv1age)"),
  paste(iss_llsurvey_ages_input, collapse = " "),
  paste("# Observed survey age compositions (proportions at age): oac_srv(1,nyrs_srv1age,1,nages*2)"),
  apply(llsurvey_ages_mat, 1, function(row) paste(row, collapse = " ")),
  paste("#: #!")
)

write.table(all_info_mat, here(dir.out, "llsurvey_ages.txt"), sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE)

pdf(here(dir.figs, "LLSurvey_AgeComps.pdf"))
ggplot(llsurvey_ages_prop, aes(x = age, y = prop, color = `sex description`)) +
  geom_line() +
  facet_wrap(~year)
dev.off()  


# Proportions within ------------------------------------------------------

# munge to get proportions across sexes
llsurvey_ages_prop = llsurvey_ages %>% 
  select(age, year, `sex description`) %>% 
  group_by(year, `sex description`) %>% 
  mutate(total = n()) %>% # get total sampled by year
  group_by(year, age, `sex description`, total) %>% 
  summarize(count_ages = n()) %>% # get total ages sampeld for a given age, year, sex
  mutate(prop = count_ages/total) %>% # get proportionsa cross sex
  select(year, age, `sex description`, prop)

# Get proportions across
fem_llsurvey_ages_prop = llsurvey_ages_prop %>% 
  filter(`sex description` == "female") %>% 
  pivot_wider(names_from = "age", values_from = "prop") %>% 
  dplyr::select(year, as.character(2:31)) # pivot wider for proporitions across

male_llsurvey_ages_prop = llsurvey_ages_prop %>% 
  filter(`sex description` == "Male") %>% 
  pivot_wider(names_from = "age", values_from = "prop") %>% 
  dplyr::select(year, as.character(2:31)) # pivot wider for proporitions across

# Construct this into a matrix (rows = years, columns = ages, females then males)
f_llsurvey_ages_mat = as.matrix(fem_llsurvey_ages_prop[,-c(1:2)])
m_llsurvey_ages_mat = as.matrix(male_llsurvey_ages_prop[,-c(1:2)])
f_llsurvey_ages_mat[is.na(f_llsurvey_ages_mat)] = 0 # replace NAs with 0s
m_llsurvey_ages_mat[is.na(m_llsurvey_ages_mat)] = 0 # replace NAs with 0s

# Now, create our text file
n_years_llsurvey_ages = nrow(llsurvey_ages_mat)
years_llsurvey_ages = unique(llsurvey_ages_prop$year)
iss_llsurvey_ages_input = round(iss_llsrv_ages$SquareRoot_Hauls)

# Write this out into a text file matching the admb format
all_info_mat <- c(
  paste("# Number of years age comps: nyrs_srv1age"),
  paste(n_years_llsurvey_ages, collapse = " "),
  paste("# Fishery age comps years"),
  paste(years_llsurvey_ages, collapse = " "),
  paste("# Number of samples: nsamples_srv1_age(1,nyrs_srv1age)"),
  paste(iss_llsurvey_ages_input, collapse = " "),
  paste("# Male obs prop:"),
  apply(m_llsurvey_ages_mat, 1, function(row) paste(row, collapse = " ")),
  paste("# Female obs prop:"),
  apply(f_llsurvey_ages_mat, 1, function(row) paste(row, collapse = " ")),
  paste("#: #!")
)

write.table(all_info_mat, here(dir.out, "propwithin_llsurvey_ages.txt"), sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE)


# LL Survey Length Compositions -------------------------------------------------
# Munging for ll survey lengths
llsurvey_lengths = ll_srv_lengths %>% 
  rename_all(tolower) %>%  # Change to lower case
  # Filtering certai things out
  filter(
    !str_detect(`geographic area name`, "Gully"),
    !is.na(length),
    !length < 40,
    `sex` != "3",
    `common name` == "Sablefish"
  ) %>% 
  # Create plus group
  mutate(length_2 = ifelse(length < 41, 41,
                           ifelse(
                             length >= max(len_bins), max(len_bins), length
                           )),
         
         length_bins = cut(length_2, breaks = seq(39.99, 99.99, 2),
                           labels = paste(seq(41, 99, 2)))) 

# get input sample sizes for ll survey
iss_llsrv_lengths = llsurvey_lengths %>% 
  group_by(year) %>%
  summarize(SquareRoot_Hauls = sqrt(n_distinct(paste(`station number`)))) %>% 
  mutate(Data = "LL Survey Lengths")

# munge to get proportions across sexes
llsurvey_lengths_prop = llsurvey_lengths %>% 
  select(length_bins, year, sex) %>% 
  group_by(year) %>% 
  mutate(total = n()) %>% # get total sampled by year
  group_by(year, length_bins, sex, total) %>% 
  summarize(count_lengths = n()) %>% # get total ages sampeld for a given age, year, sex
  mutate(prop = count_lengths/total) %>% # get proportionsa cross sex
  select(year, length_bins, sex, prop)

# Get proportions across
fem_llsurvey_lengths_prop = llsurvey_lengths_prop %>% 
  filter(sex == "2") %>% 
  pivot_wider(names_from = "length_bins", values_from = "prop") %>% 
  dplyr::select(year, as.character(seq(41,99,2))) # pivot wider for proporitions across

male_llsurvey_lengths_prop = llsurvey_lengths_prop %>% 
  filter(sex == "1") %>% 
  pivot_wider(names_from = "length_bins", values_from = "prop") %>% 
  dplyr::select(year, as.character(seq(41,99,2))) # pivot wider for proporitions across

# Construct this into a matrix (rows = years, columns = ages, females then males)
llsurvey_lengths_mat = as.matrix(cbind(fem_llsurvey_lengths_prop[,-c(1:2)], male_llsurvey_lengths_prop[,-c(1:2)]))
llsurvey_lengths_mat[is.na(llsurvey_lengths_mat)] = 0 # replace NAs with 0s

# Now, create our text file
n_years_llsurvey_lengths = nrow(llsurvey_lengths_mat)
years_llsurvey_lengths = unique(llsurvey_lengths_prop$year)
iss_llsurvey_lengths = round(iss_llsrv_lengths$SquareRoot_Hauls)

# Write this out into a text file matching the admb format
all_info_mat <- c(
  paste("# Number of years length comps: nyrs_srv1lens"),
  paste(n_years_llsurvey_lengths, collapse = " "),
  paste("# Survey length comps years"),
  paste(years_llsurvey_lengths, collapse = " "),
  paste("# Number of samples: nsamples_srv1_age(1,nyrs_srv1lens)"),
  paste(iss_llsurvey_lengths, collapse = " "),
  paste("# Observed survey length compositions (proportions at age): osc_srv(1,nyrs_srv1lens,1,nlens*2)"),
  apply(llsurvey_lengths_mat, 1, function(row) paste(row, collapse = " ")),
  paste("#: #!")
)

write.table(all_info_mat, here(dir.out, "llsurvey_lengths.txt"), sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE)

pdf(here(dir.figs, "LLSurvey_LengthComps.pdf"))
ggplot(llsurvey_lengths_prop, aes(x = length_bins, y = prop, color = factor(sex), group = factor(sex))) +
  geom_line() +
  facet_wrap(~year)
dev.off()  

# GOA Trawl Survey Length Compositions -------------------------------------------------

# Munging for ll survey lengths
goasurvey_lengths = goa_trawl_srv_lengths %>% 
  mutate(`Length (mm)` = `Length (mm)`/10) %>% # mm to cm
  rename_all(tolower) %>%  # Change to lower case
  # Filtering certai things out
  filter(
    year != 2001,
    !is.na(`length (mm)`),
    !`length (mm)` < 40,
    `sex` != "Unknown",
    `common name` == "sablefish", 
    `stratum max depth (m)` <= 500
  ) %>% 
  # Create plus group
  mutate(length_2 = ifelse(`length (mm)` < 41, 41,
                           ifelse(
                             `length (mm)` >= max(len_bins), max(len_bins), `length (mm)`
                           )),
         length_bins = cut(length_2, breaks = seq(39.99, 99.99, 2),
                           labels = paste(seq(41, 99, 2)))) 

# get input sample sizes for ll survey
iss_goasrv_lengths = goasurvey_lengths %>% 
  group_by(year) %>%
  summarize(SquareRoot_Hauls = sqrt(n_distinct(paste(`haul join id`)))) %>% 
  mutate(Data = "GOA Survey Lengths")

# munge to get proportions across sexes
goasurvey_lengths_prop = goasurvey_lengths %>% 
  select(length_bins, year, sex) %>% 
  group_by(year) %>% 
  mutate(total = n()) %>% # get total sampled by year
  group_by(year, length_bins, sex, total) %>% 
  summarize(count_lengths = n()) %>% # get total ages sampeld for a given age, year, sex
  mutate(prop = count_lengths/total) %>% # get proportionsa cross sex
  select(year, length_bins, sex, prop)

# Get proportions across
fem_goasurvey_lengths_prop = goasurvey_lengths_prop %>% 
  filter(sex == "Female") %>% 
  pivot_wider(names_from = "length_bins", values_from = "prop") %>% 
  dplyr::select(year, as.character(seq(41,99,2))) # pivot wider for proporitions across

male_goasurvey_lengths_prop = goasurvey_lengths_prop %>% 
  filter(sex == "Male") %>% 
  pivot_wider(names_from = "length_bins", values_from = "prop") %>%
  mutate(`85` = NA, `89` = NA, `95` = NA, `97` = NA, `99` = NA) %>% # make columns without data
  dplyr::select(year, as.character(seq(41,99,2))) # pivot wider for proporitions across

# Construct this into a matrix (rows = years, columns = ages, females then males)
goasurvey_lengths_mat = as.matrix(cbind(fem_goasurvey_lengths_prop[,-c(1:2)], male_goasurvey_lengths_prop[,-c(1:2)]))
goasurvey_lengths_mat[is.na(goasurvey_lengths_mat)] = 0 # replace NAs with 0s

# Now, create our text file
n_years_goasurvey_lengths = nrow(goasurvey_lengths_mat)
years_goasurvey_lengths = unique(goasurvey_lengths_prop$year)
iss_goasurvey_lengths = round(iss_goasrv_lengths$SquareRoot_Hauls)

# Write this out into a text file matching the admb format
all_info_mat <- c(
  paste("# Number of years length comps: nyrs_srv1lens"),
  paste(n_years_goasurvey_lengths, collapse = " "),
  paste("# Survey length comps years"),
  paste(years_goasurvey_lengths, collapse = " "),
  paste("# Number of samples: nsamples_srv1_age(1,nyrs_srv1lens)"),
  paste(iss_goasurvey_lengths, collapse = " "),
  paste("# Observed survey length compositions (proportions at age): osc_srv(1,nyrs_srv1lens,1,nlens*2)"),
  apply(goasurvey_lengths_mat, 1, function(row) paste(row, collapse = " ")),
  paste("#: #!")
)

write.table(all_info_mat, here(dir.out, "goasurvey_lengths.txt"), sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE)

pdf(here(dir.figs, "goasurvey_LengthComps.pdf"))
ggplot(goasurvey_lengths_prop, aes(x = length_bins, y = prop, color = factor(sex), group = factor(sex))) +
  geom_line() +
  facet_wrap(~year)
dev.off()  


# Japanese LL Survey Length Compositions ----------------------------------

# Munging for ll survey lengths
jpllsurvey_lengths = jp_ll_srv_lengths %>% 
  rename_all(tolower) %>%  # Change to lower case
  # Filtering certai things out
  filter(
    country == "Japan",
    !str_detect(`geographic area name`, "Gully"),
    !is.na(length),
    !length < 40,
    `sex` != "3",
    `common name` == "Sablefish"
  ) %>% 
  # Create plus group
  mutate(length_2 = ifelse(length < 41, 41,
                           ifelse(
                             length >= max(len_bins), max(len_bins), length
                           )),
         
         length_bins = cut(length_2, breaks = seq(39.99, 99.99, 2),
                           labels = paste(seq(41, 99, 2)))) 

# get input sample sizes for ll survey
iss_jpllsrv_lengths = jpllsurvey_lengths %>% 
  group_by(year) %>%
  summarize(SquareRoot_Hauls = sqrt(n_distinct(`station number`))) %>% 
  mutate(Data = "Japanese LL Survey Lengths")

# munge to get proportions across sexes
jpllsurvey_lengths_prop = jpllsurvey_lengths %>% 
  select(length_bins, year, sex) %>% 
  group_by(year) %>% 
  mutate(total = n()) %>% # get total sampled by year
  group_by(year, length_bins, sex, total) %>% 
  summarize(count_lengths = n()) %>% # get total ages sampeld for a given age, year, sex
  mutate(prop = count_lengths/total) %>% # get proportionsa cross sex
  select(year, length_bins, sex, prop)

# Get proportions across
fem_jpllsurvey_lengths_prop = jpllsurvey_lengths_prop %>% 
  filter(sex == "2") %>% 
  pivot_wider(names_from = "length_bins", values_from = "prop") %>% 
  dplyr::select(year, as.character(seq(41,99,2))) # pivot wider for proporitions across

male_jpllsurvey_lengths_prop = jpllsurvey_lengths_prop %>% 
  filter(sex == "1") %>% 
  pivot_wider(names_from = "length_bins", values_from = "prop") %>% 
  dplyr::select(year, as.character(seq(41,99,2))) # pivot wider for proporitions across

# Construct this into a matrix (rows = years, columns = ages, females then males)
jpllsurvey_lengths_mat = as.matrix(cbind(fem_jpllsurvey_lengths_prop[,-c(1:2)], male_jpllsurvey_lengths_prop[,-c(1:2)]))
jpllsurvey_lengths_mat[is.na(jpllsurvey_lengths_mat)] = 0 # replace NAs with 0s

# Now, create our text file
n_years_jpllsurvey_lengths = nrow(jpllsurvey_lengths_mat)
years_jpllsurvey_lengths = unique(jpllsurvey_lengths_prop$year)
iss_jpllsurvey_lengths = round(iss_jpllsrv_lengths$SquareRoot_Hauls)

# Write this out into a text file matching the admb format
all_info_mat <- c(
  paste("# Number of years length comps: nyrs_srv2lens"),
  paste(n_years_jpllsurvey_lengths, collapse = " "),
  paste("# Survey length comps years"),
  paste(years_jpllsurvey_lengths, collapse = " "),
  paste("# Number of samples: nsamples_srv1_age(1,nyrs_srv2lens)"),
  paste(iss_jpllsurvey_lengths, collapse = " "),
  paste("# Observed survey length compositions (proportions at age): osc_srv(1,nyrs_srv2lens,1,nlens*2)"),
  apply(jpllsurvey_lengths_mat, 1, function(row) paste(row, collapse = " ")),
  paste("#: #!")
)

write.table(all_info_mat, here(dir.out, "jpllsurvey_lengths.txt"), sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE)

pdf(here(dir.figs, "jpllsurvey_LengthComps.pdf"))
ggplot(jpllsurvey_lengths_prop, aes(x = length_bins, y = prop, color = factor(sex), group = factor(sex))) +
  geom_line() +
  facet_wrap(~year)
dev.off()  


# Japanese LL Survey Age Compositions -------------------------------------
# 
# Munging for ll survey ages
jp_llsurvey_ages = jp_ll_srv_ages %>%
  rename_all(tolower) %>%  # Change to lower case
  # Filtering certai things out
  filter(
    !str_detect(`geographic area name`, "Gully"),
    !is.na(age),
    age > 1,
    `sex description` != "Unknown"
  ) %>%
  # Create plus group
  mutate(age = ifelse(age > 31, 31, age))

# get input sample sizes for ll survey
iss_jpllsrv_ages = jp_llsurvey_ages %>%
  group_by(year) %>%
  summarize(SquareRoot_Hauls = sqrt(n_distinct(`station number`))) %>%
  mutate(Data = "Japanese LL Survey Ages")

# munge to get proportions across sexes
jp_llsurvey_ages_prop = jp_llsurvey_ages %>%
  select(age, year, `sex description`) %>%
  group_by(year) %>%
  mutate(total = n()) %>% # get total sampled by year
  group_by(year, age, `sex description`, total) %>%
  summarize(count_ages = n()) %>% # get total ages sampeld for a given age, year, sex
  mutate(prop = count_ages/total) %>% # get proportionsa cross sex
  select(year, age, `sex description`, prop)

# Get proportions across
fem_jp_llsurvey_ages_prop = jp_llsurvey_ages_prop %>%
  filter(`sex description` == "female") %>%
  pivot_wider(names_from = "age", values_from = "prop") %>%
  dplyr::select(year, as.character(2:31)) # pivot wider for proporitions across

male_jp_llsurvey_ages_prop = jp_llsurvey_ages_prop %>%
  filter(`sex description` == "Male") %>%
  pivot_wider(names_from = "age", values_from = "prop") %>%
  dplyr::select(year, as.character(2:31)) # pivot wider for proporitions across

# Construct this into a matrix (rows = years, columns = ages, females then males)
jp_llsurvey_ages_mat = as.matrix(cbind(fem_jp_llsurvey_ages_prop[,-c(1:2)], male_jp_llsurvey_ages_prop[,-c(1:2)]))
jp_llsurvey_ages_mat[is.na(jp_llsurvey_ages_mat)] = 0 # replace NAs with 0s

# Now, create our text file
n_years_jp_llsurvey_ages = nrow(jp_llsurvey_ages_mat)
years_jp_llsurvey_ages = unique(jp_llsurvey_ages_prop$year)
iss_jp_llsurvey_ages_input = round(iss_jpllsrv_ages$SquareRoot_Hauls)

# Write this out into a text file matching the admb format
all_info_mat <- c(
  paste("# Number of years age comps: nyrs_srv1age"),
  paste(n_years_jp_llsurvey_ages, collapse = " "),
  paste("# Fishery age comps years"),
  paste(years_jp_llsurvey_ages, collapse = " "),
  paste("# Number of samples: nsamples_srv1_age(1,nyrs_srv1age)"),
  paste(iss_jp_llsurvey_ages_input, collapse = " "),
  paste("# Observed survey age compositions (proportions at age): oac_srv(1,nyrs_srv1age,1,nages*2)"),
  apply(jp_llsurvey_ages_mat, 1, function(row) paste(row, collapse = " ")),
  paste("#: #!")
)

write.table(all_info_mat, here(dir.out, "jp_llsurvey_ages.txt"), sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE)

pdf(here(dir.figs, "jp_llsurvey_AgeComps.pdf"))
ggplot(jp_llsurvey_ages_prop, aes(x = age, y = prop, color = `sex description`)) +
  geom_line() +
  facet_wrap(~year)
dev.off()


# Proportions within ------------------------------------------------------

# munge to get proportions across sexes
jp_llsurvey_ages_prop = jp_llsurvey_ages %>%
  select(age, year, `sex description`) %>%
  group_by(year, `sex description`) %>%
  mutate(total = n()) %>% # get total sampled by year
  group_by(year, age, `sex description`, total) %>%
  summarize(count_ages = n()) %>% # get total ages sampeld for a given age, year, sex
  mutate(prop = count_ages/total) %>% # get proportionsa cross sex
  select(year, age, `sex description`, prop)

# Get proportions across
fem_jp_llsurvey_ages_prop = jp_llsurvey_ages_prop %>%
  filter(`sex description` == "female") %>%
  pivot_wider(names_from = "age", values_from = "prop") %>%
  dplyr::select(year, as.character(2:31)) # pivot wider for proporitions across

male_jp_llsurvey_ages_prop = jp_llsurvey_ages_prop %>%
  filter(`sex description` == "Male") %>%
  pivot_wider(names_from = "age", values_from = "prop") %>%
  dplyr::select(year, as.character(2:31)) # pivot wider for proporitions across

# Construct this into a matrix (rows = years, columns = ages, females then males)
f_jp_llsurvey_ages_mat = as.matrix(fem_jp_llsurvey_ages_prop[,-c(1:2)])
m_jp_llsurvey_ages_mat = as.matrix(male_jp_llsurvey_ages_prop[,-c(1:2)])
f_jp_llsurvey_ages_mat[is.na(f_jp_llsurvey_ages_mat)] = 0 # replace NAs with 0s
m_jp_llsurvey_ages_mat[is.na(m_jp_llsurvey_ages_mat)] = 0 # replace NAs with 0s

# Now, create our text file
n_years_jp_llsurvey_ages = nrow(jp_llsurvey_ages_mat)
years_jp_llsurvey_ages = unique(jp_llsurvey_ages_prop$year)
iss_jp_llsurvey_ages_input = round(iss_jpllsrv_ages$SquareRoot_Hauls)

# Write this out into a text file matching the admb format
all_info_mat <- c(
  paste("# Number of years age comps: nyrs_srv1age"),
  paste(n_years_jp_llsurvey_ages, collapse = " "),
  paste("# Fishery age comps years"),
  paste(years_jp_llsurvey_ages, collapse = " "),
  paste("# Number of samples: nsamples_srv1_age(1,nyrs_srv1age)"),
  paste(iss_jp_llsurvey_ages_input, collapse = " "),
  paste("# Male obs prop:"),
  apply(m_jp_llsurvey_ages_mat, 1, function(row) paste(row, collapse = " ")),
  paste("# Females obs prop:"),
  apply(f_jp_llsurvey_ages_mat, 1, function(row) paste(row, collapse = " ")),
  paste("#: #!")
)

write.table(all_info_mat, here(dir.out, "propwithin_jpllsurvey_ages.txt"), sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE)


# Plot ISS ----------------------------------------------------------------

iss_all = rbind(iss_fixedgear_ages, iss_fixedgear_lengths, iss_trawlgear_lengths, iss_jpllsrv_ages,
                iss_llsrv_ages, iss_llsrv_lengths,iss_goasrv_lengths, iss_jpllsrv_lengths)

pdf(here(dir.figs, "ISS_Composition.pdf"))
ggplot(iss_all, aes(x = year, y = SquareRoot_Hauls, color = Data)) +
  geom_line(size = 1.3) +
  labs(x = "Year", y = "Square Root Hauls", color = "Gear")
dev.off()


# IFQ Fishery, relative selectivity (length-based) ------------------------
yrs = 1990:1995
all_rel_sel = data.frame()
for(i in 1:length(yrs)) {
  female_rel_sel = data.frame(Sex = "Female", 
                              Rel_Selex = t(female_al_array[1,,]) %*% fixedgear_lengths_mat[i,1:30],
                              Time = yrs[i], age = 2:31) # compute relative selectivity using age length matrix
  male_rel_sel = data.frame(Sex = "Male", 
                              Rel_Selex = t(male_al_array[1,,]) %*% fixedgear_lengths_mat[i,31:60],
                              Time = yrs[i], age = 2:31)
  all_rel_sel= rbind(female_rel_sel, male_rel_sel, all_rel_sel)
} # end i loop

pdf(here(dir.figs, "IFQ_Relative_Sel.pdf"))
ggplot(all_rel_sel, aes(x = age, y = Rel_Selex, color = Sex)) +
  geom_line(size = 1.3) + 
  facet_wrap(~Time)
dev.off()
