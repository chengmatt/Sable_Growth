# Purpose: To explore growth data from the LLS for sablefish
# Creator: Matthew LH. Cheng
# Date: 11/3/23


# Set up ------------------------------------------------------------------

library(here)
library(tidyverse)

# read in data
age_dat = read.csv(here("data", "age_view.csv"), check.names = FALSE) %>% 
  filter(!is.na(Age), !is.na(`Length (cm)`), !is.na(`Weight (g)`),
         `Sex Description` != "Unknown", `Error Flag` == 0) %>%  # only keep age samples
  filter(!Age %in% c(0, 1)) %>% 
  mutate(Age = ifelse(Age >= 31, 31, Age))

# Figure out number of age samples available
n_age_samp = age_dat %>% 
  group_by(Year, `Sex Description`) %>% 
  count()

# Plot number of age samples 
ggplot(n_age_samp, aes(x = Year, y = n)) +
  geom_line() +
  facet_wrap(~factor(`Sex Description`))

# Figure out number of age samples available
count_age_samp = age_dat %>% 
  group_by(Year, `Sex Description`, Age) %>% 
  count()

ggplot(count_age_samp, aes(x = Year, y = n, color = `Sex Description`)) +
  geom_line() +
  facet_wrap(~factor(Age), scales = "free")

ggplot(count_age_samp, aes(x = Age, y = n, fill = `Sex Description`)) +
  geom_col(position = "dodge") +
  facet_wrap(~Year)

ggplot(count_age_samp %>% filter(Age %in% c(2:8)), 
       aes(x = Age, y = n, fill = `Sex Description`)) +
  geom_col(position = "dodge") +
  facet_wrap(~Year)

# distribution of ages by year
ggplot(age_dat, aes(x = Age, fill = factor(Year))) +
  geom_density() +
  facet_wrap(~Year) +
  theme(legend.position = "none")

# plot weight-at-age 
ggplot(age_dat %>% filter(Year %in% seq(2014, 2022, 1)),
       aes(x = Age, y = `Weight (g)`, color = factor(Year))) +
  geom_smooth(se = FALSE, method = "gam") +
  facet_grid(~`Sex Description`)

# plot lw relationship
ggplot(age_dat %>% filter(Year %in% seq(2014, 2022, 1)),
       aes(x = `Length (cm)`, y = `Weight (g)`, color = factor(Year))) +
  geom_smooth(se = FALSE, method = "gam") +
  facet_wrap(~`Sex Description`)

# plot laa relationship
ggplot(age_dat %>% filter(Year %in% seq(2014, 2022, 1)),
       aes(x = Age, y = `Length (cm)`, color = factor(Year))) +
  geom_smooth(se = FALSE, method = "gam") +
  facet_wrap(~`Sex Description`)

# seems like most of the variability is in the length-at-age relationship and not the lw relationship
ggplot(age_dat, aes(x = Year, y = `Weight (g)` / `Length (cm)`, group = Year)) +
  geom_boxplot() +
  facet_wrap(~`Sex Description`)

ggplot(age_dat, aes(x = Age, y = `Weight (g)`)) +
  geom_point(alpha = 0.15) +
  geom_smooth() +
  facet_grid(`Sex Description`~`NPFMC Sablefish Mgmt Area`)

ggplot(age_dat, aes(y = `End Latitude (DD)`, x = `End Longitude (DD)`, size = Age)) +
  geom_point(alpha = 0.15)+
  facet_wrap(~Year)
