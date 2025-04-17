###########################################################
## Combined clinical data
## Author: J S Huisman

## Note: what is referred to as USZ (University hospital Zurich)
# actually stems from the Institute of Medical Virology of the University of Zurich
# and includes some externally ordered tests
###########################################################

library(tidyverse)
library(readxl)
library(lubridate)

theme_set(theme_minimal() +
            theme(text = element_text(size = 20))
)
######################################################################################################################
############ USZ data ###############################
# this data comes as a line list
# dates range from April 2018 to Oct 2022
# specimen categories include "gastrointestinal"; "respiratory"; "invasive"; 
# "oral"; "eyes"; "unspecified"; "hands/feet"; "skin"; "genital/urinary tract"
USZ_raw_data <- read_xlsx("../data/USZ_data_20240531.xlsx")

##### Option 1: Counts per month #####
USZ_count_data <- USZ_raw_data %>%
  mutate(month_n = month(collection_date),
         year_n = year(collection_date),
         month_year = ym(paste0(year_n, '-', month_n)),
         week_date = ymd(floor_date(collection_date, "weeks", week_start = 1)),
         week_n = week(week_date) ) %>%
  group_by(month_year, month_n, year_n, virus, result, specimen_category) %>% 
  summarize(n_cases = n()) %>%
  mutate(hospital = 'IMV/UZH') %>%
  mutate(virus = factor(virus, levels = c('Norovirus', 'Enterovirus', 'SARS-CoV-2')) )

write_csv(USZ_count_data, '../proc_data/USZ_monthly.csv')

##### Option 2: Rolling count of cases (daily) #####

# Total tests per day
USZ_total_daily <- USZ_raw_data %>%
  group_by(collection_date, virus) %>% 
  summarize(n_tests = n()) %>%
  mutate(date = as.Date(collection_date))

# Tests per day for gastrointestinal specimens
USZ_total_gastro_daily <- USZ_raw_data %>%
  group_by(collection_date, virus, specimen_category) %>% 
  summarize(n_tests_gastro = n()) %>%
  mutate(date = as.Date(collection_date)) %>%
  filter(specimen_category == 'gastrointestinal')

# Positive gastrointestinal tests per day; includes columns for total tests and gastrointestinal tests
USZ_daily_count <- USZ_raw_data %>%
  filter(result == 'positive',
         specimen_category == 'gastrointestinal', 
         year(collection_date) != 2018) %>%
  group_by(collection_date, virus) %>% 
  summarize(n_cases = n()) %>%
  mutate(orig_data = TRUE,
         date = as.Date(collection_date)) %>%
  group_by(virus) %>%
  arrange(date) %>%
  complete(date = seq.Date(min(date), max(date), by = 'days')) %>%
  full_join(USZ_total_daily, by = c('date', 'collection_date', 'virus')) %>%
  full_join(USZ_total_gastro_daily, by = c('date', 'collection_date', 'virus')) %>%
  mutate(virus = factor(virus, levels = c('Norovirus', 'Enterovirus', 'SARS-CoV-2')),
         hospital = 'IMV/UZH') %>%
  select(-collection_date, -specimen_category) 

write_csv(USZ_daily_count, '../proc_data/USZ_daily.csv')

# Smooth the positive tests over a 21 day window
USZ_smooth_count <- USZ_daily_count %>%
  mutate(n_cases = ifelse(is.na(n_cases), 0, n_cases),
         n_tests = ifelse(is.na(n_tests), 0, n_tests),
         n_tests_gastro = ifelse(is.na(n_tests_gastro), 0, n_tests_gastro)) %>%
  group_by(virus) %>%
  arrange(virus, date) %>%
  mutate(across(where(is.numeric), ~ zoo::rollapply(.x, 21, mean, align = "center", na.rm = F, fill = NA) )) 
  #complete(date = seq.Date(min(date), max(date), by = 'days')) %>%
  #mutate(across(where(is.numeric), ~ zoo::na.approx(.x, na.rm = F) )) %>%

write_csv(USZ_smooth_count, '../proc_data/USZ_smooth_daily.csv')

######################################################################################################################
############ HUG data #############################
# Linelist data from 2011-12-30 to 2024-05-31
# Note that in this data, a single tested person will have 3 lines (Noro 1, Noro 2, Entero) associated with it
HUG_raw_data <- read_xlsx("../data/HUG_data_20240611.xlsx")

##### Option 1: Counts per month #####
HUG_count_data <- HUG_raw_data %>%
  rename(test_date = `Prise Charge Date`, analysis_method = `Libellé Analyse`, 
         virus = `Libellé Résultat`, result = 'Valeur' ) %>%
  mutate(test_date = date(test_date),
         month_n = month(test_date),
         year_n = year(test_date),
         month_year = ym(paste0(year_n, '-', month_n)),
         week_date = ymd(floor_date(test_date, "weeks", week_start = 1)),
         week_n = week(week_date) ) %>%
  rowwise() %>%
  mutate(virus = str_split_1(virus, ',')[1],
         virus = ifelse(virus == 'Norovirus 2', 'Norovirus', virus), 
         result = ifelse(result == 'POSITIF', 'positive', result)) %>%
  group_by(month_n, year_n, month_year, virus, result, analysis_method) %>% 
  summarize(n_cases = n()) %>%
  mutate(covid = ifelse(year_n %in% c(2020, 2021, 2022), T, F),
         hospital = 'HUG') %>%
  mutate(virus = factor(virus, levels = c('Norovirus', 'Enterovirus', 'SARS-CoV-2')) )  

# Subset to positive tests during the COVID years 2020-2022
HUG_covid_data <- HUG_count_data %>%
  filter(virus %in% c('Enterovirus', 'Norovirus'),
         result == 'positive',
         year_n %in% c(2020, 2021, 2022)) %>%
  full_join(crossing(month_n = seq(1, 12), year_n = c(2020, 2021, 2022), virus = c('Norovirus', 'Enterovirus')), 
            by = c('year_n', 'month_n', 'virus')) %>%
  arrange(year_n, month_n, virus) %>%
  mutate(month_year = ym(paste0(year_n, '-', month_n)),
         n_cases = ifelse(is.na(n_cases), 0, n_cases),
         hospital = 'HUG') %>%
  mutate(virus = factor(virus, levels = c('Norovirus', 'Enterovirus', 'SARS-CoV-2')) )  

# Subset to positive tests during the non-COVID years 2012-2019, 2023
HUG_clean_count_data <- HUG_count_data %>%
  filter(result == 'positive',
         covid == FALSE) %>%
  full_join(crossing(month_n = seq(1, 12), year_n = c(2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2023), 
                     virus = c('Norovirus', 'Enterovirus')), 
            by = c('year_n', 'month_n', 'virus')) %>%
  mutate(n_cases = ifelse(is.na(n_cases), 0, n_cases)) %>%
  mutate(virus = factor(virus, levels = c('Norovirus', 'Enterovirus', 'SARS-CoV-2')) )

# Calculate mean number of monthly positive tests during the non-COVID years (average 2012-2019, 2023)
HUG_monthly_averages <- HUG_clean_count_data %>%
  group_by(virus, month_n) %>%
  summarize(mean_cases = mean(n_cases),
            q25 = quantile(n_cases, probs = 0.25),
            q75 = quantile(n_cases, probs = 0.75),
            sd_cases = sd(n_cases))

# Copy of the mean number of monthly positive tests during the non-COVID years (average 2012-2019, 2023)
# but with month-year dates for the study period so they can be plotted as pink ribbons in general figure
monthly_averages_study_period <- HUG_monthly_averages %>%
  crossing(year_n = c(2020, 2021, 2022)) %>%
  mutate(month_year = ym(paste0(year_n, '-', month_n)) ) %>% 
  filter(virus != 'Norovirus 1')

##### Option 2: Rolling count of cases (daily) #####
result_mapping = c("NEGATIF" = "negative", "POSITIF" = "positive", 
  "NA" = NA, "NON INTERPRÉT." = NA, "NON DETECTÉ" = "negative", 
  "N-INT." = NA, "SANS RES." = NA, "NON REALISE" = NA)

# Note that in this data, a single tested person will have 3 lines (Noro 1, Noro 2, Entero) associated with it
# here they are each counted towards their respective virus category

# Total number of tests per day
HUG_total_daily <- HUG_raw_data %>%
  rename(collection_date = `Prise Charge Date`, analysis_method = `Libellé Analyse`, 
         virus = `Libellé Résultat`, result = 'Valeur' ) %>%
  rowwise() %>%
  mutate(virus = str_split_1(virus, ',')[1],
         virus = ifelse(virus == 'Norovirus 2', 'Norovirus', virus), 
         result = result_mapping[result],
         date = as.Date(collection_date)) %>%
  filter(virus %in% c('Norovirus', 'Enterovirus')) %>%
  group_by(date, virus) %>% 
  summarize(n_tests = n()) 

# Number of positive tests per day, total number of tests are joined on
HUG_daily_count <- HUG_raw_data %>%
  rename(collection_date = `Prise Charge Date`, analysis_method = `Libellé Analyse`, 
         virus = `Libellé Résultat`, result = 'Valeur' ) %>%
  rowwise() %>%
  mutate(virus = str_split_1(virus, ',')[1],
         virus = ifelse(virus == 'Norovirus 2', 'Norovirus', virus), 
         result = result_mapping[result],
         date = as.Date(collection_date)) %>%
  filter(result == 'positive',
         virus %in% c('Norovirus', 'Enterovirus')) %>%
  group_by(date, virus) %>% 
  summarize(n_cases = n()) %>%
  mutate(orig_data = TRUE) %>%
  group_by(virus) %>%
  arrange(date) %>%
  complete(date = seq.Date(min(date), max(date), by = 'days')) %>%
  full_join(HUG_total_daily, by = c('date', 'virus')) %>%
  mutate(virus = factor(virus, levels = c('Norovirus', 'Enterovirus', 'SARS-CoV-2')),
         hospital = 'HUG') 

# Create a 21 day moving average of the daily data
HUG_smooth_count <- HUG_daily_count %>%
  mutate(n_cases = ifelse(is.na(n_cases), 0, n_cases),
         n_tests = ifelse(is.na(n_tests), 0, n_tests)) %>%
  group_by(virus) %>%
  arrange(virus, date) %>%
  mutate(across(where(is.numeric), ~ zoo::rollapply(.x, 21, mean, align = "center", na.rm = F, fill = NA) )) 
#mutate(across(where(is.numeric), ~ zoo::na.approx(.x, na.rm = F) )) %>%
#los_test <- loess(n_cases ~ as.numeric(date), HUG_smooth_count %>% filter(virus == 'Enterovirus'), span  = 21)

# Average the smoothed daily data over all non-COVID years (2012-2019, 2023)
HUG_smooth_average <- HUG_daily_count %>%
  mutate(n_cases = ifelse(is.na(n_cases), 0, n_cases),
         n_tests = ifelse(is.na(n_tests), 0, n_tests)) %>%
  mutate(across(c('n_cases'), ~ zoo::rollapply(.x, 21, mean, align = "center", na.rm = F, fill = NA) )) %>%
  mutate(day_n = day(date),
         month_n = month(date),
         year_n = year(date)) %>%
  filter(!year_n %in% c(2020, 2021, 2022)) %>%
  group_by(virus, day_n, month_n) %>%
  summarize(mean_cases = mean(n_cases, na.rm = T),
            q25 = quantile(n_cases, probs = 0.25, na.rm = T),
            q75 = quantile(n_cases, probs = 0.75, na.rm = T),
            sd_cases = sd(n_cases, na.rm = T)) %>%
  arrange(virus, month_n, day_n) 

# Concatenate the average smoothed daily data 3 times, add dates to plot for the COVID years (2020-2022)
HUG_study_average <- HUG_smooth_average %>%
  crossing(year_n = c(2020, 2021, 2022)) %>%
  mutate(date = ymd(paste0(year_n, '-', month_n, '-', day_n)) )   %>%
  filter(!is.na(date)) # only Feb 29th

write_csv(HUG_study_average, '../proc_data/HUG_daily_average.csv')

######################################################################################################################
############ ENPEN data #############################
# ENPEN_data = data.frame(month_n = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
#                         year_n = 2021,
#                         virus = 'Enterovirus',
#                         n_cases = c(4000, 3000, 3800, 4500, 5700, 8000, 10500,
#                                     8300, 14000, 19500, 14500, 10000)) %>% # digitized by eye
#             mutate(month_year = ym(paste0(year_n, '-', month_n)))
#   
# ggplot(ENPEN_data) +
#   geom_bar(aes(x = month_n, y = n_cases), stat = 'identity')


######################################################################################################################
##### Supplementary Figure S2 ###############################

ggplot(HUG_clean_count_data %>%
         filter(virus %in% c('Enterovirus', 'Norovirus'),
                result == 'positive') ) +
  geom_boxplot(aes(x = month_n, y = n_cases, group = month_n), outliers = F) +
  geom_point(data = HUG_covid_data,
            aes(x = month_n, y = n_cases, color = factor(year_n), group = year_n)) +
  geom_line(data = HUG_covid_data,
             aes(x = month_n, y = n_cases, color = factor(year_n), group = year_n)) +
  scale_x_continuous(breaks = c(3, 6, 9, 12), labels = c('Mar', 'Jun', 'Sep', 'Dec')) +
  labs(x = 'Month', y = 'Number of cases', color = 'Covid years') +
  facet_wrap(vars(virus), ncol = 1, scales = 'free_y') +
  scale_color_manual(values = viridisLite::mako(6, end = 0.8)[c(3, 5, 6)]) +
  theme(legend.position = 'bottom')

ggsave('../figures/HUG_cases_yearly.png', width = 8, height = 6)


###### Combined clinical data for plotting ################

clinical_data <- bind_rows(USZ_count_data, HUG_count_data) %>%
  mutate(hospital = ifelse(hospital == 'USZ', 'IMV/UZH', 'HUG')) %>%
  filter(virus %in% c('Enterovirus', 'Norovirus'),
         result == 'positive', 
         specimen_category %in% c('gastrointestinal', NA),
         year_n %in% c(2020, 2021, 2022)
  ) %>%
  ungroup() %>%
  select(-analysis_method, -specimen_category, - covid, -result) %>%
  full_join(crossing(month_n = seq(1, 12), year_n = c(2020, 2021, 2022), 
                     virus = c('Norovirus', 'Enterovirus'), hospital = c('HUG', 'IMV/UZH')), 
            by = c('year_n', 'month_n', 'virus', 'hospital')) %>%
  mutate(month_year = ym(paste0(year_n, '-', month_n)),
         n_cases = ifelse(is.na(n_cases), 0, n_cases)) %>%
  mutate(virus = factor(virus, levels = c('Norovirus', 'Enterovirus', 'SARS-CoV-2')) ) 

plot_time_df <- clinical_data %>%
  group_by(virus) %>%
  summarise(min_cases = 0,
            max_cases = max(n_cases, na.rm = T))

######## Supplementary Figure S3 #########################
# Combined data

ggplot() +
  geom_rect(data = plot_time_df, aes(xmin = as_date("2020-03-21"), xmax = as_date("2020-09-21"),
                                     ymin = 0, ymax = max_cases), alpha = 0.4, fill = "lightgrey") +
  geom_rect(data = plot_time_df, aes(xmin = as_date("2021-03-21"), xmax = as_date("2021-09-21"),
                                     ymin = 0, ymax = max_cases), alpha = 0.4, fill = "lightgrey") +
  geom_rect(data = plot_time_df, aes(xmin = as_date("2022-03-21"), xmax = as_date("2022-09-21"),
                                     ymin = 0, ymax = max_cases), alpha = 0.4, fill = "lightgrey") +
  geom_ribbon(data = monthly_averages_study_period, 
              aes(x = month_year, ymin = q25, ymax = q75), fill = 'black', alpha = 0.2) + 
  #geom_bar(data = ENPEN_data, aes(x = month_year, y = n_cases/5000), stat = 'identity', position = 'dodge') +
  geom_bar(data = clinical_data,
           aes(x = month_year, y = n_cases, fill = hospital), stat = 'identity', position = 'dodge', alpha = 0.8) + 
  facet_wrap(vars(virus), ncol = 1, scales = "free_y") +
  scale_color_manual(values = viridisLite::mako(6, end = 0.8)[c(3, 5)]) +
  scale_fill_manual(values = viridisLite::mako(6, end = 0.8)[c(3, 5)]) +
  labs(x = 'Month', y = 'Number of cases', color = 'Hospital', fill = 'Hospital') +
  theme(    legend.position = 'bottom')

ggsave('../figures/HUG_USZ_cases.png', width = 12, height = 8)
