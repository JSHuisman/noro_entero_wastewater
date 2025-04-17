###########################################################
# Fig. 1 Showing raw data for 
# Norovirus (NoV GII), Enterovirus (EV), and SARS-COV-2 (COV)
#
# Author: Jana S. Huisman
###########################################################

library(tidyverse)
library(readxl)
library(lubridate)
library(patchwork)

# plotting defaults
theme_set(theme_minimal() +
            theme(text = element_text(size = 20))
)

###########################################################
## Read in raw wastewater data for NoV GII and EV ####

raw_data <- read_xlsx('../data/wastewater_data.xlsx', na = c('#VALUE!', 'NA'))

# Rename columns and filter out NA load values
core_data <- raw_data %>%
  select(date, flow = `flowrate_(kL/day)`, 
         nov_conc = `NoV GII concentration (gc/mL)`, ev_conc = `EV concentration (gc/mL)`, 
         nov_load = `NoVG II load (gc/day)`, ev_load = `EV load (gc/day)`) %>%
  mutate(date = as_date(date)) %>%
  filter(!is.na(nov_load) & !is.na(ev_load)) %>%
  mutate(flow = as.numeric(flow),
         weekday = weekdays(date)) %>%
  arrange(date)

###########################################################
## Read in EAWAG data for SARS-CoV-2 ####

## Protocol 1
# cov_data_set1 <- read_delim('https://sensors-eawag.ch/sars/__data__/processed_normed_data_zurich_v1.csv', ';') %>% 
#   filter(!is.na(`sars_cov2_rna [gc/(d*100000 capita)]`),
#          ...1 < as_date("2021-11-10"))
## Protocol 2
# cov_data_set2 <- read_delim('https://sensors-eawag.ch/sars/__data__/processed_normed_data_zurich_v2.csv', ';') %>% 
#   filter(!is.na(`sars_cov2_rna [gc/(d*100000 capita)]`))
# 
# # The data from v1-v3 is multiplied by 2.5 to arrive at a similar scale as v4
# # we subset to the same timeperiod as noro- and enterovirus
# # multiply by the pop-size of Zurich (471'000)
# cov_data <- cov_data_set1 %>%
#   mutate(`sars_cov2_rna [gc/(d*100000 capita)]` = 2.5*`sars_cov2_rna [gc/(d*100000 capita)]`) %>%
#   bind_rows(cov_data_set2) %>%
#   select(date='...1', pop_load_cov=`sars_cov2_rna [gc/(d*100000 capita)]`, flow=`flow [m^3/d]`) %>%
#   filter(date >= min(core_data$date),
#          date <= max(core_data$date)) %>%
#   mutate(virus = 'SARS-CoV-2',
#          cov_load = pop_load_cov*4.71) %>%
#   arrange(date)
# 
# write_csv(cov_data, '../data/sarscov2_data.csv')

cov_data <- read_csv('../data/sarscov2_data.csv')

###########################################################
## Add data interpolation and smooth across a [21] day window ####

# SARS-CoV-2 is rarefied to the measurement days for noro and enterovirus (left_join)
# complete all days and linearly interpolate
all_data_interpolated <- core_data %>%
  select(date, nov_load, ev_load) %>%
  left_join(cov_data[,c('date', 'cov_load')], by = 'date') %>%
  arrange(date) %>%
  mutate(orig_data = TRUE) %>%
  complete(date = seq.Date(min(date), max(date), by = 'days')) %>%
  mutate(across(where(is.numeric), ~ zoo::na.approx(.x, na.rm = F) ))

write_csv(all_data_interpolated, '../proc_data/all_load_data.csv')


## Smooth data by calculating the moving averaging across a 21 day window
# this uses the interpolated dataset
smooth_data <- all_data_interpolated %>%
  mutate(across(where(is.numeric), ~ zoo::rollapply(.x, 21, mean, align = "center", na.rm = F, fill = NA) )) 
  #mutate(across(where(is.numeric), ~ zoo::rollapply(.x, 28, mean, align = "center", na.rm = F, fill = NA) )) 

# Save the smoothed data for Re estimation in the Fig2 scripts
write_csv(smooth_data, '../proc_data/smooth_load_data.csv')

###########################################################
## Create long versions of the data for plotting ####

# smoothed, interpolated data
smooth_plot_data <- smooth_data %>%
  pivot_longer(ends_with('load'), values_to = 'load', names_to = 'virus') %>%
  mutate(virus = case_when(
    virus == 'nov_load' ~ "Norovirus",
    virus == 'ev_load' ~ "Enterovirus",
    virus == 'cov_load' ~ "SARS-CoV-2"
  ),
  virus = factor(virus, levels = c('Norovirus', 'Enterovirus', 'SARS-CoV-2'))) 

# raw, un-interpolated data
# SARS-CoV-2 is rarefied to the measurement days for noro and enterovirus (left_join)
plot_data <- core_data %>%
  left_join(cov_data[,c('date', 'cov_load')], by = 'date') %>%
  pivot_longer(ends_with('load'), values_to = 'load', names_to = 'virus') %>%
  mutate(virus = case_when(
    virus == 'nov_load' ~ "Norovirus",
    virus == 'ev_load' ~ "Enterovirus",
    virus == 'cov_load' ~ "SARS-CoV-2"
  ),
  virus = factor(virus, levels = c('Norovirus', 'Enterovirus', 'SARS-CoV-2')) 
  ) 

## Determine load ranges of smoothed and raw data ####
# Determine the min and max load (y) values for each virus
# So we can specify them separately in the final plot (this is a workaround to achieve good alpha-values)
load_range_raw_data <- plot_data %>%
  group_by(virus) %>%
  summarise(min_load = min(load, na.rm = T),
         max_load = max(load, na.rm = T))

write_csv(load_range_raw_data, '../proc_data/load_ranges.csv')

load_range_smooth_data <- smooth_plot_data %>%
  group_by(virus) %>%
  summarise(min_load = min(load, na.rm = T),
            max_load = max(load, na.rm = T))

###########################################################
## Read in  + Plot stringency data ####

all_kof_df <- read_csv('../data/kof_data_export_2023-08-10_18_33_19.csv')

kof_df <- all_kof_df %>%
  select(date, kof = 'ch.kof.stringency.zh.stringency_plus') %>%
  filter(date >= as_date('2020-12-18'),
         date <= as_date('2022-09-26')) %>%
  mutate(lag_kof = kof - lag(kof, 1),
         delta_kof = (lag_kof)/kof,
         week_date = ymd(floor_date(date, "weeks", week_start = 1)) )

# Plot the development of the stringency index over time
kof_afo_time <- ggplot(kof_df) +
  geom_line(aes(x = date, y = kof)) +
  scale_x_date(date_breaks = '3 months', date_labels = '%b %Y', 
               limits = c(date('2020-12-01'), date('2022-11-01'))) +
  labs(x = 'Date', y = 'KOF \nindex') +
  coord_cartesian(ylim = c(0, 80))
kof_afo_time

###########################################################
## Read in more data for plotting ####

# Read in changes is stringency (computed in kof script)
kof_changepoint <- read_csv('../data/kof_changepoints.csv')

HUG_daily_average <- read_csv('../proc_data/HUG_daily_average.csv')
rescaled_HUG_study_average <- HUG_daily_average %>%
  mutate(q25 = q25*5e15,
         q75 = q75*5e15) 

USZ_smooth_daily <- read_csv('../proc_data/USZ_smooth_daily.csv')
rescaled_USZ_smooth_count <- USZ_smooth_daily %>%
  mutate(n_cases = n_cases*5e15,
         n_tests_gastro = n_tests_gastro*5e15)

######### Plot Raw WW and clinical data #######################

ww_clinical_plot <- ggplot() +
  geom_rect(data = load_range_raw_data, aes(xmin = as_date("2021-04-01"), xmax = as_date("2021-10-01"),
                                     ymin = 0, ymax = max_load), alpha = 0.4, fill = "lightgrey") +
  geom_rect(data = load_range_raw_data, aes(xmin = as_date("2022-04-01"), xmax = as_date("2022-10-01"),
                                     ymin = 0, ymax = max_load), alpha = 0.4, fill = "lightgrey") +
  geom_ribbon(data = rescaled_HUG_study_average, 
              aes(x = date, ymin = q25, ymax = q75, fill = 'Clinical data'), alpha = 0.5) + 
  geom_line(data = rescaled_USZ_smooth_count ,
            aes(x = date, y = n_cases), color = 'darkred', linewidth = 1) +
  geom_point(data = plot_data , aes(x = date, y = load, colour = virus), alpha = 0.6) +
  geom_line(data = smooth_plot_data   , aes(x = date, y = load, color = virus), linewidth = 1, show.legend = F) +
  geom_vline(xintercept = kof_changepoint %>% filter(lag_kof > 0) %>% pull(date) , alpha = 0.5, linetype = 'longdash') +
  geom_vline(xintercept = kof_changepoint %>% filter(lag_kof < 0) %>% pull(date) , alpha = 0.5, linetype = 'dotted') +
  facet_wrap(vars(virus), ncol = 1, scale = 'free_y') +
  labs(y = 'Viral load in gc/day', x = 'Date', color = 'Virus', fill = '') +
  scale_x_date(date_breaks = '3 months', date_labels = '%b %Y', 
               limits = c(date('2020-12-01'), date('2022-11-01'))) +
  scale_y_continuous(sec.axis = sec_axis(~./5e15, name="Positive tests/day")) +
  scale_color_viridis_d(option = "mako", end = 0.8) +
  scale_fill_manual(values = 'pink') +
  theme(
    legend.position = 'bottom'
  )

ww_clinical_plot

ggsave('../figures/Fig1_raw_data.pdf', width = 12, height = 8)
ggsave('../figures/Fig1_raw_data.png', width = 12, height = 8)

# Combine this plot with the KOF panel to get Fig1 from the paper
ww_clinical_plot + theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + 
    kof_afo_time +
  plot_layout(ncol = 1, heights = c(4,1), guides = 'collect') +
  plot_annotation(tag_levels = c('A')) &
  theme(legend.position = 'bottom', 
        plot.tag.position  = c(.935, .96))

ggsave('../figures/kof_raw_data.png', height = 10, width = 12)
ggsave('../figures/kof_raw_data.pdf', height = 10, width = 12)


###### Supplementary Figure S4 #########################################
## Plot rescaled raw data - all in one figure pane ####
rescaled_data <- plot_data %>%
  mutate(norm_load = case_when(
    virus == 'Norovirus' ~ load/load_range_raw_data[[1,3]],
    virus == 'Enterovirus' ~ load/load_range_raw_data[[2,3]],
    virus == 'SARS-CoV-2' ~ load/load_range_raw_data[[3,3]]
  ))

rescaled_smooth_data <- smooth_plot_data %>%
  mutate(norm_load = case_when(
    virus == 'Norovirus' ~ load/load_range_smooth_data[[1,3]],
    virus == 'Enterovirus' ~ load/load_range_smooth_data[[2,3]],
    virus == 'SARS-CoV-2' ~ load/load_range_smooth_data[[3,3]]
  ))

rescaled_ww_plot <- ggplot() +
  annotate('rect', xmin = as_date("2021-03-21"), xmax = as_date("2021-09-21"),
                                     ymin = 0, ymax = 1, alpha = 0.6, fill = "lightgrey") +
  annotate('rect', xmin = as_date("2022-03-21"), xmax = as_date("2022-09-21"),
                                     ymin = 0, ymax = 1, alpha = 0.6, fill = "lightgrey") +
  geom_line(data = rescaled_smooth_data, aes(x = date, y = norm_load, colour = virus), show.legend = T) +
  geom_vline(xintercept = kof_changepoint %>% filter(lag_kof > 0) %>% pull(date) , alpha = 0.5, linetype = 'longdash') +
  geom_vline(xintercept = kof_changepoint %>% filter(lag_kof < 0) %>% pull(date) , alpha = 0.5, linetype = 'dotted') +
  labs(y = 'Rescaled viral \nload in gc/day', x = 'Date', color = 'Virus') +
  scale_x_date(date_breaks = '3 months',
               date_labels = '%b %y',
               limits = c(date('2020-12-01'), date('2022-11-01'))) +
  scale_color_viridis_d(option = "mako", end = 0.8) +
  theme(
    legend.position = 'bottom'
  )

rescaled_ww_plot

ggsave('../figures/Fig1_raw_data_rescaled.pdf', height = 5, width = 12)
ggsave('../figures/Fig1_raw_data_rescaled.png', height = 5, width = 12)
