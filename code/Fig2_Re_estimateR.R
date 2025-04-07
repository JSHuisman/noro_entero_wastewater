###########################################################
# Fig. 2 Showing Re for 
# Norovirus (NoV GII), Enterovirus (EV), and SARS-COV-2 (COV)
#
# Author: Jana S. Huisman
###########################################################

library(tidyverse)
library(readxl)
library(lubridate)
library(estimateR)

# plotting defaults
theme_set(theme_minimal() +
            theme(text = element_text(size = 20))
)

###########################################################
## Load data from Fig1 file + normalise ####
#smooth_data <- read_csv('../proc_data/smooth_load_data.csv')
all_data_interpolated <- read_csv('../proc_data/all_load_data.csv')
load_ranges <- read_csv('../proc_data/load_ranges.csv')

min_loads = setNames(load_ranges %>% pull(min_load), load_ranges %>% pull(virus))

# The estimateR pipeline includes smoothing, which is why we use the unsmoothed, interpolated
# normalized wastewater loads as input
norm_data <- all_data_interpolated %>%
  mutate(norm_nov = nov_load/min_loads[['Norovirus']],
         norm_ev = ev_load/min_loads[['Enterovirus']],
         norm_cov = cov_load/min_loads[['SARS-CoV-2']]) %>%
  filter(!is.na(norm_cov))

###########################################################
## Set all parameters for estimateR inference ####
# Same parameters as in script SuppFig_SLDs.R

# Enterovirus Shedding load distribution - Fan
sld_entero <- list(name = "gamma", 
                 shape = 13.14, 
                 scale = 0.41)

# Norovirus Shedding load distribution - Bernstein
sld_noro <- list(name = "gamma", 
                 shape = 2.72, 
                 scale = 1.49)

# SARS-CoV-2 Shedding load distribution - Benefield
# Ref: Benefield et al., medRxiv, 2020
sld_sars <- list(name = "gamma", 
                 shape = 0.929639, 
                 scale = 7.241397)

## SARS-CoV-2 serial interval (for Re estimation) in days
# Ref: Kremer et al., Emerging Infectious Diseases, 2021 and Backer et al., Eurosurveillance, 2021
mean_si_sars <- 3
std_si_sars <- 2.4

###############################
# General parameters for estimateR
estimation_window = 3 # 3-day sliding window for the Re estimation
minimum_cumul_incidence = 50 # we start estimating Re after at least 100 cases have been recorded
N_bootstrap_replicates = 100 # we take 100 replicates in the bootstrapping procedure

# We specifiy the reference date (first date of data) and the time step of data.
ref_date = min(norm_data$date)
time_step = "day"

###########################################################
## estimate Re with virus-specific SLD ####
# 
# # Currently: smoothing over 28 days
# smooth_param = 28
# 
# entero_estimates <- get_block_bootstrapped_estimate(
#   incidence_data = norm_data$norm_ev,
#   smoothing_method = "LOESS",
#   data_points_incl = smooth_param,
#   N_bootstrap_replicates = N_bootstrap_replicates,
#   delay = list(sld_entero),
#   estimation_window = estimation_window,
#   minimum_cumul_incidence = minimum_cumul_incidence,
#   mean_serial_interval = mean_si_sars,
#   std_serial_interval = std_si_sars,
#   output_Re_only = FALSE,
#   ref_date = ref_date,
#   time_step = time_step,
#   combine_bootstrap_and_estimation_uncertainties = TRUE
# )
# 
# noro_estimates <- get_block_bootstrapped_estimate(
#   incidence_data = norm_data$norm_nov,
#   smoothing_method = "LOESS",
#   data_points_incl = smooth_param,
#   N_bootstrap_replicates = N_bootstrap_replicates,
#   delay = list(sld_noro),
#   estimation_window = estimation_window,
#   minimum_cumul_incidence = minimum_cumul_incidence,
#   mean_serial_interval = mean_si_sars,
#   std_serial_interval = std_si_sars,
#   output_Re_only = FALSE,
#   ref_date = ref_date,
#   time_step = time_step,
#   combine_bootstrap_and_estimation_uncertainties = TRUE
# )
# 
# sars_estimates <- get_block_bootstrapped_estimate(
#   incidence_data = norm_data$norm_cov,
#   smoothing_method = "LOESS",
#   data_points_incl = smooth_param,
#   N_bootstrap_replicates = N_bootstrap_replicates,
#   delay = list(sld_sars),
#   estimation_window = estimation_window,
#   minimum_cumul_incidence = minimum_cumul_incidence,
#   mean_serial_interval = mean_si_sars,
#   std_serial_interval = std_si_sars,
#   output_Re_only = FALSE,
#   ref_date = ref_date,
#   time_step = time_step,
#   combine_bootstrap_and_estimation_uncertainties = TRUE
# )
# 
# all_Re_estimates <- bind_rows(
#   entero_estimates %>% mutate(virus = "Enterovirus"),
#   noro_estimates %>% mutate(virus = "Norovirus"),
#   sars_estimates %>% mutate(virus = "SARS-CoV-2")
# ) %>%
#   mutate(virus = factor(virus, levels = c('Norovirus', 'Enterovirus', 'SARS-CoV-2')))
# 
# 
# write_csv(all_Re_estimates, '../proc_data/Re_ww.csv')

###########################################################
all_Re_estimates <- read_csv('../proc_data/Re_ww.csv') %>%
  mutate(virus = factor(virus, levels = c('Norovirus', 'Enterovirus', 'SARS-CoV-2')))

###########################################################

# Date maximum Re was reached
all_Re_estimates %>% 
  select(date, Re_estimate, CI_down_Re_estimate, CI_up_Re_estimate, virus) %>% 
  group_by(virus) %>%
  summarise(max_Re = max(Re_estimate, na.rm = T),
            CI_down_max_Re = CI_down_Re_estimate[which(Re_estimate == max_Re)],
            CI_up_max_Re = CI_up_Re_estimate[which(Re_estimate == max_Re)],
            max_date = date[which(Re_estimate == max_Re)])


# Average Re over complete time period
all_Re_estimates %>% 
  select(date, Re_estimate, CI_down_Re_estimate, CI_up_Re_estimate, virus) %>% 
  group_by(virus) %>%
  summarise(mean_Re = mean(Re_estimate, na.rm = T),
            std_mean_Re = sd(Re_estimate, na.rm = T) )

###########################################################
## plot Re results
kof_changepoint <- read_csv('../data/kof_changepoints.csv')

plot_deconv_limits <- all_Re_estimates %>%
  group_by(virus) %>%
  summarise(min_deconv = min(CI_down_deconvolved_incidence, na.rm = T),
            max_deconv = max(CI_up_deconvolved_incidence, na.rm = T))

plot_Re_limits <- all_Re_estimates %>%
  group_by(virus) %>%
  summarise(min_Re = min(CI_down_Re_estimate, na.rm = T),
            max_Re = max(CI_up_Re_estimate, na.rm = T))

## Plot Infection Incidence (Deconvolution output) - Supplementary Figure S5 ####
deconv_plot <- ggplot() +
  geom_rect(data = plot_deconv_limits, aes(xmin = as_date("2021-03-21"), xmax = as_date("2021-09-21"),
                                              ymin = 0, ymax = max_deconv), alpha = 0.4, fill = "lightgrey") +
  geom_rect(data = plot_deconv_limits, aes(xmin = as_date("2022-03-21"), xmax = as_date("2022-09-21"),
                                              ymin = 0, ymax = max_deconv), alpha = 0.4, fill = "lightgrey") +
  geom_vline(xintercept = kof_changepoint %>% filter(lag_kof > 0) %>% pull(date) , alpha = 0.5, linetype = 'longdash') +
  geom_vline(xintercept = kof_changepoint %>% filter(lag_kof < 0) %>% pull(date) , alpha = 0.5, linetype = 'dotted') +
  geom_ribbon(data = all_Re_estimates, 
              aes(x=date, ymin = CI_down_deconvolved_incidence,  ymax = CI_up_deconvolved_incidence, fill = virus),
              show.legend = T) +
  labs(x = 'Date' , y='Estimated infection \nincidence per day', fill = 'Virus') +
  scale_x_date(date_breaks = '3 months', date_labels = '%b %Y', 
               limits = c(date('2020-12-01'), date('2022-11-01'))) +
  scale_fill_viridis_d(option = "mako", end = 0.8) +
  facet_wrap(vars(virus), ncol = 1, scale = 'free_y') +
  theme(legend.position = 'bottom')

deconv_plot

ggsave('../figures/deconv_smooth.pdf', width = 12, height = 8)
ggsave('../figures/deconv_smooth.png', width = 12, height = 8)

## Plot Re - Main text figure 2 ####
Re_plot <- ggplot() +
  geom_rect(data = plot_Re_limits, aes(xmin = as_date("2021-03-21"), xmax = as_date("2021-09-21"),
                                       ymin = min_Re, ymax = max_Re), alpha = 0.4, fill = "lightgrey") +
  geom_rect(data = plot_Re_limits, aes(xmin = as_date("2022-03-21"), xmax = as_date("2022-09-21"),
                                       ymin = min_Re, ymax = max_Re), alpha = 0.4, fill = "lightgrey") +
  geom_vline(xintercept = kof_changepoint %>% filter(lag_kof > 0) %>% pull(date) , alpha = 0.5, linetype = 'longdash') +
  geom_vline(xintercept = kof_changepoint %>% filter(lag_kof < 0) %>% pull(date) , alpha = 0.5, linetype = 'dotted') +
  geom_ribbon(data = all_Re_estimates, aes(x = date, ymin = CI_down_Re_estimate,
                                ymax = CI_up_Re_estimate, fill = virus),
              alpha = 0.4, show.legend = T) +
  geom_line(data = all_Re_estimates, 
            aes(x = date, y = Re_estimate, colour = virus), 
            alpha = 1, linewidth = 1.2,  show.legend = F) +
  labs(x = 'Date' , y='Re', colour = 'Virus', fill = 'Virus') +
  facet_wrap(vars(virus), ncol = 1, scale = 'free_y') +
  geom_hline(yintercept = 1) +
  scale_x_date(date_breaks = '3 months', date_labels = '%b %Y', 
               limits = c(date('2020-12-01'), date('2022-11-01'))) +
  scale_fill_viridis_d(option = "mako", end = 0.8) +
  scale_color_viridis_d(option = "mako", end = 0.8) +
  theme(legend.position = 'bottom')

Re_plot
ggsave('../figures/r_smooth.pdf', width = 12, height = 8)
ggsave('../figures/r_smooth.png', width = 12, height = 8)

###########################################################
# compute summary statistics - Table 2

all_Re_estimates %>%
  mutate(quartal = quarter(date),
         month = month(date),
         year = year(date)) %>%
  filter(year != 2020) %>%
  group_by(year, quartal, virus) %>%
  summarise(mean_Re = mean(Re_estimate),
            sd_Re = sd(Re_estimate),
            min_Re = min(Re_estimate),
            max_Re = max(Re_estimate)
            #mean_CI_bot = mean(CI_down_Re_estimate),
            #mean_CI_top = mean(CI_up_Re_estimate)
            ) %>%
  arrange(virus) %>%
  View()



