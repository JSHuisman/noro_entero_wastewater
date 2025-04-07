###########################################################
# Supplementary Fig. S6 Showing Re for 
# Norovirus (NoV GII), Enterovirus (EV), and SARS-COV-2 (COV)
# Run with EpiSewer
# 
# Author: Jana S. Huisman
###########################################################

library(tidyverse)
library(readxl)
library(lubridate)

#remotes::install_github("adrian-lison/EpiSewer", dependencies = TRUE)
#cmdstanr::check_cmdstan_toolchain()
#cmdstanr::install_cmdstan(cores = 2) # use more cores to speed up
#EpiSewer::sewer_compile()
library(EpiSewer)

library(patchwork)
library(latex2exp)

theme_set(theme_minimal() + theme(text = element_text(size = 20)))

###########################################################
## Read in Raw wastewater data for NoV GII and EV ####

raw_data <- read_xlsx('../data/WISE_wastewater_data.xlsx', sheet = 'Data', skip = 2,
                      na = c('#VALUE!', 'NA'))

cov_data <- read_csv('../data/sarscov2_data.csv') %>%
  mutate(flow = 1000*flow, # The flow data is reported in kL
         cov_conc = cov_load/flow)

# Rename columns and filter out NA load values
core_data <- raw_data %>%
  select(date, replicate, flow = `flowrate_(kL/day)`, 
         nov_conc = `NoV GII...12`, ev_conc = `EV...13`, 
         nov_load = `NoVG II`, ev_load = `EV...16`) %>%
  mutate(date = as_date(date)) %>%
  filter(!is.na(nov_load) & !is.na(ev_load)) %>%
  mutate(flow = as.numeric(flow),
         weekday = weekdays(date)) %>%
  arrange(date) %>%
  left_join(cov_data[,c('date', 'cov_load', 'cov_conc')], by = 'date') %>%
  arrange(date) %>%
  mutate(orig_data = TRUE) %>%
  complete(date = seq.Date(min(date), max(date), by = 'days'))

## The flow data can not have missing values - so we use the sars-cov-2 upload (daily measurements)
flow_data <- cov_data %>%
  select(date, flow) %>%
  mutate(orig_data = TRUE) %>%
  complete(date = seq.Date(min(date), max(date), by = 'days')) %>%
  mutate(across(where(is.numeric), ~ zoo::na.approx(.x, na.rm = F) )) 

###########################################################
## Run EpiSewer for SARS-CoV-2, Norovirus, Enterovirus data ####
###########################################################
# epi_cov_data = core_data %>% 
#   select(date, flow, concentration = cov_conc, weekday) %>%
#   filter(!is.na(concentration),
#          !is.na(flow))
# 
# ww_data_cov <- sewer_data(measurements = epi_cov_data[, c('date', 'concentration')], 
#                       flows = flow_data)
# 
# # Same values as in SuppFig_SLDs.R
# generation_dist_cov <- get_discrete_gamma_shifted(gamma_mean = 3, gamma_sd = 2.4, maxX = 12) #sars2_gtid
# incubation_dist_cov <- get_discrete_gamma(gamma_shape = 8.5, gamma_scale = 0.4, maxX = 10)
# shedding_dist_cov <- get_discrete_gamma(gamma_shape = 0.929639, gamma_scale = 7.241397, maxX = 30) #benefield
# load_per_case_cov <- 2.4e11 #1e+11
# 
# ww_assumptions_cov <- sewer_assumptions(
#   generation_dist = generation_dist_cov,
#   incubation_dist = incubation_dist_cov,
#   shedding_dist = shedding_dist_cov,
#   shedding_reference="infection",
# )
# 
# ww_shedding <- model_shedding(
#   load_per_case = load_per_case_assume(load_per_case = load_per_case_cov)
# )
# 
# ww_infections = model_infections(
#   infection_noise = infection_noise_estimate(overdispersion = TRUE)
# )
# 
# options(mc.cores = 6) # allow stan to use 4 cores, i.e. one for each chain
# ww_result_sars <- EpiSewer(
#   data = ww_data_cov,
#   assumptions = ww_assumptions_cov,
#   infections = ww_infections,
#   shedding = ww_shedding,
#   fit_opts = set_fit_opts(sampler = sampler_stan_mcmc(iter_warmup = 1000, iter_sampling = 1000, chains = 4))
# )
# saveRDS(ww_result_sars, '../proc_data/EpiSewer_sarscov2.rds')
# 
# ####Noro###################################################
# noro_data = core_data %>% 
#   select(date, flow, concentration = nov_conc, weekday) %>%
#   filter(!is.na(concentration),
#          !is.na(flow)) %>%
#   mutate(concentration = 1000*concentration)
# 
# ww_data <- sewer_data(measurements = noro_data[, c('date', 'concentration')], 
#                       flows = flow_data)
# 
# # Same values as in SuppFig_SLDs.R
# generation_dist <- get_discrete_gamma(gamma_mean = 3.1, gamma_sd = 1.8, maxX = 26)
# shedding_dist <- get_discrete_gamma(gamma_shape = 2.716096, gamma_scale = 1.492308, maxX = 30)
# load_per_case_nor <- 1.1e14
# 
# ww_assumptions <- sewer_assumptions(
#   generation_dist = generation_dist,
#   incubation_dist=c(1),
#   shedding_dist = shedding_dist,
#   shedding_reference="infection"
# )
# 
# ww_shedding <- model_shedding(
#   load_per_case = load_per_case_assume(load_per_case = load_per_case_nor)
# )
# 
# ww_infections = model_infections(
#   infection_noise = infection_noise_estimate(overdispersion = TRUE)
# )
# 
# options(mc.cores = 6) # allow stan to use 4 cores, i.e. one for each chain
# ww_result <- EpiSewer(
#   data = ww_data,
#   assumptions = ww_assumptions,
#   infections = ww_infections,
#   shedding = ww_shedding,
#   fit_opts = set_fit_opts(sampler = sampler_stan_mcmc(iter_warmup = 1000, iter_sampling = 1000, chains = 4))
# )
# saveRDS(ww_result, file = '../proc_data/EpiSewer_norovirus.rds') # after fixing the flows
# 
# ####Entero##############################################
# entero_data = core_data %>% 
#   select(date, flow, concentration = ev_conc, weekday) %>%
#   filter(!is.na(concentration),
#          !is.na(flow)) %>%
#   mutate(concentration = 1000*concentration) # The concentrations are reported per mL
# 
# ww_data_env <- sewer_data(measurements = entero_data[, c('date', 'concentration')], 
#                           flows = flow_data)
# 
# # Same values as in SuppFig_SLDs.R
# generation_dist_env <- get_discrete_gamma(gamma_mean = 3.7, gamma_sd = 2.6, maxX = 26)
# shedding_dist_env <- get_discrete_gamma(gamma_shape = 13.14, gamma_scale = 0.41, maxX = 30)
# load_per_case_env <- 8.9e12#1e13
# 
# ww_assumptions_env <- sewer_assumptions(
#   generation_dist = generation_dist_env,
#   incubation_dist=c(1),
#   shedding_dist = shedding_dist_env,
#   shedding_reference="infection"
# )
# 
# 
# ww_shedding <- model_shedding(
#   load_per_case = load_per_case_assume(load_per_case = load_per_case_env)
# )
# 
# ww_infections = model_infections(
#   infection_noise = infection_noise_estimate(overdispersion = TRUE)
# )
# 
# 
# options(mc.cores = 6) # allow stan to use 4 cores, i.e. one for each chain
# ww_result_env <- EpiSewer(
#   data = ww_data_env,
#   assumptions = ww_assumptions_env,
#   infections = ww_infections,
#   shedding = ww_shedding,
#   fit_opts = set_fit_opts(sampler = sampler_stan_mcmc(iter_warmup = 1000, iter_sampling = 1000, chains = 4))
# )
# 
# saveRDS(ww_result_env, file = '../proc_data/EpiSewer_enterovirus.rds') # max ~100 infections

###########################################################
ww_result <- readRDS('../proc_data/EpiSewer_norovirus.rds')
ww_result_env = readRDS('../proc_data/EpiSewer_enterovirus.rds')
ww_result_sars = readRDS('../proc_data/EpiSewer_sarscov2.rds')

all_ww_Re <- bind_rows(Norovirus = ww_result$summary$R, 
                       Enterovirus = ww_result_env$summary$R, 
                       'SARS-CoV-2' = ww_result_sars$summary$R, .id = 'virus') %>%
  mutate(virus = factor(virus, levels = c('Norovirus', 'Enterovirus', 'SARS-CoV-2')) ) 

kof_changepoint <- read_csv('../data/kof_changepoints.csv')


#### Plotting Re together #################################

Re_plot <- ggplot() +
  annotate(geom = 'rect', xmin = as_date("2021-03-21"), xmax = as_date("2021-09-21"),
          ymin = 0.5, ymax = 1.5, alpha = 0.4, fill = "lightgrey") +
  annotate(geom = 'rect', xmin = as_date("2022-03-21"), xmax = as_date("2022-09-21"),
           ymin = 0.5, ymax = 1.5, alpha = 0.4, fill = "lightgrey") +
  geom_vline(xintercept = kof_changepoint %>% filter(lag_kof > 0) %>% pull(date) , alpha = 0.5, linetype = 'longdash') +
  geom_vline(xintercept = kof_changepoint %>% filter(lag_kof < 0) %>% pull(date) , alpha = 0.5, linetype = 'dotted') +
  geom_ribbon(data = all_ww_Re, aes(x = date, ymin = lower_0.95,
                                ymax = upper_0.95, fill = virus),
              alpha = 0.4, show.legend = T) +
  geom_line(data = all_ww_Re, 
            aes(x = date, y = mean, colour = virus), 
            alpha = 1, linewidth = 1.2,  show.legend = F) +
  labs(x = 'Date' , y='Re', colour = 'Virus', fill = 'Virus') +
  facet_wrap(vars(virus), ncol = 1, scale = 'free_y') +
  geom_hline(yintercept = 1) +
  scale_x_date(date_breaks = '3 months', date_labels = '%b %y', 
               limits = c(date('2020-12-10'), date('2022-11-01'))) +
  scale_fill_viridis_d(option = "mako", end = 0.8) +
  scale_color_viridis_d(option = "mako", end = 0.8) +
  theme(legend.position = 'bottom')

Re_plot
ggsave('../figures/r_episewer.png', width = 12, height = 8)


###########################################################
# compute summary statistics - Supplementary Table S2
all_ww_Re %>%
  mutate(quartal = quarter(date),
         year = year(date)) %>%
  filter(year != 2020) %>%
  group_by(year, quartal, virus) %>%
  summarise(mean_Re = mean(mean),
            sd_Re = sd(mean),
            min_Re = min(mean),
            max_Re = max(mean)
            #mean_CI_bot = mean(CI_down_Re_estimate),
            #mean_CI_top = mean(CI_up_Re_estimate)
  ) %>%
  arrange(virus) %>%
  View()

