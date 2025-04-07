###########################################################
# Fig. 3 Showing % change in Re for 
# Norovirus (NoV GII), Enterovirus (EV), and SARS-COV-2 (COV)
#
# Author: Jana S. Huisman
###########################################################

library(tidyverse)
library(readxl)
library(lubridate)

library(patchwork)
library(latex2exp)

theme_set(theme_minimal() + theme(text = element_text(size = 20)))

###########################################################
## Load data for Re 
Re_ww <- read_csv('../proc_data/Re_ww.csv') %>%
  mutate(virus = factor(virus, levels = c('Norovirus', 'Enterovirus', 'SARS-CoV-2'))   ) 

###########################################################
## Load data on covid restrictions ####
#kof_link = 'https://datenservice.kof.ethz.ch/api/v1/public/sets/stringency_plus_web?mime=csv&df=Y-m-d'
all_kof_df <- read_csv('../data/kof_data_export_2023-08-10_18_33_19.csv')

kof_df <- all_kof_df %>%
  select(date, kof = 'ch.kof.stringency.zh.stringency_plus') %>%
  filter(date >= as_date('2020-12-18'),
         date <= as_date('2022-09-26')) %>%
  mutate(lag_kof = kof - lag(kof, 1),
         delta_kof = (lag_kof)/kof,
         week_date = ymd(floor_date(date, "weeks", week_start = 1)) )


## determine when kof index changed ####
kof_changepoint <- kof_df %>%
  filter(lag_kof != 0)

write_csv(kof_changepoint, '../data/kof_changepoints.csv')

# # check whether the changes in kof are "independent", i.e. occurred more than 14 days apart
# kof_changepoint <- kof_changepoint %>%
#   mutate(ind_point = ifelse(date - lag(date, 1) > 14, 1, 0),
#          ind_point = ifelse(is.na(ind_point), 1, ind_point),
#          ind_point = ifelse(ind_point == 0 & lag(ind_point == 0) & date - lag(date, 2) > 14, 1, ind_point),
#          ind_point = factor(ind_point)) %>%
#   mutate(consol_lag_kof = ifelse(lead(ind_point) == 0, lag_kof + lead(lag_kof), lag_kof))

kof_changepoint

###########################################################
# # Plot kof index over time ####
# kof_afo_time <- ggplot(kof_df) +
#   geom_line(aes(x = date, y = kof)) +
#   labs(x = 'Date', y = 'KOF Stringency \nindex') +
#   coord_cartesian(ylim = c(0, 80))
# 
# kof_afo_time
# ggsave('../figures/kof_time.pdf', width = 12, height = 6)

###########################################################
## determine change in Re ####
changes_Re <- Re_ww %>%
  group_by(virus) %>%
  mutate(R_7ago = Re_estimate - lag(Re_estimate, 7),
         R_high_7ago = CI_up_Re_estimate - lag(CI_up_Re_estimate, 7),
         R_low_7ago = CI_down_Re_estimate - lag(CI_down_Re_estimate, 7),
         #
         R7 = lead(Re_estimate, 7) - Re_estimate,
         R_high7 = lead(CI_up_Re_estimate, 7) - CI_up_Re_estimate,
         R_low7 = lead(CI_down_Re_estimate, 7) - CI_down_Re_estimate,
         #
         R14 = lead(Re_estimate, 14) - Re_estimate,
         R_high14 = lead(CI_up_Re_estimate, 14) - CI_up_Re_estimate,
         R_low14 = lead(CI_down_Re_estimate, 14) - CI_down_Re_estimate,
         #
         R21 = lead(Re_estimate, 21) - Re_estimate,
         R_high21 = lead(CI_up_Re_estimate, 21) - CI_up_Re_estimate,
         R_low21 = lead(CI_down_Re_estimate, 21) - CI_down_Re_estimate,
         #
         R28 = lead(Re_estimate, 14) - lag(Re_estimate, 14),
         R_high28 = lead(CI_up_Re_estimate, 14) - lag(CI_up_Re_estimate, 14) ,
         R_low28 = lead(CI_down_Re_estimate, 14) - lag(CI_down_Re_estimate, 14)) %>%
  right_join(kof_changepoint, by = 'date')

###########################################################
## Plot these changes ####

# Delta Re vs delta KOF - main text Fig 3A #####
deltaR_deltaKOF_plot <- ggplot(changes_Re) +
  geom_point(aes(x = lag_kof, y = R14, color = virus), size = 2, show.legend = T) +
  geom_errorbar(aes(x = lag_kof, ymin = R_low14, ymax = R_high14, color = virus), linewidth = 1, show.legend = T) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  labs(x = TeX("$\\Delta KOF$"), y = TeX("$\\Delta R_e$") ) +
  facet_wrap(vars(virus), ncol = 1, scale = 'free_y') +
  theme(legend.position = 'none') +
  scale_fill_viridis_d(option = "mako", end = 0.8) +
  scale_color_viridis_d(option = "mako", end = 0.8) +
  coord_cartesian(xlim = c(-20, 15), ylim = c(-0.45, 0.45))

deltaR_deltaKOF_plot
ggsave('../figures/Rchange_KOFchange.png', width = 6, height = 8)

## Delta Re vs time - main text Fig 3B ########
plot_changeRe_limits <- changes_Re %>%
  group_by(virus) %>%
  summarise(min_Re = min(R_low14, na.rm = T),
            max_Re = max(R_low14, na.rm = T))

deltaR_time_plot <- ggplot(changes_Re) +
  geom_rect(data = plot_changeRe_limits, aes(xmin = as_date("2021-03-01"), xmax = as_date("2021-09-01"),
                                       ymin = -0.2, ymax = 0.2), alpha = 0.4, fill = "lightgrey") +
  geom_rect(data = plot_changeRe_limits, aes(xmin = as_date("2022-03-01"), xmax = as_date("2022-09-01"),
                                       ymin = -0.2, ymax = 0.2), alpha = 0.4, fill = "lightgrey") +
  geom_vline(xintercept = kof_changepoint %>% filter(lag_kof > 0) %>% pull(date) , alpha = 0.5, linetype = 'longdash') +
  geom_vline(xintercept = kof_changepoint %>% filter(lag_kof < 0) %>% pull(date) , alpha = 0.5, linetype = 'dotted') +
  
  geom_segment(aes(x=date, y=0, xend=date, yend=R14, color = virus), 
               arrow = arrow(length=unit(.2, 'cm')), linewidth = 1, show.legend = F) +
  facet_wrap(vars(virus), ncol = 1) +
  scale_color_viridis_d(option = "mako", end = 0.8) +
  scale_x_date(date_breaks = '3 months', date_labels = '%b %Y', 
               limits = c(date('2020-12-01'), date('2022-04-01'))) +
  labs(x = 'Date', y = TeX("$\\Delta R_e$"), color = 'Virus') +
  theme(legend.position = 'bottom')

deltaR_time_plot
ggsave('../figures/Rchange_time.png', width = 12, height = 8)


deltaR_time_plot + deltaR_deltaKOF_plot +
  plot_layout(nrow = 1, widths = c(2, 1), guides = 'collect') +
  plot_annotation(tag_levels = c('A'))
ggsave('../figures/Rchange_time_KOF.png', width = 12, height = 8)

###########################################################
## plot correlation between Re changes - main text Figure 3C
# These can all be assessed for a time window of 7, 14, or 21 days

changes_Re_plotdf <- changes_Re %>%
  select(date, virus, R7, R14, R21) %>%
  pivot_wider(id_cols = c('date'), names_from = 'virus', values_from = c('R7', 'R14', 'R21'))


noro_entero = ggplot(changes_Re_plotdf) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
   geom_point(aes(x = R14_Norovirus, y = R14_Enterovirus)) +
   geom_smooth(aes(x = R14_Norovirus, y = R14_Enterovirus), method='lm', color = 'black', alpha = 0.2) +
  #geom_point(aes(x = R7_Norovirus, y = R7_Enterovirus)) +
  #geom_smooth(aes(x = R7_Norovirus, y = R7_Enterovirus), method='lm', color = 'black', alpha = 0.2) +
  labs(x = TeX("$\\Delta R_e$ Norovirus"), y = TeX('$\\Delta R_e$ Enterovirus')) +
  coord_cartesian(xlim = c(-0.4, 0.4), ylim = c(-0.4, 0.4))
  #coord_cartesian(xlim = c(-0.2, 0.2), ylim = c(-0.2, 0.2))
noro_entero
ggsave('../figures/noro_entero_rchange.pdf')

#model1 = lm(R7_Norovirus ~ R7_Enterovirus, changes_Re_plotdf)
model1 = lm(R14_Norovirus ~ R14_Enterovirus, changes_Re_plotdf)
#model1 = lm(R21_Norovirus ~ R21_Enterovirus, changes_Re_plotdf)
summary(model1)

noro_sars = ggplot(changes_Re_plotdf) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
   geom_point(aes(x = R14_Norovirus, y = `R14_SARS-CoV-2`)) +
   geom_smooth(aes(x = R14_Norovirus, y = `R14_SARS-CoV-2`), method='lm', color = 'black', alpha = 0.2) +
  #geom_point(aes(x = R7_Norovirus, y = `R7_SARS-CoV-2`)) +
  #geom_smooth(aes(x = R7_Norovirus, y = `R7_SARS-CoV-2`), method='lm', color = 'black', alpha = 0.2) +
  labs(x = TeX("$\\Delta R_e$ Norovirus"), y = TeX('$\\Delta R_e$ SARS-CoV-2')) +
  coord_cartesian(xlim = c(-0.4, 0.4), ylim = c(-0.4, 0.4))
  #coord_cartesian(xlim = c(-0.2, 0.2), ylim = c(-0.2, 0.2))
noro_sars
ggsave('../figures/noro_sars_rchange.pdf')

#model2 = lm(R7_Norovirus ~ `R7_SARS-CoV-2`, changes_Re_plotdf)
model2 = lm(R14_Norovirus ~ `R14_SARS-CoV-2`, changes_Re_plotdf)
#model2 = lm(R21_Norovirus ~ `R21_SARS-CoV-2`, changes_Re_plotdf)
summary(model2)

entero_sars = ggplot(changes_Re_plotdf) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_point(aes(x = R14_Enterovirus, y = `R14_SARS-CoV-2`)) +
  #geom_smooth(aes(x = R14_Enterovirus, y = `R14_SARS-CoV-2`), method='lm', color = 'black', alpha = 0.2) +
  #geom_point(aes(x = R7_Enterovirus, y = `R7_SARS-CoV-2`)) +
  #geom_smooth(aes(x = R7_Enterovirus, y = `R7_SARS-CoV-2`), method='lm', color = 'black', alpha = 0.2) +
  labs(x = TeX("$\\Delta R_e$ Enterovirus"), y = TeX('$\\Delta R_e$ SARS-CoV-2')) +
  coord_cartesian(xlim = c(-0.4, 0.4), ylim = c(-0.4, 0.4))
  #coord_cartesian(xlim = c(-0.2, 0.2), ylim = c(-0.2, 0.2))

entero_sars
ggsave('../figures/entero_sars_rchange.pdf')

#model3 = lm(R7_Enterovirus ~ `R7_SARS-CoV-2`, changes_Re_plotdf)
model3 = lm(R14_Enterovirus ~ `R14_SARS-CoV-2`, changes_Re_plotdf)
#model3 = lm(R21_Enterovirus ~ `R21_SARS-CoV-2`, changes_Re_plotdf)
summary(model3)

noro_sars + entero_sars +  noro_entero + 
  plot_layout(ncol = 3, guides = 'collect') +
  plot_annotation(tag_levels = c('A')) +
  theme(axis.text = element_text(size = 20))

#ggsave('../figures/rchange_viruses_d14.png', width = 12, height = 4)
ggsave('../figures/rchange_viruses_d7.png', width = 12, height = 4)


## Combine into Main text Fig 3 #######

(deltaR_time_plot + deltaR_deltaKOF_plot +
    plot_layout(widths = c(2, 1)) )/
  (noro_sars + entero_sars +  noro_entero) +
  plot_layout(nrow = 2, heights = c(2, 1), guides = 'collect') +
  plot_annotation(tag_levels = c('A'))

ggsave('../figures/Rchange_time_KOF_correlation.png', width = 12, height = 8)

###########################################################
## Investigate correlation between weekly changes in Re - supplementary Fig S7 ####
# we calculate mean weekly Re to remove some autocorrelation 

weekly_kof_df <- kof_df %>%
  group_by(week_date) %>%
  summarise(mean_kof = mean(kof, na.rm = T),
            lag_kof = sum(lag_kof, na.rm = T),
            delta_kof = sum(delta_kof, na.rm = T))

weekly_changes_Re <- Re_ww %>%
  select(date, virus, Re_estimate) %>%
  mutate(week_date = ymd(floor_date(date, "weeks", week_start = 1)) ) %>%
  group_by(virus, week_date) %>%
  summarise(mean_Re = mean(Re_estimate, na.rm = T)) %>%
  mutate(delta_Re = mean_Re - lag(mean_Re, 1),
         delta_Re14 = lead(mean_Re,1) - lag(mean_Re, 1),
         lead_Re = lead(mean_Re,1)) %>%
  full_join(weekly_kof_df, by = 'week_date') %>%
  mutate(change_week = case_when(
            is.na(lag_kof) ~ 'No change',
            lag_kof == 0 ~ 'No change',
            lag_kof > 0 ~ 'Increased',
            lag_kof < 0 ~ 'Reduced'),
         kof_class = case_when(
           mean_kof == 0 ~ 'No measures',
           mean_kof < 30 ~ 'Some measures',
           mean_kof >= 30 ~ 'Strong measures'))

## Weekly Re values are fairly normally distributed around 1, no strong effect of KOF changes in week prior
ggplot(weekly_changes_Re) +
  geom_histogram(aes(x = lead_Re, fill = change_week, color = change_week), alpha = 0.2, position = 'identity') +
  facet_wrap(vars(virus)) +
  labs(x = "Re in week t+1", y = "# Weeks", color = "KOF change \nin week t", fill = "KOF change \nin week t")

#### Plotting subpanels of Fig S7 ####
weekly_change_Re_plotdf <- weekly_changes_Re %>%
  select(week_date, virus, mean_Re, delta_Re, lead_Re, change_week) %>%
  pivot_wider(id_cols = c('week_date', 'change_week'), names_from = 'virus', values_from = c('mean_Re', 'delta_Re', 'lead_Re'))

noro_entero = ggplot(weekly_change_Re_plotdf) +
  geom_vline(xintercept = 1, linetype = 'dashed') +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  geom_point(aes(x = mean_Re_Norovirus, y = mean_Re_Enterovirus, color = change_week), size = 2, alpha = 0.8) +
  geom_smooth(aes(x = mean_Re_Norovirus, y = mean_Re_Enterovirus), method='lm', color = 'black', alpha = 0.2) +
  coord_cartesian(xlim = c(0.75, 1.25), ylim = c(0.6, 1.4))+ 
  scale_color_manual(values = viridisLite::mako(6, end = 0.8)[c(1, 4, 6)]) +
  labs(x = TeX("$R_e$ Norovirus"), y = TeX('$R_e$ Enterovirus'), color = 'Stringency')
  
noro_entero
model1 = lm(mean_Re_Enterovirus ~ mean_Re_Norovirus, weekly_change_Re_plotdf)
summary(model1)

noro_sars = ggplot(weekly_change_Re_plotdf) +
  geom_vline(xintercept = 1, linetype = 'dashed') +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  geom_point(aes(x = mean_Re_Norovirus, y = `mean_Re_SARS-CoV-2`, color = change_week), size = 2, alpha = 0.8) +
  #geom_smooth(aes(x = mean_Re_Norovirus, y = `mean_Re_SARS-CoV-2`), method='lm', color = 'black', alpha = 0.2) +
  coord_cartesian(xlim = c(0.75, 1.25), ylim = c(0.6, 1.4))+ 
  scale_color_manual(values = viridisLite::mako(6, end = 0.8)[c(1, 4, 6)]) +
  labs(x = TeX("$R_e$ Norovirus"), y = TeX('$R_e$ SARS-CoV-2'), color = 'Stringency')

noro_sars
model2 = lm(mean_Re_Norovirus ~ `mean_Re_SARS-CoV-2`, weekly_change_Re_plotdf)
summary(model2)

sars_entero = ggplot(weekly_change_Re_plotdf) +
  geom_vline(xintercept = 1, linetype = 'dashed') +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  geom_point(aes(x = `mean_Re_SARS-CoV-2`, y = mean_Re_Enterovirus, color = change_week), size = 2, alpha = 0.8) +
  #geom_smooth(aes(x = `mean_Re_SARS-CoV-2`, y = mean_Re_Enterovirus), method='lm', color = 'black', alpha = 0.2) +
  coord_cartesian(xlim = c(0.75, 1.25), ylim = c(0.6, 1.4))+ 
  scale_color_manual(values = viridisLite::mako(6, end = 0.8)[c(1, 4, 6)]) +
  labs(x = TeX("$R_e$ SARS-CoV-2"), y = TeX('$R_e$ Enterovirus'), color = 'Stringency')

sars_entero
model3 = lm(`mean_Re_SARS-CoV-2` ~ mean_Re_Enterovirus, weekly_change_Re_plotdf)
summary(model3)

## Combine into supplementary Fig S7 ####
noro_sars + sars_entero +  noro_entero + 
  plot_layout(ncol = 3, guides = 'collect') +
  plot_annotation(tag_levels = c('A')) &
  theme(axis.text = element_text(size = 20),
        legend.position = 'bottom')

ggsave('../figures/rcorr_viruses.png', width = 12, height = 4)

