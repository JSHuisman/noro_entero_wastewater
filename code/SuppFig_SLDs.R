###########################################################
# Supp Fig. Showing SLDs for
# Norovirus (NoV GII), Enterovirus (EV), and SARS-COV-2 (COV)
#
# Author: Jana S. Huisman
###########################################################

library(tidyverse)
library(readxl)
#library(lubridate)

library(patchwork)

# plotting defaults
theme_set(theme_minimal() +
            theme(text = element_text(size = 20))
)

###########################################################
# Helpful functions for the Delay/Shedding Load Distributions #####

# find gamma parameters from mean/sd of distribution
getGammaParams <- function(meanParam, sdParam){
  shapeParam <- meanParam^2 / (sdParam^2)
  scaleParam <- (sdParam^2) / meanParam
  return(list(shape = shapeParam, scale = scaleParam))
}

# fit a Gamma SLD by the first two moments
fit_gamma <- function(x, y){
  ## linear approximation, gives an unnormalized density
  d.empirical <- approxfun(x=x, y=y, yleft=max(y), yright=0)
  
  ## normalization factor (scaling for load)
  Z <- integrate(d.empirical, 0, 1000)$value
  
  ## compute the first two moments and derive sd and mean
  f.M1 <- function(x) d.empirical(x)/Z*x
  f.M2 <- function(x) d.empirical(x)/Z*x^2
  mean <- integrate(f.M1, 0, 1000)$value
  sd <- sqrt(integrate(f.M2, 0, 1000)$value - mean^2)
  
  ## calculate matching shape and scale of Gamma distribution
  gamma_params <- getGammaParams(mean, sd)
  
  return(list(norm_factor = Z, moments = c(mean, sd), params = gamma_params))
}

###########################################################
## The raw data used to fit the SLDs #####

## =================================
## Digitized data from Fig 3 of https://doi.org/10.1016/j.vaccine.2019.12.057
## Fan et al, 2020. EV71-CA16 enterovirus vaccine tested in rhesus macaques, use placebo group
## =================================
raw_fan <- data.frame(x = c(0, 2, 4, 6, 8), y = c(0, 10, 155, 225, 130), 
                      norm_factor = 910,
                      sld = 'fan', virus = 'Enterovirus')

## =================================
## Digitized data from Fig 2 of doi.org/10.1093/infdis/jiu497
## Bernstein et al, 2014. Norovirus vaccine against experimental human GII.4 Virus Illness; Placebo patients
## =================================
raw_bernstein <- data.frame(x = c(1, 2, 3, 4, 10), y = c(20000, 1200000, 1300000, 1250000, 120000),
                            norm_factor = 8544928,
                            sld = 'bernstein', virus = 'Norovirus')

# Combined raw data for SLDs#
raw_sld_points <- bind_rows(raw_fan, raw_bernstein) %>%
  mutate(norm_y = y/norm_factor) %>%
  mutate(virus = factor(virus, levels = c('Norovirus', 'Enterovirus', 'SARS-CoV-2')))

###########################################################
## Fit the SLDs #####

## Enterovirus
entero_dist = fit_gamma(raw_fan$x, raw_fan$y)

plot(raw_fan$x, raw_fan$y)
lines(0.1:50, entero_dist$norm_factor*dgamma(0.1:50, shape=entero_dist$params$shape, scale=entero_dist$params$scale), col=3)

## Norovirus
noro_dist = fit_gamma(raw_bernstein$x, raw_bernstein$y)

plot(raw_bernstein$x, raw_bernstein$y)
lines(0.1:50, noro_dist$norm_factor*dgamma(0.1:50, shape=noro_dist$params$shape, scale=noro_dist$params$scale), col=3)

###########################################################
## Summary of the SLD/Delay parameters #####

# specific distributions we use in the paper
getCountParams <- function(obs_type){
  switch(obs_type,
         incubation = getGammaParams(5.3, 3.2),
         zero = list(shape = 0, scale = 0),
         bernstein = list(shape = 2.716096, scale = 1.492308),
         fan = list(shape = 13.14, scale = 0.41),
         benefield = list(shape = 0.929639, scale = 7.241397),
         
         entero_gtid = list(shape = 2.03, scale = 1.83, dist = "Enterovirus"),
         noro_gtid = list(shape = 3.35, scale = 1.09, dist = "Norovirus"),
         sars2_gtid = list(shape = 1.56, scale = 1.92, dist = "SARS-CoV-2"))
}

###########################################################
## The fitted SLDs used in Re estimation #####

dist_for_plotting <- data.frame()
for(sld in c('bernstein', 'fan', 'benefield')){
  dist_params <- getCountParams(sld)
  x = seq(0, 30, 0.1)
  y = dgamma(x, shape = dist_params$shape, scale = dist_params$scale)
  if(sld == 'benefield'){
    # add the incubation period
    y = y + dgamma(x, shape = 2.74, scale = 1.93)
  }
  
  new_dist <- data.frame(x =x, y=y, dist = sld)  
  dist_for_plotting = bind_rows(dist_for_plotting, new_dist)
}

dist_for_plotting <- dist_for_plotting %>%
  mutate(virus = case_when(
    dist == 'bernstein' ~ "Norovirus",
    dist == 'fan' ~ "Enterovirus",
    dist == 'benefield' ~ "SARS-CoV-2"
  )) %>%
  mutate(virus = factor(virus, levels = c('Norovirus', 'Enterovirus', 'SARS-CoV-2')))

###########################################################
## Plot SLDs and raw data #####

SLD_plot <- ggplot(dist_for_plotting) +
  geom_line(aes(x, y, color = virus), linewidth = 1, show.legend = F) +
  geom_point(data = raw_sld_points, aes(x, y = norm_y, color = virus), size = 2, show.legend = F) +
  facet_wrap(vars(virus), ncol =1) +
  coord_cartesian(xlim = c(0, 15), ylim = c(0, 0.3)) +
  labs(x = 'Day after infection', y = 'Probability of shedding', color = 'Distribution') +
  scale_color_viridis_d(option = "mako", end = 0.8)

SLD_plot
ggsave('../figures/SLDs.png')

###########################################################
## Plot GTID and raw data #####


gtid_for_plotting <- data.frame()
for(gtid in list('entero_gtid', 'noro_gtid', 'sars2_gtid')){
  gtid_dist = getCountParams(gtid)
  x = seq(0, 30, 0.1)
  y = dgamma(x, shape = gtid_dist$shape, scale = gtid_dist$scale)
  new_dist <- data.frame(x =x, y=y, virus = gtid_dist$dist)  
  gtid_for_plotting = bind_rows(gtid_for_plotting, new_dist) %>%
    mutate(virus = factor(virus, levels = c('Norovirus', 'Enterovirus', 'SARS-CoV-2')))
}


GTID_plot <- ggplot(gtid_for_plotting) +
  geom_line(aes(x, y, color = virus), linewidth = 1, show.legend = F) +
  facet_wrap(vars(virus), ncol =1) +
  coord_cartesian(xlim = c(0, 15), ylim = c(0, 0.3)) +
  labs(x = 'Day after infection', y = 'Probability of infection', color = 'Distribution') +
  scale_color_viridis_d(option = "mako", end = 0.8)

GTID_plot
ggsave('../figures/GTIDs.png')

###########################################################
## Combine SLD and GTID plot #####

SLD_plot + GTID_plot +
  plot_annotation(tag_levels = 'A')

ggsave('../figures/epi_distributions.png')
