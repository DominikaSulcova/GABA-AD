# script libraries
library(readr)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(gridExtra)
library(grid)
library(gt)
library(xlsx)

# color palette - red/blue/purple(dark/light) + black/gray/light gray
col_palette = c('red3', 'tomato1','turquoise4', 'cadetblue1', '#8E3AC9', '#E9CCFD', 
                '#000000','#CAC9C9', "gray96", '#30C2D3', '#32086D')

# ***************************************************************
# TIMESERIES - ERP
# ***************************************************************

# ------------------------------------------------------
# PREPARE DATA
# ------------------------------------------------------
# import TEP datasets
TEP_timeseries <- read.csv('GABA_timeseries.csv')

# name timestamps
timestamps <- read.csv('timestamps.csv')
head(timestamps[, 1:10])
length(timestamps)

# create factors and order them
TEP_timeseries %>% 
  as_tibble() %>% 
  mutate(
    treatment = factor(treatment, levels = c('alprazolam', 'placebo'), ordered = TRUE),
    time = factor(time, levels = c('pre', 'post'), ordered = TRUE),
    stimulus = factor(stimulus, levels = c('TS', 'CS', 'ppTMS'), ordered = TRUE),
    channel = factor(channel, levels = c('Cz', 'target'), ordered = TRUE))

TEP_timeseries %>% 
  arrange(channel) %>%
  arrange(time) %>%
  arrange(treatment) %>%
  arrange(subject) %>%
  arrange(stimulus) 

# name timepoint variables as timestamps
oldnames <- sprintf("Var%s",1:ncol(timestamps))
TEP_timeseries <- TEP_timeseries %>% 
  rename_at(vars(oldnames), ~ as.character(timestamps))

TEP_timeseries_long <- TEP_timeseries_long %>% 
  mutate(subject = factor(subject))

# transform the datasets into ggplot friendly (vertical) format 
TEP_timeseries_long <- TEP_timeseries %>% 
  pivot_longer(c(-subject, -treatment, -time, -stimulus, -channel), names_to = 'timepoint', values_to = 'amplitude') 

TEP_timeseries_long$timepoint <- as.numeric(TEP_timeseries_long$timepoint)


# calculate grand average and confidence interval
# ------------------------------------------------------
# create objects for factors
Treatments <- unique(TEP_timeseries$treatment)
Times <- unique(TEP_timeseries$time)
Stimuli<- unique(TEP_timeseries$stimulus)
Channels <- unique(TEP_timeseries$channel)

# create an object ERPplot with grand averages + CI_u with upper confidence interval, CI_l with lower CI
ERPplot <- tibble()
CI_u <- tibble()
CI_l <- tibble()
for (a in 1:length(Stimuli)) {
  a0 <- 1+8*(a-1) 
  a1 <- 8*a
  ERPplot[a0:a1 ,'stimulus'] <- Stimuli[a]
  CI_u[a0:a1 ,'stimulus'] <- Stimuli[a]
  CI_l[a0:a1,'stimulus'] <- Stimuli[a]
  for (b in 1:length(Treatments)) {
    b0 <- 1+8*(a-1)+4*(b-1)
    b1 <- 4*b + 8*(a-1)
    ERPplot[b0:b1,'treatment'] <- Treatments[b]
    CI_u[b0:b1,'treatment'] <- Treatments[b]
    CI_l[b0:b1,'treatment'] <- Treatments[b]
    for (c in 1:length(Times)) {
      c0 <- 1+8*(a-1)+4*(b-1)+2*(c-1)
      c1 <- 2*c+8*(a-1)+4*(b-1)
      ERPplot[c0:c1,'time'] <- Times[c]
      CI_u[c0:c1,'time'] <- Times[c]
      CI_l[c0:c1,'time'] <- Times[c]
      for (d in 1:length(Channels)) {
        d1 <- d+8*(a-1)+4*(b-1)+2*(c-1)
        ERPplot[d1,'channel'] <- Channels[d]
        CI_u[d1,'channel'] <- Channels[d]
        CI_l[d1,'channel'] <- Channels[d]
        for (t in 1:ncol(timestamps)) {
          CI <- Rmisc::CI(unlist(TEP_timeseries %>% 
                                   filter(stimulus == Stimuli[a]) %>%
                                   filter(treatment == Treatments[b]) %>%
                                   filter(time == Times[c]) %>%
                                   filter(channel == Channels[d]) %>%
                                   select(as.character(timestamps[1,t]))), 
                          ci = 0.95)
          ERPplot[d1, as.character(timestamps[1,t])] <- CI["mean"]
          CI_u[d1, as.character(timestamps[1,t])] <- CI["upper"]
          CI_l[d1 ,as.character(timestamps[1,t])] <- CI["lower"]
        }
        
      }
      
    }
    
  }
  
}
rm(CI)

# create the final dataset
# ------------------------------------------------------
# transform the datasets into ggplot friendly (vertical) format 
ERPplot <- ERPplot %>% 
  pivot_longer(c(-treatment, -time, -stimulus, -channel), names_to = 'timepoint', values_to = 'amplitude')
CI_u <- CI_u %>%
  pivot_longer(c(-treatment, -time, -stimulus, -channel), names_to = 'timepoint', values_to = 'CI_upper')
CI_l <- CI_l %>%
  pivot_longer(c(-treatment, -time, -stimulus, -channel), names_to = 'timepoint', values_to = 'CI_lower')

# append the CI values
ERPplot <- cbind(ERPplot, CI_u$CI_upper)
ERPplot <- cbind(ERPplot, CI_l$CI_lower)
ERPplot <- ERPplot %>% 
  rename_at(vars(c('CI_u$CI_upper', 'CI_l$CI_lower')), ~ as.character(c('CI_upper', 'CI_lower')))

ERPplot <- ERPplot %>%
  mutate(
    treatment = factor(treatment, levels = c('alprazolam', 'placebo'), ordered = TRUE),
    time = factor(time, levels = c('pre', 'post'), ordered = TRUE),
    stimulus = factor(stimulus, levels = c('TS', 'CS', 'ppTMS'), ordered = TRUE),
    channel = factor(channel, levels = c('Cz', 'target'), ordered = TRUE))


# make sure all values are numeric
ERPplot$timepoint <- as.numeric(as.character(ERPplot$timepoint))


# ------------------------------------------------------
# BASELINE TEPs
# ------------------------------------------------------
# choose data 
ERP_TS_bl_Cz <- ERPplot[ERPplot$time == 'pre' & ERPplot$channel == 'Cz' & ERPplot$stimulus == 'TS' & 
                          ERPplot$timepoint >= -0.05 & ERPplot$timepoint <= 0.3, ]

# plot baseline
ERP_baseline <- ggplot(data = ERP_TS_bl_Cz,  
                    aes(x = timepoint, y = amplitude, 
                        colour = treatment)) + 
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper, fill = treatment), 
              alpha = 0.15, linetype = 0) + 
  geom_line(size = 2.25) +
  coord_cartesian() +
  labs(title = 'baseline TS-TEPs: Cz electrode', 
       x = "time (s)", 
       y = "amplitude (µV)",
       size = 20) +
  font(object = 'title', size = 24, face = 2) + 
  scale_x_continuous(minor_breaks = seq(-0.05, 0.3, 0.01), breaks = seq(-0.1, 0.3, 0.05)) +
  scale_y_continuous(breaks = seq(-20, 20, 2.5)) +
  theme(
    panel.background = element_rect(fill = col_palette[9],
                                    colour = col_palette[9],
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "white"), 
    panel.grid.minor = element_line(size = 0.15, linetype = 'solid',
                                    colour = "white"),
    axis.text = element_text(size = 18, colour = col_palette[7]),
    axis.title = element_text(size = 18),
    legend.position = 'none') +
  geom_rect(mapping = aes(xmin = -0.0005, xmax = 0.01, ymin = -10, ymax = 10), alpha = 0.5, color = 0)



# CHANGE INDUCED BY MEDICATION 
# ------------------------------------------------------
# choose data 
ERP_TS_alprazolam_Cz <- ERPplot[ERPplot$treatment == 'alprazolam' & ERPplot$channel == 'Cz' & ERPplot$stimulus == 'TS' & 
                                  ERPplot$timepoint >= -0.05 & ERPplot$timepoint <= 0.3, ]
# Cz
ERP_change <- ggplot(data = ERP_TS_alprazolam_Cz,  
                    aes(x = timepoint, y = amplitude, 
                        colour = time)) +
  scale_color_manual(values = c(col_palette[2], col_palette[1])) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper, fill = time), 
              alpha = 0.15, linetype = 0) + 
  scale_fill_manual(values = c(col_palette[2], col_palette[1])) +
  geom_line(size = 2.25) +
  coord_cartesian() +
  labs(title = 'TS-TEPs: Cz electrode - alprazolam', 
       x = "time (s)", 
       y = "amplitude (µV)") +
  font(object = 'title', size = 24, face = 2) + 
  scale_x_continuous(minor_breaks = seq(-0.05, 0.3, 0.01), breaks = seq(-0.1, 0.3, 0.05)) +
  scale_y_continuous(breaks = seq(-20, 20, 2.5)) +
  theme(
    panel.background = element_rect(fill = col_palette[9],
                                    colour = col_palette[9],
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "white"), 
    panel.grid.minor = element_line(size = 0.15, linetype = 'solid',
                                    colour = "white"),
    axis.text = element_text(size = 18, colour = col_palette[7]),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 25),
    legend.title = element_text(size = 25, face = 2),
    legend.position = c(0.8, 0.35),
    legend.justification = c('right', 'top')) +
  geom_rect(mapping = aes(xmin = -0.0005, xmax = 0.01, ymin = -10, ymax = 10), alpha = 0.5, color = 0)

# ***************************************************************
# PAIRED PLOTS
# ***************************************************************

# ------------------------------------------------------
# PREPARE DATA
# ------------------------------------------------------
load('TEP_data_tidy.RData')

# select rMT data
data_rMT <- TEP_data_tidy %>%
  select(c(1, 2, 3, 6)) %>%
  distinct()

# description stats
# ------------------------------------------------------
# append subject parameters
parameters <- read.csv('YC_parameters.csv') %>%
  slice(rep(1:nrow(parameters), each = 4)) 

parameters <- as_tibble(parameters[, 1:4])

data_all <- bind_cols(parameters, data_rMT)

data_all <- data_all[ , -5] %>%
  mutate(subject = factor(subject)) %>%
  mutate(Sex = factor(Sex)) %>%
  mutate(Sex = fct_recode(Sex,
                          'male' = '0' ,
                          'female' = '1'))

# ------------------------------------------------------
# PLOT rMT 
# ------------------------------------------------------
# plot rMT by gender
# ------------------------------------------------------
bar_rMT <- ggbarplot(data_all, x = 'time', y = 'rMT',
                     fill = 'Sex', palette = col_palette[c(3, 1)], 
                     add = "mean_sd", 
                     position = position_dodge(0.8)) + 
  scale_x_discrete(limits=c('pre', 'post')) +
  labs(x = 'time relative to medication',
       y = "rMT (%MSO)") +
  guides(fill=guide_legend(title="gender")) +
  theme(
    legend.position = 'right',
    legend.title = element_text(size = 24, face = 2),
    legend.text = element_text(size = 24),
    axis.text.x = element_text(vjust = -0.1, size = 24),
    axis.text.y = element_text(size = 24),
    axis.title = element_text(size = 24),
    panel.background = element_rect(fill = "gray96",
                                    colour = "gray96",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major.y = element_line(size = 0.5, linetype = 'solid',
                                      colour = "white"))

# compare baseline rMT between conditions
# ------------------------------------------------------
data_rMT_pre <- data_rMT[data_rMT$time == 'pre', ]

data_rMT_pre %>%
  group_by(treatment) %>%
  get_summary_stats(rMT, type = "mean_sd")

ggboxplot(data_rMT_pre, x = 'treatment', y = 'rMT', 
          fill = 'treatment', palette = col_palette[c(2, 4)],
          order = c('alprazolam', 'placebo'),
          xlab = 'treatment', ylab = 'rMT (%MSO)')

# compute the difference and test for normality
d <- with(data_rMT_pre, 
          rMT[treatment == 'alprazolam'] - rMT[treatment == 'placebo'])
ggqqplot(d)
shapiro.test(d) 

# d is not normally distributed --> compare rMT using paired Wilcoxon test
data_pre_alprazolam <- subset(data_rMT_pre,  treatment == "alprazolam", rMT,
                              drop = TRUE)
data_pre_placebo <- subset(data_rMT_pre,  treatment == 'placebo', rMT,
                           drop = TRUE)
res <- wilcox.test(data_pre_alprazolam, data_pre_placebo, paired = TRUE)
res

# paired plot
rMT_baseline <- ggpaired(data_rMT_pre, x = 'treatment', y = 'rMT',
                         fill = 'treatment', alpha = 0.25, line.color = 'grey55', line.size = 0.25, 
                         size = 1.5, point.size = 3, palette = col_palette[c(2, 4)],
                         id = 'subject') +
  stat_compare_means(paired = TRUE, label.x = 1, label.y = 64, size = 7) +
  scale_x_discrete(limits=c('alprazolam', 'placebo')) + 
  labs(title = 'rMT: baseline recording',
       x ='treatment',
       y = 'rMT (%MSO)') +
  font(object = 'title', size = 24, face = 2) +
  theme(
    legend.position = 'none',
    axis.text.x = element_text(vjust = -0.1, size = 24),
    axis.text.y = element_text(size = 24),
    axis.title = element_text(size = 24),
    panel.background = element_rect(fill = "gray96",
                                    colour = "gray96",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major.y = element_line(size = 0.5, linetype = 'solid',
                                      colour = "white")) 