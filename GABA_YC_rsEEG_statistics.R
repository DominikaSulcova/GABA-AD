library(tidyverse)
library(ggplot2)
library(MASS)
library(dplyr)
library(ggpubr)
library(corrplot) 
library(visreg) 
library(rgl)
library(knitr)
library(scatterplot3d)
library(merTools)
library(readr)
library(nlme)
library(lmerTest)
library(car)
library(lattice)

# load data
folder_dir <- choose.dir(default = getwd(), caption = "Select folder for input/output")
data_all <- read.csv(paste(folder_dir, '\\GABA_YC_medication_rsEEG.csv', sep = ""))
head(data_all, 10)
str(data_all)

# choose variable
outcome_var <- 'delta'
vars_keep <- c('subject', 'medication', outcome_var)

# make into factors
data <- as_tibble(data_all[vars_keep]) 
data$subject <- factor(data$subject)
data$medication <- factor(data$medication, levels=c('placebo', 'alprazolam'), ordered = TRUE)

# rename dependent variable
data <- rename(data, DV_change = delta)
head(data, 10)

#  descriptive statistics 
describe <- function(x){
  print(round(c(mean = mean(x, na.rm = TRUE),
                median = median(x, na.rm = TRUE),
                min = min(x, na.rm = TRUE),
                max = max(x, na.rm = TRUE),
                std = sd(x, na.rm = TRUE),
                sem = sd(x, na.rm = TRUE)/sqrt(20)), digits = 4))}

describe(data$DV_change)
cat('grouped by medication: \n')
aggregate(DV_change ~ medication, data = data, FUN = describe)

# test for normality
# placebo
data_placebo <- data$DV_change[data$medication == 'placebo']
shapiro.test(data_placebo)
qqnorm(data_placebo)
qqline(data_placebo)

# alprazolam
data_alprazolam <- data$DV_change[data$medication == 'alprazolam']
shapiro.test(data_alprazolam)
qqnorm(data_alprazolam)
qqline(data_alprazolam)

# create the empty model
model_empty <- lme(
  DV_change ~ 1,
  random = ~ 1|subject,
  data = data,
  method = 'ML'
)
summary(model_empty)
intervals(model_empty)

# random intercept
model_ri <- lme(
  DV_change ~ medication,
  random = ~ 1|subject,
  data = data,
  method = 'ML'
)
summary(model_ri)
intervals(model_ri)

# random intercept + random slope
model_rirs <- lme(
  DV_change ~ medication,
  random = ~ medication|subject,
  data = data,
  method = 'ML'
)
summary(model_rirs)
intervals(model_rirs)

# compare all models using ANOVA
anova(model_empty, model_ri, model_rirs)

# fixed effects
broom.mixed::tidy(model_rirs, effects = 'fixed')

# random effects
broom.mixed::tidy(model_rirs, effects = 'ran_pars')

# residuals
plot(model_rirs)
qqnorm(model_rirs)