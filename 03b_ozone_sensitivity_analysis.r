#' -----------------------------------------------------------------------------
#' Project: Ozone sensitivity analysis
#' Date created: July 1, 2020
#' Author: Grace Kuiper
#' Contact: grelber13@gmail.com
#' 
#' Last modified: July 21, 2020
#' 
#' 
#' Description:
#' 
#' Merge ozone and complete dataset and run models to conduct sensitivity for
#' inclusion of ozone as covariate.
#' 
#' -----------------------------------------------------------------------------

library(foreign)
library(tidyverse)
library(data.table)
library(lubridate)
library(writexl)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(lme4)
library(survival) # conditional logistic models
library(splines) # splines
library(lubridate) # works with dates
library(broom)
library(stringr)
library(dlnm)

set.seed(1234)

#' Read in ozone data

O3_data <- readRDS('ozone_data/O3_data.rds')

#' Read in full dataset

scale10 <- function(x) (x/10)

casecross_list <- readRDS("overall_results/casecross_list.rds") 
casecross_list <- casecross_list %>%
  map(~mutate_at(.,vars(contains("WFS")),scale10))

#' Filter only Fairbanks data
Fairbanks_list <- casecross_list %>%
  map(~filter(.,county=="Fairbanks North Star"))

#' Merge complete dataset with ozone data
Fairbanks_list_O3 <- Fairbanks_list %>%
  map(~left_join(.,O3_data %>%
                   select(-STATE.CODE,-COUNTY.CODE,-Site_No) %>%
                   mutate(Date=as.Date(Date)),
                 by=c("county","Date")))

##Models for dataset with ozone

outcomes <- c()
for (i in 1:length(Fairbanks_list_O3)) {
  outcomes <- c(outcomes,unique(Fairbanks_list[[i]]$out_name))
}
pm_1sd_O3_results <- Fairbanks_list_O3 %>% 
  map_dfr(. , function(df){
    mod <- clogit(outcome ~ WFS_1sd + mean_tmpf + mean_relh + ozone + strata(hfdrepisode), data = df)
    est <- cbind(WF_est <- broom::tidy(mod) %>% filter(term == "WFS_1sd") %>%
                   select(estimate, conf.low, conf.high) %>% 
                   rename(lower95 = conf.low, upper95 = conf.high) %>% 
                   setNames(paste0('WFS.', names(.))) %>%
                   mutate_all(exp),
                 temp_est <- broom::tidy(mod) %>% filter(term == "mean_tmpf") %>% 
                   select(estimate, conf.low, conf.high) %>% 
                   rename(lower95 = conf.low, upper95 = conf.high) %>% 
                   setNames(paste0('temp.', names(.))) %>%
                   mutate_all(exp),
                 relh_est <- broom::tidy(mod) %>% filter(term == "mean_relh") %>% 
                   select(estimate, conf.low, conf.high) %>% 
                   rename(lower95 = conf.low, upper95 = conf.high) %>% 
                   setNames(paste0('relh.', names(.))) %>%
                   mutate_all(exp),
                 ozone <- broom::tidy(mod) %>% filter(term == "ozone") %>%
                   select(estimate, conf.low, conf.high) %>%
                   rename(lower95 = conf.low, upper95 = conf.high) %>%
                   setNames(paste0('O3.', names(.))) %>%
                   mutate_all(exp))
  }) %>% # end map
  cbind(outcomes, .)
# print results
print(pm_1sd_O3_results)
saveRDS(pm_1sd_O3_results,"ozone_data/Fairbanks_pm_1sd_O3_results.rds")

pm_2sd_O3_results <- Fairbanks_list_O3 %>% 
  map_dfr(. , function(df){
    mod <- clogit(outcome ~ WFS_2sd + mean_tmpf + mean_relh + ozone + strata(hfdrepisode), data = df)
    est <- cbind(WF_est <- broom::tidy(mod) %>% filter(term == "WFS_2sd") %>%
                   select(estimate, conf.low, conf.high) %>% 
                   rename(lower95 = conf.low, upper95 = conf.high) %>% 
                   setNames(paste0('WFS.', names(.))) %>%
                   mutate_all(exp),
                 temp_est <- broom::tidy(mod) %>% filter(term == "mean_tmpf") %>% 
                   select(estimate, conf.low, conf.high) %>% 
                   rename(lower95 = conf.low, upper95 = conf.high) %>% 
                   setNames(paste0('temp.', names(.))) %>%
                   mutate_all(exp),
                 relh_est <- broom::tidy(mod) %>% filter(term == "mean_relh") %>% 
                   select(estimate, conf.low, conf.high) %>% 
                   rename(lower95 = conf.low, upper95 = conf.high) %>% 
                   setNames(paste0('relh.', names(.))) %>%
                   mutate_all(exp),
                 ozone <- broom::tidy(mod) %>% filter(term == "ozone") %>%
                   select(estimate, conf.low, conf.high) %>%
                   rename(lower95 = conf.low, upper95 = conf.high) %>%
                   setNames(paste0('O3.', names(.))) %>%
                   mutate_all(exp))
  }) %>% # end map
  cbind(outcomes, .)
# print results
print(pm_2sd_O3_results)
saveRDS(pm_2sd_O3_results,"ozone_data/Fairbanks_pm_2sd_O3_results.rds")

library(extrafont)
plot <- pm_1sd_O3_results %>%
  mutate(group=ifelse(outcomes=="cardiorespiratory","All",
                      ifelse(grepl("asthma|copd|pneumonia|bronchitis|respiratory",outcomes),
                             "Respiratory","Cardiovascular"))) %>%
  filter(!grepl("heatillness|forearmfracture",outcomes)) %>%
  mutate(outcomes=ifelse(outcomes=="respiratory","all respiratory",
                         ifelse(outcomes=="cardiovascular","all cardiovascular",as.character(outcomes)))) %>%
  arrange(desc(group),outcomes) %>%
  ggplot(aes(x=reorder(outcomes, desc(group)), y = WFS.estimate, colour = group)) +
  geom_point() +
  geom_errorbar(aes(ymin=WFS.lower95, ymax=WFS.upper95), width = 0.3) +
  scale_color_manual("Cardiorespiratory", values = c("grey", "darkblue","#ff00cc")) +
  guides(color = guide_legend(reverse = TRUE)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "red") +
  ylab(expression(paste("OR: Increase in Smoke PM"[2.5]," by 10 ", mu, "g/m"^3))) +
  xlab("Primary diagnosis") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(linetype = "dotted"),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=0.95, vjust = 0.9),
        text=element_text(family="ArialMT")) +
  ggtitle(expression(paste("Odds Ratio for ED Visit Using 1-Day Lagged 1SD Smoke Day")))

print(plot)

## 1 day lag Associations 

# 1 day lag association with WFS PM2.5 and ED visits during wildfire season
pm_1sd_lag1_O3_results <- Fairbanks_list_O3 %>% 
  map_dfr(. , function(df){
    mod <- clogit(outcome ~ WFS_1sd_lag1 + mean_tmpf_lag1 + mean_relh_lag1 + ozone_lag1 + strata(hfdrepisode), data = df)
    est <- cbind(WF_est <- broom::tidy(mod) %>% filter(term == "WFS_1sd_lag1") %>%
                   select(estimate, conf.low, conf.high) %>% 
                   rename(lower95 = conf.low, upper95 = conf.high) %>% 
                   setNames(paste0('WFS.', names(.))) %>%
                   mutate_all(exp),
                 temp_est <- broom::tidy(mod) %>% filter(term == "mean_tmpf_lag1") %>% 
                   select(estimate, conf.low, conf.high) %>% 
                   rename(lower95 = conf.low, upper95 = conf.high) %>% 
                   setNames(paste0('temp.', names(.))) %>%
                   mutate_all(exp),
                 relh_lag1_est <- broom::tidy(mod) %>% filter(term == "mean_relh_lag1") %>% 
                   select(estimate, conf.low, conf.high) %>% 
                   rename(lower95 = conf.low, upper95 = conf.high) %>% 
                   setNames(paste0('relh.', names(.))) %>%
                   mutate_all(exp),
                 ozone_lag1 <- broom::tidy(mod) %>% filter(term == "ozone_lag1") %>%
                   select(estimate, conf.low, conf.high) %>%
                   rename(lower95 = conf.low, upper95 = conf.high) %>%
                   setNames(paste0('O3.', names(.))) %>%
                   mutate_all(exp))
  }) %>% # end map
  cbind(outcomes, .)
# print results
print(pm_1sd_lag1_O3_results)
saveRDS(pm_1sd_lag1_O3_results,"ozone_data/Fairbanks_pm_1sd_lag1_O3_results.rds")

pm_2sd_lag1_O3_results <- Fairbanks_list_O3 %>% 
  map_dfr(. , function(df){
    mod <- clogit(outcome ~ WFS_2sd_lag1 + mean_tmpf_lag1 + mean_relh_lag1 + ozone_lag1 + strata(hfdrepisode), data = df)
    est <- cbind(WF_est <- broom::tidy(mod) %>% filter(term == "WFS_2sd_lag1") %>%
                   select(estimate, conf.low, conf.high) %>% 
                   rename(lower95 = conf.low, upper95 = conf.high) %>% 
                   setNames(paste0('WFS.', names(.))) %>%
                   mutate_all(exp),
                 temp_est <- broom::tidy(mod) %>% filter(term == "mean_tmpf_lag1") %>% 
                   select(estimate, conf.low, conf.high) %>% 
                   rename(lower95 = conf.low, upper95 = conf.high) %>% 
                   setNames(paste0('temp.', names(.))) %>%
                   mutate_all(exp),
                 relh_lag1_est <- broom::tidy(mod) %>% filter(term == "mean_relh_lag1") %>% 
                   select(estimate, conf.low, conf.high) %>% 
                   rename(lower95 = conf.low, upper95 = conf.high) %>% 
                   setNames(paste0('relh.', names(.))) %>%
                   mutate_all(exp),
                 ozone_lag1 <- broom::tidy(mod) %>% filter(term == "ozone_lag1") %>%
                   select(estimate, conf.low, conf.high) %>%
                   rename(lower95 = conf.low, upper95 = conf.high) %>%
                   setNames(paste0('O3.', names(.))) %>%
                   mutate_all(exp))
  }) %>% # end map
  cbind(outcomes, .)
# print results
print(pm_2sd_lag1_O3_results)
saveRDS(pm_2sd_lag1_O3_results,"ozone_data/Fairbanks_pm_2sd_lag1_O3_results.rds")

## 2 day lag Associations 

# 2 day lag association with WFS PM2.5 and ED visits during wildfire season
pm_1sd_lag2_O3_results <- Fairbanks_list_O3 %>% 
  map_dfr(. , function(df){
    mod <- clogit(outcome ~ WFS_1sd_lag2 + mean_tmpf_lag2 + mean_relh_lag2 + ozone_lag2 + strata(hfdrepisode), data = df)
    est <- cbind(WF_est <- broom::tidy(mod) %>% filter(term == "WFS_1sd_lag2") %>%
                   select(estimate, conf.low, conf.high) %>% 
                   rename(lower95 = conf.low, upper95 = conf.high) %>% 
                   setNames(paste0('WFS.', names(.))) %>%
                   mutate_all(exp),
                 temp_est <- broom::tidy(mod) %>% filter(term == "mean_tmpf_lag2") %>% 
                   select(estimate, conf.low, conf.high) %>% 
                   rename(lower95 = conf.low, upper95 = conf.high) %>% 
                   setNames(paste0('temp.', names(.))) %>%
                   mutate_all(exp),
                 relh_est <- broom::tidy(mod) %>% filter(term == "mean_relh_lag2") %>% 
                   select(estimate, conf.low, conf.high) %>% 
                   rename(lower95 = conf.low, upper95 = conf.high) %>% 
                   setNames(paste0('relh.', names(.))) %>%
                   mutate_all(exp),
                 ozone_lag2 <- broom::tidy(mod) %>% filter(term == "ozone_lag2") %>%
                   select(estimate, conf.low, conf.high) %>%
                   rename(lower95 = conf.low, upper95 = conf.high) %>%
                   setNames(paste0('O3.', names(.))) %>%
                   mutate_all(exp))
  }) %>% # end map
  cbind(outcomes, .)
# print results
print(pm_1sd_lag2_O3_results)
saveRDS(pm_1sd_lag2_O3_results,"ozone_data/Fairbanks_pm_1sd_lag2_O3_results.rds")

pm_2sd_lag2_O3_results <- Fairbanks_list_O3 %>% 
  map_dfr(. , function(df){
    mod <- clogit(outcome ~ WFS_2sd_lag2 + mean_tmpf_lag2 + mean_relh_lag2 + ozone_lag2 + strata(hfdrepisode), data = df)
    est <- cbind(WF_est <- broom::tidy(mod) %>% filter(term == "WFS_2sd_lag2") %>%
                   select(estimate, conf.low, conf.high) %>% 
                   rename(lower95 = conf.low, upper95 = conf.high) %>% 
                   setNames(paste0('WFS.', names(.))) %>%
                   mutate_all(exp),
                 temp_est <- broom::tidy(mod) %>% filter(term == "mean_tmpf_lag2") %>% 
                   select(estimate, conf.low, conf.high) %>% 
                   rename(lower95 = conf.low, upper95 = conf.high) %>% 
                   setNames(paste0('temp.', names(.))) %>%
                   mutate_all(exp),
                 relh_est <- broom::tidy(mod) %>% filter(term == "mean_relh_lag2") %>% 
                   select(estimate, conf.low, conf.high) %>% 
                   rename(lower95 = conf.low, upper95 = conf.high) %>% 
                   setNames(paste0('relh.', names(.))) %>%
                   mutate_all(exp),
                 ozone_lag2 <- broom::tidy(mod) %>% filter(term == "ozone_lag2") %>%
                   select(estimate, conf.low, conf.high) %>%
                   rename(lower95 = conf.low, upper95 = conf.high) %>%
                   setNames(paste0('O3.',names(.))) %>%
                   mutate_all(exp))
  }) %>% # end map
  cbind(outcomes, .)
# print results
print(pm_2sd_lag2_O3_results)
saveRDS(pm_2sd_lag2_O3_results,"ozone_data/Fairbanks_pm_2sd_lag2_O3_results.rds")

## 3 day lag Associations 

]# 3 day lag association with WFS PM2.5 and ED visits during wildfire season
pm_1sd_lag3_O3_results <- Fairbanks_list_O3 %>% 
  map_dfr(. , function(df){
    mod <- clogit(outcome ~ WFS_1sd_lag3 + mean_tmpf_lag3 + mean_relh_lag3 + ozone_lag3 + strata(hfdrepisode), data = df)
    est <- cbind(WF_est <- broom::tidy(mod) %>% filter(term == "WFS_1sd_lag3") %>%
                   select(estimate, conf.low, conf.high) %>% 
                   rename(lower95 = conf.low, upper95 = conf.high) %>% 
                   setNames(paste0('WFS.', names(.))) %>%
                   mutate_all(exp),
                 temp_est <- broom::tidy(mod) %>% filter(term == "mean_tmpf_lag3") %>% 
                   select(estimate, conf.low, conf.high) %>% 
                   rename(lower95 = conf.low, upper95 = conf.high) %>% 
                   setNames(paste0('temp.', names(.))) %>%
                   mutate_all(exp),
                 relh_est <- broom::tidy(mod) %>% filter(term == "mean_relh_lag3") %>% 
                   select(estimate, conf.low, conf.high) %>% 
                   rename(lower95 = conf.low, upper95 = conf.high) %>% 
                   setNames(paste0('relh.', names(.))) %>%
                   mutate_all(exp),
                 ozone_lag3 <- broom::tidy(mod) %>% filter(term == "ozone_lag3") %>%
                   select(estimate, conf.low, conf.high) %>%
                   rename(lower95 = conf.low, upper95 = conf.high) %>%
                   setNames(paste0('O3.', names(.))) %>%
                   mutate_all(exp))
  }) %>% # end map
  cbind(outcomes, .)
# print results
print(pm_1sd_lag3_O3_results)
saveRDS(pm_1sd_lag3_O3_results,"ozone_data/Fairbanks_pm_1sd_lag3_O3_results.rds")

pm_2sd_lag3_O3_results <- Fairbanks_list_O3 %>% 
  map_dfr(. , function(df){
    mod <- clogit(outcome ~ WFS_2sd_lag3 + mean_tmpf_lag3 + mean_relh_lag3 + ozone_lag3 + strata(hfdrepisode), data = df)
    est <- cbind(WF_est <- broom::tidy(mod) %>% filter(term == "WFS_2sd_lag3") %>%
                   select(estimate, conf.low, conf.high) %>% 
                   rename(lower95 = conf.low, upper95 = conf.high) %>% 
                   setNames(paste0('WFS.', names(.))) %>%
                   mutate_all(exp),
                 temp_est <- broom::tidy(mod) %>% filter(term == "mean_tmpf_lag3") %>% 
                   select(estimate, conf.low, conf.high) %>% 
                   rename(lower95 = conf.low, upper95 = conf.high) %>% 
                   setNames(paste0('temp.', names(.))) %>%
                   mutate_all(exp),
                 relh_est <- broom::tidy(mod) %>% filter(term == "mean_relh_lag3") %>% 
                   select(estimate, conf.low, conf.high) %>% 
                   rename(lower95 = conf.low, upper95 = conf.high) %>% 
                   setNames(paste0('relh.', names(.))) %>%
                   mutate_all(exp),
                 ozone_lag3 <- broom::tidy(mod) %>% filter(term == "ozone_lag3") %>%
                   select(estimate, conf.low, conf.high) %>%
                   rename(lower95 = conf.low, upper95 = conf.high) %>%
                   setNames(paste0('O3.', names(.))) %>%
                   mutate_all(exp))
  }) %>% # end map
  cbind(outcomes, .)
# print results
print(pm_2sd_lag3_O3_results)
saveRDS(pm_2sd_lag3_O3_results,"ozone_data/Fairbanks_pm_2sd_lag3_O3_results.rds")

## 4 day lag Associations 

# 4 day lag association with WFS PM2.5 and ED visits during wildfire season
pm_1sd_lag4_O3_results <- Fairbanks_list_O3 %>% 
  map_dfr(. , function(df){
    mod <- clogit(outcome ~ WFS_1sd_lag4 + mean_tmpf_lag4 + mean_relh_lag4 + ozone_lag4 + strata(hfdrepisode), data = df)
    est <- cbind(WF_est <- broom::tidy(mod) %>% filter(term == "WFS_1sd_lag4") %>%
                   select(estimate, conf.low, conf.high) %>% 
                   rename(lower95 = conf.low, upper95 = conf.high) %>% 
                   setNames(paste0('WFS.', names(.))) %>%
                   mutate_all(exp),
                 temp_est <- broom::tidy(mod) %>% filter(term == "mean_tmpf_lag4") %>% 
                   select(estimate, conf.low, conf.high) %>% 
                   rename(lower95 = conf.low, upper95 = conf.high) %>% 
                   setNames(paste0('temp.', names(.))) %>%
                   mutate_all(exp),
                 relh_est <- broom::tidy(mod) %>% filter(term == "mean_relh_lag4") %>% 
                   select(estimate, conf.low, conf.high) %>% 
                   rename(lower95 = conf.low, upper95 = conf.high) %>% 
                   setNames(paste0('relh.', names(.))) %>%
                   mutate_all(exp),
                 ozone_lag4 <- broom::tidy(mod) %>% filter(term == "ozone_lag4") %>%
                   select(estimate, conf.low, conf.high) %>%
                   rename(lower95 = conf.low, upper95 = conf.high) %>%
                   setNames(paste0('O3.', names(.))) %>%
                   mutate_all(exp))
  }) %>% # end map
  cbind(outcomes, .)
# print results
print(pm_1sd_lag4_O3_results)
saveRDS(pm_1sd_lag4_O3_results,"ozone_data/Fairbanks_pm_1sd_lag4_O3_results.rds")

pm_2sd_lag4_O3_results <- Fairbanks_list_O3 %>% 
  map_dfr(. , function(df){
    mod <- clogit(outcome ~ WFS_2sd_lag4 + mean_tmpf_lag4 + mean_relh_lag4 + ozone_lag4 + strata(hfdrepisode), data = df)
    est <- cbind(WF_est <- broom::tidy(mod) %>% filter(term == "WFS_2sd_lag4") %>%
                   select(estimate, conf.low, conf.high) %>% 
                   rename(lower95 = conf.low, upper95 = conf.high) %>% 
                   setNames(paste0('WFS.', names(.))) %>%
                   mutate_all(exp),
                 temp_est <- broom::tidy(mod) %>% filter(term == "mean_tmpf_lag4") %>% 
                   select(estimate, conf.low, conf.high) %>% 
                   rename(lower95 = conf.low, upper95 = conf.high) %>% 
                   setNames(paste0('temp.', names(.))) %>%
                   mutate_all(exp),
                 relh_est <- broom::tidy(mod) %>% filter(term == "mean_relh_lag4") %>% 
                   select(estimate, conf.low, conf.high) %>% 
                   rename(lower95 = conf.low, upper95 = conf.high) %>% 
                   setNames(paste0('relh.', names(.))) %>%
                   mutate_all(exp),
                 ozone_lag4 <- broom::tidy(mod) %>% filter(term == "ozone_lag4") %>%
                   select(estimate, conf.low, conf.high) %>%
                   rename(lower95 = conf.low, upper95 = conf.high) %>%
                   setNames(paste0('O3.',names(.))) %>%
                   mutate_all(exp))
  }) %>% # end map
  cbind(outcomes, .)
# print results
print(pm_2sd_lag4_O3_results)
saveRDS(pm_2sd_lag4_O3_results,"ozone_data/Fairbanks_pm_2sd_lag4_O3_results.rds")

## 5 day lag Associations 

# 5 day lag association with WFS PM2.5 and ED visits during wildfire season
pm_1sd_lag5_O3_results <- Fairbanks_list_O3 %>% 
  map_dfr(. , function(df){
    mod <- clogit(outcome ~ WFS_1sd_lag5 + mean_tmpf_lag5 + mean_relh_lag5 + ozone_lag5 + strata(hfdrepisode), data = df)
    est <- cbind(WF_est <- broom::tidy(mod) %>% filter(term == "WFS_1sd_lag5") %>%
                   select(estimate, conf.low, conf.high) %>% 
                   rename(lower95 = conf.low, upper95 = conf.high) %>% 
                   setNames(paste0('WFS.', names(.))) %>%
                   mutate_all(exp),
                 temp_est <- broom::tidy(mod) %>% filter(term == "mean_tmpf_lag5") %>% 
                   select(estimate, conf.low, conf.high) %>% 
                   rename(lower95 = conf.low, upper95 = conf.high) %>% 
                   setNames(paste0('temp.', names(.))) %>%
                   mutate_all(exp),
                 relh_est <- broom::tidy(mod) %>% filter(term == "mean_relh_lag5") %>% 
                   select(estimate, conf.low, conf.high) %>% 
                   rename(lower95 = conf.low, upper95 = conf.high) %>% 
                   setNames(paste0('relh.', names(.))) %>%
                   mutate_all(exp),
                 ozone_lag5 <- broom::tidy(mod) %>% filter(term == "ozone_lag5") %>%
                   select(estimate, conf.low, conf.high) %>%
                   rename(lower95 = conf.low, upper95 = conf.high) %>%
                   setNames(paste0('O3.', names(.))) %>%
                   mutate_all(exp))
  }) %>% # end map
  cbind(outcomes, .)
# print results
print(pm_1sd_lag5_O3_results)
saveRDS(pm_1sd_lag5_O3_results,"ozone_data/Fairbanks_pm_1sd_lag5_O3_results.rds")

pm_2sd_lag5_O3_results <- Fairbanks_list_O3 %>% 
  map_dfr(. , function(df){
    mod <- clogit(outcome ~ WFS_2sd_lag5 + mean_tmpf_lag5 + mean_relh_lag5 + ozone_lag5 + strata(hfdrepisode), data = df)
    est <- cbind(WF_est <- broom::tidy(mod) %>% filter(term == "WFS_2sd_lag5") %>%
                   select(estimate, conf.low, conf.high) %>% 
                   rename(lower95 = conf.low, upper95 = conf.high) %>% 
                   setNames(paste0('WFS.', names(.))) %>%
                   mutate_all(exp),
                 temp_est <- broom::tidy(mod) %>% filter(term == "mean_tmpf_lag5") %>% 
                   select(estimate, conf.low, conf.high) %>% 
                   rename(lower95 = conf.low, upper95 = conf.high) %>% 
                   setNames(paste0('temp.', names(.))) %>%
                   mutate_all(exp),
                 relh_est <- broom::tidy(mod) %>% filter(term == "mean_relh_lag5") %>% 
                   select(estimate, conf.low, conf.high) %>% 
                   rename(lower95 = conf.low, upper95 = conf.high) %>% 
                   setNames(paste0('relh.', names(.))) %>%
                   mutate_all(exp),
                 ozone_lag5 <- broom::tidy(mod) %>% filter(term == "ozone_lag5") %>%
                   select(estimate, conf.low, conf.high) %>%
                   rename(lower95 = conf.low, upper95 = conf.high) %>%
                   setNames(paste0('O3.', names(.))) %>%
                   mutate_all(exp))
  }) %>% # end map
  cbind(outcomes, .)
# print results
print(pm_2sd_lag5_O3_results)
saveRDS(pm_2sd_lag5_O3_results,"ozone_data/Fairbanks_pm_2sd_lag5_O3_results.rds")

#' merge all results
pm_1sd_O3_results <- rename_at(pm_1sd_O3_results,vars(-outcomes,-O3.estimate,-O3.lower95,-O3.upper95),function(x) paste0(x,"_ozone")) %>%
  mutate(lag=0,sd=1)
pm_2sd_O3_results <- rename_at(pm_2sd_O3_results,vars(-outcomes,-O3.estimate,-O3.lower95,-O3.upper95),function(x) paste0(x,"_ozone")) %>%
  mutate(lag=0,sd=2)
pm_1sd_lag1_O3_results <- rename_at(pm_1sd_lag1_O3_results,vars(-outcomes,-O3.estimate,-O3.lower95,-O3.upper95),function(x) paste0(x,"_ozone")) %>%
  mutate(lag=1,sd=1)
pm_2sd_lag1_O3_results <- rename_at(pm_2sd_lag1_O3_results,vars(-outcomes,-O3.estimate,-O3.lower95,-O3.upper95),function(x) paste0(x,"_ozone")) %>%
  mutate(lag=1,sd=2)
pm_1sd_lag2_O3_results <- rename_at(pm_1sd_lag2_O3_results,vars(-outcomes,-O3.estimate,-O3.lower95,-O3.upper95),function(x) paste0(x,"_ozone")) %>%
  mutate(lag=2,sd=1)
pm_2sd_lag2_O3_results <- rename_at(pm_2sd_lag2_O3_results,vars(-outcomes,-O3.estimate,-O3.lower95,-O3.upper95),function(x) paste0(x,"_ozone")) %>%
  mutate(lag=2,sd=2)
pm_1sd_lag3_O3_results <- rename_at(pm_1sd_lag3_O3_results,vars(-outcomes,-O3.estimate,-O3.lower95,-O3.upper95),function(x) paste0(x,"_ozone")) %>%
  mutate(lag=3,sd=1)
pm_2sd_lag3_O3_results <- rename_at(pm_2sd_lag3_O3_results,vars(-outcomes,-O3.estimate,-O3.lower95,-O3.upper95),function(x) paste0(x,"_ozone")) %>%
  mutate(lag=3,sd=2)
pm_1sd_lag4_O3_results <- rename_at(pm_1sd_lag4_O3_results,vars(-outcomes,-O3.estimate,-O3.lower95,-O3.upper95),function(x) paste0(x,"_ozone")) %>%
  mutate(lag=4,sd=1)
pm_2sd_lag4_O3_results <- rename_at(pm_2sd_lag4_O3_results,vars(-outcomes,-O3.estimate,-O3.lower95,-O3.upper95),function(x) paste0(x,"_ozone")) %>%
  mutate(lag=4,sd=2)
pm_1sd_lag5_O3_results <- rename_at(pm_1sd_lag5_O3_results,vars(-outcomes,-O3.estimate,-O3.lower95,-O3.upper95),function(x) paste0(x,"_ozone")) %>%
  mutate(lag=5,sd=1)
pm_2sd_lag5_O3_results <- rename_at(pm_2sd_lag5_O3_results,vars(-outcomes,-O3.estimate,-O3.lower95,-O3.upper95),function(x) paste0(x,"_ozone")) %>%
  mutate(lag=5,sd=2)


pm_1sd_results <- readRDS('stratified/county/Fairbanks North Star/pm_1sd_results.rds')
pm_2sd_results <- readRDS('stratified/county/Fairbanks North Star/pm_2sd_results.rds')
pm_1sd_lag1_results <- readRDS('stratified/county/Fairbanks North Star/pm_1sd_lag1_results.rds')
pm_2sd_lag1_results <- readRDS('stratified/county/Fairbanks North Star/pm_2sd_lag1_results.rds')
pm_1sd_lag2_results <- readRDS('stratified/county/Fairbanks North Star/pm_1sd_lag2_results.rds')
pm_2sd_lag2_results <- readRDS('stratified/county/Fairbanks North Star/pm_2sd_lag2_results.rds')
pm_1sd_lag3_results <- readRDS('stratified/county/Fairbanks North Star/pm_1sd_lag3_results.rds')
pm_2sd_lag3_results <- readRDS('stratified/county/Fairbanks North Star/pm_2sd_lag4_results.rds')
pm_1sd_lag4_results <- readRDS('stratified/county/Fairbanks North Star/pm_1sd_lag4_results.rds')
pm_2sd_lag4_results <- readRDS('stratified/county/Fairbanks North Star/pm_2sd_lag3_results.rds')
pm_1sd_lag5_results <- readRDS('stratified/county/Fairbanks North Star/pm_1sd_lag3_results.rds')
pm_2sd_lag5_results <- readRDS('stratified/county/Fairbanks North Star/pm_2sd_lag3_results.rds')

pm_1sd_results <- mutate(pm_1sd_results,lag=0,sd=1)
pm_2sd_results <- mutate(pm_2sd_results,lag=0,sd=2)
pm_1sd_lag1_results <- mutate(pm_1sd_lag1_results,lag=1,sd=1)
pm_2sd_lag1_results <- mutate(pm_2sd_lag1_results,lag=1,sd=2)
pm_1sd_lag2_results <- mutate(pm_1sd_lag2_results,lag=2,sd=1)
pm_2sd_lag2_results <- mutate(pm_2sd_lag2_results,lag=2,sd=2)
pm_1sd_lag3_results <- mutate(pm_1sd_lag3_results,lag=3,sd=1)
pm_2sd_lag3_results <- mutate(pm_2sd_lag4_results,lag=3,sd=2)
pm_1sd_lag4_results <- mutate(pm_1sd_lag4_results,lag=4,sd=1)
pm_2sd_lag4_results <- mutate(pm_2sd_lag3_results,lag=4,sd=2)
pm_1sd_lag5_results <- mutate(pm_1sd_lag3_results,lag=5,sd=1)
pm_2sd_lag5_results <- mutate(pm_2sd_lag3_results,lag=5,sd=2)


pm_1sd_comb_results <- full_join(pm_1sd_results,pm_1sd_O3_results,by=c("outcomes","lag","sd"))
pm_2sd_comb_results <- full_join(pm_2sd_results,pm_2sd_O3_results,by=c("outcomes","lag","sd"))
pm_1sd_lag1_comb_results <- full_join(pm_1sd_lag1_results,pm_1sd_lag1_O3_results,by=c("outcomes","lag","sd"))
pm_2sd_lag1_comb_results <- full_join(pm_2sd_lag1_results,pm_2sd_lag1_O3_results,by=c("outcomes","lag","sd"))
pm_1sd_lag2_comb_results <- full_join(pm_1sd_lag2_results,pm_1sd_lag2_O3_results,by=c("outcomes","lag","sd"))
pm_2sd_lag2_comb_results <- full_join(pm_2sd_lag2_results,pm_2sd_lag2_O3_results,by=c("outcomes","lag","sd"))
pm_1sd_lag3_comb_results <- full_join(pm_1sd_lag3_results,pm_1sd_lag3_O3_results,by=c("outcomes","lag","sd"))
pm_2sd_lag3_comb_results <- full_join(pm_2sd_lag3_results,pm_2sd_lag3_O3_results,by=c("outcomes","lag","sd"))
pm_1sd_lag4_comb_results <- full_join(pm_1sd_lag4_results,pm_1sd_lag4_O3_results,by=c("outcomes","lag","sd"))
pm_2sd_lag4_comb_results <- full_join(pm_2sd_lag4_results,pm_2sd_lag4_O3_results,by=c("outcomes","lag","sd"))
pm_1sd_lag5_comb_results <- full_join(pm_1sd_lag5_results,pm_1sd_lag5_O3_results,by=c("outcomes","lag","sd"))
pm_2sd_lag5_comb_results <- full_join(pm_2sd_lag5_results,pm_2sd_lag5_O3_results,by=c("outcomes","lag","sd"))



pm_comb_results <- rbind(pm_1sd_comb_results,pm_2sd_comb_results,
                         pm_1sd_lag1_comb_results,pm_2sd_lag1_comb_results,
                         pm_1sd_lag2_comb_results,pm_2sd_lag2_comb_results,
                         pm_1sd_lag3_comb_results,pm_2sd_lag3_comb_results,
                         pm_1sd_lag4_comb_results,pm_2sd_lag4_comb_results,
                         pm_1sd_lag5_comb_results,pm_2sd_lag5_comb_results)

#' Identify confounding
pm_comb_results <- pm_comb_results %>%
  mutate(perc_change=(abs(WFS.estimate-WFS.estimate_ozone)/WFS.estimate)*100)
conf_results <- pm_comb_results %>%
  filter(perc_change>=10)
