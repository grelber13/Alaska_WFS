#' -----------------------------------------------------------------------------
#' Project: Stratified analysis
#' Date created: July 1, 2020
#' Author: Grace Kuiper
#' Contact: grelber13@gmail.com
#' 
#' Last modified: August 11, 2020
#' 
#' 
#' Description:
#' 
#' conduct single-day and DL models using stratified data
#' 
#' -----------------------------------------------------------------------------
setwd('~/Alaska')
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

set.seed(1234)

#' Read in full dataset

scale10 <- function(x) (x/10)

casecross_list <- readRDS("overall_results/casecross_list.rds") 
casecross_list <- casecross_list %>%
  map(~mutate_at(.,vars(contains("WFS")),scale10)) %>%
  map(~mutate(.,age=ifelse(age<15,"<15",
                          ifelse(age<66,"15-65",
                                 ifelse(age>65,">65",NA)))))
casecross_list <- casecross_list %>%
  map(~mutate(.,race=ifelse(race==3,"Alaskan Native",
                            ifelse(race<9 & race!=3,"Non-Alaskan Native",NA))))

strat_var <- "age"
strat_obs <- ">65"

#' Filter only Anchorage data
temp_list <- casecross_list %>%
  map(~rename(.,strat_var=!!strat_var)) %>%
  map(~filter(.,strat_var==strat_obs))

if (!dir.exists(paste(getwd(),'stratified',strat_var,sep="/"))) {
  dir.create(paste(getwd(),'stratified',strat_var,sep="/")) }

if (!dir.exists(paste(getwd(),'stratified',strat_var,strat_obs,sep="/"))) {
  dir.create(paste(getwd(),'stratified',strat_var,strat_obs,sep="/")) }
## Same-Day Associations 

# Evaluating a 10 ug/m^3^ increase in wildfire smoke-related PM~2.5~ and risk for emergency department visit for certain cardiopulmonary reasons.
# same-day association with WFS PM2.5 and ED visits during wildfire season
full_list <- temp_list
temp_list <- list()
outcomes <- c()
for (i in 1:length(full_list)) {
  if (nrow(full_list[[i]])>0 & length(unique(full_list[[i]]$hfdrepisode))>1) {
    outcomes <- c(outcomes,unique(full_list[[i]]$out_name))
    temp_list[[length(temp_list)+1]] <- full_list[[i]]
  }
}
pm_1sd_results <- temp_list %>% 
  map_dfr(. , function(df){
    mod <- clogit(outcome ~ WFS_1sd + mean_tmpf + mean_relh + strata(hfdrepisode), data = df)
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
                   mutate_all(exp))
  }) %>% # end map
  cbind(outcomes, .)
# print results
print(pm_1sd_results)
saveRDS(pm_1sd_results,paste('stratified',strat_var,strat_obs,"pm_1sd_results.rds",sep="/"))

pm_2sd_results <- temp_list %>% 
  map_dfr(. , function(df){
    mod <- clogit(outcome ~ WFS_2sd + mean_tmpf + mean_relh + strata(hfdrepisode), data = df)
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
                   mutate_all(exp))
  }) %>% # end map
  cbind(outcomes, .)
# print results
print(pm_2sd_results)
saveRDS(pm_2sd_results,paste('stratified',strat_var,strat_obs,"pm_2sd_results.rds",sep="/"))

# library(extrafont)
# plot <- pm_1sd_lag1_results %>%
#   mutate(group=ifelse(outcomes=="cardiorespiratory","All",
#                       ifelse(grepl("asthma|copd|pneumonia|bronchitis|respiratory",outcomes),
#                              "Respiratory","Cardiovascular"))) %>%
#   filter(!grepl("heatillness|forearmfracture",outcomes)) %>%
#   mutate(outcomes=ifelse(outcomes=="respiratory","all respiratory",
#                          ifelse(outcomes=="cardiovascular","all cardiovascular",as.character(outcomes)))) %>%
#   arrange(desc(group),outcomes) %>%
#   ggplot(aes(x=reorder(outcomes, desc(group)), y = WFS.estimate, colour = group)) +
#   geom_point() +
#   geom_errorbar(aes(ymin=WFS.lower95, ymax=WFS.upper95), width = 0.3) +
#   scale_color_manual("Cardiorespiratory", values = c("grey", "darkblue","#ff00cc")) +
#   guides(color = guide_legend(reverse = TRUE)) +
#   geom_hline(yintercept = 1, linetype = "dashed", colour = "red") +
#   ylab(expression(paste("OR: Increase in Smoke PM"[2.5]," by 10 ", mu, "g/m"^3))) +
#   xlab("Primary diagnosis") +
#   theme_bw() +
#   theme(panel.grid.major.x = element_blank(),
#         panel.grid.major.y = element_line(linetype = "dotted"),
#         panel.grid.minor = element_blank(),
#         axis.text.x = element_text(angle = 45, hjust=0.95, vjust = 0.9),
#         text=element_text(family="ArialMT")) +
#   ggtitle(expression(paste("Odds Ratio for ED Visit Using 1-Day Lagged 1SD Smoke Day")))
# 
# print(plot)

## 1 day lag Associations 

# 1 day lag association with WFS PM2.5 and ED visits during wildfire season
pm_1sd_lag1_results <- temp_list %>% 
  map_dfr(. , function(df){
    mod <- clogit(outcome ~ WFS_1sd_lag1 + mean_tmpf_lag1 + mean_relh_lag1 + strata(hfdrepisode), data = df)
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
                   mutate_all(exp))
  }) %>% # end map
  cbind(outcomes, .)
# print results
print(pm_1sd_lag1_results)
saveRDS(pm_1sd_lag1_results,paste('stratified',strat_var,strat_obs,"pm_1sd_lag1_results.rds",sep="/"))

pm_2sd_lag1_results <- temp_list %>% 
  map_dfr(. , function(df){
    mod <- clogit(outcome ~ WFS_2sd_lag1 + mean_tmpf_lag1 + mean_relh_lag1 + strata(hfdrepisode), data = df)
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
                   mutate_all(exp))
  }) %>% # end map
  cbind(outcomes, .)
# print results
print(pm_2sd_lag1_results)
saveRDS(pm_2sd_lag1_results,paste('stratified',strat_var,strat_obs,"pm_2sd_lag1_results.rds",sep="/"))

## 2 day lag Associations 

# 2 day lag association with WFS PM2.5 and ED visits during wildfire season
pm_1sd_lag2_results <- temp_list %>% 
  map_dfr(. , function(df){
    mod <- clogit(outcome ~ WFS_1sd_lag2 + mean_tmpf_lag2 + mean_relh_lag2 + strata(hfdrepisode), data = df)
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
                   mutate_all(exp))
  }) %>% # end map
  cbind(outcomes, .)
# print results
print(pm_1sd_lag2_results)
saveRDS(pm_1sd_lag2_results,paste('stratified',strat_var,strat_obs,"pm_1sd_lag2_results.rds",sep="/"))

pm_2sd_lag2_results <- temp_list %>% 
  map_dfr(. , function(df){
    mod <- clogit(outcome ~ WFS_2sd_lag2 + mean_tmpf_lag2 + mean_relh_lag2 + strata(hfdrepisode), data = df)
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
                   mutate_all(exp))
  }) %>% # end map
  cbind(outcomes, .)
# print results
print(pm_2sd_lag2_results)
saveRDS(pm_2sd_lag2_results,paste('stratified',strat_var,strat_obs,"pm_2sd_lag2_results.rds",sep="/"))

## 3 day lag Associations 

# 3 day lag association with WFS PM2.5 and ED visits during wildfire season
pm_1sd_lag3_results <- temp_list %>% 
  map_dfr(. , function(df){
    mod <- clogit(outcome ~ WFS_1sd_lag3 + mean_tmpf_lag3 + mean_relh_lag3 + strata(hfdrepisode), data = df)
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
                   mutate_all(exp))
  }) %>% # end map
  cbind(outcomes, .)
# print results
print(pm_1sd_lag3_results)
saveRDS(pm_1sd_lag3_results,paste('stratified',strat_var,strat_obs,"pm_1sd_lag3_results.rds",sep="/"))

pm_2sd_lag3_results <- temp_list %>% 
  map_dfr(. , function(df){
    mod <- clogit(outcome ~ WFS_2sd_lag3 + mean_tmpf_lag3 + mean_relh_lag3 + strata(hfdrepisode), data = df)
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
                   mutate_all(exp))
  }) %>% # end map
  cbind(outcomes, .)
# print results
print(pm_2sd_lag3_results)
saveRDS(pm_2sd_lag3_results,paste('stratified',strat_var,strat_obs,"pm_2sd_lag3_results.rds",sep="/"))

## 4 day lag Associations 

# 4 day lag association with WFS PM2.5 and ED visits during wildfire season
pm_1sd_lag4_results <- temp_list %>% 
  map_dfr(. , function(df){
    mod <- clogit(outcome ~ WFS_1sd_lag4 + mean_tmpf_lag4 + mean_relh_lag4 + strata(hfdrepisode), data = df)
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
                   mutate_all(exp))
  }) %>% # end map
  cbind(outcomes, .)
# print results
print(pm_1sd_lag4_results)
saveRDS(pm_1sd_lag4_results,paste('stratified',strat_var,strat_obs,"pm_1sd_lag4_results.rds",sep="/"))

pm_2sd_lag4_results <- temp_list %>% 
  map_dfr(. , function(df){
    mod <- clogit(outcome ~ WFS_2sd_lag4 + mean_tmpf_lag4 + mean_relh_lag4 + strata(hfdrepisode), data = df)
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
                   mutate_all(exp))
  }) %>% # end map
  cbind(outcomes, .)
# print results
print(pm_2sd_lag4_results)
saveRDS(pm_2sd_lag4_results,paste('stratified',strat_var,strat_obs,"pm_2sd_lag4_results.rds",sep="/"))

## 5 day lag Associations 

# 5 day lag association with WFS PM2.5 and ED visits during wildfire season
pm_1sd_lag5_results <- temp_list %>% 
  map_dfr(. , function(df){
    mod <- clogit(outcome ~ WFS_1sd_lag5 + mean_tmpf_lag5 + mean_relh_lag5 + strata(hfdrepisode), data = df)
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
                   mutate_all(exp))
  }) %>% # end map
  cbind(outcomes, .)
# print results
print(pm_1sd_lag5_results)
saveRDS(pm_1sd_lag5_results,paste('stratified',strat_var,strat_obs,"pm_1sd_lag5_results.rds",sep="/"))

pm_2sd_lag5_results <- temp_list %>% 
  map_dfr(. , function(df){
    mod <- clogit(outcome ~ WFS_2sd_lag5 + mean_tmpf_lag5 + mean_relh_lag5 + strata(hfdrepisode), data = df)
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
                   mutate_all(exp))
  }) %>% # end map
  cbind(outcomes, .)
# print results
print(pm_2sd_lag5_results)
saveRDS(pm_2sd_lag5_results,paste('stratified',strat_var,strat_obs,"pm_2sd_lag5_results.rds",sep="/"))


#' table of health outcomes
ED_visits_table <- data.frame()
for (i in 1:length(temp_list)) {
  npatients <- nrow(temp_list[[i]][which(!is.na(temp_list[[i]]$enctype)),])
  ED_visits_table[i,"Diagnosis"] <- unique(temp_list[[i]]$out_name)
  ED_visits_table[i,"All (n)"] <- sum(temp_list[[i]]$outcome)
  ED_visits_table[i,"F (%)"] <- round(nrow(temp_list[[i]][which(temp_list[[i]]$sex=="F"),])/npatients*100,1)
  ED_visits_table[i,"M (%)"] <- round(nrow(temp_list[[i]][which(temp_list[[i]]$sex=="M"),])/npatients*100,1)
  ED_visits_table[i,"<15 years (%)"] <- round(nrow(temp_list[[i]][which(temp_list[[i]]$age<15),])/npatients*100,1)
  ED_visits_table[i,"15-65 years (%)"] <- round(nrow(temp_list[[i]][which(temp_list[[i]]$age>=15 &
                                                                                 temp_list[[i]]$age<66),])/npatients*100,1)
  ED_visits_table[i,">65 years (%)"] <- round(nrow(temp_list[[i]][which(temp_list[[i]]$age>65),])/npatients*100,1)
}
write.csv(ED_visits_table,'tables_figures/ED_visits_table.csv')
## Distributed Lag Associations
# 
# I believe Rish might use the dlnm package created by Antonio Gasparrini, and I think Rish would like to include him as a co-author, which would be great! I have some trouble using his package for models like the conditional logistic regression or mixed models, so I use steps Ander Wilson has shown me to estimate distributed lag for linear relationships. I think keeping the term as linear is fine in this case since it seems like most outcomes between PM~2.5~ and cardiopulmonary are linear in other studies that have specifically tried to address this question.
# 
# First step is to define a function (called distribut_that_lag) that is really long (sorry!) that will take a model with a distributed lag basis function in it and extract the lagged and cumulative results from it. The lag_mod term is where you reference the model, and strata is the strata you'd like to estimate it for. This is used to get out terms for PM~2.5~ smoke and non-smoke. Right now this function only handles lagged matrices of 0-6 days, and degrees of freedom, and will need to be modified if you want to vary this. I can try and make this in to a general function/package for use by Rish if need be.

distribute_that_lag <- function(lag_mod, strata, exposure_basis) {
  # output pm basis estimates
  parms <- broom::tidy(lag_mod) %>% 
    filter(stringr::str_detect(term, strata)) %>% 
    select(estimate) %>% 
    as_vector()
  # output estimate names for cov matrix
  names <- stringr::str_subset(names(lag_mod$coefficients), strata)
  # estimate associations
  est <- exposure_basis %*% parms
  # estimate standard error for each interval
  # time variable
  time <- ((rep(1:length(est))-1))
  # covariance matrix for knots 
  cov_mat <- as.matrix(vcov(lag_mod))[names, names]
  # estimate variance of spline
  var <- exp_b %*% cov_mat %*% t(exp_b)
  # estimate lag ----
  # estimate standard error for each lag day for smoke
  l_se <- sqrt(diag(var))
  # calculate lower and upper bound for smoke
  l_est_l95 <- est + (l_se*qnorm(1-0.975))
  l_est_u95 <- est + (l_se*qnorm(0.975))
  l_type <- "lag"
  # lag dataframe
  l_df <- data.frame(strata, l_type, time, 
                     exp(est), exp(l_est_l95), exp(l_est_u95), 
                     row.names = NULL) 
  # assign column names
  colnames(l_df) <- c("strata", "type", "time", 
                      "odds_ratio", "lower_95", "upper_95")
  # cumulative estimates
  c_est <- sapply(seq_along(est), function(x){
    sum(est[1:x])
  })
  # stderr cumulative effect smk
  c_se <- sapply(seq_along(c_est), function(y){
    sqrt(sum(var[1:y,1:y]))
  })
  # estimate 95% CI
  c_l95 <- c_est+(c_se*qnorm(1-0.975))
  c_u95 <- c_est+(c_se*qnorm(0.975))
  # type
  c_type <- "cumulative"
  # return dataframe
  c_df <- data.frame(strata, c_type, time, exp(c_est), 
                     exp(c_l95), exp(c_u95), row.names = NULL) 
  # assign column names
  colnames(c_df) <- c("strata", "type", "time", 
                      "odds_ratio", "lower_95", "upper_95")
  # bind lagged and cumulative 
  lag_est <- rbind(l_df, c_df) %>% 
    mutate(strata = as.character(strata),
           type = as.character(type))
  # return lagged estimate
  return(lag_est)
} # end lag estimate function


# distributed lag function
lag_est <- data.frame()
#temp_list <- temp_list[-12]
for (i in 1:length(temp_list)) {
  # output dataframe from list
  data <- temp_list[[i]] %>% 
    ungroup() %>%
    mutate(outcome = as.numeric(as.character(outcome))) %>%
    # remove missing lagged data
    filter(!is.na(WFS_1sd_lag5))
  # output outcome name
  out_name <- as.character(unique(data$out_name))
  # print(out_name) # track which outcome dataframe it's on
  # create lagged matrix
  pm_mat <- as.matrix(select(data, WFS_1sd, contains("1sd_lag")))
  temp_mat <- as.matrix(select(data, contains("mean_tmpf")))
  relh_mat <- as.matrix(select(data, contains("mean_relh")))
  # define lagged basis spline
  exp_b <- ns(0:(ncol(pm_mat)-1), df = 3, intercept = T)
  # pm basis
  pm_basis <- pm_mat %*% exp_b
  # temp basis
  temp_basis <- temp_mat %*% exp_b
  relh_basis <- relh_mat %*% exp_b
  # run lagged model
  lag_mod <- clogit(outcome ~ pm_basis + temp_basis + relh_basis + strata(hfdrepisode), data = data)
  
  # estimate lag estimate
  lag_est_temp <- distribute_that_lag(lag_mod, strata = "pm", 
                                      exposure_basis = exp_b) %>% 
    mutate(outcome = out_name) %>% select(outcome, strata:upper_95)
  lag_est <- rbind(lag_est,lag_est_temp)
}
mort_dl_pm_results <- lag_est

lagged_results <- mort_dl_pm_results %>% filter(type == "lag")
cumulative_results <- mort_dl_pm_results %>% filter(type == "cumulative")

saveRDS(lagged_results,paste('stratified',strat_var,strat_obs,'lagged_results.rds',sep="/"))
saveRDS(cumulative_results,paste('stratified',strat_var,strat_obs,'cumulative_results.rds',sep="/"))
