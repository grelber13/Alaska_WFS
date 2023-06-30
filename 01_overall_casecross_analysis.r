# ---
#   title: "Distributed Lag Model"
# author: "Grace Kuiper"
# date: '2018-03-22'
# output:
#   html_document:
#   df_print: paged
# editor_options:
#   chunk_output_type: inline

## Purpose: For the purpose of AK WFS analysis using wildfire-related PM2.5: 
##			identify threshold values for cardiovascular and respiratory 
##          outcomes using conditional logistic regression models for a case-crossover 
##          study design.

## Contributor: Ryan Gan for code to run distributed lag models

## Email: grace.kuiper@colostate.edu
## ---------------------------
  
library(tidyverse) # data wrangle/plot
library(survival) # conditional logistic models
library(splines) # splines
library(lubridate) # works with dates
library(broom)
library(stringr)
library(dlnm)

# Read in case-crossover dataset from `00d_HFDR_data_cleaning.R` 
scale10 <- function(x) (x/10)

casecross_list <- readRDS("health_data/casecross_list.rds") 
casecross_list <- casecross_list %>%
  map(~mutate_at(.,vars(contains("WFS")),scale10))

#### Same-Day Associations ####

# Evaluating a 10 ug/m^3^ increase in wildfire smoke-related PM~2.5~ and risk for emergency department visit for certain cardiopulmonary reasons.
# Same-day association with WFS PM2.5 and ED visits during wildfire season
outcomes <- c()
for (i in 1:length(casecross_list)) {
  outcomes <- c(outcomes,unique(casecross_list[[i]]$out_name))
}
pm_1sd_results <- casecross_list %>% 
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
saveRDS(pm_1sd_results,"pm_1sd_results.rds")

pm_2sd_results <- casecross_list %>% 
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
saveRDS(pm_2sd_results,"pm_2sd_results.rds")
														
## 1 day lag Associations 

# 1 day lag association with WFS PM2.5 and ED visits during wildfire season
pm_1sd_lag1_results <- casecross_list %>% 
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
saveRDS(pm_1sd_lag1_results,"pm_1sd_lag1_results.rds")

pm_2sd_lag1_results <- casecross_list %>% 
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
saveRDS(pm_2sd_lag1_results,"pm_2sd_lag1_results.rds")

## 2 day lag Associations 

# 2 day lag association with WFS PM2.5 and ED visits during wildfire season
pm_1sd_lag2_results <- casecross_list %>% 
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
saveRDS(pm_1sd_lag2_results,"pm_1sd_lag2_results.rds")

pm_2sd_lag2_results <- casecross_list %>% 
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
saveRDS(pm_2sd_lag2_results,"pm_2sd_lag2_results.rds")

## 3 day lag Associations 

# 3 day lag association with WFS PM2.5 and ED visits during wildfire season
pm_1sd_lag3_results <- casecross_list %>% 
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
saveRDS(pm_1sd_lag3_results,"pm_1sd_lag3_results.rds")

pm_2sd_lag3_results <- casecross_list %>% 
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
saveRDS(pm_2sd_lag3_results,"pm_2sd_lag3_results.rds")

## 4 day lag Associations 

# 4 day lag association with WFS PM2.5 and ED visits during wildfire season
pm_1sd_lag4_results <- casecross_list %>% 
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
saveRDS(pm_1sd_lag4_results,"pm_1sd_lag4_results.rds")

pm_2sd_lag4_results <- casecross_list %>% 
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
saveRDS(pm_2sd_lag4_results,"pm_2sd_lag4_results.rds")

## 5 day lag Associations 

# 5 day lag association with WFS PM2.5 and ED visits during wildfire season
pm_1sd_lag5_results <- casecross_list %>% 
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
saveRDS(pm_1sd_lag5_results,"pm_1sd_lag5_results.rds")

pm_2sd_lag5_results <- casecross_list %>% 
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
saveRDS(pm_2sd_lag5_results,"pm_2sd_lag5_results.rds")


## Distributed Lag Associations
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
for (i in 1:length(casecross_list)) {
  # output dataframe from list
  data <- casecross_list[[i]] %>% 
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

saveRDS(lagged_results,'lagged_results.rds')
saveRDS(cumulative_results,'cumulative_results.rds')

# plot results
plot <- lagged_results %>%
  mutate(group=ifelse(outcome=="cardiorespiratory","All",
                      ifelse(grepl("asthma|copd|pneumonia|bronchitis|respiratory",outcome),
                             "Respiratory","Cardiovascular"))) %>%
  filter(!grepl("heatillness|forearmfracture",outcome)) %>%
  mutate(outcome=ifelse(outcome=="respiratory","all respiratory",
                         ifelse(outcome=="cardiovascular","all cardiovascular",as.character(outcome)))) %>%
  mutate(outcome_f = factor(outcome, levels=c('all respiratory','asthma','bronchitis',
                                              'copd','pneumonia','all cardiovascular',
                                              'arrhythmia','cerebrovascular',
                                              'heartfailure','ischemic','myocardial',
                                              'cardiorespiratory'))) %>%
  ggplot(aes(x=time, y=odds_ratio,color=group)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95,fill=group), 
               alpha = 0.3) + 
  scale_x_continuous(breaks = c(seq(0,7, by=1))) +
  geom_hline(yintercept = 1, linetype = 2, colour = "red") +
  scale_color_manual("Cardiorespiratory", values = c("grey", "darkblue","#ff00cc")) +
  scale_fill_manual("Cardiorespiratory", values = c("grey", "darkblue","#ff00cc")) +
  guides(fill = guide_legend(reverse = TRUE),color=FALSE) +
  # adding facet wrap to estimate for each outcome
  facet_wrap(~outcome_f, scales = "free_y") +
  ylab(expression("Odds Ratio: 10 ug/m^3 increase PM2.5")) +
  xlab("Lagged Days") +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(linetype = "dotted"),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=9),
        axis.text.x=element_text(size=7),
        axis.text.y=element_text(size=7),
        axis.ticks = element_line(color="black"),
        panel.background = element_rect(fill = "white", colour = "black"))
print(plot)
