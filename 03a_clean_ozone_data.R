#' -----------------------------------------------------------------------------
#' Project: Pull DEC ozone data
#' Date created:July 1, 2020
#' Author: Grace Kuiper
#' Contact: grelber13@gmail.com
#' 
#' Last modified: July 1, 2020
#' 
#' 
#' Description:
#' 
#' Pull and clean ozone data for Fairbanks to conduct sensitivity analysis.
#' 
#' -----------------------------------------------------------------------------
setwd('~/Alaska/ozone_data')
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

set.seed(1234)

#' read in annual data and filter to only include study sites

O3_2015 <- read.csv("hourly_44201_2015.csv") %>% 
  filter(State.Code==2 & County.Code==90) %>%
  group_by(`Date.Local`,`State.Code`,`County.Code`,`Site.Num`,`Parameter.Code`,POC,
           Latitude,Longitude,Datum,`Parameter.Name`,`Units.of.Measure`,Qualifier,
           MDL,Method.Type,Method.Code,Method.Name,`County.Name`,`State.Name`,
           Date.of.Last.Change) %>%
  select(-Date.GMT,-Time.Local,-Time.GMT) %>%
  add_tally() %>%
  rename(Observation.Count=n) %>%
  mutate(Observation.Percent=(Observation.Count/24)*100) %>%
  summarise_all(mean,na.rm=T) %>%
  filter(Observation.Percent>=75)
  
O3_2016 <- read.csv("hourly_44201_2016.csv") %>%
  filter(State.Code==2 & County.Code==90) %>%
  group_by(`Date.Local`,`State.Code`,`County.Code`,`Site.Num`,`Parameter.Code`,POC,
           Latitude,Longitude,Datum,`Parameter.Name`,`Units.of.Measure`,Qualifier,
           MDL,Method.Type,Method.Code,Method.Name,`County.Name`,`State.Name`,
           Date.of.Last.Change) %>%
  select(-Date.GMT,-Time.Local,-Time.GMT) %>%
  add_tally() %>%
  rename(Observation.Count=n) %>%
  mutate(Observation.Percent=(Observation.Count/24)*100) %>%
  summarise_all(mean,na.rm=T) %>%
  filter(Observation.Percent>=75)
O3_2017 <- read.csv("hourly_44201_2017.csv") %>% 
  filter(State.Code==2 & County.Code==90) %>%
  group_by(`Date.Local`,`State.Code`,`County.Code`,`Site.Num`,`Parameter.Code`,POC,
           Latitude,Longitude,Datum,`Parameter.Name`,`Units.of.Measure`,Qualifier,
           MDL,Method.Type,Method.Code,Method.Name,`County.Name`,`State.Name`,
           Date.of.Last.Change) %>%
  select(-Date.GMT,-Time.Local,-Time.GMT) %>%
  add_tally() %>%
  rename(Observation.Count=n) %>%
  mutate(Observation.Percent=(Observation.Count/24)*100) %>%
  summarise_all(mean,na.rm=T) %>%
  filter(Observation.Percent>=75)
O3_2018 <- read.csv("hourly_44201_2018.csv") %>%
  filter(State.Code==2 & County.Code==90) %>%
  group_by(`Date.Local`,`State.Code`,`County.Code`,`Site.Num`,`Parameter.Code`,POC,
           Latitude,Longitude,Datum,`Parameter.Name`,`Units.of.Measure`,Qualifier,
           MDL,Method.Type,Method.Code,Method.Name,`County.Name`,`State.Name`,
           Date.of.Last.Change) %>%
  select(-Date.GMT,-Time.Local,-Time.GMT) %>%
  add_tally() %>%
  rename(Observation.Count=n) %>%
  mutate(Observation.Percent=(Observation.Count/24)*100) %>%
  summarise_all(mean,na.rm=T) %>%
  filter(Observation.Percent>=75)
O3_2019 <- read.csv("hourly_44201_2019.csv") %>% 
  filter(State.Code==2 & County.Code==90) %>%
  group_by(`Date.Local`,`State.Code`,`County.Code`,`Site.Num`,`Parameter.Code`,POC,
           Latitude,Longitude,Datum,`Parameter.Name`,`Units.of.Measure`,Qualifier,
           MDL,Method.Type,Method.Code,Method.Name,`County.Name`,`State.Name`,
           Date.of.Last.Change) %>%
  select(-Date.GMT,-Time.Local,-Time.GMT) %>%
  add_tally() %>%
  rename(Observation.Count=n) %>%
  mutate(Observation.Percent=(Observation.Count/24)*100) %>%
  summarise_all(mean,na.rm=T) %>%
  filter(Observation.Percent>=75)

O3_data <- rbind(O3_2015,O3_2016,O3_2017,O3_2018,O3_2019) %>%
  rename(ozone=Sample.Measurement)
#### Impute for PM2.5 that is missing for one or two days ####
every.day <- seq(min(as.Date(O3_data$Date.Local)), max(as.Date(O3_data$Date.Local)), by="1 day")
every.day <- data.frame(Date.Local=every.day)
O3_data <- O3_data %>%
  mutate(Date.Local=as.Date(Date.Local)) %>%
  arrange(Date.Local) %>%
  full_join(every.day,by=c("Date.Local"))

for (i in 1:nrow(O3_data)) {
  if (i==1) {
    O3_data[i,"group"]=0
  } else {
    if (is.na(O3_data[i,"ozone"])) {
      if (is.na(O3_data[i-1,"ozone"])) {
        O3_data[i,"group"]=O3_data[i-1,"group"]
      } else if (!is.na(O3_data[i-1,"ozone"])) {
        O3_data[i,"group"]=O3_data[i-1,"group"]+1
      }
    } else if (!is.na(O3_data[i,"ozone"])) {
      O3_data[i,"group"]=O3_data[i-1,"group"]+1
    }
  }
}
O3_data <- O3_data %>%
  group_by(group) %>%
  add_tally() %>%
  ungroup() %>%
  mutate(group=ifelse(is.na(ozone),group,NA),
         n=ifelse(is.na(ozone),n,NA))

O3_data <- O3_data %>%
  mutate(ozone=ifelse(!is.na(n) & n==1,(lag(ozone)+lead(ozone))/2,ozone))
O3_data <- O3_data %>%
  mutate(ozone=ifelse(!is.na(n) & n==2 & !is.na(lag(ozone)),(lag(ozone)+lead(ozone,2))/2,
                                ifelse(!is.na(n) & n==2 & is.na(lag(ozone)),(lag(ozone,2)+lead(ozone))/2,
                                       ozone))) #59 days missing during study period
O3_data <- O3_data %>%
  select(County.Name,ozone,State.Code,Date.Local,County.Code,Site.Num) %>%
  mutate(County.Name="Fairbanks North Star",
         State.Code=2,
         County.Code=90,
         Site.Num=34,
         ozone_lag1=lag(ozone),
         ozone_lag2=lag(ozone,2),
         ozone_lag3=lag(ozone,3),
         ozone_lag4=lag(ozone,4),
         ozone_lag5=lag(ozone,5)) %>%
  rename(county=County.Name,
         STATE.CODE=State.Code,
         COUNTY.CODE=County.Code,
         Site_No=Site.Num,
         Date=Date.Local)
saveRDS(O3_data,"O3_data.rds")
