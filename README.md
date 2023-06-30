# Wildfire Smoke Is Associated With an Increased Risk of Cardiorespiratory Emergency Department Visits in Alaska
### Created by: Grace Kuiper
### Version: 6/8/2023

## Background

This repository provides R code that was used for analysis of Alaska ambient air pollution and emergency department visit data in order to understand the impact of wildfire smoke (WFS) in Alaska on cardiorespiratory-related health outcomes. The finidngs of these analyses are available here: https://doi.org/10.1029/2020GH000349

## Methods

This was a case-crossover study to detect associations between WFS (as measured using satellite imagery of smoke plumes and regulatory ambient air pollution monitoring) and cardiorespiratory-related emergency department visits. The three largest population centers in Alaska were included in this study: Anchorage, Matanuska-Susitna Valley, and Fairbanks. Emergency department data were obtained from the Alaska Health Facilities Data Reporting Program for 2015-2019. Monitored PM2.5 data were obtained from the Alaska Department of Environmental Conservation (DEC), which maintains ground-based monitoring stations for  Environmental Protection Agency (EPA) regulatory monitoring. Smoke plume locations identified using satellite imagery were obtained from NOAA's Hazard Mapping System. Daily temperature and relative humidity for each study site were included as a covariates in all models; these data were obtained from the NOAA ASOS network. 

![plot](https://github.com/grelber13/Alaska_heat_2022/blob/main/AK_map.jpeg?raw=true)

For each emergency department visit that occurred within the study sites during April-September of 2015-2019 that had a primary diagnosis indicative of a cardiorespiratory event, patient age, sex, and race was available. The day on which the health event occurred was a case day; control days were sampled by selecting the same day of the week from every other week during the wildfire season (April-September) of the same year.

### Wilfire smoke exposure

To isolate WFS PM2.5, first each day was evaluated to determine if it was a wildfire day (WFD). To be considered a WFD, two criteria needed to be met: 1) the daily PM2.5 concentration exceeded one standard deviation above the long-term (2008-2019) monthly mean for daily PM2.5, and 2) the DEC monitor site had to be located within 50 km of a smoke plume. On days that were determined to be a WFD, WFS PM2.5 was estimated by subtracting the long-term monthly mean PM2.5 from the PM2.5 measured that day. Non-WFDs were assigned a WFS PM2.5 concentration of zero.

### Case-crossover analysis

Conditional logistic regression models were fit for same-day and lagged WFS PM2.5 to estimate the assocation between WFS PM2.5 and cardiorespiratory ED visits on days 0-5 after exposure. In all models, daily relative humidity and maximum temperature were included as covaraites. In the primary analysis, models were fit using the entire dataset. Additional stratified analyses were also conducted by race, sex, and age group.

## Data

Available in the **/raw_data/** folder are the following files:
1. **daily_88101_0819_NCORE_Garden_Parkgate_Butte.xlsx** - this Excel spreadsheet contains daily PM2.5 data from the Alaska Department of Environmental Conservation.

Data requests can be submitted to the Alaska Health Facilities Data Reporting Program [here](https://dhss.alaska.gov/dph/VitalStats/Pages/HFDR/default.aspx#:~:text=The%20Alaska%20Health%20Facilities%20Data,and%20the%20Alaska%20Outpatient%20Database.).

## Scripts

Available in the **/scripts/** folder are the following files:
1. **00a_DEC_data_cleaning.R** - cleans PM2.5 covariate data
2. **00b_pull_NOAA_data.R** - pulls and cleans relative humidity and temperature data from the NOAA ASOS network using the `riem` and `rnoaa` packages
3. **00c_finalize_exposure_data.R** - merges and cleans PM2.5 and weather data to create dataset with single daily observations for each study site
4. **00d_HFDR_data_cleaning.R** - cleans HFDR emergency department visits data and creates case-crossover dataset
6. **01_overall_casecross_analysis.R** - conducts the case-crossover single day and distributed lag model analyses for the entire dataset using conditional logistic regression
7. **02_stratified_casecross_analysis.R** - conducts stratified case-crossover single day and distributed lag model analyses using conditional logistic regression
8. **03a_clean_ozone_data.R** - clean DEC ozone data
9. **03b_ozone_sensitivity_analysis.R** - conduct sensitivity analysis with ozone included as covariate
