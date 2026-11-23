# ---------------------------
## load in packages that are required
library(activity)
library(DHARMa)
library(emmeans)
library(glmmTMB)
library(ggpubr)
library(viridis)
library(tidyverse)

#set random seed
set.seed(180)



#Load and format data ####
#read csv of bat data
bat_metrics <-  read_csv("output/bat_data_fil_class_join_red_night_met_weather_today_daily.csv")
bat_metrics <- distinct(bat_metrics, bat_name, date, .keep_all = TRUE)

bat_tags <-  read_csv("output/bat_tags.csv")
bat_tags <- distinct(bat_tags, bat_name, .keep_all = TRUE)


#check for NA values
colSums(is.na(bat_metrics))

#drop na values from moon columns
bat_metrics <- drop_na(bat_metrics, meanMoonlightIntensity)


# define factors
bat_metrics <- mutate(arrange(bat_metrics, time),
                      time = as.factor(time),
                      study_year = as.factor(study_year),
                      season = as.factor(season),
                      month = as.factor(month),
                      week = week,
                      week_integer = as.factor(week_integer), 
                      two_week = two_week,
                      date = as.Date(date),
                      bat_name = as.factor(bat_name),
                      tag_id = as.factor(tag_id), 
                      location = as.factor(location), 
                      survey_type = as.factor(survey_type), 
                      ring_no = as.factor(ring_no),
                      forearm = forearm, 
                      sex = as.factor(sex),
                      repro = as.factor(repro),
                      age = as.factor(age),
                      weight = weight, 
                      bat_emerged = bat_emerged, 
                      emergence_time = emergence_time, 
                      emergence_time_minutes = emergence_time_minutes, 
                      reentry_time = reentry_time, 
                      reentry_time_minutes = reentry_time_minutes, 
                      total_activity_minutes = total_activity_minutes, 
                      total_activity_proportion = total_activity_proportion, 
                      absolute_activity_minutes = absolute_activity_minutes, 
                      absolute_activity_proportion = absolute_activity_proportion, 
                      sunsetStart = as.POSIXct(sunsetStart),
                      sunriseEnd = as.POSIXct(sunriseEnd),
                      night_length = night_length, 
                      night_length_rounded = night_length_rounded, 
                      av_rel_humidity = av_rel_humidity, 
                      av_wind_speed = av_wind_speed, 
                      av_cloud_okta = av_cloud_okta, 
                      av_air_pressure = av_air_pressure, 
                      av_air_temp = av_air_temp, 
                      daily_rain = daily_rain, 
                      meanMoonlightIntensity = meanMoonlightIntensity,
                      meanMoonPhase = meanMoonPhase)


#make column of tagging date
bat_tags <- select(mutate(bat_tags,
                          tagging_date = as.Date(time)),
                  bat_name, tagging_date)

#join to main dataframe then filter to remove first night of data
bat_metrics <- left_join(bat_metrics, bat_tags, by = "bat_name")
bat_metrics <- filter(bat_metrics, date > tagging_date) #42 rows removed


#remove bats that had little data
bat_metrics <- filter(bat_metrics, !bat_name %in% c('Gerty'))

#remove bats tech error
bat_metrics <- filter(bat_metrics, !bat_name %in% c('Conchita', 'Jennet'))

#remove spring records
bat_metrics <- filter(bat_metrics, survey_type != "hibernation")

# remove weeks with only data from Vincent
bat_metrics <- filter(bat_metrics, week_integer != 17 & week_integer != 18 & week_integer != 19)

#keep only bats that emerged
bat_metrics <- filter(bat_metrics, bat_emerged == 1)                                                                          

bat_metrics <- droplevels(bat_metrics)



#Prepare variables for modelling ####
#standardise variables
bat_metrics$av_rel_humidity_std <- bat_metrics$av_rel_humidity / max(bat_metrics$av_rel_humidity)
bat_metrics$av_wind_speed_std <- bat_metrics$av_wind_speed / max(bat_metrics$av_wind_speed)
bat_metrics$av_cloud_okta_std <- bat_metrics$av_cloud_okta / max(bat_metrics$av_cloud_okta)
bat_metrics$av_air_pressure_std <- (bat_metrics$av_air_pressure - mean(bat_metrics$av_air_pressure)) / sd(bat_metrics$av_air_pressure)
bat_metrics$av_air_temp_std <- (bat_metrics$av_air_temp - mean(bat_metrics$av_air_temp)) / sd(bat_metrics$av_air_temp)
bat_metrics$daily_rain_std <- bat_metrics$daily_rain / max(bat_metrics$daily_rain)
bat_metrics$forearm_std <- (bat_metrics$forearm - mean(bat_metrics$forearm)) / sd(bat_metrics$forearm)
bat_metrics$weight_std <- (bat_metrics$weight - mean(bat_metrics$weight)) / sd(bat_metrics$weight)



# fitting model ####
model_1 <- glmmTMB(absolute_activity_proportion ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + weight_std + sex + age + study_year + week_integer + sex*age + weight_std*av_air_temp_std + av_wind_speed_std*av_air_temp_std + (1 | ring_no),
                   data = bat_metrics,
                   weights = night_length_rounded,
                   family = binomial(link = "logit"))

summary(model_1)
fixef(model_1)$cond
check_collinearity(model_1)

#DHARMa simulation
model_1_simulationOutput <- simulateResiduals(fittedModel = model_1, plot = F, n = 1000)
plotQQunif(model_1_simulationOutput, testUniformity = TRUE, 
           testOutliers = TRUE, testDispersion = TRUE)
plotResiduals(model_1_simulationOutput, quantreg = TRUE, smoothScatter = FALSE) # deviations in 0.25 quantile
testDispersion(model_1_simulationOutput)
testOutliers(model_1_simulationOutput)
testZeroInflation(model_1_simulationOutput)
mod_1_res = recalculateResiduals(model_1_simulationOutput, group = bat_metrics$time)
testTemporalAutocorrelation(mod_1_res, time = unique(bat_metrics$time))


#create inactive column
bat_metrics <- mutate(bat_metrics,
                      inactive = night_length_rounded - absolute_activity_minutes,
                      inactive = as.numeric(inactive),
                      active = as.numeric(absolute_activity_minutes))

#refit with betabinomial
model_1 <- glmmTMB(cbind(active, inactive) ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + weight_std + sex + age + study_year + week_integer + sex*age + weight_std*av_air_temp_std + av_wind_speed_std*av_air_temp_std + (1 | ring_no),
                   data = bat_metrics,
                   family = betabinomial(link = "logit"))

summary(model_1)
fixef(model_1)$cond
check_collinearity(model_1)

#DHARMa simulation
model_1_simulationOutput <- simulateResiduals(fittedModel = model_1, plot = F, n = 1000)
plotQQunif(model_1_simulationOutput, testUniformity = TRUE, 
           testOutliers = TRUE, testDispersion = TRUE)
plotResiduals(model_1_simulationOutput, quantreg = TRUE, smoothScatter = FALSE) # deviations in 0.25 quantile
testDispersion(model_1_simulationOutput)
testOutliers(model_1_simulationOutput)
testZeroInflation(model_1_simulationOutput)
mod_1_res = recalculateResiduals(model_1_simulationOutput, group = bat_metrics$time)
testTemporalAutocorrelation(mod_1_res, time = unique(bat_metrics$time))

plotResiduals(model_1_simulationOutput, form = bat_metrics$av_wind_speed_std, xlab = "Rank-transformed Standardised Wind Speed")
plotResiduals(model_1_simulationOutput, form = bat_metrics$av_air_temp_std, xlab = "Rank-transformed Standardised Air Temperature")
plotResiduals(model_1_simulationOutput, form = bat_metrics$meanMoonlightIntensity, xlab = "Rank-transformed Standardised Moonlight Intensity") #gotten worse
plotResiduals(model_1_simulationOutput, form = bat_metrics$av_cloud_okta_std, xlab = "Rank-transformed Standardised Cloud Cover") 
plotResiduals(model_1_simulationOutput, form = bat_metrics$daily_rain_std, xlab = "Rank-transformed Standardised Rainfall") # slight deviation in 0.25 quartile
plotResiduals(model_1_simulationOutput, form = bat_metrics$av_rel_humidity_std, xlab = "Rank-transformed Standardised Humidity")
plotResiduals(model_1_simulationOutput, form = bat_metrics$weight_std, xlab = "Rank-transformed Standardised Weight")
plotResiduals(model_1_simulationOutput, form = bat_metrics$age, xlab = "Rank-transformed Age")
plotResiduals(model_1_simulationOutput, form = bat_metrics$sex, xlab = "Rank-transformed Sex")
plotResiduals(model_1_simulationOutput, form = bat_metrics$study_year, xlab = "Rank-transformed Study Year")
plotResiduals(model_1_simulationOutput, form = bat_metrics$week_integer, xlab = "Rank-transformed Week")


#remove outliers (3)
bat_metrics <- filter(bat_metrics, absolute_activity_proportion < 0.5)
bat_metrics <- mutate(arrange(bat_metrics, time), 
                      tim_fac = numFactor(time))
bat_metrics <- droplevels(bat_metrics)


#refit betabinomial without outliers
model_1 <- glmmTMB(cbind(active, inactive) ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + weight_std + sex + age + study_year + week_integer + sex*age + weight_std*av_air_temp_std + av_wind_speed_std*av_air_temp_std + (1 | ring_no),
                   data = bat_metrics,
                   family = betabinomial(link = "logit"))

summary(model_1)
fixef(model_1)$cond
check_collinearity(model_1)

#DHARMa simulation
model_1_simulationOutput <- simulateResiduals(fittedModel = model_1, plot = F, n = 1000)
plotQQunif(model_1_simulationOutput, testUniformity = TRUE, 
           testOutliers = TRUE, testDispersion = TRUE)
plotResiduals(model_1_simulationOutput, quantreg = TRUE, smoothScatter = FALSE) # deviations in 0.25 quantile
testDispersion(model_1_simulationOutput)
testOutliers(model_1_simulationOutput)
testZeroInflation(model_1_simulationOutput)
mod_1_res = recalculateResiduals(model_1_simulationOutput, group = bat_metrics$time)
testTemporalAutocorrelation(mod_1_res, time = unique(bat_metrics$time))

plotResiduals(model_1_simulationOutput, form = bat_metrics$av_wind_speed_std, xlab = "Rank-transformed Standardised Wind Speed")
plotResiduals(model_1_simulationOutput, form = bat_metrics$av_air_temp_std, xlab = "Rank-transformed Standardised Air Temperature")
plotResiduals(model_1_simulationOutput, form = bat_metrics$meanMoonlightIntensity, xlab = "Rank-transformed Standardised Moonlight Intensity") #gotten worse
plotResiduals(model_1_simulationOutput, form = bat_metrics$av_cloud_okta_std, xlab = "Rank-transformed Standardised Cloud Cover") 
plotResiduals(model_1_simulationOutput, form = bat_metrics$daily_rain_std, xlab = "Rank-transformed Standardised Rainfall") # slight deviation in 0.25 quartile
plotResiduals(model_1_simulationOutput, form = bat_metrics$av_rel_humidity_std, xlab = "Rank-transformed Standardised Humidity")
plotResiduals(model_1_simulationOutput, form = bat_metrics$weight_std, xlab = "Rank-transformed Standardised Weight")
plotResiduals(model_1_simulationOutput, form = bat_metrics$age, xlab = "Rank-transformed Age")
plotResiduals(model_1_simulationOutput, form = bat_metrics$sex, xlab = "Rank-transformed Sex")
plotResiduals(model_1_simulationOutput, form = bat_metrics$study_year, xlab = "Rank-transformed Study Year")
plotResiduals(model_1_simulationOutput, form = bat_metrics$week_integer, xlab = "Rank-transformed Week")


#change dispersion
model_1_disp_week <- update(model_1, dispformula = ~ week_integer)

anova(model_1, model_1_disp_week)
#                   Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# model_1           33 6149.0 6298.2 -3041.5   6083.0
# model_1_disp_week 48 6035.3 6252.3 -2969.6   5939.3 143.73     15  < 2.2e-16 ***


model_1_disp_year <- update(model_1_disp_week, dispformula = ~ week_integer + study_year)

anova(model_1_disp_week, model_1_disp_year)
#                   Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# model_1_disp_week 48 6035.3 6252.3 -2969.6   5939.3
# model_1_disp_year 49 6034.6 6256.1 -2968.3   5936.6 2.693      1     0.1008


model_1_disp_sex <- update(model_1_disp_week, dispformula = ~ week_integer + sex)

anova(model_1_disp_week, model_1_disp_sex)
#                   Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# model_1_disp_week 48 6035.3 6252.3 -2969.6   5939.3
# model_1_disp_sex  49 6036.0 6257.5 -2969.0   5938.0 1.2522      1     0.2631


model_1_disp_age <- update(model_1_disp_week, dispformula = ~ week_integer + age)

anova(model_1_disp_week, model_1_disp_age)
#                   Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# model_1_disp_week 48 6035.3 6252.3 -2969.6   5939.3
# model_1_disp_age  50 6019.5 6245.5 -2959.7   5919.5 19.808      2  4.997e-055 ***


model_1_disp_moon <- update(model_1_disp_age, dispformula = ~ week_integer + age + meanMoonlightIntensity)

anova(model_1_disp_age, model_1_disp_moon)
#                   Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# model_1_disp_age  50 6019.5 6245.5 -2959.7   5919.5
# model_1_disp_moon 51 6018.1 6248.6 -2958.0   5916.1 3.3881      1    0.06567 .


model_1_disp_cloud <- update(model_1_disp_age, dispformula = ~ week_integer + age + av_cloud_okta_std)

anova(model_1_disp_age, model_1_disp_cloud)
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# model_1_disp_age   50 6019.5 6245.5 -2959.7   5919.5
# model_1_disp_cloud 51 6021.4 6251.9 -2959.7   5919.4 0.0801      1     0.7771


model_1_disp_temp <- update(model_1_disp_age, dispformula = ~ week_integer + age + av_air_temp_std)

anova(model_1_disp_age, model_1_disp_temp)
#                   Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)   
# model_1_disp_age  50 6019.5 6245.5 -2959.7   5919.5
# model_1_disp_temp 51 6015.6 6246.1 -2956.8   5913.6 5.8776      1    0.01533 *


model_1_disp_wind <- update(model_1_disp_temp, dispformula = ~ week_integer + age + av_air_temp_std + av_wind_speed_std)

anova(model_1_disp_temp, model_1_disp_wind)
#                   Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)  
# model_1_disp_temp 51 6015.6 6246.1 -2956.8   5913.6
# model_1_disp_wind 52 6015.5 6250.6 -2955.8   5911.5 2.0478      1     0.1524


model_1_disp_rain <- update(model_1_disp_temp, dispformula = ~ week_integer + age + av_air_temp_std + daily_rain_std)

anova(model_1_disp_temp, model_1_disp_rain)
#                   Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)   
# model_1_disp_temp 51 6015.6 6246.1 -2956.8   5913.6
# model_1_disp_rain 52 6011.6 6246.6 -2953.8   5907.6 6.0406      1    0.01398 *


model_1_disp_humid <- update(model_1_disp_rain, dispformula = ~ week_integer + age + av_air_temp_std + daily_rain_std + av_rel_humidity_std)

anova(model_1_disp_rain, model_1_disp_humid)
#                    Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# model_1_disp_rain  52 6011.6 6246.6 -2953.8   5907.6
# model_1_disp_humid 53 6013.5 6253.1 -2953.8   5907.5 0.0257      1     0.8726


model_1_disp_weight <- update(model_1_disp_rain, dispformula = ~ week_integer + age + av_air_temp_std + daily_rain_std + weight_std)

anova(model_1_disp_rain, model_1_disp_weight)
#                     Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# model_1_disp_rain   52 6011.6 6246.6 -2953.8   5907.6
# model_1_disp_weight 53 6010.2 6249.8 -2952.1   5904.2 3.3423      1    0.06752 .


# final dispersion model
summary(model_1_disp_rain)
#  Family: betabinomial  ( logit )
# Formula:
# cbind(active, inactive) ~ av_air_temp_std + av_wind_speed_std +
#     meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std +
#     av_rel_humidity_std + weight_std + sex + age + study_year +
#     week_integer + sex * age + weight_std * av_air_temp_std +
#     av_wind_speed_std * av_air_temp_std + (1 | ring_no)
# Dispersion:
# ~week_integer + age + av_air_temp_std + daily_rain_std
# Data: bat_metrics

#       AIC       BIC    logLik -2*log(L)  df.resid
#    6011.6    6246.6   -2953.8    5907.6       627

# Random effects:

# Conditional model:
#  Groups  Name        Variance Std.Dev.
#  ring_no (Intercept) 0.1063   0.3261
# Number of obs: 679, groups:  ring_no, 36

# Conditional model:
#                                    Estimate Std. Error z value Pr(>|z|)
# (Intercept)                       -0.798503   0.865346  -0.923 0.356135
# av_air_temp_std                    0.259358   0.081541   3.181 0.001469 **
# av_wind_speed_std                 -0.335628   0.343734  -0.976 0.328857
# meanMoonlightIntensity            -0.157983   0.393108  -0.402 0.687771
# av_cloud_okta_std                  0.004890   0.147742   0.033 0.973599    
# daily_rain_std                    -0.102878   0.180827  -0.569 0.569405
# av_rel_humidity_std               -1.280471   0.810065  -1.581 0.113946
# weight_std                        -0.284013   0.104818  -2.710 0.006737 **
# sexM                               0.066990   0.205737   0.326 0.744720
# ageJ                              -0.038284   0.306031  -0.125 0.900446
# ageYA                             -0.474140   0.277697  -1.707 0.087748 .
# study_yearSecond Winter           -0.086068   0.142281  -0.605 0.545236
# week_integer2                      0.229481   0.263053   0.872 0.383004
# week_integer3                     -0.110406   0.267773  -0.412 0.680112
# week_integer4                     -0.509920   0.349020  -1.461 0.144014
# week_integer5                     -0.460337   0.351772  -1.309 0.190662    
# week_integer6                     -0.005152   0.333357  -0.015 0.987670
# week_integer7                     -0.565175   0.297917  -1.897 0.057816 .
# week_integer8                     -0.509302   0.284370  -1.791 0.073296 .
# week_integer9                     -0.953866   0.273016  -3.494 0.000476 ***
# week_integer10                    -0.937640   0.278622  -3.365 0.000765 ***
# week_integer11                    -1.275753   0.303107  -4.209 2.57e-05 ***
# week_integer12                    -1.908018   0.300713  -6.345 2.22e-10 ***
# week_integer13                    -2.031016   0.296368  -6.853 7.23e-12 ***
# week_integer14                    -2.239789   0.314318  -7.126 1.03e-12 ***
# week_integer15                    -2.655941   0.352542  -7.534 4.93e-14 ***
# week_integer16                    -1.901444   0.394251  -4.823 1.41e-06 ***
# sexM:ageJ                         -0.535686   0.373092  -1.436 0.151059
# sexM:ageYA                        -0.401797   0.417612  -0.962 0.335985
# av_air_temp_std:weight_std         0.065240   0.032041   2.036 0.041738 *
# av_air_temp_std:av_wind_speed_std -0.290890   0.199619  -1.457 0.145053
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Dispersion model:
#                 Estimate Std. Error z value Pr(>|z|)    
# (Intercept)      2.48900    0.40664   6.121 9.31e-10 ***
# week_integer2   -0.14694    0.45915  -0.320 0.748947
# week_integer3    0.73363    0.48345   1.517 0.129143
# week_integer4    0.43192    0.58557   0.738 0.460753
# week_integer5    1.02424    0.60723   1.687 0.091652 .
# week_integer6   -0.18051    0.49853  -0.362 0.717290
# week_integer7    0.68534    0.49993   1.371 0.170416    
# week_integer8    0.38153    0.44513   0.857 0.391372
# week_integer9    1.28600    0.44539   2.887 0.003885 **
# week_integer10   0.85327    0.43475   1.963 0.049688 *
# week_integer11   1.40111    0.52337   2.677 0.007427 **
# week_integer12   1.68258    0.47652   3.531 0.000414 ***
# week_integer13   1.69139    0.47302   3.576 0.000349 ***
# week_integer14   1.49149    0.46715   3.193 0.001409 **
# week_integer15   2.32761    0.56109   4.148 3.35e-05 ***
# week_integer16   1.17910    0.58763   2.007 0.044799 *
# ageJ             0.49441    0.15762   3.137 0.001708 **
# ageYA            0.81429    0.23865   3.412 0.000645 ***
# av_air_temp_std -0.25554    0.09305  -2.746 0.006026 **
# daily_rain_std   1.11831    0.44821   2.495 0.012594 *
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

fixef(model_1_disp_rain)$cond
check_collinearity(model_1_disp_rain)

#DHARMa simulation
model_1_simulationOutput <- simulateResiduals(fittedModel = model_1_disp_rain, plot = F, n = 1000)
plotQQunif(model_1_simulationOutput, testUniformity = TRUE, 
           testOutliers = TRUE, testDispersion = TRUE)
plotResiduals(model_1_simulationOutput, quantreg = TRUE, smoothScatter = FALSE) # deviations in 0.25 quantile
testDispersion(model_1_simulationOutput)
testOutliers(model_1_simulationOutput)
testZeroInflation(model_1_simulationOutput)
mod_1_res = recalculateResiduals(model_1_simulationOutput, group = bat_metrics$time)
testTemporalAutocorrelation(mod_1_res, time = unique(bat_metrics$time))

plotResiduals(model_1_simulationOutput, form = bat_metrics$av_wind_speed_std, xlab = "Rank-transformed Standardised Wind Speed")
plotResiduals(model_1_simulationOutput, form = bat_metrics$av_air_temp_std, xlab = "Rank-transformed Standardised Air Temperature")
plotResiduals(model_1_simulationOutput, form = bat_metrics$meanMoonlightIntensity, xlab = "Rank-transformed Standardised Moonlight Intensity") #gotten worse
plotResiduals(model_1_simulationOutput, form = bat_metrics$av_cloud_okta_std, xlab = "Rank-transformed Standardised Cloud Cover") 
plotResiduals(model_1_simulationOutput, form = bat_metrics$daily_rain_std, xlab = "Rank-transformed Standardised Rainfall") # slight deviation in 0.25 quartile
plotResiduals(model_1_simulationOutput, form = bat_metrics$av_rel_humidity_std, xlab = "Rank-transformed Standardised Humidity")
plotResiduals(model_1_simulationOutput, form = bat_metrics$weight_std, xlab = "Rank-transformed Standardised Weight")
plotResiduals(model_1_simulationOutput, form = bat_metrics$age, xlab = "Rank-transformed Age")
plotResiduals(model_1_simulationOutput, form = bat_metrics$sex, xlab = "Rank-transformed Sex")
plotResiduals(model_1_simulationOutput, form = bat_metrics$study_year, xlab = "Rank-transformed Study Year")
plotResiduals(model_1_simulationOutput, form = bat_metrics$week_integer, xlab = "Rank-transformed Week")


#make simpler models
model_1.b <- glmmTMB(cbind(active, inactive) ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + weight_std + sex + age + study_year + week_integer + weight_std*av_air_temp_std + av_wind_speed_std*av_air_temp_std + (1 | ring_no),
                   data = bat_metrics,
                   dispformula = ~ week_integer + age + av_air_temp_std + daily_rain_std,
                   family = betabinomial(link = "logit"))

model_1.c <- glmmTMB(cbind(active, inactive) ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + weight_std + sex + age + study_year + week_integer + sex*age + av_wind_speed_std*av_air_temp_std + (1 | ring_no),
                   data = bat_metrics,
                   dispformula = ~ week_integer + age + av_air_temp_std + daily_rain_std,
                   family = betabinomial(link = "logit"))

model_1.d <- glmmTMB(cbind(active, inactive) ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + weight_std + sex + age + study_year + week_integer + sex*age + weight_std*av_air_temp_std + (1 | ring_no),
                   data = bat_metrics,
                   dispformula = ~ week_integer + age + av_air_temp_std + daily_rain_std,
                   family = betabinomial(link = "logit"))
        

# test interactions
anova(model_1_disp_rain, model_1.b) #sex*age
#              Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# model_1.b         50 6010.0 6236.0 -2955.0   5910.0
# model_1_disp_rain 52 6011.6 6246.6 -2953.8   5907.6 2.4444      2     0.2946

anova(model_1_disp_rain, model_1.c) #weight_std*av_air_temp_std
#              Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# model_1.c         51 6013.7 6244.2 -2955.8   5911.7
# model_1_disp_rain 52 6011.6 6246.6 -2953.8   5907.6 4.1141      1    0.04253 *

anova(model_1_disp_rain, model_1.d) #av_wind_speed_std*av_air_temp_std
#              Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# model_1.d         51 6011.7 6242.2 -2954.8   5909.7
# model_1_disp_rain 52 6011.6 6246.6 -2953.8   5907.6 2.1324      1     0.1442



#model selection drop1 ####
print("Dropping")
#model selection
drop1(model_1_disp_rain, test = "Chi")

#                                   Df    AIC     LRT Pr(>Chi)
# <none>                               6011.6
# meanMoonlightIntensity             1 6009.7   0.162  0.68694
# av_cloud_okta_std                  1 6009.6   0.001  0.97364
# daily_rain_std                     1 6009.9   0.326  0.56783
# av_rel_humidity_std                1 6012.0   2.484  0.11499    
# study_year                         1 6009.9   0.364  0.54652
# week_integer                      15 6152.9 171.363  < 2e-16 ***
# sex:age                            2 6010.0   2.444  0.29458
# av_air_temp_std:weight_std         1 6013.7   4.114  0.04253 *
# av_air_temp_std:av_wind_speed_std  1 6011.7   2.132  0.14422


# drop sex:age interaction
drop1(model_1.b, test = "Chi")

#                                   Df    AIC     LRT Pr(>Chi)    
# <none>                               6010.0
# meanMoonlightIntensity             1 6008.2   0.205  0.65102
# av_cloud_okta_std                  1 6008.0   0.004  0.94897
# daily_rain_std                     1 6008.3   0.307  0.57957
# av_rel_humidity_std                1 6010.6   2.625  0.10517
# sex                                1 6008.4   0.363  0.54660
# age                                2 6013.1   7.119  0.02846 *
# study_year                         1 6008.6   0.635  0.42567
# week_integer                      15 6150.0 169.971  < 2e-16 ***
# av_air_temp_std:weight_std         1 6011.9   3.912  0.04794 *
# av_air_temp_std:av_wind_speed_std  1 6010.1   2.053  0.15186    


# drop av_air_temp_std:av_wind_speed_std interaction
model_2 <- glmmTMB(cbind(active, inactive) ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + weight_std + sex + age + study_year + week_integer + weight_std*av_air_temp_std + (1 | ring_no),
                   data = bat_metrics,
                   dispformula = ~ week_integer + age + av_air_temp_std + daily_rain_std,
                   family = betabinomial(link = "logit"))

#model selection
drop1(model_2, test = "Chi")

#                            Df    AIC     LRT Pr(>Chi)    
# <none>                        6010.1
# av_wind_speed_std           1 6010.1   2.013  0.15599
# meanMoonlightIntensity      1 6008.3   0.232  0.62976
# av_cloud_okta_std           1 6008.1   0.001  0.97419
# daily_rain_std              1 6008.1   0.056  0.81291
# av_rel_humidity_std         1 6011.9   3.816  0.05075 .
# sex                         1 6008.5   0.494  0.48216    
# age                         2 6013.4   7.382  0.02495 *
# study_year                  1 6008.7   0.625  0.42922
# week_integer               15 6148.3 168.202  < 2e-16 ***
# av_air_temp_std:weight_std  1 6011.0   2.941  0.08637 .


# drop av_air_temp_std:weight_std interaction
model_3 <- glmmTMB(cbind(active, inactive) ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + weight_std + sex + age + study_year + week_integer + (1 | ring_no),
                   data = bat_metrics,
                   dispformula = ~ week_integer + age + av_air_temp_std + daily_rain_std,
                   family = betabinomial(link = "logit"))

#model selection
drop1(model_3, test = "Chi")

#                        Df    AIC     LRT Pr(>Chi)
# <none>                    6011.0
# av_air_temp_std         1 6016.5   7.497  0.00618 **
# av_wind_speed_std       1 6011.4   2.381  0.12280
# meanMoonlightIntensity  1 6009.1   0.141  0.70719
# av_cloud_okta_std       1 6009.0   0.005  0.94462
# daily_rain_std          1 6009.0   0.032  0.85710
# av_rel_humidity_std     1 6013.0   4.048  0.04421 *
# weight_std              1 6014.6   5.637  0.01758 *
# sex                     1 6009.7   0.677  0.41069    
# age                     2 6014.6   7.639  0.02194 *
# study_year              1 6009.6   0.605  0.43675
# week_integer           15 6148.9 167.958  < 2e-16 ***



#check model fits ####
summary(model_3)
#  Family: betabinomial  ( logit )
# Formula:
# cbind(active, inactive) ~ av_air_temp_std + av_wind_speed_std +
#     meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std +
#     av_rel_humidity_std + weight_std + sex + age + study_year +  
#     week_integer + (1 | ring_no)
# Dispersion:
# ~week_integer + age + av_air_temp_std + daily_rain_std
# Data: bat_metrics

#       AIC       BIC    logLik -2*log(L)  df.resid
#    6011.0    6228.0   -2957.5    5915.0       631

# Random effects:

# Conditional model:
#  Groups  Name        Variance Std.Dev.
#  ring_no (Intercept) 0.1162   0.3408
# Number of obs: 679, groups:  ring_no, 36

# Conditional model:
#                         Estimate Std. Error z value Pr(>|z|)
# (Intercept)             -0.43373    0.84422  -0.514 0.607416    
# av_air_temp_std          0.17718    0.06464   2.741 0.006124 **
# av_wind_speed_std       -0.51395    0.33511  -1.534 0.125112
# meanMoonlightIntensity  -0.14605    0.38994  -0.375 0.707992
# av_cloud_okta_std       -0.01032    0.14852  -0.069 0.944592
# daily_rain_std          -0.03172    0.17656  -0.180 0.857414
# av_rel_humidity_std     -1.58188    0.78479  -2.016 0.043835 *
# weight_std              -0.25795    0.10444  -2.470 0.013517 *  
# sexM                    -0.14551    0.17635  -0.825 0.409296
# ageJ                    -0.36849    0.21486  -1.715 0.086342 .
# ageYA                   -0.65934    0.23411  -2.816 0.004857 **
# study_yearSecond Winter -0.11269    0.14425  -0.781 0.434689
# week_integer2            0.25537    0.26827   0.952 0.341133
# week_integer3           -0.06340    0.27287  -0.232 0.816260
# week_integer4           -0.50181    0.37070  -1.354 0.175832    
# week_integer5           -0.39775    0.36003  -1.105 0.269263
# week_integer6            0.09669    0.33986   0.285 0.776025
# week_integer7           -0.46220    0.30237  -1.529 0.126372
# week_integer8           -0.38497    0.28801  -1.337 0.181332
# week_integer9           -0.85393    0.27743  -3.078 0.002084 **
# week_integer10          -0.80682    0.28184  -2.863 0.004200 **
# week_integer11          -1.14266    0.30570  -3.738 0.000186 ***
# week_integer12          -1.75886    0.30165  -5.831 5.52e-09 ***
# week_integer13          -1.90687    0.30044  -6.347 2.20e-10 ***
# week_integer14          -2.12006    0.31750  -6.677 2.43e-11 ***
# week_integer15          -2.50232    0.35290  -7.091 1.33e-12 ***
# week_integer16          -1.76561    0.39538  -4.466 7.99e-06 ***
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Dispersion model:
#                  Estimate Std. Error z value Pr(>|z|)    
# (Intercept)      2.371190   0.403345   5.879 4.13e-09 ***
# week_integer2    0.003727   0.452881   0.008 0.993434
# week_integer3    0.832260   0.482350   1.725 0.084450 .
# week_integer4    0.361546   0.585553   0.617 0.536943
# week_integer5    1.078806   0.606997   1.777 0.075522 .
# week_integer6   -0.106047   0.497801  -0.213 0.831302
# week_integer7    0.802057   0.496188   1.616 0.105999
# week_integer8    0.516010   0.442182   1.167 0.243226
# week_integer9    1.425407   0.441159   3.231 0.001233 **
# week_integer10   0.960319   0.432271   2.222 0.026313 *
# week_integer11   1.607891   0.514608   3.124 0.001781 **
# week_integer12   1.819374   0.469853   3.872 0.000108 ***
# week_integer13   1.792500   0.469069   3.821 0.000133 ***
# week_integer14   1.621891   0.463844   3.497 0.000471 ***
# week_integer15   2.446370   0.556892   4.393 1.12e-05 ***
# week_integer16   1.310743   0.583868   2.245 0.024773 *
# ageJ             0.476841   0.157747   3.023 0.002504 **
# ageYA            0.811094   0.239801   3.382 0.000719 ***
# av_air_temp_std -0.225192   0.093275  -2.414 0.015767 *
# daily_rain_std   1.040166   0.447057   2.327 0.019982 *
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


check_collinearity(model_3)


#DHARMa simulation
model_3_simulationOutput <- simulateResiduals(fittedModel = model_3, plot = F, n = 1000)
plotQQunif(model_3_simulationOutput, testUniformity = TRUE, 
           testOutliers = TRUE, testDispersion = TRUE)
plotResiduals(model_3_simulationOutput, quantreg = TRUE, smoothScatter = FALSE) # deviations in 0.25 quantile
testDispersion(model_3_simulationOutput)
testOutliers(model_3_simulationOutput)
testZeroInflation(model_3_simulationOutput)
mod_3_res = recalculateResiduals(model_3_simulationOutput, group = bat_metrics$time)
testTemporalAutocorrelation(mod_3_res, time = unique(bat_metrics$time))

plotResiduals(model_3_simulationOutput, form = bat_metrics$av_wind_speed_std, xlab = "Rank-transformed Standardised Wind Speed")
plotResiduals(model_3_simulationOutput, form = bat_metrics$av_air_temp_std, xlab = "Rank-transformed Standardised Air Temperature")
plotResiduals(model_3_simulationOutput, form = bat_metrics$meanMoonlightIntensity, xlab = "Rank-transformed Standardised Moonlight Intensity") #gotten worse
plotResiduals(model_3_simulationOutput, form = bat_metrics$av_cloud_okta_std, xlab = "Rank-transformed Standardised Cloud Cover") 
plotResiduals(model_3_simulationOutput, form = bat_metrics$daily_rain_std, xlab = "Rank-transformed Standardised Rainfall") # slight deviation in 0.25 quartile
plotResiduals(model_3_simulationOutput, form = bat_metrics$av_rel_humidity_std, xlab = "Rank-transformed Standardised Humidity")
plotResiduals(model_3_simulationOutput, form = bat_metrics$weight_std, xlab = "Rank-transformed Standardised Weight")
plotResiduals(model_3_simulationOutput, form = bat_metrics$age, xlab = "Rank-transformed Age")
plotResiduals(model_3_simulationOutput, form = bat_metrics$sex, xlab = "Rank-transformed Sex")
plotResiduals(model_3_simulationOutput, form = bat_metrics$study_year, xlab = "Rank-transformed Study Year")
plotResiduals(model_3_simulationOutput, form = bat_metrics$week_integer, xlab = "Rank-transformed Week")



#Test significance of main effects ####
model_3.a <- glmmTMB(cbind(active, inactive) ~ av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + weight_std + sex + age + study_year + week_integer + (1 | ring_no),
                   data = bat_metrics,
                   dispformula = ~ week_integer + age + av_air_temp_std + daily_rain_std,
                   family = betabinomial(link = "logit"))

model_3.b <- glmmTMB(cbind(active, inactive) ~ av_air_temp_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + weight_std + sex + age + study_year + week_integer + (1 | ring_no),
                   data = bat_metrics,
                   dispformula = ~ week_integer + age + av_air_temp_std + daily_rain_std,
                   family = betabinomial(link = "logit"))

model_3.c <- glmmTMB(cbind(active, inactive) ~ av_air_temp_std + av_wind_speed_std + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + weight_std + sex + age + study_year + week_integer + (1 | ring_no),
                   data = bat_metrics,
                   dispformula = ~ week_integer + age + av_air_temp_std + daily_rain_std,
                   family = betabinomial(link = "logit"))

model_3.d <- glmmTMB(cbind(active, inactive) ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + daily_rain_std + av_rel_humidity_std + weight_std + sex + age + study_year + week_integer + (1 | ring_no),
                   data = bat_metrics,
                   dispformula = ~ week_integer + age + av_air_temp_std + daily_rain_std,
                   family = betabinomial(link = "logit"))

model_3.e <- glmmTMB(cbind(active, inactive) ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + av_rel_humidity_std + weight_std + sex + age + study_year + week_integer + (1 | ring_no),
                   data = bat_metrics,
                   dispformula = ~ week_integer + age + av_air_temp_std + daily_rain_std,
                   family = betabinomial(link = "logit"))

model_3.f <- glmmTMB(cbind(active, inactive) ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + weight_std + sex + age + study_year + week_integer + (1 | ring_no),
                   data = bat_metrics,
                   dispformula = ~ week_integer + age + av_air_temp_std + daily_rain_std,
                   family = betabinomial(link = "logit"))

model_3.g <- glmmTMB(cbind(active, inactive) ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + sex + age + study_year + week_integer + (1 | ring_no),
                   data = bat_metrics,
                   dispformula = ~ week_integer + age + av_air_temp_std + daily_rain_std,
                   family = betabinomial(link = "logit"))

model_3.h <- glmmTMB(cbind(active, inactive) ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + weight_std + age + study_year + week_integer + (1 | ring_no),
                   data = bat_metrics,
                   dispformula = ~ week_integer + age + av_air_temp_std + daily_rain_std,
                   family = betabinomial(link = "logit"))

model_3.i <- glmmTMB(cbind(active, inactive) ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + weight_std + sex + study_year + week_integer + (1 | ring_no),
                   data = bat_metrics,
                   dispformula = ~ week_integer + age + av_air_temp_std + daily_rain_std,
                   family = betabinomial(link = "logit"))

model_3.j <- glmmTMB(cbind(active, inactive) ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + weight_std + sex + age + week_integer + (1 | ring_no),
                   data = bat_metrics,
                   dispformula = ~ week_integer + age + av_air_temp_std + daily_rain_std,
                   family = betabinomial(link = "logit"))

model_3.k <- glmmTMB(cbind(active, inactive) ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + weight_std + sex + age + study_year + (1 | ring_no),
                   data = bat_metrics,
                   dispformula = ~ week_integer + age + av_air_temp_std + daily_rain_std,
                   family = betabinomial(link = "logit"))

#main effect significance
anova(model_3, model_3.a) #av_air_temp_std
#           Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# model_3.a 47 6016.5 6229 -2961.2   5922.5
# model_3   48 6011.0 6228 -2957.5   5915.0 7.4971      1    0.00618 **

anova(model_3, model_3.b) #av_wind_speed_std
#           Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# model_3.b 47 6011.4 6223.8 -2958.7   5917.4
# model_3   48 6011.0 6228.0 -2957.5   5915.0 2.3812      1     0.1228

anova(model_3, model_3.c) #meanMoonlightIntensity
#           Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# model_3.c 47 6009.1 6221.6 -2957.6   5915.1
# model_3   48 6011.0 6228.0 -2957.5   5915.0 0.1411      1     0.7072

anova(model_3, model_3.d) #av_cloud_okta_std
#           Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# model_3.d 47 6009 6221.5 -2957.5     5915
# model_3   48 6011 6228.0 -2957.5     5915 0.0048      1     0.9446

anova(model_3, model_3.e) #daily_rain_std
#           Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# model_3.e 47 6009 6221.5 -2957.5     5915
# model_3   48 6011 6228.0 -2957.5     5915 0.0324      1     0.8571

anova(model_3, model_3.f) #av_rel_humidity_std
#           Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# model_3.f 47 6013 6225.5 -2959.5     5919
# model_3   48 6011 6228.0 -2957.5     5915 4.0484      1    0.04421 *

anova(model_3, model_3.g) #weight_std
#           Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)   
# model_3.g 47 6014.6 6227.1 -2960.3   5920.6
# model_3   48 6011.0 6228.0 -2957.5   5915.0 5.6371      1    0.01758 *

anova(model_3, model_3.h) #sex
#           Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# model_3.h 47 6009.7 6222.1 -2957.8   5915.7
# model_3   48 6011.0 6228.0 -2957.5   5915.0 0.6768      1     0.4107

anova(model_3, model_3.i) #age
#           Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)   
# model_3.i 46 6014.6 6222.6 -2961.3   5922.6
# model_3   48 6011.0 6228.0 -2957.5   5915.0 7.6392      2    0.02194 *

anova(model_3, model_3.j) #study_year
#           Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# model_3.j 47 6009.6 6222.1 -2957.8   5915.6
# model_3   48 6011.0 6228.0 -2957.5   5915.0 0.6048      1     0.4368

anova(model_3, model_3.k) #week_integer
#           Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# model_3.k 33 6148.9 6298.1 -3041.5   6082.9
# model_3   48 6011.0 6228.0 -2957.5   5915.0 167.96     15  < 2.2e-16 ***



# plotting ####
#av_air_temp_std
MyData_temp <- crossing(age = "A",
                             sex = "F",
                             study_year = "Second Winter",
                             week_integer = 2,
                             weight_std = mean(bat_metrics$weight_std),
                             av_wind_speed_std = mean(bat_metrics$av_wind_speed_std),
                             daily_rain_std = mean(bat_metrics$daily_rain_std),
                             av_air_temp_std = seq(from = min(bat_metrics$av_air_temp_std),
                                                   to = max(bat_metrics$av_air_temp_std),
                                                   length = 100),
                             av_rel_humidity_std = mean(bat_metrics$av_rel_humidity_std),
                             av_cloud_okta_std = mean(bat_metrics$av_cloud_okta_std),
                             meanMoonlightIntensity = mean(bat_metrics$meanMoonlightIntensity))

#' #' To compute population-level predictions for a given grouping
MyData_temp$ring_no <- NA
MyData_temp$ntrials <- round(mean(bat_metrics$inactive + bat_metrics$active))

#' Calculate the predicted  values.
Fit_temp <- predict(model_3,
                           newdata = MyData_temp,
                           re.form = NULL,  #~0, Keep it to NULL..gives odd results with ~0
                           type = "link",
                           se = TRUE)

attr(Fit_temp$fit, "names") <- NULL
attr(Fit_temp$se.fit, "names") <- NULL


#' Extract the predicted values and standard errors, and calculate
#' standard errors (SE) for predicted values.
MyData_temp$mu   <- plogis(Fit_temp$fit) 
MyData_temp$SeUp <- plogis(Fit_temp$fit + 1.96 * Fit_temp$se.fit)
MyData_temp$SeLo <- plogis(Fit_temp$fit - 1.96 * Fit_temp$se.fit)
MyData_temp$se   <- Fit_temp$se.fit

#' dispersion parameter (σ)
Fit_mod <- predict(model_3,
                   newdata = MyData_temp,
                   re.form = NA,
                   type = "disp")

attr(Fit_mod, "names") <- NULL

MyData_temp$phi   <- Fit_mod

#calculate variance and back standardize
MyData_temp <- mutate(MyData_temp,
                      con_var = mu * (1-mu) * (ntrials * ((phi + ntrials) / (phi + 1))),
                      prop_var = con_var/ntrials^2,
                      sd_prop = sqrt(prop_var),
                      upper_var = mu + 2 * sd_prop,
                      lower_var = mu - 2 * sd_prop,
                      av_air_temp = av_air_temp_std * sd(bat_metrics$av_air_temp) + mean(bat_metrics$av_air_temp))

#plot the result
bat_plot <- ggplot(MyData_temp, aes(x = av_air_temp, y = mu)) +
  geom_ribbon(aes(ymax = SeUp,
                  ymin = SeLo),
              alpha = 0.1,
              show.legend = FALSE) +
  geom_line(aes(y = upper_var), linetype = "dashed", colour = "blue", linewidth = 0.8) +
  geom_line(aes(y = lower_var), linetype = "dashed", colour = "blue", linewidth = 0.8) +
  geom_line(linewidth = 1) +
  geom_jitter(data = bat_metrics,
              aes(y = absolute_activity_proportion),
              alpha = 0.3,
              height = 0.001,
              width = 0.01) +
  ylim(c(0, (max(bat_metrics$absolute_activity_proportion)+0.05))) +
  labs(y = "Proportion of Activity",
       x = "Average Air Temperature (DegC)") +
  theme_bw()
bat_plot

ggsave(filename = paste0("figs/bat_daily_activity_model_no_hibernation/activity_model_trapping_temp.png"), plot = bat_plot, width = 10, height = 6, dpi = 300)


#daily_rain_std
MyData_rain <- crossing(age = "A",
                             sex = "F",
                             study_year = "Second Winter",
                             week_integer = 2,
                             weight_std = mean(bat_metrics$weight_std),
                             av_wind_speed_std = mean(bat_metrics$av_wind_speed_std),
                             daily_rain_std = seq(from = min(bat_metrics$daily_rain_std),
                                                  to = max(bat_metrics$daily_rain_std),
                                                  length = 100),
                             av_air_temp_std = mean(bat_metrics$av_air_temp_std),
                             av_rel_humidity_std = mean(bat_metrics$av_rel_humidity_std),
                             av_cloud_okta_std = mean(bat_metrics$av_cloud_okta_std),
                             meanMoonlightIntensity = mean(bat_metrics$meanMoonlightIntensity))

#' #' To compute population-level predictions for a given grouping
MyData_rain$ring_no <- NA
MyData_rain$ntrials <- round(mean(bat_metrics$inactive + bat_metrics$active))

#' Calculate the predicted  values.
Fit_rain <- predict(model_3,
                         newdata = MyData_rain,
                         re.form = NULL,  #~0, Keep it to NULL..gives odd results with ~0
                         type = "link",
                         se = TRUE)

attr(Fit_rain$fit, "names") <- NULL
attr(Fit_rain$se.fit, "names") <- NULL

#' Extract the predicted values and standard errors, and calculate
#' standard errors (SE) for predicted values.
MyData_rain$mu   <- plogis(Fit_rain$fit)
MyData_rain$SeUp <- plogis(Fit_rain$fit + 1.96 * Fit_rain$se.fit)
MyData_rain$SeLo <- plogis(Fit_rain$fit - 1.96 * Fit_rain$se.fit)
MyData_rain$se   <- Fit_rain$se.fit

#' dispersion parameter (phi)
Fit_mod <- predict(model_3,
                   newdata = MyData_rain,
                   re.form = NA,
                   type = "disp")

MyData_rain$phi   <- Fit_mod

#calculate variance and back standardize
MyData_rain <- mutate(MyData_rain,
                      con_var = mu * (1-mu) * (ntrials * ((phi + ntrials) / (phi + 1))),
                      prop_var = con_var/ntrials^2,
                      sd_prop = sqrt(prop_var),
                      upper_var = mu + 2 * sd_prop,
                      lower_var = mu - 2 * sd_prop,
                      daily_rain = daily_rain_std * max(bat_metrics$daily_rain))

#plot the result
bat_plot <- ggplot(MyData_rain, aes(x = daily_rain, y = mu)) +
  geom_ribbon(aes(ymax = SeUp,
                  ymin = SeLo),
              alpha = 0.1,
              show.legend = FALSE) +
  geom_line(aes(y = upper_var), linetype = "dashed", colour = "blue", linewidth = 0.8) +
  geom_line(aes(y = lower_var), linetype = "dashed", colour = "blue", linewidth = 0.8) +
  geom_line(linewidth = 1) +
  geom_jitter(data = bat_metrics,
              aes(y = absolute_activity_proportion),
              alpha = 0.3,
              height = 0.001,
              width = 0.01) +
  ylim(c(0, (max(bat_metrics$absolute_activity_proportion)+0.05))) +
  labs(y = "Proportion of Activity",
       x = "Total Rainfall (mm)") +
  theme_bw()
bat_plot

ggsave(filename = paste0("figs/bat_daily_activity_model_no_hibernation/activity_model_trapping_rain.png"), plot = bat_plot, width = 10, height = 6, dpi = 300)


#av_wind_speed_std
MyData_wind <- crossing(age = "A",
                        sex = "F",
                        study_year = "Second Winter",
                        week_integer = 2,
                        weight_std = mean(bat_metrics$weight_std),
                        av_wind_speed_std = seq(from = min(bat_metrics$av_wind_speed_std),
                                                to = max(bat_metrics$av_wind_speed_std),
                                                length = 100),
                        daily_rain_std = mean(bat_metrics$daily_rain_std),
                        av_air_temp_std = mean(bat_metrics$av_air_temp_std),
                        av_rel_humidity_std = mean(bat_metrics$av_rel_humidity_std),
                        av_cloud_okta_std = mean(bat_metrics$av_cloud_okta_std),
                        meanMoonlightIntensity = mean(bat_metrics$meanMoonlightIntensity))

#' #' To compute population-level predictions for a given grouping
MyData_wind$ring_no <- NA
MyData_wind$ntrials <- round(mean(bat_metrics$inactive + bat_metrics$active))

#' Calculate the predicted  values.
Fit_wind <- predict(model_3,
                    newdata = MyData_wind,
                    re.form = NULL,  #~0, Keep it to NULL..gives odd results with ~0
                    type = "link",
                    se = TRUE)

attr(Fit_wind$fit, "names") <- NULL
attr(Fit_wind$se.fit, "names") <- NULL


#' Extract the predicted values and standard errors, and calculate
#' standard errors (SE) for predicted values.
MyData_wind$mu   <- plogis(Fit_wind$fit)
MyData_wind$SeUp <- plogis(Fit_wind$fit + 1.96 * Fit_wind$se.fit)
MyData_wind$SeLo <- plogis(Fit_wind$fit - 1.96 * Fit_wind$se.fit)
MyData_wind$se   <- Fit_wind$se.fit

#back standardise
MyData_wind$av_wind_speed <- MyData_wind$av_wind_speed_std  * max(bat_metrics$av_wind_speed)

#plot the result
bat_plot <- ggplot(MyData_wind, aes(x = av_wind_speed, y = mu)) +
  geom_ribbon(aes(ymax = SeUp,
                  ymin = SeLo),
              alpha = 0.03,
              show.legend = FALSE) +
  geom_line(linewidth = 1) +
  geom_jitter(data = bat_metrics,
              aes(y = absolute_activity_proportion),
              alpha = 0.3,
              height = 0.001,
              width = 0.01) +
  labs(y = "Proportion of Activity",
       x = "Average Wind Speed (kn)") +
  ylim(c(0, (max(bat_metrics$absolute_activity_proportion)+0.05))) +
  theme_bw()
bat_plot

ggsave(filename = paste0("figs/bat_daily_activity_model_no_hibernation/activity_model_trapping_wind.png"), plot = bat_plot, width = 10, height = 6, dpi = 300)


#meanMoonlightIntensity
MyData_moon <- crossing(age = "A",
                        sex = "F",
                        study_year = "Second Winter",
                        week_integer = 2,
                        weight_std = mean(bat_metrics$weight_std),
                        av_wind_speed_std = mean(bat_metrics$av_wind_speed_std),
                        daily_rain_std = mean(bat_metrics$daily_rain_std),
                        av_air_temp_std = mean(bat_metrics$av_air_temp_std),
                        av_rel_humidity_std = mean(bat_metrics$av_rel_humidity_std),
                        av_cloud_okta_std = mean(bat_metrics$av_cloud_okta_std),
                        meanMoonlightIntensity = seq(from = min(bat_metrics$meanMoonlightIntensity),
                                                     to = max(bat_metrics$meanMoonlightIntensity),
                                                     length = 100))

#' #' To compute population-level predictions for a given grouping
MyData_moon$ring_no <- NA
MyData_moon$ntrials <- round(mean(bat_metrics$inactive + bat_metrics$active))

#' Calculate the predicted  values.
Fit_moon <- predict(model_3,
                    newdata = MyData_moon,
                    re.form = NULL,  #~0, Keep it to NULL..gives odd results with ~0
                    type = "link",
                    se = TRUE)

attr(Fit_moon$fit, "names") <- NULL
attr(Fit_moon$se.fit, "names") <- NULL

#' Extract the predicted values and standard errors, and calculate
#' standard errors (SE) for predicted values.
MyData_moon$mu   <- plogis(Fit_moon$fit)
MyData_moon$SeUp <- plogis(Fit_moon$fit + 1.96 * Fit_moon$se.fit)
MyData_moon$SeLo <- plogis(Fit_moon$fit - 1.96 * Fit_moon$se.fit)
MyData_moon$se   <- Fit_moon$se.fit

#back standardise
MyData_moon$meanMoonlightIntensity_lux <- MyData_moon$meanMoonlightIntensity * 0.3
bat_metrics$meanMoonlightIntensity_lux <- bat_metrics$meanMoonlightIntensity * 0.3

#plot the result
bat_plot <- ggplot(MyData_moon, aes(x = meanMoonlightIntensity_lux, y = mu)) +
  geom_ribbon(aes(ymax = SeUp,
                  ymin = SeLo),
              alpha = 0.03,
              show.legend = FALSE) +
  geom_line(linewidth = 1) +
  geom_jitter(data = bat_metrics,
              aes(y = absolute_activity_proportion),
              alpha = 0.3,
              height = 0.001,
              width = 0.01) +
  labs(y = "Proportion of Activity",
       x = "Average Moonlight Intensity (lux)") +
  ylim(c(0, (max(bat_metrics$absolute_activity_proportion)+0.05))) +
  theme_bw()
bat_plot

ggsave(filename = paste0("figs/bat_daily_activity_model_no_hibernation/activity_model_trapping_moon.png"), plot = bat_plot, width = 10, height = 6, dpi = 300)


#av_cloud_cover_Std
MyData_cloud <- crossing(age = "A",
                        sex = "F",
                        study_year = "Second Winter",
                        week_integer = 2,
                        weight_std = mean(bat_metrics$weight_std),
                        av_wind_speed_std = mean(bat_metrics$av_wind_speed_std),
                        daily_rain_std = mean(bat_metrics$daily_rain_std),
                        av_air_temp_std = mean(bat_metrics$av_air_temp_std),
                        av_rel_humidity_std = mean(bat_metrics$av_rel_humidity_std),
                        av_cloud_okta_std =seq(from = min(bat_metrics$av_cloud_okta_std),
                                               to = max(bat_metrics$av_cloud_okta_std),
                                               length = 100),
                        meanMoonlightIntensity = mean(bat_metrics$meanMoonlightIntensity))

#' #' To compute population-level predictions for a given grouping
MyData_cloud$ring_no <- NA
MyData_cloud$ntrials <- round(mean(bat_metrics$inactive + bat_metrics$active))

#' Calculate the predicted  values.
Fit_cloud <- predict(model_3,
                    newdata = MyData_cloud,
                    re.form = NULL,  #~0, Keep it to NULL..gives odd results with ~0
                    type = "link",
                    se = TRUE)

attr(Fit_cloud$fit, "names") <- NULL
attr(Fit_cloud$se.fit, "names") <- NULL

#' Extract the predicted values and standard errors, and calculate
#' standard errors (SE) for predicted values.
MyData_cloud$mu   <- plogis(Fit_cloud$fit)
MyData_cloud$SeUp <- plogis(Fit_cloud$fit + 1.96 * Fit_cloud$se.fit)
MyData_cloud$SeLo <- plogis(Fit_cloud$fit - 1.96 * Fit_cloud$se.fit)
MyData_cloud$se   <- Fit_cloud$se.fit

#back standardise
MyData_cloud$av_cloud_okta <- MyData_cloud$av_cloud_okta_std * max(bat_metrics$av_cloud_okta)

#plot the result
bat_plot <- ggplot(MyData_cloud, aes(x = av_cloud_okta, y = mu)) +
  geom_ribbon(aes(ymax = SeUp,
                  ymin = SeLo),
              alpha = 0.03,
              show.legend = FALSE) +
  geom_line(linewidth = 1) +
  geom_jitter(data = bat_metrics,
              aes(y = absolute_activity_proportion),
              alpha = 0.3,
              height = 0.001,
              width = 0.01) +
  labs(y = "Proportion of Activity",
       x = "Average Cloud Cover (Oktas)") +
  ylim(c(0, (max(bat_metrics$absolute_activity_proportion)+0.05))) +
  theme_bw()
bat_plot

ggsave(filename = paste0("figs/bat_daily_activity_model_no_hibernation/activity_model_trapping_cloud.png"), plot = bat_plot, width = 10, height = 6, dpi = 300)


#av_rel_humidity_std
MyData_humid <- crossing(age = "A",
                         sex = "F",
                         study_year = "Second Winter",
                         week_integer = 2,
                         weight_std = mean(bat_metrics$weight_std),
                         av_wind_speed_std = mean(bat_metrics$av_wind_speed_std),
                         daily_rain_std = mean(bat_metrics$daily_rain_std),
                         av_air_temp_std = mean(bat_metrics$av_air_temp_std),
                         av_rel_humidity_std = seq(from = min(bat_metrics$av_rel_humidity_std),
                                                   to = max(bat_metrics$av_rel_humidity_std),
                                                   length = 100),
                         av_cloud_okta_std = mean(bat_metrics$av_cloud_okta_std),
                         meanMoonlightIntensity = mean(bat_metrics$meanMoonlightIntensity))

#' #' To compute population-level predictions for a given grouping
MyData_humid$ring_no <- NA
MyData_humid$ntrials <- round(mean(bat_metrics$inactive + bat_metrics$active))

#' Calculate the predicted  values.
Fit_humid <- predict(model_3,
                     newdata = MyData_humid,
                     re.form = NULL,  #~0, Keep it to NULL..gives odd results with ~0
                     type = "link",
                     se = TRUE)

attr(Fit_humid$fit, "names") <- NULL
attr(Fit_humid$se.fit, "names") <- NULL

#' Extract the predicted values and standard errors, and calculate
#' standard errors (SE) for predicted values.
MyData_humid$mu   <- plogis(Fit_humid$fit)
MyData_humid$SeUp <- plogis(Fit_humid$fit + 1.96 * Fit_humid$se.fit)
MyData_humid$SeLo <- plogis(Fit_humid$fit - 1.96 * Fit_humid$se.fit)
MyData_humid$se   <- Fit_humid$se.fit

#back standardise
MyData_humid$av_rel_humidity <- MyData_humid$av_rel_humidity_std * max(bat_metrics$av_rel_humidity)

#plot the result
bat_plot <- ggplot(MyData_humid, aes(x = av_rel_humidity, y = mu)) +
  geom_ribbon(aes(ymax = SeUp,
                  ymin = SeLo),
              alpha = 0.03,
              show.legend = FALSE) +
  geom_line(linewidth = 1) +
  geom_jitter(data = bat_metrics,
              aes(y = absolute_activity_proportion),
              alpha = 0.3,
              height = 0.001,
              width = 0.01) +
  labs(y = "Proportion of Activity",
       x = "Average Relative Humidity (%)") +
  ylim(c(0, (max(bat_metrics$absolute_activity_proportion)+0.05))) +
  theme_bw()
bat_plot

ggsave(filename = paste0("figs/bat_daily_activity_model_no_hibernation/activity_model_trapping_humid.png"), plot = bat_plot, width = 10, height = 6, dpi = 300)


#weight_std
MyData_weight <- crossing(age = "A",
                         sex = "F",
                         study_year = "Second Winter",
                         week_integer = 2,
                         weight_std = seq(from = min(bat_metrics$weight_std),
                                          to = max(bat_metrics$weight_std),
                                          length = 100),
                         av_wind_speed_std = mean(bat_metrics$av_wind_speed_std),
                         daily_rain_std = mean(bat_metrics$daily_rain_std),
                         av_air_temp_std = mean(bat_metrics$av_air_temp_std),
                         av_rel_humidity_std = mean(bat_metrics$av_rel_humidity_std),
                         av_cloud_okta_std = mean(bat_metrics$av_cloud_okta_std),
                         meanMoonlightIntensity = mean(bat_metrics$meanMoonlightIntensity))

#' #' To compute population-level predictions for a given grouping
MyData_weight$ring_no <- NA
MyData_weight$ntrials <- round(mean(bat_metrics$inactive + bat_metrics$active))

#' Calculate the predicted  values.
Fit_weight <- predict(model_3,
                     newdata = MyData_weight,
                     re.form = NULL,  #~0, Keep it to NULL..gives odd results with ~0
                     type = "link",
                     se = TRUE)

attr(Fit_weight$fit, "names") <- NULL
attr(Fit_weight$se.fit, "names") <- NULL

#' Extract the predicted values and standard errors, and calculate
#' standard errors (SE) for predicted values.
MyData_weight$mu   <- plogis(Fit_weight$fit)
MyData_weight$SeUp <- plogis(Fit_weight$fit + 1.96 * Fit_weight$se.fit)
MyData_weight$SeLo <- plogis(Fit_weight$fit - 1.96 * Fit_weight$se.fit)
MyData_weight$se   <- Fit_weight$se.fit

#back standardise
MyData_weight$weight <- MyData_weight$weight_std  * sd(bat_metrics$weight) +
  mean(bat_metrics$weight)

#plot the result
bat_plot <- ggplot(MyData_weight, aes(x = weight, y = mu)) +
  geom_ribbon(aes(ymax = SeUp,
                  ymin = SeLo),
              alpha = 0.03,
              show.legend = FALSE) +
  geom_line(linewidth = 1) +
  geom_jitter(data = bat_metrics,
              aes(y = absolute_activity_proportion),
              alpha = 0.3,
              height = 0.001,
              width = 0.01) +
  labs(y = "Proportion of Activity",
       x = "Weight (g)") +
  ylim(c(0, (max(bat_metrics$absolute_activity_proportion)+0.05))) +
  theme_bw()
bat_plot

ggsave(filename = paste0("figs/bat_daily_activity_model_no_hibernation/activity_model_trapping_weight.png"), plot = bat_plot, width = 10, height = 6, dpi = 300)



#eemeans for categorical variable ####
#age effect
m3_emm <- emmeans(model_3, ~ age, component = "cmean")

### pairwise comparisons
m3_contrast <- summary(contrast(m3_emm,
                                  method = "pairwise"),
                         infer = TRUE)

#Make into a dataframe
m3_emm_df <- data.frame(m3_emm)
m3_emm_contrast_df <- data.frame(m3_contrast)

# > m3_emm_df
#   age     emmean          SE  df  asymp.LCL  asymp.UCL
# 1   A 0.05639461 0.006775997 Inf 0.04311389 0.06967532
# 2   J 0.04000257 0.006582228 Inf 0.02710164 0.05290350
# 3  YA 0.03034629 0.005820845 Inf 0.01893764 0.04175493

# > m3_emm_contrast_df
#   contrast    estimate          SE  df    asymp.LCL  asymp.UCL  z.ratio     p.value
# 1    A - J 0.016392034 0.009302539 Inf -0.005410332 0.03819440 1.762103 0.182484668
# 2   A - YA 0.026048320 0.008707657 Inf  0.005640178 0.04645646 2.991427 0.007819646
# 3   J - YA 0.009656286 0.008019867 Inf -0.009139881 0.02845245 1.204046 0.450736908

#age
MyData <- crossing(age = levels(bat_metrics$age),
                         sex = "F",
                         study_year = "Second Winter",
                         week_integer = 2,
                         weight_std = mean(bat_metrics$weight_std),
                         av_wind_speed_std = mean(bat_metrics$av_wind_speed_std),
                         daily_rain_std = mean(bat_metrics$daily_rain_std),
                         av_air_temp_std = mean(bat_metrics$av_air_temp_std),
                         av_rel_humidity_std = mean(bat_metrics$av_rel_humidity_std),
                         av_cloud_okta_std = mean(bat_metrics$av_cloud_okta_std),
                         meanMoonlightIntensity = mean(bat_metrics$meanMoonlightIntensity))

#' #' To compute population-level predictions for a given grouping
MyData$ring_no <- NA
MyData$ntrials <- round(mean(bat_metrics$inactive + bat_metrics$active))

#' Calculate the predicted  values.
Fit_mod <- predict(model_3,
                     newdata = MyData,
                     re.form = NULL,  #~0, Keep it to NULL..gives odd results with ~0
                     type = "disp")

attr(Fit_mod, "names") <- NULL

MyData$phi   <- Fit_mod

# calculate variance ± SD
m3_emm_df <- m3_emm_df %>%
  left_join(MyData %>% select(age, phi, ntrials), by = "age")

# calculate variance and convert to proportion
m3_emm_df <- m3_emm_df %>% 
  mutate(
    variance_counts = emmean * (1-emmean) * (ntrials * ((phi + ntrials) / (phi + 1))),
    variance_props = variance_counts / (ntrials^2),
    sd_prop = sqrt(variance_props),
    SD_plus = emmean + 2*sd_prop,
    SD_minus = emmean - 2*sd_prop
  )


#Change names
m3_emm_df$age <- factor(m3_emm_df$age , levels=c("J", "YA", "A"))
m3_emm_df$age <- plyr::revalue(m3_emm_df$age, c("J" = "Juvenile", "YA" = "Young Adult", "A" = "Adult"))

bat_metrics$age <- factor(bat_metrics$age , levels=c("J", "YA", "A"))
bat_metrics$age <- plyr::revalue(bat_metrics$age, c("J" = "Juvenile", "YA" = "Young Adult", "A" = "Adult"))

#Make contrasting groups match variable names
m3_emm_contrast_df$group1 <- NA
m3_emm_contrast_df$group2 <- NA
m3_emm_contrast_df$contrast_pretty <- NA

m3_emm_contrast_df$group1 <- c("Adult", "Adult", "Juvenile")
m3_emm_contrast_df$group2 <- c("Juvenile", "Young Adult", "Young Adult")
m3_emm_contrast_df$contrast_pretty <- c("Adult - Juvenile", "Adult - Young Adult", "Juvenile - Young Adult")

m3_emm_contrast_df$p_rounded <- rstatix::p_round(m3_emm_contrast_df$p.value, digits = 2)
m3_emm_contrast_df$p_pretty <- rstatix::p_format(rstatix::p_round(m3_emm_contrast_df$p.value, digits = 2),
                                              accuracy = 1e-04, add.p = TRUE)


#plot the result
bat_plot <- ggplot(bat_metrics, aes(x = age, y = absolute_activity_proportion)) +
  geom_jitter(alpha = 0.3,
              height = 0.001,
              width = 0.3) +
  geom_point(data = m3_emm_df,
             aes(y = emmean,
                 x = age), 
             size = 4,
             colour = "red") +
  geom_errorbar(data = m3_emm_df, 
                aes(y = emmean,
                    x = age,
                    ymin = asymp.LCL, 
                    ymax = asymp.UCL), 
                width = 0.05,
                colour = "red") +
  geom_segment(data = m3_emm_df,
             aes(x = as.numeric(age) -  0.3,
                 xend = as.numeric(age) +  0.3,
                 y = SD_minus, yend = SD_minus),
             colour = "blue",
             linetype = "dashed",
             linewidth = 0.8,
             inherit.aes = FALSE) +
geom_segment(data = m3_emm_df,
             aes(x = as.numeric(age) -  0.3,
                 xend = as.numeric(age) +  0.3,
                 y = SD_plus, yend = SD_plus),
             colour = "blue",
             linetype = "dashed",
             linewidth = 0.8,
             inherit.aes = FALSE) +
  ylim(c(0, (max(bat_metrics$absolute_activity_proportion)+0.05))) +
  labs(y = "Proportion of Activity",
       x = "Age") +
  theme_pubr()
bat_plot

bat_effects <- ggplot(data = m3_emm_contrast_df, aes(x = estimate, y = contrast_pretty)) +
  geom_point(size = 4) +
  xlab("Effects (Proportion)") +
  geom_errorbar(aes(x = estimate,
                    xmin = asymp.LCL,
                    xmax = asymp.UCL),
                width = 0.1) +
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "blue") +
  annotate(geom = "text",
           x = m3_emm_contrast_df$estimate,
           y = 1:3 + 0.2,
           label = m3_emm_contrast_df$p_pretty,
           size = 3) +
  theme_pubr() +
  rremove("ylab")
bat_effects


# make int one plot
m3_bottom <- bat_plot +
  coord_flip() # rotate to horizontal

m3_top <- bat_effects + # assign to new plot object
  scale_x_continuous(position="top") # move x-axis to "top"

m3_plot <- cowplot::plot_grid(m3_top,
                     m3_bottom,
                     nrow = 2,
                     align = "v",
                     axis = "lr",
                     rel_heights = c(1, 1))
m3_plot

ggsave(filename = paste0("figs/bat_daily_activity_model_no_hibernation/activity_model_trapping_age.png"), plot = bat_plot, width = 10, height = 6, dpi = 300)

ggsave(filename = paste0("figs/bat_daily_activity_model_no_hibernation/activity_model_trapping_age_effect.png"), plot = m3_plot, width = 10, height = 6, dpi = 300)


#sex effect
m3_emm <- emmeans(model_3, ~ sex)

### pairwise comparisons
m3_contrast <- summary(contrast(m3_emm,
                                  method = "pairwise"),
                         infer = TRUE)

#Make into a dataframe
m3_emm_df <- data.frame(m3_emm)
m3_emm_contrast_df <- data.frame(m3_contrast)


# > m3_emm_df
#   sex    emmean        SE  df asymp.LCL asymp.UCL
# 1   F -3.362090 0.1286353 Inf -3.614210 -3.109969
# 2   M -3.507599 0.1314020 Inf -3.765142 -3.250055

# > m3_emm_contrast_df
#   contrast  estimate       SE  df  asymp.LCL asymp.UCL   z.ratio   p.value
# 1    F - M 0.1455089 0.176346 Inf -0.2001228 0.4911406 0.8251333 0.4092959


#Change names
m3_emm_df$sex <- factor(m3_emm_df$sex , levels=c("F", "M"))
m3_emm_df$sex <- plyr::revalue(m3_emm_df$sex, c("F" = "Female", "M" = "Male"))

bat_metrics$sex <- factor(bat_metrics$sex , levels=c("F", "M"))
bat_metrics$sex <- plyr::revalue(bat_metrics$sex, c("F" = "Female", "M" = "Male"))

m3_emm_df$prob <- exp(m3_emm_df$emmean)/(1 + exp(m3_emm_df$emmean))
m3_emm_df$LCL <- exp(m3_emm_df$asymp.LCL)/(1 + exp(m3_emm_df$asymp.LCL))
m3_emm_df$UCL <- exp(m3_emm_df$asymp.UCL)/(1 + exp(m3_emm_df$asymp.UCL))

#Make contrasting groups match variable names
m3_emm_contrast_df$group1 <- NA
m3_emm_contrast_df$group2 <- NA
m3_emm_contrast_df$contrast_pretty <- NA

m3_emm_contrast_df$group1 <- c("Female")
m3_emm_contrast_df$group2 <- c("Male")
m3_emm_contrast_df$contrast_pretty <- c("Female - Male")

m3_emm_contrast_df$p_rounded <- rstatix::p_round(m3_emm_contrast_df$p.value, digits = 2)
m3_emm_contrast_df$p_pretty <- rstatix::p_format(rstatix::p_round(m3_emm_contrast_df$p.value, digits = 2),
                                              accuracy = 1e-04, add.p = TRUE)


#plot the result
bat_plot <- ggplot(bat_metrics, aes(x = sex, y = absolute_activity_proportion)) +
  geom_jitter(alpha = 0.3,
              height = 0.001,
              width = 0.3) +
  geom_point(data = m3_emm_df,
             aes(y = prob,
                 x = sex), 
             size = 4,
             colour = "red") +
  geom_errorbar(data = m3_emm_df, 
                aes(y = prob,
                    x = sex,
                    ymin = LCL, 
                    ymax = UCL), 
                width = 0.05,
                colour = "red") +
  ylim(c(0, (max(bat_metrics$absolute_activity_proportion)+0.05))) +
  labs(y = "Proportion of Activity",
       x = "Sex") +
  theme_pubr()
bat_plot

bat_effects <- ggplot(data = m3_emm_contrast_df, aes(x = estimate, y = contrast_pretty)) +
  geom_point(size = 4) +
  xlab("Effects (Proportion)") +
  geom_errorbar(aes(x = estimate,
                    xmin = asymp.LCL,
                    xmax = asymp.UCL),
                width = 0.1) +
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "blue") +
  annotate(geom = "text",
           x = m3_emm_contrast_df$estimate,
           y = 1 + 0.2,
           label = m3_emm_contrast_df$p_pretty,
           size = 3) +
  theme_pubr() +
  rremove("ylab")
bat_effects


# make into one plot
m3_bottom <- bat_plot +
  coord_flip() # rotate to horizontal

m3_top <- bat_effects + # assign to new plot object
  scale_x_continuous(position="top") # move x-axis to "top"

m3_plot <- cowplot::plot_grid(m3_top,
                     m3_bottom,
                     nrow = 2,
                     align = "v",
                     axis = "lr",
                     rel_heights = c(1, 1))
m3_plot

ggsave(filename = paste0("figs/bat_daily_activity_model_no_hibernation/activity_model_trapping_sex.png"), plot = bat_plot, width = 10, height = 6, dpi = 300)

ggsave(filename = paste0("figs/bat_daily_activity_model_no_hibernation/activity_model_trapping_sex_effect.png"), plot = m3_plot, width = 10, height = 6, dpi = 300)


#year effect
m3_emm <- emmeans(model_3, ~ study_year)

### pairwise comparisons
m3_contrast <- summary(contrast(m3_emm,
                                  method = "pairwise"),
                         infer = TRUE)

#Make into a dataframe
m3_emm_df <- data.frame(m3_emm)
m3_emm_contrast_df <- data.frame(m3_contrast)


# > m3_emm_df
#      study_year    emmean         SE  df asymp.LCL asymp.UCL
# 1  First Winter -3.378501 0.13714433 Inf -3.647299 -3.109703
# 2 Second Winter -3.491188 0.09929743 Inf -3.685807 -3.296568

# > m3_emm_contrast_df
#                       contrast  estimate        SE  df  asymp.LCL asymp.UCL   z.ratio   p.value
# 1 First Winter - Second Winter 0.1126868 0.1442496 Inf -0.1700372 0.3954108 0.7811932 0.4346889


#Change names
m3_emm_df$study_year <- factor(m3_emm_df$study_year , levels=c("First Winter", "Second Winter"))

bat_metrics$study_year <- factor(bat_metrics$study_year , levels=c("First Winter", "Second Winter"))

m3_emm_df$prob <- exp(m3_emm_df$emmean)/(1 + exp(m3_emm_df$emmean))
m3_emm_df$LCL <- exp(m3_emm_df$asymp.LCL)/(1 + exp(m3_emm_df$asymp.LCL))
m3_emm_df$UCL <- exp(m3_emm_df$asymp.UCL)/(1 + exp(m3_emm_df$asymp.UCL))

#Make contrasting groups match variable names
m3_emm_contrast_df$group1 <- NA
m3_emm_contrast_df$group2 <- NA
m3_emm_contrast_df$contrast_pretty <- NA

m3_emm_contrast_df$group1 <- c("First Winter")
m3_emm_contrast_df$group2 <- c("Second Winter")
m3_emm_contrast_df$contrast_pretty <- c("First - Second")

m3_emm_contrast_df$p_rounded <- rstatix::p_round(m3_emm_contrast_df$p.value, digits = 2)
m3_emm_contrast_df$p_pretty <- rstatix::p_format(rstatix::p_round(m3_emm_contrast_df$p.value, digits = 2),
                                              accuracy = 1e-04, add.p = TRUE)


#plot the result
bat_plot <- ggplot(bat_metrics, aes(x = study_year, y = absolute_activity_proportion)) +
  geom_jitter(alpha = 0.3,
              height = 0.001,
              width = 0.3) +
  geom_point(data = m3_emm_df,
             aes(y = prob,
                 x = study_year), 
             size = 4,
             colour = "red") +
  geom_errorbar(data = m3_emm_df, 
                aes(y = prob,
                    x = study_year,
                    ymin = LCL, 
                    ymax = UCL), 
                width = 0.05,
                colour = "red") +
  ylim(c(0, (max(bat_metrics$absolute_activity_proportion)+0.05))) +
  labs(y = "Proportion of Activity",
       x = "Study Year") +
  theme_pubr()
bat_plot

bat_effects <- ggplot(data = m3_emm_contrast_df, aes(x = estimate, y = contrast_pretty)) +
  geom_point(size = 4) +
  xlab("Effects (Proportion)") +
  geom_errorbar(aes(x = estimate,
                    xmin = asymp.LCL,
                    xmax = asymp.UCL),
                width = 0.1) +
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "blue") +
  annotate(geom = "text",
           x = m3_emm_contrast_df$estimate,
           y = 1 + 0.2,
           label = m3_emm_contrast_df$p_pretty,
           size = 3) +
  theme_pubr() +
  rremove("ylab")
bat_effects


# make into one plot
m3_bottom <- bat_plot +
  coord_flip() # rotate to horizontal

m3_top <- bat_effects + # assign to new plot object
  scale_x_continuous(position="top") # move x-axis to "top"

m3_plot <- cowplot::plot_grid(m3_top,
                     m3_bottom,
                     nrow = 2,
                     align = "v",
                     axis = "lr",
                     rel_heights = c(1, 1))
m3_plot

ggsave(filename = paste0("figs/bat_daily_activity_model_no_hibernation/activity_model_trapping_year.png"), plot = bat_plot, width = 10, height = 6, dpi = 300)

ggsave(filename = paste0("figs/bat_daily_activity_model_no_hibernation/activity_model_trapping_year_effect.png"), plot = m3_plot, width = 10, height = 6, dpi = 300)


#week effect
m3_emm <- emmeans(model_3, ~ week_integer, component = "cmean")

### pairwise comparisons
m3_contrast <- summary(contrast(m3_emm,
                                method = "pairwise"),
                       infer = TRUE)

#Make into a dataframe
m3_emm_df <- data.frame(m3_emm)
m3_emm_contrast_df <- data.frame(m3_contrast)


# > m3_emm_df
#    week_integer      emmean          SE  df   asymp.LCL   asymp.UCL
# 1             1 0.075465297 0.017842975 Inf 0.040493708 0.110436886
# 2             2 0.095188209 0.017620626 Inf 0.060652416 0.129724001
# 3             3 0.071182310 0.012762196 Inf 0.046168866 0.096195754
# 4             4 0.047196123 0.013711772 Inf 0.020321544 0.074070703
# 5             5 0.052082529 0.013850190 Inf 0.024936656 0.079228401
# 6             6 0.082450918 0.018398836 Inf 0.046389862 0.118511974
# 7             7 0.049002993 0.007602254 Inf 0.034102849 0.063903137
# 8             8 0.052714506 0.007247419 Inf 0.038509825 0.066919187
# 9             9 0.033695866 0.003894410 Inf 0.026062963 0.041328769
# 10           10 0.035259581 0.004207658 Inf 0.027012722 0.043506439
# 11           11 0.025477442 0.003728240 Inf 0.018170226 0.032784659
# 12           12 0.013933868 0.002083534 Inf 0.009850216 0.018017520
# 13           13 0.012041850 0.001957349 Inf 0.008205516 0.015878184
# 14           14 0.009754221 0.001798497 Inf 0.006229232 0.013279211
# 15           15 0.006677810 0.001556381 Inf 0.003627360 0.009728261
# 16           16 0.013841611 0.004066164 Inf 0.005872077 0.021811146

# > m3_emm_contrast_df
#                            contrast      estimate          SE  df     asymp.LCL   asymp.UCL     z.ratio      p.value
# 1     week_integer1 - week_integer2 -1.972291e-02 0.020284743 Inf -8.921928e-02 0.049773457 -0.97230277 9.998650e-01
# 2     week_integer1 - week_integer3  4.282987e-03 0.018635683 Inf -5.956363e-02 0.068129608  0.22982719 1.000000e+00
# 3     week_integer1 - week_integer4  2.826917e-02 0.020895421 Inf -4.331940e-02 0.099857751  1.35288845 9.938598e-01
# 4     week_integer1 - week_integer5  2.338277e-02 0.021354437 Inf -4.977842e-02 0.096543955  1.09498406 9.994248e-01
# 5     week_integer1 - week_integer6 -6.985621e-03 0.024534912 Inf -9.104324e-02 0.077072002 -0.28472166 1.000000e+00
# 6     week_integer1 - week_integer7  2.646230e-02 0.019263364 Inf -3.953478e-02 0.092459386  1.37371145 9.928060e-01
# 7     week_integer1 - week_integer8  2.275079e-02 0.018854919 Inf -4.184694e-02 0.087348523  1.20662366 9.982287e-01
# 8     week_integer1 - week_integer9  4.176943e-02 0.018039997 Inf -2.003635e-02 0.103575209  2.31537900 6.121626e-01
# 9    week_integer1 - week_integer10  4.020572e-02 0.018194932 Inf -2.213087e-02 0.102542307  2.20972057 6.904383e-01
# 10   week_integer1 - week_integer11  4.998785e-02 0.018445452 Inf -1.320703e-02 0.113182736  2.71003688 3.237205e-01
# 11   week_integer1 - week_integer12  6.153143e-02 0.018011372 Inf -1.762766e-04 0.123239134  3.41625442 5.157978e-02
# 12   week_integer1 - week_integer13  6.342345e-02 0.017881584 Inf  2.160401e-03 0.124686493  3.54685850 3.366303e-02
# 13   week_integer1 - week_integer14  6.571108e-02 0.017931544 Inf  4.276863e-03 0.127145289  3.66455194 2.244020e-02
# 14   week_integer1 - week_integer15  6.878749e-02 0.017953056 Inf  7.279575e-03 0.130295398  3.83151972 1.221585e-02
# 15   week_integer1 - week_integer16  6.162369e-02 0.018333686 Inf -1.188281e-03 0.124435653  3.36122729 6.127356e-02
# 16    week_integer2 - week_integer3  2.400590e-02 0.016653355 Inf -3.304918e-02 0.081060982  1.44150524 9.882947e-01
# 17    week_integer2 - week_integer4  4.799209e-02 0.019416564 Inf -1.852987e-02 0.114514036  2.47170850 4.928116e-01
# 18    week_integer2 - week_integer5  4.310568e-02 0.020536280 Inf -2.725247e-02 0.113463826  2.09900136 7.661636e-01
# 19    week_integer2 - week_integer6  1.273729e-02 0.023361068 Inf -6.729869e-02 0.092773276  0.54523582 9.999999e-01
# 20    week_integer2 - week_integer7  4.618522e-02 0.018854541 Inf -1.841122e-02 0.110781654  2.44955392 5.096398e-01
# 21    week_integer2 - week_integer8  4.247370e-02 0.018475680 Inf -2.082474e-02 0.105772148  2.29889791 6.246228e-01
# 22    week_integer2 - week_integer9  6.149234e-02 0.017924864 Inf  8.101847e-05 0.122903667  3.43056127 4.928446e-02
# 23   week_integer2 - week_integer10  5.992863e-02 0.017988012 Inf -1.699045e-03 0.121556301  3.33158706 6.710174e-02
# 24   week_integer2 - week_integer11  6.971077e-02 0.018349682 Inf  6.843997e-03 0.132577536  3.79901773 1.379153e-02
# 25   week_integer2 - week_integer12  8.125434e-02 0.017974527 Inf  1.967287e-02 0.142835813  4.52052742 6.830230e-04
# 26   week_integer2 - week_integer13  8.314636e-02 0.017666486 Inf  2.262025e-02 0.143672473  4.70644569 2.850120e-04
# 27   week_integer2 - week_integer14  8.543399e-02 0.017725823 Inf  2.470458e-02 0.146163391  4.81974723 1.641984e-04
# 28   week_integer2 - week_integer15  8.851040e-02 0.017708053 Inf  2.784188e-02 0.149178921  4.99831335 6.693074e-05
# 29   week_integer2 - week_integer16  8.134660e-02 0.018217169 Inf  1.893382e-02 0.143759373  4.46538080 8.785315e-04
# 30    week_integer3 - week_integer4  2.398619e-02 0.015457604 Inf -2.897220e-02 0.076944578  1.55174029 9.763159e-01
# 31    week_integer3 - week_integer5  1.909978e-02 0.016388299 Inf -3.704721e-02 0.075246771  1.16545239 9.988072e-01
# 32    week_integer3 - week_integer6 -1.126861e-02 0.020162288 Inf -8.034544e-02 0.057808224 -0.55889529 9.999999e-01
# 33    week_integer3 - week_integer7  2.217932e-02 0.014231273 Inf -2.657761e-02 0.070936247  1.55849144 9.753511e-01
# 34    week_integer3 - week_integer8  1.846780e-02 0.013834241 Inf -2.892888e-02 0.065864486  1.33493444 9.946623e-01
# 35    week_integer3 - week_integer9  3.748644e-02 0.012848700 Inf -6.533734e-03 0.081506623  2.91752817 2.060876e-01
# 36   week_integer3 - week_integer10  3.592273e-02 0.013114276 Inf -9.007323e-03 0.080852782  2.73920795 3.052470e-01
# 37   week_integer3 - week_integer11  4.570487e-02 0.013557117 Inf -7.423768e-04 0.092152113  3.37128218 5.939566e-02
# 38   week_integer3 - week_integer12  5.724844e-02 0.013067406 Inf  1.247897e-02 0.102017917  4.38101036 1.282528e-03
# 39   week_integer3 - week_integer13  5.914046e-02 0.012779827 Inf  1.535624e-02 0.102924677  4.62764169 4.147683e-04
# 40   week_integer3 - week_integer14  6.142809e-02 0.012872444 Inf  1.732656e-02 0.105529614  4.77206115 2.074494e-04
# 41   week_integer3 - week_integer15  6.450450e-02 0.012907539 Inf  2.028274e-02 0.108726262  4.99742836 6.723477e-05
# 42   week_integer3 - week_integer16  5.734070e-02 0.013505328 Inf  1.107089e-02 0.103610512  4.24578347 2.311328e-03
# 43    week_integer4 - week_integer5 -4.886405e-03 0.017244620 Inf -6.396719e-02 0.054194377 -0.28335825 1.000000e+00
# 44    week_integer4 - week_integer6 -3.525479e-02 0.021247958 Inf -1.080512e-01 0.037541590 -1.65920857 9.570608e-01
# 45    week_integer4 - week_integer7 -1.806870e-03 0.015123436 Inf -5.362039e-02 0.050006648 -0.11947484 1.000000e+00
# 46    week_integer4 - week_integer8 -5.518383e-03 0.014785934 Inf -5.617560e-02 0.045138838 -0.37321841 1.000000e+00
# 47    week_integer4 - week_integer9  1.350026e-02 0.013871735 Inf -3.402488e-02 0.061025395  0.97322055 9.998634e-01
# 48   week_integer4 - week_integer10  1.193654e-02 0.014055940 Inf -3.621969e-02 0.060092775  0.84921693 9.999761e-01
# 49   week_integer4 - week_integer11  2.171868e-02 0.014258279 Inf -2.713077e-02 0.070568135  1.52323298 9.800700e-01
# 50   week_integer4 - week_integer12  3.326225e-02 0.013931689 Inf -1.446829e-02 0.080992799  2.38752487 5.570836e-01
# 51   week_integer4 - week_integer13  3.515427e-02 0.013724021 Inf -1.186479e-02 0.082173338  2.56151407 4.259863e-01
# 52   week_integer4 - week_integer14  3.744190e-02 0.013765331 Inf -9.718691e-03 0.084602495  2.72001467 3.173354e-01
# 53   week_integer4 - week_integer15  4.051831e-02 0.013814422 Inf -6.810467e-03 0.087847093  2.93304447 1.986294e-01
# 54   week_integer4 - week_integer16  3.335451e-02 0.014282907 Inf -1.557932e-02 0.082288341  2.33527482 5.970417e-01
# 55    week_integer5 - week_integer6 -3.036839e-02 0.021044512 Inf -1.024678e-01 0.041730979 -1.44305507 9.881696e-01
# 56    week_integer5 - week_integer7  3.079535e-03 0.014900050 Inf -4.796865e-02 0.054127721  0.20667952 1.000000e+00
# 57    week_integer5 - week_integer8 -6.319774e-04 0.014765065 Inf -5.121770e-02 0.049953746 -0.04280221 1.000000e+00
# 58    week_integer5 - week_integer9  1.838666e-02 0.013923874 Inf -2.931710e-02 0.066090430  1.32051347 9.952418e-01
# 59   week_integer5 - week_integer10  1.682295e-02 0.014111111 Inf -3.152230e-02 0.065168198  1.19217741 9.984544e-01
# 60   week_integer5 - week_integer11  2.660509e-02 0.014294755 Inf -2.236934e-02 0.075579509  1.86117819 8.934137e-01
# 61   week_integer5 - week_integer12  3.814866e-02 0.013982412 Inf -9.755660e-03 0.086052981  2.72833196 3.120653e-01
# 62   week_integer5 - week_integer13  4.004068e-02 0.013831751 Inf -7.347473e-03 0.087428831  2.89483798 2.173356e-01
# 63   week_integer5 - week_integer14  4.232831e-02 0.013881774 Inf -5.231224e-03 0.089887839  3.04920018 1.487527e-01
# 64   week_integer5 - week_integer15  4.540472e-02 0.013927736 Inf -2.312281e-03 0.093121717  3.26002153 8.309669e-02
# 65   week_integer5 - week_integer16  3.824092e-02 0.014410932 Inf -1.113153e-02 0.087613365  2.65360481 3.610686e-01
# 66    week_integer6 - week_integer7  3.344792e-02 0.018664739 Inf -3.049824e-02 0.097394093  1.79203815 9.197064e-01
# 67    week_integer6 - week_integer8  2.973641e-02 0.018475198 Inf -3.356038e-02 0.093033205  1.60953143 9.670395e-01
# 68    week_integer6 - week_integer9  4.875505e-02 0.018198695 Inf -1.359443e-02 0.111104532  2.67904118 3.439813e-01
# 69   week_integer6 - week_integer10  4.719134e-02 0.018309275 Inf -1.553700e-02 0.109919671  2.57745525 4.144427e-01
# 70   week_integer6 - week_integer11  5.697348e-02 0.018781184 Inf -7.371638e-03 0.121318589  3.03354014 1.548742e-01
# 71   week_integer6 - week_integer12  6.851705e-02 0.018460254 Inf  5.271455e-03 0.131762644  3.71159842 1.897904e-02
# 72   week_integer6 - week_integer13  7.040907e-02 0.018322080 Inf  7.636864e-03 0.133181272  3.84285342 1.170598e-02
# 73   week_integer6 - week_integer14  7.269670e-02 0.018352462 Inf  9.820402e-03 0.135572991  3.96114137 7.426944e-03
# 74   week_integer6 - week_integer15  7.577311e-02 0.018427903 Inf  1.263835e-02 0.138907867  4.11186807 4.053083e-03
# 75   week_integer6 - week_integer16  6.860931e-02 0.018757648 Inf  4.344828e-03 0.132873785  3.65767108 2.299091e-02
# 76    week_integer7 - week_integer8 -3.711513e-03 0.008787073 Inf -3.381639e-02 0.026393362 -0.42238327 1.000000e+00
# 77    week_integer7 - week_integer9  1.530713e-02 0.007197634 Inf -9.352265e-03 0.039966521  2.12668871 7.480204e-01
# 78   week_integer7 - week_integer10  1.374341e-02 0.007527859 Inf -1.204735e-02 0.039534171  1.82567338 9.075281e-01
# 79   week_integer7 - week_integer11  2.352555e-02 0.007318902 Inf -1.549309e-03 0.048600411  3.21435542 9.484377e-02
# 80   week_integer7 - week_integer12  3.506913e-02 0.007358767 Inf  9.857686e-03 0.060280565  4.76562528 2.140593e-04
# 81   week_integer7 - week_integer13  3.696114e-02 0.007475783 Inf  1.134880e-02 0.062573487  4.94411648 8.820463e-05
# 82   week_integer7 - week_integer14  3.924877e-02 0.007476413 Inf  1.363427e-02 0.064863273  5.24967926 1.786624e-05
# 83   week_integer7 - week_integer15  4.232518e-02 0.007519236 Inf  1.656397e-02 0.068086395  5.62892107 2.151783e-06
# 84   week_integer7 - week_integer16  3.516138e-02 0.008229815 Inf  6.965694e-03 0.063357070  4.27243874 2.061514e-03
# 85    week_integer8 - week_integer9  1.901864e-02 0.006579752 Inf -3.523863e-03 0.041561143  2.89047973 2.195425e-01
# 86   week_integer8 - week_integer10  1.745493e-02 0.006756915 Inf -5.694544e-03 0.040604395  2.58326852 4.102621e-01
# 87   week_integer8 - week_integer11  2.723706e-02 0.007059626 Inf  3.050494e-03 0.051423634  3.85814557 1.104863e-02
# 88   week_integer8 - week_integer12  3.878064e-02 0.006878936 Inf  1.521312e-02 0.062348157  5.63759257 2.046607e-06
# 89   week_integer8 - week_integer13  4.067266e-02 0.006921228 Inf  1.696024e-02 0.064385069  5.87650876 4.993911e-07
# 90   week_integer8 - week_integer14  4.296028e-02 0.006995989 Inf  1.899174e-02 0.066928831  6.14070237 9.820341e-08
# 91   week_integer8 - week_integer15  4.603670e-02 0.007095950 Inf  2.172568e-02 0.070347714  6.48774233 1.043687e-08
# 92   week_integer8 - week_integer16  3.887289e-02 0.007731326 Inf  1.238505e-02 0.065360737  5.02797254 5.747158e-05
# 93   week_integer9 - week_integer10 -1.563715e-03 0.003960696 Inf -1.513322e-02 0.012005793 -0.39480814 1.000000e+00
# 94   week_integer9 - week_integer11  8.218424e-03 0.003846609 Inf -4.960218e-03 0.021397065  2.13653735 7.414284e-01
# 95   week_integer9 - week_integer12  1.976200e-02 0.003336382 Inf  8.331414e-03 0.031192581  5.92318163 3.765846e-07
# 96   week_integer9 - week_integer13  2.165402e-02 0.003478328 Inf  9.737122e-03 0.033570910  6.22541028 5.744800e-08
# 97   week_integer9 - week_integer14  2.394164e-02 0.003600159 Inf  1.160735e-02 0.036275937  6.65016395 3.508854e-09
# 98   week_integer9 - week_integer15  2.701806e-02 0.003790872 Inf  1.403037e-02 0.040005741  7.12713429 1.230240e-10
# 99   week_integer9 - week_integer16  1.985425e-02 0.004989169 Inf  2.761155e-03 0.036947355  3.97947110 6.910175e-03
# 100 week_integer10 - week_integer11  9.782139e-03 0.004042123 Inf -4.066341e-03 0.023630618  2.42004990 5.321688e-01
# 101 week_integer10 - week_integer12  2.132571e-02 0.003700702 Inf  8.646955e-03 0.034004470  5.76261304 9.853429e-07
# 102 week_integer10 - week_integer13  2.321773e-02 0.003837744 Inf  1.006946e-02 0.036365999  6.04983914 1.731633e-07
# 103 week_integer10 - week_integer14  2.550536e-02 0.003889697 Inf  1.217910e-02 0.038831621  6.55715894 6.570928e-09
# 104 week_integer10 - week_integer15  2.858177e-02 0.003994992 Inf  1.489476e-02 0.042268777  7.15440072 1.009069e-10
# 105 week_integer10 - week_integer16  2.141797e-02 0.005054120 Inf  4.102344e-03 0.038733595  4.23772447 2.392256e-03
# 106 week_integer11 - week_integer12  1.154357e-02 0.003219568 Inf  5.131998e-04 0.022573948  3.58544158 2.953593e-02
# 107 week_integer11 - week_integer13  1.343559e-02 0.003496116 Inf  1.457754e-03 0.025413430  3.84300534 1.169928e-02
# 108 week_integer11 - week_integer14  1.572322e-02 0.003561392 Inf  3.521744e-03 0.027924698  4.41490856 1.102762e-03
# 109 week_integer11 - week_integer15  1.879963e-02 0.003548021 Inf  6.643965e-03 0.030955299  5.29862452 1.370900e-05
# 110 week_integer11 - week_integer16  1.163583e-02 0.004819192 Inf -4.874921e-03 0.028146583  2.41447744 5.364340e-01
# 111 week_integer12 - week_integer13  1.892018e-03 0.002226623 Inf -5.736485e-03 0.009520522  0.84972543 9.999759e-01
# 112 week_integer12 - week_integer14  4.179647e-03 0.002216343 Inf -3.413635e-03 0.011772929  1.88583045 8.828453e-01
# 113 week_integer12 - week_integer15  7.256058e-03 0.002187287 Inf -2.376794e-04 0.014749795  3.31737736 7.005578e-02
# 114 week_integer12 - week_integer16  9.225702e-05 0.004101708 Inf -1.396036e-02 0.014144878  0.02249234 1.000000e+00
# 115 week_integer13 - week_integer14  2.287628e-03 0.002234628 Inf -5.368300e-03 0.009943557  1.02371767 9.997446e-01
# 116 week_integer13 - week_integer15  5.364040e-03 0.002210564 Inf -2.209444e-03 0.012937523  2.42654791 5.271984e-01
# 117 week_integer13 - week_integer16 -1.799761e-03 0.004174012 Inf -1.610010e-02 0.012500578 -0.43118258 1.000000e+00
# 118 week_integer14 - week_integer15  3.076411e-03 0.002065515 Inf -4.000129e-03 0.010152952  1.48941588 9.838983e-01
# 119 week_integer14 - week_integer16 -4.087390e-03 0.004081972 Inf -1.807239e-02 0.009897615 -1.00132730 9.998055e-01
# 120 week_integer15 - week_integer16 -7.163801e-03 0.004076080 Inf -2.112862e-02 0.006801016 -1.75752238 9.310218e-01


#week
MyData <- crossing(age = "A",
                         sex = "F",
                         study_year = "Second Winter",
                         week_integer = levels(bat_metrics$week_integer),
                         weight_std = mean(bat_metrics$weight_std),
                         av_wind_speed_std = mean(bat_metrics$av_wind_speed_std),
                         daily_rain_std = mean(bat_metrics$daily_rain_std),
                         av_air_temp_std = mean(bat_metrics$av_air_temp_std),
                         av_rel_humidity_std = mean(bat_metrics$av_rel_humidity_std),
                         av_cloud_okta_std = mean(bat_metrics$av_cloud_okta_std),
                         meanMoonlightIntensity = mean(bat_metrics$meanMoonlightIntensity))

#' #' To compute population-level predictions for a given grouping
MyData$ring_no <- NA
MyData$ntrials <- round(mean(bat_metrics$inactive + bat_metrics$active))

#' Calculate the predicted  values.
Fit_mod <- predict(model_3,
                     newdata = MyData,
                     re.form = NULL,  #~0, Keep it to NULL..gives odd results with ~0
                     type = "disp")

attr(Fit_mod, "names") <- NULL

MyData$phi   <- Fit_mod

# calculate variance ± SD
m3_emm_df <- m3_emm_df %>%
  left_join(MyData %>% select(week_integer, phi, ntrials), by = "week_integer")

# calculate variance and convert to proportion
m3_emm_df <- m3_emm_df %>% 
  mutate(
    variance_counts = emmean * (1-emmean) * (ntrials * ((phi + ntrials) / (phi + 1))),
    variance_props = variance_counts / (ntrials^2),
    sd_prop = sqrt(variance_props),
    SD_plus = emmean + 2*sd_prop,
    SD_minus = emmean - 2*sd_prop
  )


#Change names
m3_emm_df$week_integer <- factor(m3_emm_df$week_integer , levels=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16"))
m3_emm_df$week_integer <- plyr::revalue(m3_emm_df$week_integer, c("1" = "38", "2" = "39", "3" = "40", "4" = "41", "5" = "42", "6" = "43", "7" = "44", "8" = "45", "9" = "46", "10" = "47", "11" = "48", "12" = "49", "13" = "50", "14" = "51", "15" = "52", "16" = "1"))

bat_metrics$week_integer <- factor(bat_metrics$week_integer , levels=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16"))
bat_metrics$week_integer_rename <- plyr::revalue(bat_metrics$week_integer, c("1" = "38", "2" = "39", "3" = "40", "4" = "41", "5" = "42", "6" = "43", "7" = "44", "8" = "45", "9" = "46", "10" = "47", "11" = "48", "12" = "49", "13" = "50", "14" = "51", "15" = "52", "16" = "1"))

#Make contrasting groups match variable names
m3_emm_contrast_df$group1 <- NA
m3_emm_contrast_df$group2 <- NA
m3_emm_contrast_df$contrast_pretty <- NA

m3_emm_contrast_df$p_rounded <- rstatix::p_round(m3_emm_contrast_df$p.value, digits = 2)
m3_emm_contrast_df$p_pretty <- rstatix::p_format(rstatix::p_round(m3_emm_contrast_df$p.value, digits = 2),
                                                 accuracy = 1e-04, add.p = TRUE)

#remove non-significant differences
m3_emm_contrast_df <- filter(m3_emm_contrast_df, p_rounded <= 0.05)


m3_emm_contrast_df$group1 <- c("38", "38", "38", "38", "39", "39", "39", "39", "39", "39", "39", "40", "40", "40", "40", "40", "43", "43", "43", "43", "43", "44", "44", "44", "44", "44", "45", "45", "45", "45", "45", "45", "46", "46", "46", "46", "46", "47", "47", "47", "47", "47", "48", "48", "48", "48")
m3_emm_contrast_df$group2 <- c("49", "50", "51", "52", "46", "48", "49", "50", "51", "52", "1", "49", "50", "51", "52", "1", "49", "50", "51", "52", "1", "49", "50", "51", "52", "1", "48", "49", "50", "51", "52", "1", "49", "50", "51", "52", "1", "49", "50", "51", "52", "1", "49", "50", "51", "52")
m3_emm_contrast_df$contrast_pretty <- c("38 - 49", "38 - 50", "38 - 51", "38 - 52", "39 - 46", "39 - 48", "39 - 49", "39 - 50", "39 - 51", "39 - 52", "39 - 1", "40 - 49", "40 - 50", "40 - 51", "40 - 52", "40 - 1", "43 - 49", "43 - 50", "43 - 51", "43 - 52", "43 - 1", "44 - 49", "44 - 50", "44 - 51", "44 - 52", "44 - 1", "45 - 48", "45 - 49", "45 - 50", "45 - 51", "45 - 52", "45 - 1", "46 - 49", "46 - 50", "46 - 51", "46 - 52", "46 - 1", "47 - 49", "47 - 50", "47 - 51", "47 - 52", "47 - 1", "48 - 49", "48 - 50", "48 - 51", "48 - 52")


#plot the result
bat_plot <- ggplot(bat_metrics, aes(x = week_integer_rename, y = absolute_activity_proportion)) +
  geom_jitter(alpha = 0.3,
              height = 0.001,
              width = 0.3) +
  geom_point(data = m3_emm_df,
             aes(y = emmean,
                 x = week_integer), 
             size = 4,
             colour = "red") +
  geom_errorbar(data = m3_emm_df, 
                aes(y = emmean,
                    x = week_integer,
                    ymin = asymp.LCL, 
                    ymax = asymp.UCL), 
                width = 0.05,
                colour = "red") +
    geom_segment(data = m3_emm_df,
             aes(x = as.numeric(week_integer) -  0.3,
                 xend = as.numeric(week_integer) +  0.3,
                 y = SD_minus, yend = SD_minus),
             colour = "blue",
             linetype = "dashed",
             linewidth = 0.8,
             inherit.aes = FALSE) +
geom_segment(data = m3_emm_df,
             aes(x = as.numeric(week_integer) -  0.3,
                 xend = as.numeric(week_integer) +  0.3,
                 y = SD_plus, yend = SD_plus),
             colour = "blue",
             linetype = "dashed",
             linewidth = 0.8,
             inherit.aes = FALSE) +
  ylim(c(0, (max(bat_metrics$absolute_activity_proportion)+0.05))) +
  labs(y = "Proportion of Activity",
       x = "Week in Year") +
  theme_pubr()
bat_plot

bat_effects <- ggplot(data = m3_emm_contrast_df, aes(x = estimate, y = contrast_pretty)) +
  geom_point(size = 4) +
  xlab("Effects (Emergence)") +
  geom_errorbar(aes(x = estimate,
                    xmin = asymp.LCL,
                    xmax = asymp.UCL),
                width = 0.1) +
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "blue") +
  annotate(geom = "text",
           x = 0.1,
           y = 40,
           label = "p <= 0.05",
           size = 5) +
  theme_pubr() +
  rremove("ylab")
bat_effects

ggsave(filename = paste0("figs/bat_daily_activity_model_no_hibernation/activity_model_trapping_week.png"), plot = bat_plot, width = 10, height = 6, dpi = 300)

ggsave(filename = paste0("figs/bat_daily_activity_model_no_hibernation/activity_model_trapping_week_effect.png"), plot = bat_effects, width = 10, height = 6, dpi = 300)
