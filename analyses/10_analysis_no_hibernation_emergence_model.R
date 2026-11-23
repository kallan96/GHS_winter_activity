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


#check for NA values
colSums(is.na(bat_metrics))

#drop na values from moon columns
bat_metrics <- drop_na(bat_metrics, meanMoonlightIntensity)


# define factors
bat_metrics <- mutate(arrange(bat_metrics, time),
                      time = numFactor(time),
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


#remove bats that had little data
bat_metrics <- filter(bat_metrics, !bat_name %in% c('Gerty'))

#remove bats tech error
bat_metrics <- filter(bat_metrics, !bat_name %in% c('Conchita', 'Jennet'))

#remove spring records
bat_metrics <- filter(bat_metrics, survey_type != "hibernation")

# remove weeks with only data from Vincent
bat_metrics <- filter(bat_metrics, week_integer != 17 & week_integer != 18 & week_integer != 19)

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
model_1 <- glmmTMB(bat_emerged ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + weight_std + sex + age + study_year + week_integer + sex*age + weight_std*av_air_temp_std + av_wind_speed_std*av_air_temp_std  + (1 | ring_no),
                           data = bat_metrics,
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

plotResiduals(model_1_simulationOutput, form = bat_metrics$av_wind_speed_std, xlab = "Rank-transformed Standardised Wind Speed")
plotResiduals(model_1_simulationOutput, form = bat_metrics$av_air_temp_std, xlab = "Rank-transformed Standardised Air Temperature")
plotResiduals(model_1_simulationOutput, form = bat_metrics$meanMoonlightIntensity, xlab = "Rank-transformed Standardised Moonlight Intensity") #gotten worse
plotResiduals(model_1_simulationOutput, form = bat_metrics$av_cloud_okta_std, xlab = "Rank-transformed Standardised Cloud Cover") 
plotResiduals(model_1_simulationOutput, form = bat_metrics$daily_rain_std, xlab = "Rank-transformed Standardised Rainfall") # slight deviation in 0.25 quartile
plotResiduals(model_1_simulationOutput, form = bat_metrics$av_rel_humidity_std, xlab = "Rank-transformed Standardised Humidity")
plotResiduals(model_1_simulationOutput, form = bat_metrics$weight_std, xlab = "Rank-transformed Standardised Weight")
plotResiduals(model_1_simulationOutput, form = bat_metrics$age, xlab = "Rank-transformed Age")
plotResiduals(model_1_simulationOutput, form = bat_metrics$sex, xlab = "Rank-transformed Sex")
plotResiduals(model_1_simulationOutput, form = bat_metrics$week_integer, xlab = "Rank-transformed Week")
plotResiduals(model_1_simulationOutput, form = bat_metrics$study_year, xlab = "Rank-transformed Study Year")


#make simpler mdoels without interactions
# sex*age
model_1.c <- glmmTMB(bat_emerged ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + weight_std + sex + age + study_year + week_integer + weight_std*av_air_temp_std + av_wind_speed_std*av_air_temp_std + (1 | ring_no),
                           data = bat_metrics,
                           family = binomial(link = "logit"))

# weight_std*av_air_temp_std
model_1.d <- glmmTMB(bat_emerged ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + weight_std + sex + age + study_year + week_integer + sex*age + av_wind_speed_std*av_air_temp_std + (1 | ring_no),
                           data = bat_metrics,
                           family = binomial(link = "logit"))

# av_wind_speed_std*av_air_temp_std
model_1.e <- glmmTMB(bat_emerged ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + weight_std + sex + age + study_year + week_integer + sex*age + weight_std*av_air_temp_std + (1 | ring_no),
                           data = bat_metrics,
                           family = binomial(link = "logit"))


# test interactions
anova(model_1, model_1.c) #sex*age
#           Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# model_1.c 30 1647.8 1804.7 -793.91   1587.8
# model_1   32 1648.6 1815.9 -792.31   1584.6 3.1976      2     0.2021

anova(model_1, model_1.d) #weight_std*av_air_temp_std
#           Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# model_1.d 31 1656.7 1818.7 -797.34   1594.7
# model_1   32 1648.6 1815.9 -792.31   1584.6 10.046      1   0.001527 **

anova(model_1, model_1.e) #av_wind_speed_std*av_air_temp_std
#           Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# model_1.e 31 1647.9 1810.0 -792.95   1585.9
# model_1   32 1648.6 1815.9 -792.31   1584.6 1.2669      1     0.2603



#drop1 model selection ####
print("Dropping")
#model selection
drop1(model_1, test = "Chi")

#                                   Df    AIC    LRT  Pr(>Chi)
# <none>                               1648.6
# meanMoonlightIntensity             1 1646.7  0.025  0.875469    
# av_cloud_okta_std                  1 1648.2  1.626  0.202232
# daily_rain_std                     1 1647.3  0.650  0.419986
# av_rel_humidity_std                1 1649.5  2.914  0.087799 .
# study_year                         1 1646.9  0.251  0.616638
# week_integer                      15 1692.4 73.766 9.451e-10 ***
# sex:age                            2 1647.8  3.198  0.202140
# av_air_temp_std:weight_std         1 1656.7 10.046  0.001527 **
# av_air_temp_std:av_wind_speed_std  1 1647.9  1.267  0.260345  


# drop sex:age interaction
drop1(model_1.c, test = "Chi")

#                                   Df    AIC    LRT  Pr(>Chi)
# <none>                               1647.8
# meanMoonlightIntensity             1 1645.8  0.016  0.900143
# av_cloud_okta_std                  1 1647.4  1.571  0.210111
# daily_rain_std                     1 1646.5  0.659  0.416985
# av_rel_humidity_std                1 1648.8  3.020  0.082266 .
# sex                                1 1646.5  0.680  0.409603
# age                                2 1644.9  1.044  0.593284
# study_year                         1 1646.3  0.490  0.484059    
# week_integer                      15 1693.1 75.241 5.121e-10 ***
# av_air_temp_std:weight_std         1 1655.8  9.931  0.001625 **
# av_air_temp_std:av_wind_speed_std  1 1647.0  1.227  0.268008


# drop av_air_temp_std:av_wind_speed_std interaction
model_2 <- glmmTMB(bat_emerged ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + weight_std + sex + age + study_year + week_integer + weight_std*av_air_temp_std + (1 | ring_no),
                           data = bat_metrics,
                           family = binomial(link = "logit"))

#model selection               
drop1(model_2, test = "Chi")

#                            Df    AIC    LRT  Pr(>Chi)
# <none>                        1647.0
# av_wind_speed_std           1 1694.1 49.049 2.497e-12 ***
# meanMoonlightIntensity      1 1645.1  0.019  0.889060
# av_cloud_okta_std           1 1646.5  1.445  0.229292
# daily_rain_std              1 1645.9  0.865  0.352476    
# av_rel_humidity_std         1 1648.9  3.881  0.048822 *
# sex                         1 1645.7  0.636  0.425274
# age                         2 1644.1  1.067  0.586469    
# study_year                  1 1645.5  0.463  0.496180
# week_integer               15 1691.1 74.015 8.524e-10 ***
# av_air_temp_std:weight_std  1 1654.6  9.579  0.001969 **



#check model fits ####
summary(model_2)
#  Family: binomial  ( logit )
# Formula:
# bat_emerged ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity +  
#     av_cloud_okta_std + daily_rain_std + av_rel_humidity_std +
#     weight_std + sex + age + study_year + week_integer + weight_std *
#     av_air_temp_std + (1 | ring_no)
# Data: bat_metrics

#       AIC       BIC    logLik -2*log(L)  df.resid
#    1647.1    1798.7    -794.5    1589.1      1348

# Random effects:

# Conditional model:
#  Groups  Name        Variance Std.Dev.
#  ring_no (Intercept) 0.5986   0.7737
# Number of obs: 1377, groups:  ring_no, 36

# Conditional model:
#                             Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                 5.637198   1.783597   3.161  0.00157 **
# av_air_temp_std             1.284204   0.149940   8.565  < 2e-16 ***
# av_wind_speed_std          -4.251779   0.627566  -6.775 1.24e-11 ***
# meanMoonlightIntensity     -0.118644   0.850458  -0.140  0.88905    
# av_cloud_okta_std          -0.409515   0.341225  -1.200  0.23009    
# daily_rain_std              0.429854   0.463826   0.927  0.35405
# av_rel_humidity_std        -3.335285   1.698663  -1.963  0.04959 *
# weight_std                 -0.103634   0.222856  -0.465  0.64191    
# sexM                        0.301323   0.380953   0.791  0.42896
# ageJ                        0.190830   0.473700   0.403  0.68706
# ageYA                      -0.334442   0.478821  -0.698  0.48488    
# study_yearSecond Winter    -0.203907   0.298773  -0.682  0.49493
# week_integer2              -2.198968   0.533723  -4.120 3.79e-05 ***
# week_integer3              -1.606773   0.564632  -2.846  0.00443 **
# week_integer4              -1.865991   0.700938  -2.662  0.00776 ** 
# week_integer5              -0.866717   0.913738  -0.949  0.34285
# week_integer6              -0.671920   0.710113  -0.946  0.34404    
# week_integer7               0.111968   0.665218   0.168  0.86633
# week_integer8              -0.007067   0.610089  -0.012  0.99076
# week_integer9              -0.440716   0.587727  -0.750  0.45334    
# week_integer10             -0.231966   0.596964  -0.389  0.69759
# week_integer11             -0.566828   0.648622  -0.874  0.38218
# week_integer12              0.071797   0.628208   0.114  0.90901    
# week_integer13             -1.294168   0.614846  -2.105  0.03530 *
# week_integer14             -0.742322   0.621614  -1.194  0.23241    
# week_integer15             -1.555668   0.659391  -2.359  0.01831 *
# week_integer16             -0.812534   0.695398  -1.168  0.24263    
# av_air_temp_std:weight_std  0.249591   0.081372   3.067  0.00216 **
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


check_collinearity(model_2)


#DHARMa simulation
model_2_simulationOutput <- simulateResiduals(fittedModel = model_2, plot = F, n = 1000)
plotQQunif(model_2_simulationOutput, testUniformity = TRUE, 
           testOutliers = TRUE, testDispersion = TRUE)
plotResiduals(model_2_simulationOutput, quantreg = TRUE, smoothScatter = FALSE) # deviations in 0.25 quantile
testDispersion(model_2_simulationOutput)
testOutliers(model_2_simulationOutput)
testZeroInflation(model_2_simulationOutput)
mod_3_res = recalculateResiduals(model_2_simulationOutput, group = bat_metrics$time)
testTemporalAutocorrelation(mod_3_res, time = unique(bat_metrics$time))

plotResiduals(model_2_simulationOutput, form = bat_metrics$av_wind_speed_std, xlab = "Rank-transformed Standardised Wind Speed")
plotResiduals(model_2_simulationOutput, form = bat_metrics$av_air_temp_std, xlab = "Rank-transformed Standardised Air Temperature")
plotResiduals(model_2_simulationOutput, form = bat_metrics$meanMoonlightIntensity, xlab = "Rank-transformed Standardised Moonlight Intensity") #gotten worse
plotResiduals(model_2_simulationOutput, form = bat_metrics$av_cloud_okta_std, xlab = "Rank-transformed Standardised Cloud Cover") 
plotResiduals(model_2_simulationOutput, form = bat_metrics$daily_rain_std, xlab = "Rank-transformed Standardised Rainfall") # slight deviation in 0.25 quartile
plotResiduals(model_2_simulationOutput, form = bat_metrics$av_rel_humidity_std, xlab = "Rank-transformed Standardised Humidity")
plotResiduals(model_2_simulationOutput, form = bat_metrics$weight_std, xlab = "Rank-transformed Standardised Weight")
plotResiduals(model_2_simulationOutput, form = bat_metrics$age, xlab = "Rank-transformed Age")
plotResiduals(model_2_simulationOutput, form = bat_metrics$sex, xlab = "Rank-transformed Sex")
plotResiduals(model_2_simulationOutput, form = bat_metrics$study_year, xlab = "Rank-transformed Study Year")



#Test significance of main effects ####
model_2.a <- glmmTMB(bat_emerged ~ av_air_temp_std + av_wind_speed_std + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + weight_std + sex + age + study_year + week_integer + weight_std*av_air_temp_std + (1 | ring_no),
                           data = bat_metrics,
                           family = binomial(link = "logit"))

model_2.b <- glmmTMB(bat_emerged ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + daily_rain_std + av_rel_humidity_std + weight_std + sex + age + study_year + week_integer + weight_std*av_air_temp_std + (1 | ring_no),
                           data = bat_metrics,
                           family = binomial(link = "logit"))

model_2.c <- glmmTMB(bat_emerged ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + av_rel_humidity_std + weight_std + sex + age + study_year + week_integer + weight_std*av_air_temp_std + (1 | ring_no),
                           data = bat_metrics,
                           family = binomial(link = "logit"))

model_2.d <- glmmTMB(bat_emerged ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + weight_std + sex + age + study_year + week_integer + weight_std*av_air_temp_std + (1 | ring_no),
                           data = bat_metrics,
                           family = binomial(link = "logit"))

model_2.e <- glmmTMB(bat_emerged ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + weight_std + age + study_year + week_integer + weight_std*av_air_temp_std + (1 | ring_no),
                           data = bat_metrics,
                           family = binomial(link = "logit"))

model_2.f <- glmmTMB(bat_emerged ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + weight_std + sex + study_year + week_integer + weight_std*av_air_temp_std + (1 | ring_no),
                           data = bat_metrics,
                           family = binomial(link = "logit"))

model_2.g <- glmmTMB(bat_emerged ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + weight_std + sex + age + week_integer + weight_std*av_air_temp_std + (1 | ring_no),
                           data = bat_metrics,
                           family = binomial(link = "logit"))

model_2.h <- glmmTMB(bat_emerged ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + weight_std + sex + age + study_year + weight_std*av_air_temp_std + (1 | ring_no),
                           data = bat_metrics,
                           family = binomial(link = "logit"))
                    
model_2.i <- glmmTMB(bat_emerged ~ av_air_temp_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + weight_std + sex + age + study_year + week_integer + weight_std*av_air_temp_std + (1 | ring_no),
                           data = bat_metrics,
                           family = binomial(link = "logit"))


anova(model_2, model_2.a) #moon
#           Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# model_2.a 28 1645.1 1791.5 -794.54   1589.1
# model_2   29 1647.0 1798.7 -794.53   1589.0 0.0195      1     0.8891

anova(model_2, model_2.b) #cloud
#           Df    AIC    BIC  logLik deviance Chisq Chi Df Pr(>Chisq)
# model_2.b 28 1646.5 1792.9 -795.25   1590.5
# model_2   29 1647.0 1798.7 -794.53   1589.0 1.4452      1     0.2293

anova(model_2, model_2.c) #rain
#           Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# model_2.c 28 1645.9 1792.3 -794.96   1589.9
# model_2   29 1647.0 1798.7 -794.53   1589.0 0.8645      1     0.3525

anova(model_2, model_2.d) #humidity
#           Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# model_2.d 28 1648.9 1795.3 -796.47   1592.9
# model_2   29 1647.0 1798.7 -794.53   1589.0 3.8814      1    0.04882 *

anova(model_2, model_2.e) #sex
#           Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# model_2.e 28 1645.7 1792.1 -794.84   1589.7
# model_2   29 1647.0 1798.7 -794.53   1589.0 0.6357      1     0.4253

anova(model_2, model_2.f) #age
#           Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# model_2.f 27 1644.1 1785.3 -795.06   1590.1
# model_2   29 1647.0 1798.7 -794.53   1589.0 1.0673      2     0.5865

anova(model_2, model_2.g) #year
#           Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# model_2.g 28 1645.5 1791.9 -794.76   1589.5
# model_2   29 1647.0 1798.7 -794.53   1589.0 0.4631      1     0.4962

anova(model_2, model_2.h) #week
#           Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# model_2.h 14 1691.1 1764.3 -831.53   1663.1
# model_2   29 1647.0 1798.7 -794.53   1589.0 74.015     15  8.524e-10 ***

anova(model_2, model_2.i) #wind
#           Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# model_2.i 28 1694.1 1840.5 -819.05   1638.1
# model_2   29 1647.0 1798.7 -794.53   1589.0 49.049      1  2.497e-12 ***



# plotting - model_2 #####
#weight_std * av_air_temp_std
MyData_weight_temp <- crossing(age = "A",
                               sex = "F",
                               study_year = "Second Winter",
                               week_integer = 14,
                               weight_std = seq(from = min(bat_metrics$weight_std),
                                                to = max(bat_metrics$weight_std),
                                                length = 20),
                               av_wind_speed_std = mean(bat_metrics$av_wind_speed_std),
                               daily_rain_std = mean(bat_metrics$daily_rain_std),
                               av_air_temp_std = seq(from = min(bat_metrics$av_air_temp_std),
                                                     to = max(bat_metrics$av_air_temp_std),
                                                     length = 100),
                               av_rel_humidity_std = mean(bat_metrics$av_rel_humidity_std),
                               av_cloud_okta_std = mean(bat_metrics$av_cloud_okta_std),
                               meanMoonlightIntensity = mean(bat_metrics$meanMoonlightIntensity))

#' #' To compute population-level predictions for a given grouping
MyData_weight_temp$ring_no <- NA
MyData_weight_temp$time <- NA

#' Calculate the predicted  values.
Fit_weight_temp <- predict(model_2,
                           newdata = MyData_weight_temp,
                           re.form = NULL,  #~0, Keep it to NULL..gives odd results with ~0
                           type = "link",
                           se = TRUE)

#' Extract the predicted values and standard errors, and calculate
#' standard errors (SE) for predicted values.
MyData_weight_temp$Pi   <- exp(Fit_weight_temp$fit) / (1 +  exp(Fit_weight_temp$fit))
MyData_weight_temp$SeUp <- exp(Fit_weight_temp$fit + 1.96 * Fit_weight_temp$se.fit) / (1 + exp(Fit_weight_temp$fit + 1.96 * Fit_weight_temp$se.fit))
MyData_weight_temp$SeLo <- exp(Fit_weight_temp$fit - 1.96 * Fit_weight_temp$se.fit) / (1 +  exp(Fit_weight_temp$fit - 1.96 * Fit_weight_temp$se.fit))
MyData_weight_temp$se   <- Fit_weight_temp$se.fit
head(MyData_weight_temp)

#back standardise
MyData_weight_temp$av_air_temp <- MyData_weight_temp$av_air_temp_std * sd(bat_metrics$av_air_temp) +
  mean(bat_metrics$av_air_temp)
MyData_weight_temp$weight <- MyData_weight_temp$weight_std  * sd(bat_metrics$weight) +
  mean(bat_metrics$weight)


#plot the result
bat_plot <- ggplot(MyData_weight_temp, aes(x = av_air_temp, y = Pi)) +
  geom_ribbon(aes(ymax = SeUp,
                  ymin = SeLo,
                  fill = as.factor(round(weight, 1))),
              alpha = 0.03,
              show.legend = FALSE) +
  geom_line(aes(colour = as.factor(round(weight, 1))),
            linewidth = 1) +
  geom_jitter(data = bat_metrics,
              aes(y = bat_emerged),
              alpha = 0.3,
              height = 0.01,
              width = 0.01) +
  labs(y = "Probability of Emergence",
       x = "Average Air Temperature (DegC)",
       colour = "Weight (g)") +
  scale_color_viridis(discrete=TRUE) +
  scale_fill_viridis(discrete=TRUE) +
  theme_bw()
  

ggsave(filename = paste0("figs/bat_daily_emergence_model_no_hibernation/emergence_model_trapping_temp_weight.png"), plot = bat_plot, width = 10, height = 6, dpi = 300)


#av_wind_speed_std
MyData_wind <- crossing(age = "A",
                             sex = "F",
                             study_year = "Second Winter",
                             week_integer = 14,
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
MyData_wind$time <- NA

#' Calculate the predicted  values.
Fit_wind <- predict(model_2,
                         newdata = MyData_wind,
                         re.form = NULL,  #~0, Keep it to NULL..gives odd results with ~0
                         type = "link",
                         se = TRUE)

#' Extract the predicted values and standard errors, and calculate
#' standard errors (SE) for predicted values.
MyData_wind$Pi   <- exp(Fit_wind$fit) / (1 +  exp(Fit_wind$fit))
MyData_wind$SeUp <- exp(Fit_wind$fit + 1.96 * Fit_wind$se.fit) / (1 + exp(Fit_wind$fit + 1.96 * Fit_wind$se.fit))
MyData_wind$SeLo <- exp(Fit_wind$fit - 1.96 * Fit_wind$se.fit) / (1 +  exp(Fit_wind$fit - 1.96 * Fit_wind$se.fit))
MyData_wind$se   <- Fit_wind$se.fit
head(MyData_wind)

#back standardise
MyData_wind$av_wind_speed <- MyData_wind$av_wind_speed_std * max(bat_metrics$av_wind_speed)


#plot the result
bat_plot <- ggplot(MyData_wind, aes(x = av_wind_speed, y = Pi)) +
  geom_ribbon(aes(ymax = SeUp,
                  ymin = SeLo),
              alpha = 0.03) +
  geom_line(linewidth = 1) +
  geom_jitter(data = bat_metrics,
              aes(y = bat_emerged),
              alpha = 0.3,
              height = 0.01,
              width = 0.01) +
  labs(y = "Probability of Emergence",
       x = "Average Wind \n Speed (kn)") +
  theme_bw()

ggsave(filename = paste0("figs/bat_daily_emergence_model_no_hibernation/emergence_model_trapping_wind.png"), plot = bat_plot, width = 10, height = 6, dpi = 300)


#cloud cover
MyData_cloud <- crossing(age = "A",
                         sex = "F",
                         study_year = "Second Winter",
                         week_integer = 14,
                         weight_std = mean(bat_metrics$weight_std),
                         av_wind_speed_std = mean(bat_metrics$av_wind_speed_std),
                         daily_rain_std = mean(bat_metrics$daily_rain_std),
                         av_air_temp_std = mean(bat_metrics$av_air_temp_std),
                         av_rel_humidity_std = mean(bat_metrics$av_rel_humidity_std),
                         av_cloud_okta_std = seq(from = min(bat_metrics$av_cloud_okta_std),
                                                 to = max(bat_metrics$av_cloud_okta_std),
                                                 length = 100),
                         meanMoonlightIntensity = mean(bat_metrics$meanMoonlightIntensity))

#To compute population-level predictions for a given grouping
MyData_cloud$ring_no <- NA
MyData_cloud$time <- NA

#Calculate the predicted  values.
Fit_cloud <- predict(model_2,
                     newdata = MyData_cloud,
                     re.form = NULL,  #~0, Keep it to NULL..gives odd results with ~0
                     type = "link",
                     se = TRUE)

#' Extract the predicted values and standard errors, and calculate
#' standard errors (SE) for predicted values.
MyData_cloud$Pi   <- exp(Fit_cloud$fit) / (1 +  exp(Fit_cloud$fit))
MyData_cloud$SeUp <- exp(Fit_cloud$fit + 1.96 * Fit_cloud$se.fit) / (1 + exp(Fit_cloud$fit + 1.96 * Fit_cloud$se.fit))
MyData_cloud$SeLo <- exp(Fit_cloud$fit - 1.96 * Fit_cloud$se.fit) / (1 +  exp(Fit_cloud$fit - 1.96 * Fit_cloud$se.fit))
MyData_cloud$se   <- Fit_cloud$se.fit
head(MyData_cloud)

#back standardise
MyData_cloud$av_cloud_okta <- MyData_cloud$av_cloud_okta_std * max(bat_metrics$av_cloud_okta)


#plot the result
bat_plot <- ggplot(MyData_cloud, aes(x = av_cloud_okta, y = Pi)) +
  geom_ribbon(aes(ymax = SeUp,
                  ymin = SeLo),
              alpha = 0.03,
              show.legend = FALSE) +
  geom_line(linewidth = 1) +
  geom_jitter(data = bat_metrics,
              aes(y = bat_emerged),
              alpha = 0.3,
              height = 0.01,
              width = 0.01) +
  scale_color_viridis(discrete=TRUE) +
  scale_fill_viridis(discrete=TRUE) +
  labs(y = "Probability of Emergence",
       x = "Average Cloud Cover (Oktas)") +
  theme_bw()

ggsave(filename = paste0("figs/bat_daily_emergence_model_no_hibernation/emergence_model_trapping_cloud.png"), plot = bat_plot, width = 10, height = 6, dpi = 300)


#meanMoonlightIntensity
MyData_moonlight <- crossing(age = "A",
                             sex = "F",
                             study_year = "Second Winter",
                             week_integer = 14,
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
MyData_moonlight$ring_no <- NA
MyData_moonlight$time <- NA

#' Calculate the predicted  values.
Fit_moonlight <- predict(model_2,
                         newdata = MyData_moonlight,
                         re.form = NULL,  #~0, Keep it to NULL..gives odd results with ~0
                         type = "link",
                         se = TRUE)

#' Extract the predicted values and standard errors, and calculate
#' standard errors (SE) for predicted values.
MyData_moonlight$Pi   <- exp(Fit_moonlight$fit) / (1 +  exp(Fit_moonlight$fit))
MyData_moonlight$SeUp <- exp(Fit_moonlight$fit + 1.96 * Fit_moonlight$se.fit) / (1 + exp(Fit_moonlight$fit + 1.96 * Fit_moonlight$se.fit))
MyData_moonlight$SeLo <- exp(Fit_moonlight$fit - 1.96 * Fit_moonlight$se.fit) / (1 +  exp(Fit_moonlight$fit - 1.96 * Fit_moonlight$se.fit))
MyData_moonlight$se   <- Fit_moonlight$se.fit
head(MyData_moonlight)

#back standardise
MyData_moonlight$meanMoonlightIntensity_lux <- MyData_moonlight$meanMoonlightIntensity * 0.3


#plot the result
bat_plot <- ggplot(MyData_moonlight, aes(x = meanMoonlightIntensity, y = Pi)) +
  geom_ribbon(aes(ymax = SeUp,
                  ymin = SeLo),
              alpha = 0.03,
              show.legend = FALSE) +
  geom_line(linewidth = 1) +
  geom_jitter(data = bat_metrics,
              aes(y = bat_emerged),
              alpha = 0.3,
              height = 0.01,
              width = 0.01) +
  labs(y = "Probability of Emergence",
       x = "Average Moonlight Intensity (lux)") +
  theme_bw()

ggsave(filename = paste0("figs/bat_daily_emergence_model_no_hibernation/emergence_model_trapping_moonlight.png"), plot = bat_plot, width = 10, height = 6, dpi = 300)


#daily rain
MyData_rain <- crossing(age = "A",
                        sex = "F",
                        study_year = "Second Winter",
                        week_integer = 14,
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
MyData_rain$time <- NA

#' Calculate the predicted  values.
Fit_rain <- predict(model_2,
                    newdata = MyData_rain,
                    re.form = NULL,  #~0, Keep it to NULL..gives odd results with ~0
                    type = "link",
                    se = TRUE)

#' Extract the predicted values and standard errors, and calculate
#' standard errors (SE) for predicted values.
MyData_rain$Pi   <- exp(Fit_rain$fit) / (1 +  exp(Fit_rain$fit))
MyData_rain$SeUp <- exp(Fit_rain$fit + 1.96 * Fit_rain$se.fit) / (1 + exp(Fit_rain$fit + 1.96 * Fit_rain$se.fit))
MyData_rain$SeLo <- exp(Fit_rain$fit - 1.96 * Fit_rain$se.fit) / (1 +  exp(Fit_rain$fit - 1.96 * Fit_rain$se.fit))
MyData_rain$se   <- Fit_rain$se.fit
head(MyData_rain)

#back standardise
MyData_rain$daily_rain <- MyData_rain$daily_rain_std * max(bat_metrics$daily_rain)


#plot the result
bat_plot <- ggplot(MyData_rain, aes(x = daily_rain, y = Pi)) +
  geom_ribbon(aes(ymax = SeUp,
                  ymin = SeLo),
              alpha = 0.03,
              show.legend = FALSE) +
  geom_line(linewidth = 1) +
  geom_jitter(data = bat_metrics,
              aes(y = bat_emerged),
              alpha = 0.3,
              height = 0.01,
              width = 0.01) +
  labs(y = "Probability of Emergence",
       x = "Total Rainfall (mm)") +
  theme_bw()

ggsave(filename = paste0("figs/bat_daily_emergence_model_no_hibernation/emergence_model_trapping_rain.png"), plot = bat_plot, width = 10, height = 6, dpi = 300)


#av_rel_humidity
MyData_humid <- crossing(age = "A",
                         sex = "F",
                         study_year = "Second Winter",
                         week_integer = 14,
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
MyData_humid$time <- NA

#' Calculate the predicted  values.
Fit_humid <- predict(model_2,
                     newdata = MyData_humid,
                     re.form = NULL,  #~0, Keep it to NULL..gives odd results with ~0
                     type = "link",
                     se = TRUE)

#' Extract the predicted values and standard errors, and calculate
#' standard errors (SE) for predicted values.
MyData_humid$Pi   <- exp(Fit_humid$fit) / (1 +  exp(Fit_humid$fit))
MyData_humid$SeUp <- exp(Fit_humid$fit + 1.96 * Fit_humid$se.fit) / (1 + exp(Fit_humid$fit + 1.96 * Fit_humid$se.fit))
MyData_humid$SeLo <- exp(Fit_humid$fit - 1.96 * Fit_humid$se.fit) / (1 +  exp(Fit_humid$fit - 1.96 * Fit_humid$se.fit))
MyData_humid$se   <- Fit_humid$se.fit
head(MyData_humid)

#back standardise
MyData_humid$av_rel_humidity <- MyData_humid$av_rel_humidity_std * max(bat_metrics$av_rel_humidity)


#plot the result
bat_plot <- ggplot(MyData_humid, aes(x = av_rel_humidity, y = Pi)) +
  geom_ribbon(aes(ymax = SeUp,
                  ymin = SeLo),
              alpha = 0.03,
              show.legend = FALSE) +
  geom_line(linewidth = 1) +
  geom_jitter(data = bat_metrics,
              aes(y = bat_emerged),
              alpha = 0.3,
              height = 0.01,
              width = 0.01) +
  labs(y = "Probability of Emergence",
       x = "Average Relative Humidity (%)") +
  theme_bw()

ggsave(filename = paste0("figs/bat_daily_emergence_model_no_hibernation/emergence_model_trapping_humidity.png"), plot = bat_plot, width = 10, height = 6, dpi = 300)



#eemeans for categorical variable ####
#age effect
m3_emm <- emmeans(model_2, ~ age, component = "cmean")

### pairwise comparisons
m3_contrast <- summary(contrast(m3_emm,
                                  method = "pairwise"),
                         infer = TRUE)

#Make into a dataframe
m3_emm_df <- data.frame(m3_emm)
m3_emm_contrast_df <- data.frame(m3_contrast)


# > m3_emm_df
#   age    emmean         SE  df asymp.LCL asymp.UCL
# 1   A 0.4854078 0.05501551 Inf 0.3775793 0.5932362
# 2   J 0.5277555 0.08304191 Inf 0.3649963 0.6905146
# 3  YA 0.4117047 0.08883957 Inf 0.2375824 0.5858271

# > m3_emm_contrast_df
#   contrast    estimate        SE  df  asymp.LCL asymp.UCL    z.ratio   p.value
# 1    A - J -0.04234770 0.1048060 Inf -0.2879815 0.2032861 -0.4040581 0.9139489
# 2   A - YA  0.07370302 0.1049652 Inf -0.1723040 0.3197100  0.7021663 0.7622015
# 3   J - YA  0.11605073 0.1129764 Inf -0.1487322 0.3808337  1.0272117 0.5596397

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
bat_plot <- ggplot(bat_metrics, aes(x = age, y = bat_emerged)) +
  geom_jitter(alpha = 0.3,
              height = 0.1,
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
  labs(y = "Probability of Emergence",
       x = "Age") +
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
           x = m3_emm_contrast_df$estimate,
           y = 1:3 + 0.2,
           label = m3_emm_contrast_df$p_pretty,
           size = 3) +
  theme_pubr() +
  rremove("ylab")
bat_effects

#make into one plot
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

ggsave(filename = paste0("figs/bat_daily_emergence_model_no_hibernation/emergence_model_trapping_age.png"), plot = bat_plot, width = 10, height = 6, dpi = 300)

ggsave(filename = paste0("figs/bat_daily_emergence_model_no_hibernation/emergence_model_trapping_age_effect.png"), plot = m3_plot, width = 10, height = 6, dpi = 300)


#sex effect
m3_emm <- emmeans(model_2, ~sex, component = "cmean")

### pairwise comparisons
m3_contrast <- summary(contrast(m3_emm,
                                method = "pairwise"),
                       infer = TRUE)

#Make into a dataframe
m3_emm_df <- data.frame(m3_emm)
m3_emm_contrast_df <- data.frame(m3_contrast)


# > m3_emm_df
#   sex    emmean         SE  df asymp.LCL asymp.UCL
# 1   F 0.4416991 0.06220995 Inf 0.3197698 0.5636283
# 2   M 0.5082129 0.06185893 Inf 0.3869716 0.6294542

# > m3_emm_contrast_df
#   contrast    estimate         SE  df  asymp.LCL  asymp.UCL    z.ratio   p.value
# 1    F - M -0.06651381 0.08420868 Inf -0.2315598 0.09853217 -0.7898688 0.4296044


#Change names
m3_emm_df$sex <- factor(m3_emm_df$sex , levels=c("F", "M"))
m3_emm_df$sex <- plyr::revalue(m3_emm_df$sex, c("F" = "Female", "M" = "Male"))

bat_metrics$sex <- factor(bat_metrics$sex , levels=c("F", "M"))
bat_metrics$sex <- plyr::revalue(bat_metrics$sex, c("F" = "Female", "M" = "Male"))


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
bat_plot <- ggplot(bat_metrics, aes(x = sex, y = bat_emerged)) +
  geom_jitter(alpha = 0.3,
              height = 0.1,
              width = 0.3) +
  geom_point(data = m3_emm_df,
             aes(y = emmean,
                 x = sex), 
             size = 4,
             colour = "red") +
  geom_errorbar(data = m3_emm_df, 
                aes(y = emmean,
                    x = sex,
                    ymin = asymp.LCL, 
                    ymax = asymp.UCL), 
                width = 0.05,
                colour = "red") +
  labs(y = "Probability of Emergence",
       x = "Sex") +
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
           x = m3_emm_contrast_df$estimate,
           y = 1 + 0.2,
           label = m3_emm_contrast_df$p_pretty,
           size = 3) +
  theme_pubr() +
  rremove("ylab")
bat_effects

#make into one plot
m5_bottom <- bat_plot +
  coord_flip() # rotate to horizontal

m5_top <- bat_effects + # assign to new plot object
  scale_x_continuous(position="top") # move x-axis to "top"

m3_plot <- cowplot::plot_grid(m5_top,
                              m5_bottom,
                              nrow = 2,
                              align = "v",
                              axis = "lr",
                              rel_heights = c(1, 1))
m3_plot

ggsave(filename = paste0("figs/bat_daily_emergence_model_no_hibernation/emergence_model_trapping_sex.png"), plot = bat_plot, width = 10, height = 6, dpi = 300)

ggsave(filename = paste0("figs/bat_daily_emergence_model_no_hibernation/emergence_model_trapping_sex_effect.png"), plot = m3_plot, width = 10, height = 6, dpi = 300)


#year effect
m3_emm <- emmeans(model_2, ~ study_year, component = "cmean")

### pairwise comparisons
m3_contrast <- summary(contrast(m3_emm,
                                method = "pairwise"),
                       infer = TRUE)

#Make into a dataframe
m3_emm_df <- data.frame(m3_emm)
m3_emm_contrast_df <- data.frame(m3_contrast)


# > m3_emm_df
#      study_year    emmean         SE  df asymp.LCL asymp.UCL
# 1  First Winter 0.4974285 0.06419083 Inf 0.3716168 0.6232402
# 2 Second Winter 0.4524835 0.04694337 Inf 0.3604761 0.5444908

# > m3_emm_contrast_df
#                       contrast   estimate         SE  df   asymp.LCL asymp.UCL   z.ratio   p.value
# 1 First Winter - Second Winter 0.04494503 0.06592554 Inf -0.08426666 0.1741567 0.6817543 0.4953943


#Change names
m3_emm_df$study_year <- factor(m3_emm_df$study_year , levels=c("First Winter", "Second Winter"))

bat_metrics$study_year <- factor(bat_metrics$study_year , levels=c("First Winter", "Second Winter"))


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
bat_plot <- ggplot(bat_metrics, aes(x = study_year, y = bat_emerged)) +
  geom_jitter(alpha = 0.3,
              height = 0.1,
              width = 0.3) +
  geom_point(data = m3_emm_df,
             aes(y = emmean,
                 x = study_year), 
             size = 4,
             colour = "red") +
  geom_errorbar(data = m3_emm_df, 
                aes(y = emmean,
                    x = study_year,
                    ymin = asymp.LCL, 
                    ymax = asymp.UCL), 
                width = 0.05,
                colour = "red") +
  labs(y = "Probability of Emergence",
       x = "Study Year") +
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

ggsave(filename = paste0("figs/bat_daily_emergence_model_no_hibernation/emergence_model_trapping_year.png"), plot = bat_plot, width = 10, height = 6, dpi = 300)

ggsave(filename = paste0("figs/bat_daily_emergence_model_no_hibernation/emergence_model_trapping_year_effect.png"), plot = m3_plot, width = 10, height = 6, dpi = 300)


#week effect
m3_emm <- emmeans(model_2, ~ week_integer, component = "cmean")

### pairwise comparisons
m3_contrast <- summary(contrast(m3_emm,
                                method = "pairwise"),
                       infer = TRUE)

#Make into a dataframe
m3_emm_df <- data.frame(m3_emm)
m3_emm_contrast_df <- data.frame(m3_contrast)


# > m3_emm_df
#    week_integer    emmean         SE  df  asymp.LCL asymp.UCL
# 1             1 0.6580810 0.12018687 Inf 0.42251909 0.8936429
# 2             2 0.1814505 0.06508787 Inf 0.05388067 0.3090204
# 3             3 0.2844692 0.09340798 Inf 0.10139292 0.4675455
# 4             4 0.2354506 0.10199385 Inf 0.03554635 0.4353549
# 5             5 0.4512886 0.19132344 Inf 0.07630151 0.8262756
# 6             6 0.4988900 0.12700967 Inf 0.24995560 0.7478244
# 7             7 0.6823868 0.08917055 Inf 0.50761571 0.8571579
# 8             8 0.6565175 0.07274022 Inf 0.51394924 0.7990857
# 9             9 0.5553403 0.06856651 Inf 0.42095242 0.6897282
# 10           10 0.6051663 0.06494581 Inf 0.47787483 0.7324577
# 11           11 0.5246316 0.07843246 Inf 0.37090683 0.6783564
# 12           12 0.6737701 0.06729742 Inf 0.54186963 0.8056707
# 13           13 0.3509183 0.07026810 Inf 0.21319530 0.4886412
# 14           14 0.4816399 0.07711141 Inf 0.33050429 0.6327754
# 15           15 0.2948167 0.07419354 Inf 0.14940001 0.4402333
# 16           16 0.4644784 0.10409786 Inf 0.26045034 0.6685064

# > m3_emm_contrast_df
#                            contrast     estimate         SE  df    asymp.LCL   asymp.UCL     z.ratio      p.value
# 1     week_integer1 - week_integer2  0.476630478 0.10937701 Inf  0.101900321  0.85136063  4.35768435 1.421868e-03
# 2     week_integer1 - week_integer3  0.373611826 0.12020073 Inf -0.038200861  0.78542451  3.10823249 1.273029e-01
# 3     week_integer1 - week_integer4  0.422630389 0.13899909 Inf -0.053586241  0.89884702  3.04052632 1.521205e-01
# 4     week_integer1 - week_integer5  0.206792467 0.21552594 Inf -0.531608308  0.94519324  0.95947834 9.998857e-01
# 5     week_integer1 - week_integer6  0.159191040 0.16471569 Inf -0.405131738  0.72351382  0.96645946 9.998748e-01
# 6     week_integer1 - week_integer7 -0.024305763 0.14520158 Inf -0.521772370  0.47316084 -0.16739324 1.000000e+00
# 7     week_integer1 - week_integer8  0.001563564 0.13489865 Inf -0.460604803  0.46373193  0.01159066 1.000000e+00
# 8     week_integer1 - week_integer9  0.102740701 0.13191943 Inf -0.349220722  0.55470213  0.77881402 9.999924e-01
# 9    week_integer1 - week_integer10  0.052914754 0.13330667 Inf -0.403799426  0.50962893  0.39694002 1.000000e+00
# 10   week_integer1 - week_integer11  0.133449393 0.14712743 Inf -0.370615268  0.63751405  0.90703272 9.999440e-01
# 11   week_integer1 - week_integer12 -0.015689124 0.13811625 Inf -0.488881095  0.45750285 -0.11359361 1.000000e+00
# 12   week_integer1 - week_integer13  0.307162766 0.13530258 Inf -0.156389463  0.77071500  2.27019153 6.461426e-01
# 13   week_integer1 - week_integer14  0.176441154 0.14052379 Inf -0.304999155  0.65788146  1.25559635 9.972411e-01
# 14   week_integer1 - week_integer15  0.363264348 0.14143461 Inf -0.121296494  0.84782519  2.56842605 4.209670e-01
# 15   week_integer1 - week_integer16  0.193602631 0.15931867 Inf -0.352229741  0.73943500  1.21519107 9.980819e-01
# 16    week_integer2 - week_integer3 -0.103018652 0.08120873 Inf -0.381243118  0.17520581 -1.26856623 9.969122e-01
# 17    week_integer2 - week_integer4 -0.054000089 0.10155783 Inf -0.401941411  0.29394123 -0.53171765 1.000000e+00
# 18    week_integer2 - week_integer5 -0.269838011 0.19312716 Inf -0.931499645  0.39182362 -1.39720386 9.914422e-01
# 19    week_integer2 - week_integer6 -0.317439438 0.13137514 Inf -0.767536103  0.13265723 -2.41628240 5.350522e-01
# 20    week_integer2 - week_integer7 -0.500936241 0.10606714 Inf -0.864326648 -0.13754583 -4.72282221 2.634061e-04
# 21    week_integer2 - week_integer8 -0.475066914 0.09202661 Inf -0.790353892 -0.15977994 -5.16227759 2.849342e-05
# 22    week_integer2 - week_integer9 -0.373889776 0.09026266 Inf -0.683133373 -0.06464618 -4.14224211 3.575139e-03
# 23   week_integer2 - week_integer10 -0.423715724 0.08833327 Inf -0.726349169 -0.12108228 -4.79678510 1.838228e-04
# 24   week_integer2 - week_integer11 -0.343181085 0.10663273 Inf -0.708509227  0.02214706 -3.21834664 9.376639e-02
# 25   week_integer2 - week_integer12 -0.492319601 0.09784960 Inf -0.827556369 -0.15708283 -5.03139121 5.646764e-05
# 26   week_integer2 - week_integer13 -0.169467711 0.09225267 Inf -0.485529177  0.14659375 -1.83699519 9.031687e-01
# 27   week_integer2 - week_integer14 -0.300189324 0.09887378 Inf -0.638934992  0.03855634 -3.03608619 1.538664e-01
# 28   week_integer2 - week_integer15 -0.113366129 0.09734646 Inf -0.446879145  0.22014689 -1.16456340 9.988176e-01
# 29   week_integer2 - week_integer16 -0.283027847 0.12442696 Inf -0.709319762  0.14326407 -2.27465049 6.428173e-01
# 30    week_integer3 - week_integer4  0.049018562 0.11048535 Inf -0.329508802  0.42754593  0.44366574 1.000000e+00
# 31    week_integer3 - week_integer5 -0.166819359 0.19745666 Inf -0.843314048  0.50967533 -0.84484038 9.999777e-01
# 32    week_integer3 - week_integer6 -0.214420786 0.14223054 Inf -0.701708486  0.27286691 -1.50755803 9.819249e-01
# 33    week_integer3 - week_integer7 -0.397917589 0.12193361 Inf -0.815667199  0.01983202 -3.26339533 8.227813e-02
# 34    week_integer3 - week_integer8 -0.372048262 0.11022181 Inf -0.749672754  0.00557623 -3.37545040 5.863152e-02
# 35    week_integer3 - week_integer9 -0.270871125 0.10779200 Inf -0.640170993  0.09842874 -2.51290553 4.618277e-01
# 36   week_integer3 - week_integer10 -0.320697072 0.10831687 Inf -0.691795155  0.05040101 -2.96073057 1.857919e-01
# 37   week_integer3 - week_integer11 -0.240162433 0.12578626 Inf -0.671111363  0.19078650 -1.90928989 8.722027e-01
# 38   week_integer3 - week_integer12 -0.389300950 0.11748516 Inf -0.791809957  0.01320806 -3.31361818 7.085511e-02
# 39   week_integer3 - week_integer13 -0.066449060 0.11181341 Inf -0.449526414  0.31662829 -0.59428527 9.999998e-01
# 40   week_integer3 - week_integer14 -0.197170672 0.11905744 Inf -0.605066403  0.21072506 -1.65609697 9.577447e-01
# 41   week_integer3 - week_integer15 -0.010347478 0.11854367 Inf -0.416482987  0.39578803 -0.08728832 1.000000e+00
# 42   week_integer3 - week_integer16 -0.180009195 0.14124394 Inf -0.663916767  0.30389838 -1.27445609 9.967522e-01
# 43    week_integer4 - week_integer5 -0.215837921 0.19571752 Inf -0.886374230  0.45469839 -1.10280329 9.993741e-01
# 44    week_integer4 - week_integer6 -0.263439349 0.14669448 Inf -0.766020718  0.23914202 -1.79583678 9.183885e-01
# 45    week_integer4 - week_integer7 -0.446936151 0.12692192 Inf -0.881775905 -0.01209640 -3.52134720 3.666086e-02
# 46    week_integer4 - week_integer8 -0.421066824 0.11669582 Inf -0.820871541 -0.02126211 -3.60824248 2.731184e-02
# 47    week_integer4 - week_integer9 -0.319889687 0.11449720 Inf -0.712161844  0.07238247 -2.79386463 2.722764e-01
# 48   week_integer4 - week_integer10 -0.369715634 0.11492165 Inf -0.763441966  0.02401070 -3.21711037 9.409905e-02
# 49   week_integer4 - week_integer11 -0.289180996 0.12920756 Inf -0.731851448  0.15348946 -2.23811202 6.698428e-01
# 50   week_integer4 - week_integer12 -0.438319512 0.12245967 Inf -0.857871414 -0.01876761 -3.57929681 3.016190e-02
# 51   week_integer4 - week_integer13 -0.115467622 0.11769106 Inf -0.518682070  0.28774683 -0.98110783 9.998490e-01
# 52   week_integer4 - week_integer14 -0.246189235 0.12491507 Inf -0.674153435  0.18177497 -1.97085295 8.415885e-01
# 53   week_integer4 - week_integer15 -0.059366040 0.12485908 Inf -0.487138401  0.36840632 -0.47546436 1.000000e+00
# 54   week_integer4 - week_integer16 -0.229027758 0.14583277 Inf -0.728656846  0.27060133 -1.57048217 9.735627e-01
# 55    week_integer5 - week_integer6 -0.047601427 0.20964880 Inf -0.765866901  0.67066405 -0.22705318 1.000000e+00
# 56    week_integer5 - week_integer7 -0.231098230 0.20139457 Inf -0.921084350  0.45888789 -1.14748989 9.990031e-01
# 57    week_integer5 - week_integer8 -0.205228903 0.19523826 Inf -0.874123257  0.46366545 -1.05117155 9.996475e-01
# 58    week_integer5 - week_integer9 -0.104051766 0.19431551 Inf -0.769784732  0.56168120 -0.53547845 1.000000e+00
# 59   week_integer5 - week_integer10 -0.153877713 0.19528528 Inf -0.822933155  0.51517773 -0.78796372 9.999911e-01
# 60   week_integer5 - week_integer11 -0.073343074 0.20494901 Inf -0.775506847  0.62882070 -0.35786011 1.000000e+00
# 61   week_integer5 - week_integer12 -0.222481591 0.19971204 Inf -0.906703294  0.46174011 -1.11401192 9.992948e-01
# 62   week_integer5 - week_integer13  0.100370299 0.19665006 Inf -0.573360948  0.77410155  0.51040055 1.000000e+00
# 63   week_integer5 - week_integer14 -0.030351313 0.20038446 Inf -0.716876772  0.65617415 -0.15146540 1.000000e+00
# 64   week_integer5 - week_integer15  0.156471881 0.20261453 Inf -0.537693878  0.85063764  0.77226388 9.999932e-01
# 65   week_integer5 - week_integer16 -0.013189836 0.21456768 Inf -0.748307590  0.72192792 -0.06147168 1.000000e+00
# 66    week_integer6 - week_integer7 -0.183496803 0.13674802 Inf -0.652001163  0.28500756 -1.34186081 9.943637e-01
# 67    week_integer6 - week_integer8 -0.157627475 0.12991161 Inf -0.602710021  0.28745507 -1.21334404 9.981144e-01
# 68    week_integer6 - week_integer9 -0.056450338 0.12914643 Inf -0.498911343  0.38601067 -0.43710337 1.000000e+00
# 69   week_integer6 - week_integer10 -0.106276286 0.12907681 Inf -0.548498776  0.33594620 -0.82335693 9.999841e-01
# 70   week_integer6 - week_integer11 -0.025741647 0.14402580 Inf -0.519179983  0.46769669 -0.17872942 1.000000e+00
# 71   week_integer6 - week_integer12 -0.174880163 0.13770155 Inf -0.646651380  0.29689105 -1.26999413 9.968740e-01
# 72   week_integer6 - week_integer13  0.147971727 0.13363063 Inf -0.309852346  0.60579580  1.10731892 9.993432e-01
# 73   week_integer6 - week_integer14  0.017250114 0.13711197 Inf -0.452501173  0.48700140  0.12581042 1.000000e+00
# 74   week_integer6 - week_integer15  0.204073309 0.13979191 Inf -0.274859566  0.68300618  1.45983631 9.867444e-01
# 75   week_integer6 - week_integer16  0.034411591 0.15655943 Inf -0.501967483  0.57079067  0.21979891 1.000000e+00
# 76    week_integer7 - week_integer8  0.025869327 0.09518755 Inf -0.300247150  0.35198580  0.27177218 1.000000e+00
# 77    week_integer7 - week_integer9  0.127046464 0.09121709 Inf -0.185467060  0.43955999  1.39279234 9.917133e-01
# 78   week_integer7 - week_integer10  0.077220517 0.09169884 Inf -0.236943496  0.39138453  0.84211009 9.999786e-01
# 79   week_integer7 - week_integer11  0.157755156 0.10220802 Inf -0.192413738  0.50792405  1.54347145 9.774573e-01
# 80   week_integer7 - week_integer12  0.008616639 0.09678256 Inf -0.322964409  0.34019769  0.08903091 1.000000e+00
# 81   week_integer7 - week_integer13  0.331468529 0.09876669 Inf -0.006910222  0.66984728  3.35607627 6.225477e-02
# 82   week_integer7 - week_integer14  0.200746917 0.10297760 Inf -0.152058617  0.55355245  1.94942306 8.526816e-01
# 83   week_integer7 - week_integer15  0.387570111 0.10298128 Inf  0.034751971  0.74038825  3.76350047 1.572132e-02
# 84   week_integer7 - week_integer16  0.217908394 0.12420410 Inf -0.207619999  0.64343679  1.75443798 9.319759e-01
# 85    week_integer8 - week_integer9  0.101177137 0.07169599 Inf -0.144456290  0.34681056  1.41119660 9.905337e-01
# 86   week_integer8 - week_integer10  0.051351190 0.07110562 Inf -0.192259596  0.29496198  0.72218190 9.999972e-01
# 87   week_integer8 - week_integer11  0.131885829 0.08474063 Inf -0.158439077  0.42221073  1.55634703 9.756608e-01
# 88   week_integer8 - week_integer12 -0.017252688 0.07631962 Inf -0.278726873  0.24422150 -0.22605835 1.000000e+00
# 89   week_integer8 - week_integer13  0.305599202 0.07721922 Inf  0.041042951  0.57015545  3.95755348 7.532153e-03
# 90   week_integer8 - week_integer14  0.174877590 0.08164482 Inf -0.104840925  0.45459610  2.14193136 7.377885e-01
# 91   week_integer8 - week_integer15  0.361700784 0.08453436 Inf  0.072082567  0.65131900  4.27874278 2.006245e-03
# 92   week_integer8 - week_integer16  0.192039067 0.10857517 Inf -0.179943944  0.56402208  1.76872000 9.274800e-01
# 93   week_integer9 - week_integer10 -0.049825947 0.06335358 Inf -0.266877923  0.16722603 -0.78647410 9.999913e-01
# 94   week_integer9 - week_integer11  0.030708691 0.07736708 Inf -0.234354136  0.29577152  0.39692192 1.000000e+00
# 95   week_integer9 - week_integer12 -0.118429825 0.06625023 Inf -0.345405869  0.10854622 -1.78761368 9.212231e-01
# 96   week_integer9 - week_integer13  0.204422065 0.06725618 Inf -0.026000379  0.43484451  3.03945415 1.525408e-01
# 97   week_integer9 - week_integer14  0.073700452 0.07316923 Inf -0.176980359  0.32438126  1.00726018 9.997907e-01
# 98   week_integer9 - week_integer15  0.260523647 0.07797384 Inf -0.006617951  0.52766525  3.34116739 6.516931e-02
# 99   week_integer9 - week_integer16  0.090861930 0.10253195 Inf -0.260416771  0.44214063  0.88618163 9.999585e-01
# 100 week_integer10 - week_integer11  0.080534639 0.07018758 Inf -0.159930904  0.32100018  1.14742013 9.990038e-01
# 101 week_integer10 - week_integer12 -0.068603878 0.06328528 Inf -0.285421864  0.14821411 -1.08404163 9.994896e-01
# 102 week_integer10 - week_integer13  0.254248012 0.06531041 Inf  0.030491858  0.47800417  3.89291733 9.677006e-03
# 103 week_integer10 - week_integer14  0.123526400 0.06997230 Inf -0.116201590  0.36325439  1.76536147 9.285552e-01
# 104 week_integer10 - week_integer15  0.310349594 0.07295453 Inf  0.060404347  0.56029484  4.25401388 2.231323e-03
# 105 week_integer10 - week_integer16  0.140687877 0.09889391 Inf -0.198126744  0.47950250  1.42261419 9.897353e-01
# 106 week_integer11 - week_integer12 -0.149138516 0.06831076 Inf -0.383174001  0.08489697 -2.18323614 7.092702e-01
# 107 week_integer11 - week_integer13  0.173713374 0.07591422 Inf -0.086371878  0.43379863  2.28828510 6.326083e-01
# 108 week_integer11 - week_integer14  0.042991761 0.08127260 Inf -0.235451524  0.32143505  0.52898224 1.000000e+00
# 109 week_integer11 - week_integer15  0.229814956 0.07950786 Inf -0.042582248  0.50221216  2.89046854 2.195482e-01
# 110 week_integer11 - week_integer16  0.060153238 0.10239943 Inf -0.290671447  0.41097792  0.58743723 9.999998e-01
# 111 week_integer12 - week_integer13  0.322851890 0.06825825 Inf  0.088996307  0.55670747  4.72985900 2.546097e-04
# 112 week_integer12 - week_integer14  0.192130277 0.06909604 Inf -0.044595610  0.42885617  2.78062652 2.800599e-01
# 113 week_integer12 - week_integer15  0.378953472 0.07284380 Inf  0.129387589  0.62851936  5.20227468 2.303609e-05
# 114 week_integer12 - week_integer16  0.209291754 0.09647366 Inf -0.121230994  0.53981450  2.16941864 7.189324e-01
# 115 week_integer13 - week_integer14 -0.130721612 0.07388004 Inf -0.383837672  0.12239445 -1.76937668 9.272685e-01
# 116 week_integer13 - week_integer15  0.056101582 0.07997635 Inf -0.217900700  0.33010386  0.70147716 9.999981e-01
# 117 week_integer13 - week_integer16 -0.113560135 0.10413297 Inf -0.470323994  0.24320372 -1.09053009 9.994520e-01
# 118 week_integer14 - week_integer15  0.186823194 0.07656541 Inf -0.075493057  0.44913945  2.44004705 5.168871e-01
# 119 week_integer14 - week_integer16  0.017161477 0.10019487 Inf -0.326110286  0.36043324  0.17128100 1.000000e+00
# 120 week_integer15 - week_integer16 -0.169661717 0.10136023 Inf -0.516926071  0.17760264 -1.67384892 9.537331e-01


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


m3_emm_contrast_df$group1 <- c("38", "39", "39", "39", "39", "39", "41", "41", "41", "44", "45", "45", "47", "47", "49", "49")
m3_emm_contrast_df$group2 <- c("2", "7", "8", "9", "10", "11", "12", "14", "7", "8", "10", "12", "7", "8", "9", "10", "12", "15", "15", "15", "13", "15""39", "44", "45", "46", "47", "49", "44", "45", "49", "52", "50", "52", "50", "52", "50", "52")
m3_emm_contrast_df$contrast_pretty <- c("38 - 39", "39 - 44", "39 - 45", "39 - 46", "39 - 47", "39 - 49", "41 - 44", "41 - 45", "41 - 49", "44 - 52", "45 - 50", "45 - 52", "47 - 50", "47 - 52", "49 - 50", "49 - 52")


#plot the result
bat_plot <- ggplot(filter(bat_metrics, week_integer_rename %in% c(38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 1)), aes(x = week_integer_rename, y = bat_emerged)) +
  geom_jitter(alpha = 0.3,
              height = 0.1,
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
  labs(y = "Probability of Emergence",
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
           x = -0.5,
           y = 14,
           label = "p <= 0.05",
           size = 4) +
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

ggsave(filename = paste0("figs/bat_daily_emergence_model_no_hibernation/emergence_model_trapping_week.png"), plot = bat_plot, width = 10, height = 6, dpi = 300)

ggsave(filename = paste0("figs/bat_daily_emergence_model_no_hibernation/emergence_model_trapping_week_effect.png"), plot = m3_plot, width = 10, height = 6, dpi = 300)
