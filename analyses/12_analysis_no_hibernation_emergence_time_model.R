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
                      night_length_rounded = as.integer(round((sunriseEnd - sunsetStart)*60, 0)), 
                      av_rel_humidity = av_rel_humidity, 
                      av_wind_speed = av_wind_speed, 
                      av_cloud_okta = av_cloud_okta, 
                      av_air_pressure = av_air_pressure, 
                      av_air_temp = av_air_temp, 
                      daily_rain = daily_rain, 
                      meanMoonlightIntensity = meanMoonlightIntensity,
                      meanMoonPhase = meanMoonPhase,
                      time_not_active = as.integer(night_length_rounded - absolute_activity_minutes),
                      absolute_activity_minutes = as.integer(absolute_activity_minutes))


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

bat_metrics <- droplevels(bat_metrics)



# create anchored emergence times relative to sunset ####
#create solar time column
tmp <- solartime (dat = bat_metrics$emergence_time, 
                  lat = 50.590305, 
                  lon = -3.612610, 
                  tz = 0)

bat_metrics$clock_emergence_time <- tmp$clock
bat_metrics$solar_emergence_time <- tmp$solar

bat_metrics$warp_emergence_time <- bat_metrics$solar_emergence_time*24/(2*pi)

tmp <- solartime (dat = bat_metrics$sunsetStart, 
                  lat = 50.590305, 
                  lon = -3.612610, 
                  tz = 0)

bat_metrics$clock_sunset_time <- tmp$clock
bat_metrics$solar_sunset_time <- tmp$solar

bat_metrics$warp_sunset_time <- bat_metrics$solar_sunset_time*24/(2*pi)

tmp <- solartime (dat = bat_metrics$sunriseEnd, 
                  lat = 50.590305, 
                  lon = -3.612610, 
                  tz = 0)

bat_metrics$clock_sunrise_time <- tmp$clock
bat_metrics$solar_sunrise_time <- tmp$solar

bat_metrics$warp_sunrise_time <- bat_metrics$solar_sunrise_time*24/(2*pi)



#keep only bats that emerged
bat_active_data <- filter(bat_metrics, bat_emerged == 1) 
bat_active_data <- droplevels(bat_active_data)


#make relative emergence time ans time period variables
bat_active_data <- mutate(bat_active_data, 
                          tim_fac = numFactor(time),
                          emergence_period = case_when(warp_emergence_time < 12 ~ "Late",
                                                       TRUE ~ "Early"),
                          relative_emergence_time = case_when(emergence_period == "Early" ~ warp_emergence_time - warp_sunset_time,
                                                              TRUE ~ (24-warp_sunset_time) + warp_emergence_time),
                          relative_emergence_time = relative_emergence_time*60,
                          emergence_period = as.factor(emergence_period))



#Prepare variables for modelling ####
#standardise variables
bat_active_data$av_rel_humidity_std <- bat_active_data$av_rel_humidity / max(bat_active_data$av_rel_humidity)
bat_active_data$av_wind_speed_std <- bat_active_data$av_wind_speed / max(bat_active_data$av_wind_speed)
bat_active_data$av_cloud_okta_std <- bat_active_data$av_cloud_okta / max(bat_active_data$av_cloud_okta)
bat_active_data$av_air_pressure_std <- (bat_active_data$av_air_pressure - mean(bat_active_data$av_air_pressure)) / sd(bat_active_data$av_air_pressure)
bat_active_data$av_air_temp_std <- (bat_active_data$av_air_temp - mean(bat_active_data$av_air_temp)) / sd(bat_active_data$av_air_temp)
bat_active_data$daily_rain_std <- bat_active_data$daily_rain / max(bat_active_data$daily_rain)
bat_active_data$forearm_std <- (bat_active_data$forearm - mean(bat_active_data$forearm)) / sd(bat_active_data$forearm)
bat_active_data$weight_std <- (bat_active_data$weight - mean(bat_active_data$weight)) / sd(bat_active_data$weight)



#fit first model ####
gamma_log_model_1 <- glmmTMB(relative_emergence_time ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + weight_std + sex + age + study_year + week_integer + sex*age + weight_std*av_air_temp_std + av_wind_speed_std*av_air_temp_std + (1 | ring_no),
                             data = bat_active_data,
                             family = Gamma("log"),
                             control = glmmTMBControl(parallel = 10))

summary(gamma_log_model_1)
fixef(gamma_log_model_1)$cond
check_collinearity(gamma_log_model_1)


#DHARMa simulation
gamma_log_model_1_simulationOutput <- simulateResiduals(fittedModel = gamma_log_model_1, plot = F, n = 1000)
plotQQunif(gamma_log_model_1_simulationOutput, testUniformity = TRUE, 
           testOutliers = TRUE, testDispersion = TRUE)
plotResiduals(gamma_log_model_1_simulationOutput, quantreg = TRUE, smoothScatter = FALSE) # deviations in 0.25 quantile
testDispersion(gamma_log_model_1_simulationOutput)
testOutliers(gamma_log_model_1_simulationOutput)
mod_1_res = recalculateResiduals(gamma_log_model_1_simulationOutput, group = bat_active_data$time)
testTemporalAutocorrelation(mod_1_res, time = unique(bat_active_data$time))



# #reduce to early only data try refitting #####
bat_active_data <- filter(bat_active_data, warp_emergence_time > min(bat_active_data$warp_sunset_time) & warp_emergence_time < (max(bat_active_data$warp_sunset_time) + 2))

bat_active_data <- droplevels(bat_active_data) 


#restandardize variables
bat_active_data$av_rel_humidity_std <- bat_active_data$av_rel_humidity / max(bat_active_data$av_rel_humidity)
bat_active_data$av_wind_speed_std <- bat_active_data$av_wind_speed / max(bat_active_data$av_wind_speed)
bat_active_data$av_cloud_okta_std <- bat_active_data$av_cloud_okta / max(bat_active_data$av_cloud_okta)
bat_active_data$av_air_pressure_std <- (bat_active_data$av_air_pressure - mean(bat_active_data$av_air_pressure)) / sd(bat_active_data$av_air_pressure)
bat_active_data$av_air_temp_std <- (bat_active_data$av_air_temp - mean(bat_active_data$av_air_temp)) / sd(bat_active_data$av_air_temp)
bat_active_data$daily_rain_std <- bat_active_data$daily_rain / max(bat_active_data$daily_rain)
bat_active_data$forearm_std <- (bat_active_data$forearm - mean(bat_active_data$forearm)) / sd(bat_active_data$forearm)
bat_active_data$weight_std <- (bat_active_data$weight - mean(bat_active_data$weight)) / sd(bat_active_data$weight)

#fit first model ####
gamma_log_model_1 <- glmmTMB(relative_emergence_time ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + weight_std + sex + age + study_year + week_integer + sex*age + weight_std*av_air_temp_std + av_wind_speed_std*av_air_temp_std + (1 | ring_no),
                             data = bat_active_data,
                             family = Gamma("log"),
                             control = glmmTMBControl(parallel = 10))

summary(gamma_log_model_1)
fixef(gamma_log_model_1)$cond
check_collinearity(gamma_log_model_1)


#DHARMa simulation
gamma_log_model_1_simulationOutput <- simulateResiduals(fittedModel = gamma_log_model_1, plot = F, n = 1000)
plotQQunif(gamma_log_model_1_simulationOutput, testUniformity = TRUE, 
           testOutliers = TRUE, testDispersion = TRUE)
plotResiduals(gamma_log_model_1_simulationOutput, quantreg = TRUE, smoothScatter = FALSE) # deviations in 0.25 quantile
testDispersion(gamma_log_model_1_simulationOutput)
testOutliers(gamma_log_model_1_simulationOutput)
mod_1_res = recalculateResiduals(gamma_log_model_1_simulationOutput, group = bat_active_data$time)
testTemporalAutocorrelation(mod_1_res, time = unique(bat_active_data$time))

plotResiduals(gamma_log_model_1_simulationOutput, form = bat_active_data$av_wind_speed_std, xlab = "Rank-transformed Standardised Wind Speed")
plotResiduals(gamma_log_model_1_simulationOutput, form = bat_active_data$av_air_temp_std, xlab = "Rank-transformed Standardised Air Temperature")
plotResiduals(gamma_log_model_1_simulationOutput, form = bat_active_data$meanMoonlightIntensity, xlab = "Rank-transformed Standardised Moonlight Intensity") #gotten worse
plotResiduals(gamma_log_model_1_simulationOutput, form = bat_active_data$av_cloud_okta_std, xlab = "Rank-transformed Standardised Cloud Cover") 
plotResiduals(gamma_log_model_1_simulationOutput, form = bat_active_data$daily_rain_std, xlab = "Rank-transformed Standardised Rainfall") # slight deviation in 0.25 quartile
plotResiduals(gamma_log_model_1_simulationOutput, form = bat_active_data$av_rel_humidity_std, xlab = "Rank-transformed Standardised Humidity")
plotResiduals(gamma_log_model_1_simulationOutput, form = bat_active_data$weight_std, xlab = "Rank-transformed Standardised Weight")
plotResiduals(gamma_log_model_1_simulationOutput, form = bat_active_data$age, xlab = "Rank-transformed Age")
plotResiduals(gamma_log_model_1_simulationOutput, form = bat_active_data$sex, xlab = "Rank-transformed Sex")
plotResiduals(gamma_log_model_1_simulationOutput, form = bat_active_data$study_year, xlab = "Rank-transformed Study Year")
plotResiduals(gamma_log_model_1_simulationOutput, form = bat_active_data$week_integer, xlab = "Rank-transformed Week")


#chanmge dispersion
gamma_log_model_1_disp_week <- update(gamma_log_model_1, dispformula = ~week_integer)

anova(gamma_log_model_1, gamma_log_model_1_disp_week)
#                             Df    AIC    BIC  logLik deviance  Chisq Chi Df  Pr(>Chisq)
# gamma_log_model_1           33 4348.2 4487.8 -2141.1   4282.2
# gamma_log_model_1_disp_week 48 4337.4 4540.4 -2120.7   4241.4 40.781     15  0.0003454 ***


gamma_log_model_1_disp_year <- update(gamma_log_model_1_disp_week, dispformula = ~week_integer + study_year)

anova(gamma_log_model_1_disp_week, gamma_log_model_1_disp_year)
#                             Df    AIC    BIC  logLik deviance  Chisq Chi Df  Pr(>Chisq)
# gamma_log_model_1_disp_week 48 4337.4 4540.4 -2120.7   4241.4
# gamma_log_model_1_disp_year 49 4322.8 4530.1 -2112.4   4224.8 16.539      1  4.766e-05 ***


gamma_log_model_1_disp_sex <- update(gamma_log_model_1_disp_year, dispformula = ~week_integer + study_year + sex)

anova(gamma_log_model_1_disp_year, gamma_log_model_1_disp_sex)
#                             Df    AIC    BIC  logLik deviance  Chisq Chi Df  Pr(>Chisq)
# gamma_log_model_1_disp_year 49 4322.8 4530.1 -2112.4   4224.8
# gamma_log_model_1_disp_sex  50 4323.2 4534.7 -2111.6   4223.2 1.6231      1      0.2027


gamma_log_model_1_disp_age <- update(gamma_log_model_1_disp_year, dispformula = ~week_integer + study_year + age)

anova(gamma_log_model_1_disp_year, gamma_log_model_1_disp_age)
#                             Df    AIC    BIC  logLik deviance  Chisq Chi Df  Pr(>Chisq)
# gamma_log_model_1_disp_year 49 4322.8 4530.1 -2112.4   4224.8
# gamma_log_model_1_disp_age  51 4321.9 4537.7 -2110.0   4219.9 4.8893      2     0.08676 .


gamma_log_model_1_disp_moon <- update(gamma_log_model_1_disp_year, dispformula = ~week_integer + study_year + meanMoonlightIntensity)

anova(gamma_log_model_1_disp_year, gamma_log_model_1_disp_moon)
#                             Df    AIC    BIC  logLik deviance  Chisq Chi Df  Pr(>Chisq)
# gamma_log_model_1_disp_year 49 4322.8 4530.1 -2112.4   4224.8
# gamma_log_model_1_disp_moon 50 4320.3 4531.8 -2110.1   4220.3 4.5785      1    0.03237 *


gamma_log_model_1_disp_cloud <- update(gamma_log_model_1_disp_moon, dispformula = ~week_integer + study_year + meanMoonlightIntensity + av_cloud_okta_std)

anova(gamma_log_model_1_disp_moon, gamma_log_model_1_disp_cloud)
#                             Df    AIC    BIC  logLik deviance  Chisq Chi Df  Pr(>Chisq)
# gamma_log_model_1_disp_moon  50 4320.3 4531.8 -2110.1   4220.3
# gamma_log_model_1_disp_cloud 51 4318.8 4534.5 -2108.4   4216.8 3.4864      1    0.06188 .


gamma_log_model_1_disp_temp <- update(gamma_log_model_1_disp_moon, dispformula = ~week_integer + study_year + meanMoonlightIntensity + av_air_temp_std)

anova(gamma_log_model_1_disp_moon, gamma_log_model_1_disp_temp)
#                             Df    AIC    BIC  logLik deviance  Chisq Chi Df  Pr(>Chisq)
# gamma_log_model_1_disp_moon 50 4320.3 4531.8 -2110.1   4220.3
# gamma_log_model_1_disp_temp 51 4318.6 4534.4 -2108.3   4216.6 3.6059      1    0.05757 .


gamma_log_model_1_disp_wind <- update(gamma_log_model_1_disp_moon, dispformula = ~week_integer + study_year + meanMoonlightIntensity + av_wind_speed_std)

anova(gamma_log_model_1_disp_moon, gamma_log_model_1_disp_wind)
#                             Df    AIC    BIC  logLik deviance  Chisq Chi Df  Pr(>Chisq)
# gamma_log_model_1_disp_moon 50 4320.3 4531.8 -2110.1   4220.3
# gamma_log_model_1_disp_wind 51 4315.3 4531.1 -2106.7   4213.3 6.949      1   0.008387 **


gamma_log_model_1_disp_rain <- update(gamma_log_model_1_disp_wind, dispformula = ~week_integer + study_year + meanMoonlightIntensity + av_wind_speed_std + daily_rain_std)

anova(gamma_log_model_1_disp_wind, gamma_log_model_1_disp_rain)
#                             Df    AIC    BIC  logLik deviance  Chisq Chi Df  Pr(>Chisq)
# gamma_log_model_1_disp_wind 51 4315.3 4531.1 -2106.7   4213.3
# gamma_log_model_1_disp_rain 52 4317.0 4537.0 -2106.5   4213.0 0.318      1     0.5728


gamma_log_model_1_disp_humid <- update(gamma_log_model_1_disp_wind, dispformula = ~week_integer + study_year + meanMoonlightIntensity + av_wind_speed_std + av_rel_humidity_std)

anova(gamma_log_model_1_disp_wind, gamma_log_model_1_disp_humid)
#                             Df    AIC    BIC  logLik deviance  Chisq Chi Df  Pr(>Chisq)
# gamma_log_model_1_disp_wind  51 4315.3 4531.1 -2106.7   4213.3
# gamma_log_model_1_disp_humid 52 4317.1 4537.1 -2106.5   4213.1 0.2248      1     0.6354


gamma_log_model_1_disp_weight <- update(gamma_log_model_1_disp_wind, dispformula = ~week_integer + study_year + meanMoonlightIntensity + av_wind_speed_std + weight_std)

anova(gamma_log_model_1_disp_wind, gamma_log_model_1_disp_weight)
#                             Df    AIC    BIC  logLik deviance  Chisq Chi Df  Pr(>Chisq)
# gamma_log_model_1_disp_wind   51 4315.3 4531.1 -2106.7   4213.3
# gamma_log_model_1_disp_weight 52 4313.8 4533.8 -2104.9   4209.8 3.4878      1    0.06182 .


model_1_disp <- update(gamma_log_model_1_disp_wind, dispformula = ~week_integer + study_year + av_wind_speed_std)

anova(gamma_log_model_1_disp_wind, model_1_disp)
#                             Df    AIC    BIC  logLik deviance  Chisq Chi Df  Pr(>Chisq)
# model_1_disp                50 4317.0 4528.5 -2108.5   4217.0
# gamma_log_model_1_disp_wind 51 4315.3 4531.1 -2106.7   4213.3 3.6474      1    0.05616 .


summary(model_1_disp)
#  Family: Gamma  ( log )
# Formula:
# relative_emergence_time ~ av_air_temp_std + av_wind_speed_std +
#     meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std +
#     av_rel_humidity_std + weight_std + sex + age + study_year +
#     week_integer + sex * age + weight_std * av_air_temp_std +
#     av_wind_speed_std * av_air_temp_std + (1 | ring_no)
# Dispersion:
# ~week_integer + study_year + av_wind_speed_std
# Data: bat_active_data

#       AIC       BIC    logLik -2*log(L)  df.resid
#    4317.0    4528.5   -2108.5    4217.0       458

# Random effects:

# Conditional model:
#  Groups  Name        Variance Std.Dev.
#  ring_no (Intercept) 0.04093  0.2023
# Number of obs: 508, groups:  ring_no, 34

# Conditional model:
#                                    Estimate Std. Error z value Pr(>|z|)
# (Intercept)                        4.382023   0.640593   6.841 7.89e-12 ***
# av_air_temp_std                   -0.041846   0.047759  -0.876   0.3809
# av_wind_speed_std                  0.007083   0.234904   0.030   0.9759
# meanMoonlightIntensity            -0.391157   0.264415  -1.479   0.1391
# av_cloud_okta_std                 -0.249426   0.100944  -2.471   0.0135 *
# daily_rain_std                     0.217063   0.164003   1.324   0.1857
# av_rel_humidity_std               -0.505339   0.576145  -0.877   0.3804
# weight_std                         0.029972   0.068829   0.435   0.6632
# sexM                               0.050608   0.131201   0.386   0.6997
# ageJ                              -0.297327   0.193963  -1.533   0.1253
# ageYA                              0.269552   0.172028   1.567   0.1171
# study_yearSecond Winter            0.040138   0.088877   0.452   0.6515
# week_integer2                     -0.070204   0.268123  -0.262   0.7934
# week_integer3                     -0.091612   0.288782  -0.317   0.7511
# week_integer4                     -0.336904   0.301679  -1.117   0.2641
# week_integer5                     -0.310326   0.282659  -1.098   0.2723
# week_integer6                      0.051472   0.292738   0.176   0.8604
# week_integer7                     -0.374562   0.292131  -1.282   0.1998
# week_integer8                      0.036475   0.281959   0.129   0.8971
# week_integer9                     -0.322816   0.281749  -1.146   0.2519
# week_integer10                    -0.215936   0.280954  -0.769   0.4421
# week_integer11                    -0.288443   0.294625  -0.979   0.3276
# week_integer12                    -0.143153   0.289462  -0.495   0.6209
# week_integer13                     0.100345   0.284562   0.353   0.7244
# week_integer14                     0.013740   0.289506   0.047   0.9621
# week_integer15                     0.160357   0.298319   0.538   0.5909
# week_integer16                    -0.001245   0.312750  -0.004   0.9968
# sexM:ageJ                          0.277327   0.239300   1.159   0.2465
# sexM:ageYA                         0.302799   0.295911   1.023   0.3062
# av_air_temp_std:weight_std        -0.011892   0.020310  -0.586   0.5582
# av_air_temp_std:av_wind_speed_std  0.021783   0.127034   0.171   0.8639
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Dispersion model:
#                         Estimate Std. Error z value Pr(>|z|)
# (Intercept)               1.5263     0.4538   3.363 0.000770 ***
# week_integer2             1.7150     0.5220   3.285 0.001019 **
# week_integer3             0.7058     0.5272   1.339 0.180677
# week_integer4             0.8493     0.5687   1.494 0.135282
# week_integer5             2.3358     0.6165   3.789 0.000151 ***
# week_integer6             0.8034     0.5437   1.478 0.139513
# week_integer7             0.9126     0.4829   1.890 0.058751 .
# week_integer8             1.1754     0.4644   2.531 0.011375 *
# week_integer9             1.0526     0.4264   2.468 0.013569 *
# week_integer10            1.5153     0.4319   3.509 0.000450 ***
# week_integer11            0.8003     0.4660   1.717 0.085902 .
# week_integer12            1.6903     0.4472   3.780 0.000157 ***
# week_integer13            1.6973     0.4693   3.616 0.000299 ***
# week_integer14            1.2808     0.4962   2.581 0.009848 **
# week_integer15            2.1052     0.5916   3.558 0.000373 ***
# week_integer16            1.4536     0.6815   2.133 0.032937 *
# study_yearSecond Winter  -0.6755     0.1871  -3.610 0.000306 ***
# av_wind_speed_std        -1.3285     0.4708  -2.822 0.004779 **
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

fixef(model_1_disp)$cond
check_collinearity(model_1_disp)

#DHARMa simulation
gamma_log_model_1_disp_simulationOutput <- simulateResiduals(fittedModel = model_1_disp, plot = F, n = 1000)
plotQQunif(gamma_log_model_1_disp_simulationOutput, testUniformity = TRUE, 
           testOutliers = TRUE, testDispersion = TRUE)
plotResiduals(gamma_log_model_1_disp_simulationOutput, quantreg = TRUE, smoothScatter = FALSE) # deviations in 0.25 quantile
testDispersion(gamma_log_model_1_disp_simulationOutput)
testOutliers(gamma_log_model_1_disp_simulationOutput)
mod_1_res = recalculateResiduals(gamma_log_model_1_disp_simulationOutput, group = bat_active_data$time)
testTemporalAutocorrelation(mod_1_res, time = unique(bat_active_data$time))

plotResiduals(gamma_log_model_1_disp_simulationOutput, form = bat_active_data$av_wind_speed_std, xlab = "Rank-transformed Standardised Wind Speed")
plotResiduals(gamma_log_model_1_disp_simulationOutput, form = bat_active_data$av_air_temp_std, xlab = "Rank-transformed Standardised Air Temperature")
plotResiduals(gamma_log_model_1_disp_simulationOutput, form = bat_active_data$meanMoonlightIntensity, xlab = "Rank-transformed Standardised Moonlight Intensity") #gotten worse
plotResiduals(gamma_log_model_1_disp_simulationOutput, form = bat_active_data$av_cloud_okta_std, xlab = "Rank-transformed Standardised Cloud Cover") 
plotResiduals(gamma_log_model_1_disp_simulationOutput, form = bat_active_data$daily_rain_std, xlab = "Rank-transformed Standardised Rainfall") # slight deviation in 0.25 quartile
plotResiduals(gamma_log_model_1_disp_simulationOutput, form = bat_active_data$av_rel_humidity_std, xlab = "Rank-transformed Standardised Humidity")
plotResiduals(gamma_log_model_1_disp_simulationOutput, form = bat_active_data$weight_std, xlab = "Rank-transformed Standardised Weight")
plotResiduals(gamma_log_model_1_disp_simulationOutput, form = bat_active_data$age, xlab = "Rank-transformed Age")
plotResiduals(gamma_log_model_1_disp_simulationOutput, form = bat_active_data$sex, xlab = "Rank-transformed Sex")
plotResiduals(gamma_log_model_1_disp_simulationOutput, form = bat_active_data$study_year, xlab = "Rank-transformed Study Year")
plotResiduals(gamma_log_model_1_disp_simulationOutput, form = bat_active_data$week_integer, xlab = "Rank-transformed Week")



#continue with model ####
#sex*age
model_1_disp.b <- glmmTMB(relative_emergence_time ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + weight_std + sex + age + study_year + week_integer + weight_std*av_air_temp_std + av_wind_speed_std*av_air_temp_std + (1 | ring_no),
                             dispformula = ~week_integer + study_year + av_wind_speed_std,
                             data = bat_active_data,
                             family = Gamma("log"),
                             control = glmmTMBControl(parallel = 10))

#weight_std*av_air_temp_std
model_1_disp.c <- glmmTMB(relative_emergence_time ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + weight_std + sex + age + study_year + week_integer + sex*age + av_wind_speed_std*av_air_temp_std + (1 | ring_no),
                             dispformula = ~week_integer + study_year + av_wind_speed_std,
                             data = bat_active_data,
                             family = Gamma("log"),
                             control = glmmTMBControl(parallel = 10))

#av_wind_speed_std*av_air_temp_std
model_1_disp.d <- glmmTMB(relative_emergence_time ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + weight_std + sex + age + study_year + week_integer + sex*age + weight_std*av_air_temp_std + (1 | ring_no),
                             dispformula = ~week_integer + study_year + av_wind_speed_std,
                             data = bat_active_data,
                             family = Gamma("log"),
                             control = glmmTMBControl(parallel = 10))


# test interactions
anova(model_1_disp, model_1_disp.b) #sex*age
#                               Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# model_1_disp.b 48 4314.9 4518.0 -2109.4   4218.9
# model_1_disp   50 4317.0 4528.5 -2108.5   4217.0 1.9525      2     0.3767

anova(model_1_disp, model_1_disp.c) #weight_std*av_air_temp_std
#                               Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# model_1_disp.c 49 4315.3 4522.6 -2108.7   4217.3
# model_1_disp   50 4317.0 4528.5 -2108.5   4217.0 0.3394      1     0.5602

anova(model_1_disp, model_1_disp.d) #av_wind_speed_std*av_air_temp_std
#                               Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# model_1_disp.d 49 4315 4522.3 -2108.5     4217
# model_1_disp   50 4317 4528.5 -2108.5     4217 0.0294      1     0.8639



#drop1 model selection ####
print("Dropping")
#model selection
drop1(model_1_disp, test = "Chi")

#                                   Df    AIC    LRT  Pr(>Chi)
# <none>                               4317.0
# meanMoonlightIntensity             1 4317.1  2.159   0.14170
# av_cloud_okta_std                  1 4321.1  6.108   0.01345 *
# daily_rain_std                     1 4316.7  1.788   0.18111
# av_rel_humidity_std                1 4315.7  0.770   0.38019
# study_year                         1 4315.2  0.203   0.65202
# week_integer                      15 4346.5 59.588 2.969e-07 ***
# sex:age                            2 4314.9  1.952   0.37673
# av_air_temp_std:weight_std         1 4315.3  0.339   0.56020
# av_air_temp_std:av_wind_speed_std  1 4315.0  0.029   0.86394


# drop sex:age interaction
drop1(model_1_disp.b, test = "Chi")


#                                   Df    AIC    LRT  Pr(>Chi)
# <none>                               4314.9
# meanMoonlightIntensity             1 4315.2  2.262   0.13255
# av_cloud_okta_std                  1 4319.1  6.207   0.01272 *
# daily_rain_std                     1 4314.6  1.730   0.18841
# av_rel_humidity_std                1 4313.6  0.719   0.39662
# sex                                1 4314.6  1.672   0.19593
# age                                2 4319.4  8.477   0.01443 *
# study_year                         1 4313.3  0.386   0.53438
# week_integer                      15 4344.2 59.261 3.377e-07 ***
# av_air_temp_std:weight_std         1 4313.2  0.287   0.59200
# av_air_temp_std:av_wind_speed_std  1 4312.9  0.024   0.87785


# drop av_air_temp_std:av_wind_speed_std interaction
model_2 <- glmmTMB(relative_emergence_time ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + weight_std + sex + age + study_year + week_integer + weight_std*av_air_temp_std + (1 | ring_no),
                             dispformula = ~week_integer + study_year + av_wind_speed_std,
                             data = bat_active_data,
                             family = Gamma("log"),
                             control = glmmTMBControl(parallel = 10))

 drop1(model_2, test = "Chi")

#                            Df    AIC    LRT  Pr(>Chi)
# <none>                        4312.9
# av_wind_speed_std           1 4310.9  0.009   0.92472
# meanMoonlightIntensity      1 4313.2  2.245   0.13406
# av_cloud_okta_std           1 4317.2  6.261   0.01234 *
# daily_rain_std              1 4312.6  1.713   0.19058
# av_rel_humidity_std         1 4311.6  0.696   0.40420
# sex                         1 4312.6  1.686   0.19407
# age                         2 4317.4  8.499   0.01427 *
# study_year                  1 4311.3  0.395   0.52981
# week_integer               15 4342.3 59.339 3.275e-07 ***
# av_air_temp_std:weight_std  1 4311.2  0.265   0.60663


# drop av_air_temp_std:weight_std interaction
model_3 <- glmmTMB(relative_emergence_time ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + weight_std + sex + age + study_year + week_integer + (1 | ring_no),
                             dispformula = ~week_integer + study_year + av_wind_speed_std,
                             data = bat_active_data,
                             family = Gamma("log"),
                             control = glmmTMBControl(parallel = 10))

drop1(model_3, test = "Chi")   

#                        Df    AIC    LRT  Pr(>Chi)
# <none>                    4311.2
# av_air_temp_std         1 4310.1  0.908   0.34066
# av_wind_speed_std       1 4309.2  0.009   0.92238
# meanMoonlightIntensity  1 4311.4  2.252   0.13340
# av_cloud_okta_std       1 4315.6  6.360   0.01167 *
# daily_rain_std          1 4310.8  1.603   0.20545
# av_rel_humidity_std     1 4310.0  0.755   0.38484
# weight_std              1 4309.3  0.109   0.74153
# sex                     1 4310.9  1.746   0.18641
# age                     2 4315.9  8.709   0.01285 *
# study_year              1 4309.6  0.416   0.51877
# week_integer           15 4340.3 59.076 3.634e-07 ***



summary(model_3)
#  Family: Gamma  ( log )
# Formula:
# relative_emergence_time ~ av_air_temp_std + av_wind_speed_std +
#     meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std +
#     av_rel_humidity_std + weight_std + sex + age + study_year +
#     week_integer + (1 | ring_no)
# Dispersion:
# ~week_integer + study_year + av_wind_speed_std
# Data: bat_active_data

#       AIC       BIC    logLik -2*log(L)  df.resid
#    4311.2    4505.8   -2109.6    4219.2       462

# Random effects:

# Conditional model:
#  Groups  Name        Variance Std.Dev.
#  ring_no (Intercept) 0.0474   0.2177
# Number of obs: 508, groups:  ring_no, 34

# Conditional model:
#                          Estimate Std. Error z value Pr(>|z|)
# (Intercept)              4.315909   0.630539   6.845 7.66e-12 ***
# av_air_temp_std         -0.037963   0.039955  -0.950   0.3420
# av_wind_speed_std        0.022366   0.229740   0.097   0.9224
# meanMoonlightIntensity  -0.400259   0.264938  -1.511   0.1308
# av_cloud_okta_std       -0.256077   0.101318  -2.527   0.0115 *
# daily_rain_std           0.204802   0.163462   1.253   0.2102
# av_rel_humidity_std     -0.487043   0.560339  -0.869   0.3847
# weight_std               0.022844   0.069222   0.330   0.7414
# sexM                     0.155158   0.117058   1.325   0.1850
# ageJ                    -0.123265   0.138339  -0.891   0.3729
# ageYA                    0.371124   0.150346   2.468   0.0136 *
# study_yearSecond Winter  0.059220   0.091455   0.648   0.5173
# week_integer2           -0.064115   0.269211  -0.238   0.8118
# week_integer3           -0.090412   0.289011  -0.313   0.7544
# week_integer4           -0.324196   0.304085  -1.066   0.2864
# week_integer5           -0.299758   0.284313  -1.054   0.2917
# week_integer6            0.024306   0.292098   0.083   0.9337
# week_integer7           -0.401402   0.291591  -1.377   0.1686
# week_integer8            0.006586   0.282232   0.023   0.9814
# week_integer9           -0.354255   0.281767  -1.257   0.2087
# week_integer10          -0.247893   0.280677  -0.883   0.3771
# week_integer11          -0.323530   0.294162  -1.100   0.2714
# week_integer12          -0.180432   0.289037  -0.624   0.5325
# week_integer13           0.067578   0.284664   0.237   0.8124
# week_integer14          -0.017905   0.289617  -0.062   0.9507
# week_integer15           0.127360   0.297369   0.428   0.6684
# week_integer16          -0.023525   0.313396  -0.075   0.9402
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Dispersion model:
#                         Estimate Std. Error z value Pr(>|z|)
# (Intercept)               1.4970     0.4524   3.309 0.000936 ***
# week_integer2             1.7906     0.5141   3.483 0.000497 ***
# week_integer3             0.7502     0.5236   1.433 0.151875
# week_integer4             0.8276     0.5682   1.457 0.145239
# week_integer5             2.3444     0.6094   3.847 0.000120 ***
# week_integer6             0.8686     0.5416   1.604 0.108803
# week_integer7             0.9234     0.4812   1.919 0.054966 .
# week_integer8             1.1906     0.4636   2.568 0.010218 *
# week_integer9             1.0609     0.4265   2.487 0.012880 *
# week_integer10            1.5260     0.4313   3.538 0.000402 ***
# week_integer11            0.8217     0.4642   1.770 0.076741 .
# week_integer12            1.7068     0.4459   3.828 0.000129 ***
# week_integer13            1.7204     0.4669   3.685 0.000229 ***
# week_integer14            1.3044     0.4950   2.635 0.008409 **
# week_integer15            2.1407     0.5872   3.645 0.000267 ***
# week_integer16            1.4266     0.6796   2.099 0.035799 *
# study_yearSecond Winter  -0.6594     0.1873  -3.521 0.000429 ***
# av_wind_speed_std        -1.3282     0.4679  -2.839 0.004530 **
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


fixef(model_3)$cond
check_collinearity(model_3)


model_3_simulationOutput <- simulateResiduals(fittedModel = model_3, plot = F, n = 1000)
plotQQunif(model_3_simulationOutput, testUniformity = TRUE, 
           testOutliers = TRUE, testDispersion = TRUE)
plotResiduals(model_3_simulationOutput, quantreg = TRUE, smoothScatter = FALSE) # deviations in 0.25 quantile
testDispersion(model_3_simulationOutput)
testOutliers(model_3_simulationOutput)
mod_3_res = recalculateResiduals(model_3_simulationOutput, group = bat_active_data$time)
testTemporalAutocorrelation(mod_3_res, time = unique(bat_active_data$time))

plotResiduals(model_3_simulationOutput, form = bat_active_data$av_wind_speed_std, xlab = "Rank-transformed Standardised Wind Speed")
plotResiduals(model_3_simulationOutput, form = bat_active_data$av_air_temp_std, xlab = "Rank-transformed Standardised Air Temperature")
plotResiduals(model_3_simulationOutput, form = bat_active_data$meanMoonlightIntensity, xlab = "Rank-transformed Standardised Moonlight Intensity") #gotten worse
plotResiduals(model_3_simulationOutput, form = bat_active_data$av_cloud_okta_std, xlab = "Rank-transformed Standardised Cloud Cover") 
plotResiduals(model_3_simulationOutput, form = bat_active_data$daily_rain_std, xlab = "Rank-transformed Standardised Rainfall") # slight deviation in 0.25 quartile
plotResiduals(model_3_simulationOutput, form = bat_active_data$av_rel_humidity_std, xlab = "Rank-transformed Standardised Humidity")
plotResiduals(model_3_simulationOutput, form = bat_active_data$weight_std, xlab = "Rank-transformed Standardised Weight")
plotResiduals(model_3_simulationOutput, form = bat_active_data$age, xlab = "Rank-transformed Age")
plotResiduals(model_3_simulationOutput, form = bat_active_data$sex, xlab = "Rank-transformed Sex")
plotResiduals(model_3_simulationOutput, form = bat_active_data$study_year, xlab = "Rank-transformed Study Year")
plotResiduals(model_3_simulationOutput, form = bat_active_data$week_integer, xlab = "Rank-transformed Week")



#Test significance of main effects ####
model_3.a <- glmmTMB(relative_emergence_time ~ av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + weight_std + sex + age + study_year + week_integer + (1 | ring_no),
                             dispformula = ~week_integer + study_year + av_wind_speed_std,
                             data = bat_active_data,
                             family = Gamma("log"),
                             control = glmmTMBControl(parallel = 10))  

model_3.b <- glmmTMB(relative_emergence_time ~ av_air_temp_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + weight_std + sex + age + study_year + week_integer + (1 | ring_no),
                             dispformula = ~week_integer + study_year + av_wind_speed_std,
                             data = bat_active_data,
                             family = Gamma("log"),
                             control = glmmTMBControl(parallel = 10))  

model_3.c <- glmmTMB(relative_emergence_time ~ av_air_temp_std + av_wind_speed_std + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + weight_std + sex + age + study_year + week_integer + (1 | ring_no),
                             dispformula = ~week_integer + study_year + av_wind_speed_std,
                             data = bat_active_data,
                             family = Gamma("log"),
                             control = glmmTMBControl(parallel = 10))  

model_3.d <- glmmTMB(relative_emergence_time ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + daily_rain_std + av_rel_humidity_std + weight_std + sex + age + study_year + week_integer + (1 | ring_no),
                             dispformula = ~week_integer + study_year + av_wind_speed_std,
                             data = bat_active_data,
                             family = Gamma("log"),
                             control = glmmTMBControl(parallel = 10))  

model_3.e <- glmmTMB(relative_emergence_time ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + av_rel_humidity_std + weight_std + sex + age + study_year + week_integer + (1 | ring_no),
                             dispformula = ~week_integer + study_year + av_wind_speed_std,
                             data = bat_active_data,
                             family = Gamma("log"),
                             control = glmmTMBControl(parallel = 10))  

model_3.f <- glmmTMB(relative_emergence_time ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + weight_std + sex + age + study_year + week_integer + (1 | ring_no),
                             dispformula = ~week_integer + study_year + av_wind_speed_std,
                             data = bat_active_data,
                             family = Gamma("log"),
                             control = glmmTMBControl(parallel = 10))  

model_3.g <- glmmTMB(relative_emergence_time ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + sex + age + study_year + week_integer + (1 | ring_no),
                             dispformula = ~week_integer + study_year + av_wind_speed_std,
                             data = bat_active_data,
                             family = Gamma("log"),
                             control = glmmTMBControl(parallel = 10))  

model_3.h <- glmmTMB(relative_emergence_time ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + weight_std + age + study_year + week_integer + (1 | ring_no),
                             dispformula = ~week_integer + study_year + av_wind_speed_std,
                             data = bat_active_data,
                             family = Gamma("log"),
                             control = glmmTMBControl(parallel = 10))   

model_3.i <- glmmTMB(relative_emergence_time ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + weight_std + sex + study_year + week_integer + (1 | ring_no),
                             dispformula = ~week_integer + study_year + av_wind_speed_std,
                             data = bat_active_data,
                             family = Gamma("log"),
                             control = glmmTMBControl(parallel = 10))  

model_3.j <- glmmTMB(relative_emergence_time ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + weight_std + sex + age + week_integer + (1 | ring_no),
                             dispformula = ~week_integer + study_year + av_wind_speed_std,
                             data = bat_active_data,
                             family = Gamma("log"),
                             control = glmmTMBControl(parallel = 10))  

model_3.k <- glmmTMB(relative_emergence_time ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + weight_std + sex + age + study_year + (1 | ring_no),
                             dispformula = ~week_integer + study_year + av_wind_speed_std,
                             data = bat_active_data,
                             family = Gamma("log"),
                             control = glmmTMBControl(parallel = 10))  


anova(model_3, model_3.a) #temp
#               Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# model_3.a 45 4310.1 4500.5 -2110.1   4220.1
# model_3   46 4311.2 4505.8 -2109.6   4219.2 0.908      1     0.3407

anova(model_3, model_3.b) #wind 
#           Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# model_3.b 45 4309.2 4499.6 -2109.6   4219.2
# model_3   46 4311.2 4505.8 -2109.6   4219.2 0.0095      1     0.9224

anova(model_3, model_3.c) #moonlight 
#           Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# model_3.c 45 4311.4 4501.8 -2110.7   4221.4
# model_3   46 4311.2 4505.8 -2109.6   4219.2 2.2525      1     0.1334

anova(model_3, model_3.d) #cloud
#           Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# model_3.d 45 4315.6 4505.9 -2112.8   4225.6
# model_3   46 4311.2 4505.8 -2109.6   4219.2 6.3599      1    0.01167 *

anova(model_3, model_3.e) #rain
#           Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# model_3.e 45 4310.8 4501.2 -2110.4   4220.8
# model_3   46 4311.2 4505.8 -2109.6   4219.2 1.6032      1     0.2054

anova(model_3, model_3.f) #humidity
#           Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# model_3.f 45 4310.0 4500.3 -2110.0   4220.0
# model_3   46 4311.2 4505.8 -2109.6   4219.2 0.7552      1     0.3848

anova(model_3, model_3.g) #weight
#           Df  AIC    BIC logLik deviance Chisq Chi Df Pr(>Chisq)
# model_3.g 45 4309.3 4499.7 -2109.7   4219.3
# model_3   46 4311.2 4505.8 -2109.6   4219.2 0.1088      1     0.7415

anova(model_3, model_3.h) #sex
#           Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# model_3.h 45 4310.9 4501.3 -2110.5   4220.9
# model_3   46 4311.2 4505.8 -2109.6   4219.2 1.7457      1     0.1864

anova(model_3, model_3.i) #age
#           Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# model_3.i 44 4315.9 4502.0 -2113.9   4227.9
# model_3   46 4311.2 4505.8 -2109.6   4219.2 8.7093      2    0.01285 *

anova(model_3, model_3.j) #year
#           Df  AIC    BIC logLik deviance Chisq Chi Df Pr(>Chisq)
# model_3.j 45 4309.6 4500.0 -2109.8   4219.6
# model_3   46 4311.2 4505.8 -2109.6   4219.2 0.4163      1     0.5188

anova(model_3, model_3.k) #week
#           Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# model_3.k 31 4340.3 4471.4 -2139.1   4278.3
# model_3   46 4311.2 4505.8 -2109.6   4219.2 59.076     15  3.634e-07 ***



# plotting - model_3 #####
# av_air_temp_std
MyData <- crossing(age = "A",
                   sex = "F",
                   study_year = "Second Winter",
                   week_integer = 14,
                   weight_std = mean(bat_active_data$weight_std),
                   av_wind_speed_std = mean(bat_active_data$av_wind_speed_std),
                   daily_rain_std = mean(bat_active_data$daily_rain_std),
                   av_air_temp_std = seq(from = min(bat_active_data$av_air_temp_std),
                                           to = max(bat_active_data$av_air_temp_std),
                                       length = 100),
                   av_rel_humidity_std = mean(bat_active_data$av_rel_humidity_std),
                   av_cloud_okta_std = mean(bat_active_data$av_cloud_okta_std),
                   meanMoonlightIntensity = mean(bat_active_data$meanMoonlightIntensity))

#' #' To compute population-level predictions for a given grouping
MyData$ring_no <- NA
MyData$time <- NA

#' Calculate the predicted  values.
Fit_mod <- predict(model_3,
                   newdata = MyData,
                   re.form = NULL,  #~0, Keep it to NULL..gives odd results with ~0
                   type = "link",
                   se = TRUE)

attr(Fit_mod$fit, "names") <- NULL
attr(Fit_mod$se.fit, "names") <- NULL

#' Extract the predicted values and standard errors, and calculate
#' standard errors (SE) for predicted values.
MyData$log_mu   <- Fit_mod$fit
MyData$mu <- exp(Fit_mod$fit)
MyData$SeUp <- exp(Fit_mod$fit + 1.96 * Fit_mod$se.fit)
MyData$SeLo <- exp(Fit_mod$fit - 1.96 * Fit_mod$se.fit)
MyData$se   <- Fit_mod$se.fit
head(MyData)

#back standardise
MyData$av_air_temp <- MyData$av_air_temp_std * sd(bat_active_data$av_air_temp) +
  mean(bat_active_data$av_air_temp)


#plot the result
bat_plot <- ggplot(MyData, aes(x = av_air_temp, y = mu)) +
  geom_ribbon(aes(ymax = SeUp,
                  ymin = SeLo),
              alpha = 0.03,
              show.legend = FALSE) +
  geom_line(linewidth = 1) +
  geom_jitter(data = bat_active_data,
              aes(y = relative_emergence_time),
              alpha = 0.3,
              height = 0.01,
              width = 0.01) +
  labs(y = "Relative Emergence Time \n (minutes from sunset)",
       x = "Average Air Temperature (DegC)") +
  theme_bw()
bat_plot

ggsave(filename = paste0("figs/bat_daily_emergence_time_model_no_hibernation/emergence_time_model_trapping_temp.png"), plot = bat_plot, width = 10, height = 6, dpi = 300)


# av_wind_speed_std
MyData <- crossing(age = "A",
                   sex = "F",
                   study_year = "Second Winter",
                   week_integer = 14,
                   weight_std = mean(bat_active_data$weight_std),
                   av_wind_speed_std = seq(from = min(bat_active_data$av_wind_speed_std),
                                           to = max(bat_active_data$av_wind_speed_std),
                                       length = 100),
                   daily_rain_std = mean(bat_active_data$daily_rain_std),
                   av_air_temp_std = mean(bat_active_data$av_air_temp_std),
                   av_rel_humidity_std = mean(bat_active_data$av_rel_humidity_std),
                   av_cloud_okta_std = mean(bat_active_data$av_cloud_okta_std),
                   meanMoonlightIntensity = mean(bat_active_data$meanMoonlightIntensity))

#' #' To compute population-level predictions for a given grouping
MyData$ring_no <- NA
MyData$time <- NA

#' Calculate the predicted  values.
Fit_mod <- predict(model_3,
                   newdata = MyData,
                   re.form = NULL,  #~0, Keep it to NULL..gives odd results with ~0
                   type = "link",
                   se = TRUE)

attr(Fit_mod$fit, "names") <- NULL
attr(Fit_mod$se.fit, "names") <- NULL

#' Extract the predicted values and standard errors, and calculate
#' standard errors (SE) for predicted values.
MyData$log_mu   <- Fit_mod$fit
MyData$mu <- exp(Fit_mod$fit)
MyData$SeUp <- exp(Fit_mod$fit + 1.96 * Fit_mod$se.fit)
MyData$SeLo <- exp(Fit_mod$fit - 1.96 * Fit_mod$se.fit)
MyData$se   <- Fit_mod$se.fit

#' dispersion parameter (σ)
Fit_mod <- predict(model_3,
                   newdata = MyData,
                   re.form = NA,
                   type = "disp")

attr(Fit_mod, "names") <- NULL

MyData$sigma   <- Fit_mod


MyData$variance <- (MyData$mu^2) * (MyData$sigma^2)

MyData$upper_var <- MyData$mu + 2*sqrt(MyData$variance)
MyData$lower_var <- MyData$mu - 2*sqrt(MyData$variance)

head(MyData)

#back standardise
MyData$av_wind_speed <- MyData$av_wind_speed_std * max(bat_active_data$av_wind_speed)


#plot the result
bat_plot <- ggplot(MyData, aes(x = av_wind_speed, y = mu)) +
  geom_ribbon(aes(ymax = SeUp,
                  ymin = SeLo),
              alpha = 0.1,
              show.legend = FALSE) +
  geom_line(aes(y = upper_var), linetype = "dashed", colour = "blue", linewidth = 0.8) +
  geom_line(aes(y = lower_var), linetype = "dashed", colour = "blue", linewidth = 0.8) +
  geom_line(linewidth = 1) +
  geom_jitter(data = bat_active_data,
              aes(y = relative_emergence_time),
              alpha = 0.3,
              height = 0.01,
              width = 0.01) +
  ylim(c(0, (max(bat_active_data$relative_emergence_time)+5))) +
  labs(y = "Relative Emergence Time \n (minutes from sunset)",
       x = "Average Wind Speed (kn)") +
  theme_bw()
bat_plot

ggsave(filename = paste0("figs/bat_daily_emergence_time_model_no_hibernation/emergence_time_model_trapping_wind.png"), plot = bat_plot, width = 10, height = 6, dpi = 300)


# meanMoonlightIntensity
MyData <- crossing(age = "A",
                   sex = "F",
                   study_year = "Second Winter",
                   week_integer = 14,
                   weight_std = mean(bat_active_data$weight_std),
                   av_wind_speed_std = mean(bat_active_data$av_wind_speed_std),
                   daily_rain_std = mean(bat_active_data$daily_rain_std),
                   av_air_temp_std = mean(bat_active_data$av_air_temp_std),
                   av_rel_humidity_std = mean(bat_active_data$av_rel_humidity_std),
                   av_cloud_okta_std = mean(bat_active_data$av_cloud_okta_std),
                   meanMoonlightIntensity = seq(from = min(bat_active_data$meanMoonlightIntensity),
                                           to = max(bat_active_data$meanMoonlightIntensity),
                                       length = 100))

#' #' To compute population-level predictions for a given grouping
MyData$ring_no <- NA
MyData$time <- NA

#' Calculate the predicted  values.
Fit_mod <- predict(model_3,
                   newdata = MyData,
                   re.form = NULL,  #~0, Keep it to NULL..gives odd results with ~0
                   type = "link",
                   se = TRUE)

attr(Fit_mod$fit, "names") <- NULL
attr(Fit_mod$se.fit, "names") <- NULL

#' Extract the predicted values and standard errors, and calculate
#' standard errors (SE) for predicted values.
MyData$log_mu   <- Fit_mod$fit
MyData$mu <- exp(Fit_mod$fit)
MyData$SeUp <- exp(Fit_mod$fit + 1.96 * Fit_mod$se.fit)
MyData$SeLo <- exp(Fit_mod$fit - 1.96 * Fit_mod$se.fit)
MyData$se   <- Fit_mod$se.fit
head(MyData)

#back standardise
MyData$meanMoonlightIntensity_lux <- MyData$meanMoonlightIntensity * 0.3


#plot the result
bat_plot <- ggplot(MyData, aes(x = meanMoonlightIntensity, y = mu)) +
  geom_ribbon(aes(ymax = SeUp,
                  ymin = SeLo),
              alpha = 0.03,
              show.legend = FALSE) +
  geom_line(linewidth = 1) +
  geom_jitter(data = bat_active_data,
              aes(y = relative_emergence_time),
              alpha = 0.3,
              height = 0.01,
              width = 0.01) +
  labs(y = "Relative Emergence Time \n (minutes from sunset)",
       x = "Average Moonlight Intensity (lux)") +
  theme_bw()
bat_plot

ggsave(filename = paste0("figs/bat_daily_emergence_time_model_no_hibernation/emergence_time_model_trapping_moon.png"), plot = bat_plot, width = 10, height = 6, dpi = 300)


# av_cloud_okta_std
MyData <- crossing(age = "A",
                   sex = "F",
                   study_year = "Second Winter",
                   week_integer = 14,
                   weight_std = mean(bat_active_data$weight_std),
                   av_wind_speed_std = mean(bat_active_data$av_wind_speed_std),
                   daily_rain_std = mean(bat_active_data$daily_rain_std),
                   av_air_temp_std = mean(bat_active_data$av_air_temp_std),
                   av_rel_humidity_std = mean(bat_active_data$av_rel_humidity_std),
                   av_cloud_okta_std = seq(from = min(bat_active_data$av_cloud_okta_std),
                                           to = max(bat_active_data$av_cloud_okta_std),
                                       length = 100),
                   meanMoonlightIntensity = mean(bat_active_data$meanMoonlightIntensity))

#' #' To compute population-level predictions for a given grouping
MyData$ring_no <- NA
MyData$time <- NA

#' Calculate the predicted  values.
Fit_mod <- predict(model_3,
                   newdata = MyData,
                   re.form = NULL,  #~0, Keep it to NULL..gives odd results with ~0
                   type = "link",
                   se = TRUE)

attr(Fit_mod$fit, "names") <- NULL
attr(Fit_mod$se.fit, "names") <- NULL

#' Extract the predicted values and standard errors, and calculate
#' standard errors (SE) for predicted values.
MyData$log_mu   <- Fit_mod$fit
MyData$mu <- exp(Fit_mod$fit)
MyData$SeUp <- exp(Fit_mod$fit + 1.96 * Fit_mod$se.fit)
MyData$SeLo <- exp(Fit_mod$fit - 1.96 * Fit_mod$se.fit)
MyData$se   <- Fit_mod$se.fit
head(MyData)

#back standardise
MyData$av_cloud_okta <- MyData$av_cloud_okta_std * max(bat_active_data$av_cloud_okta)


#plot the result
bat_plot <- ggplot(MyData, aes(x = av_cloud_okta, y = mu)) +
  geom_ribbon(aes(ymax = SeUp,
                  ymin = SeLo),
              alpha = 0.03,
              show.legend = FALSE) +
  geom_line(linewidth = 1) +
  geom_jitter(data = bat_active_data,
              aes(y = relative_emergence_time),
              alpha = 0.3,
              height = 0.01,
              width = 0.01) +
  labs(y = "Relative Emergence Time \n (minutes from sunset)",
       x = "Average Cloud Cover (Oktas)") +
  theme_bw()
bat_plot

ggsave(filename = paste0("figs/bat_daily_emergence_time_model_no_hibernation/emergence_time_model_trapping_cloud.png"), plot = bat_plot, width = 10, height = 6, dpi = 300)


# daily_rain_std 
MyData <- crossing(age = "A",
                   sex = "F",
                   study_year = "Second Winter",
                   week_integer = 14,
                   weight_std = mean(bat_active_data$weight_std),
                   av_wind_speed_std = mean(bat_active_data$av_wind_speed_std),
                   daily_rain_std = seq(from = min(bat_active_data$daily_rain_std),
                                           to = max(bat_active_data$daily_rain_std),
                                       length = 100),
                   av_air_temp_std = mean(bat_active_data$av_air_temp_std),
                   av_rel_humidity_std = mean(bat_active_data$av_rel_humidity_std),
                   av_cloud_okta_std = mean(bat_active_data$av_cloud_okta_std),
                   meanMoonlightIntensity = mean(bat_active_data$meanMoonlightIntensity))

#' #' To compute population-level predictions for a given grouping
MyData$ring_no <- NA
MyData$time <- NA

#' Calculate the predicted  values.
Fit_mod <- predict(model_3,
                   newdata = MyData,
                   re.form = NULL,  #~0, Keep it to NULL..gives odd results with ~0
                   type = "link",
                   se = TRUE)

attr(Fit_mod$fit, "names") <- NULL
attr(Fit_mod$se.fit, "names") <- NULL

#' Extract the predicted values and standard errors, and calculate
#' standard errors (SE) for predicted values.
MyData$log_mu   <- Fit_mod$fit
MyData$mu <- exp(Fit_mod$fit)
MyData$SeUp <- exp(Fit_mod$fit + 1.96 * Fit_mod$se.fit)
MyData$SeLo <- exp(Fit_mod$fit - 1.96 * Fit_mod$se.fit)
MyData$se   <- Fit_mod$se.fit
head(MyData)

#back standardise
MyData$daily_rain <- MyData$daily_rain_std * max(bat_active_data$daily_rain)


#plot the result
bat_plot <- ggplot(MyData, aes(x = daily_rain, y = mu)) +
  geom_ribbon(aes(ymax = SeUp,
                  ymin = SeLo),
              alpha = 0.03,
              show.legend = FALSE) +
  geom_line(linewidth = 1) +
  geom_jitter(data = bat_active_data,
              aes(y = relative_emergence_time),
              alpha = 0.3,
              height = 0.01,
              width = 0.01) +
  labs(y = "Relative Emergence Time \n (minutes from sunset)",
       x = "Total Rainfall (mm)") +
  theme_bw()
bat_plot

ggsave(filename = paste0("figs/bat_daily_emergence_time_model_no_hibernation/emergence_time_model_trapping_rain.png"), plot = bat_plot, width = 10, height = 6, dpi = 300)


# av_rel_humidity_std
MyData <- crossing(age = "A",
                   sex = "F",
                   study_year = "Second Winter",
                   week_integer = 14,
                   weight_std = mean(bat_active_data$weight_std),
                   av_wind_speed_std = mean(bat_active_data$av_wind_speed_std),
                   daily_rain_std = mean(bat_active_data$daily_rain_std),
                   av_air_temp_std = mean(bat_active_data$av_air_temp_std),
                   av_rel_humidity_std = seq(from = min(bat_active_data$av_rel_humidity_std),
                                           to = max(bat_active_data$av_rel_humidity_std),
                                       length = 100),
                   av_cloud_okta_std = mean(bat_active_data$av_cloud_okta_std),
                   meanMoonlightIntensity = mean(bat_active_data$meanMoonlightIntensity))

#' #' To compute population-level predictions for a given grouping
MyData$ring_no <- NA
MyData$time <- NA

#' Calculate the predicted  values.
Fit_mod <- predict(model_3,
                   newdata = MyData,
                   re.form = NULL,  #~0, Keep it to NULL..gives odd results with ~0
                   type = "link",
                   se = TRUE)

attr(Fit_mod$fit, "names") <- NULL
attr(Fit_mod$se.fit, "names") <- NULL

#' Extract the predicted values and standard errors, and calculate
#' standard errors (SE) for predicted values.
MyData$log_mu   <- Fit_mod$fit
MyData$mu <- exp(Fit_mod$fit)
MyData$SeUp <- exp(Fit_mod$fit + 1.96 * Fit_mod$se.fit)
MyData$SeLo <- exp(Fit_mod$fit - 1.96 * Fit_mod$se.fit)
MyData$se   <- Fit_mod$se.fit
head(MyData)

#back standardise
MyData$av_rel_humidity <- MyData$av_rel_humidity_std * max(bat_active_data$av_rel_humidity)


#plot the result
bat_plot <- ggplot(MyData, aes(x = av_rel_humidity, y = mu)) +
  geom_ribbon(aes(ymax = SeUp,
                  ymin = SeLo),
              alpha = 0.03,
              show.legend = FALSE) +
  geom_line(linewidth = 1) +
  geom_jitter(data = bat_active_data,
              aes(y = relative_emergence_time),
              alpha = 0.3,
              height = 0.01,
              width = 0.01) +
  labs(y = "Relative Emergence Time \n (minutes from sunset)",
       x = "Average Relative Humidity (%)") +
  theme_bw()
bat_plot

ggsave(filename = paste0("figs/bat_daily_emergence_time_model_no_hibernation/emergence_time_model_trapping_humid.png"), plot = bat_plot, width = 10, height = 6, dpi = 300)


# weight_std
MyData <- crossing(age = "A",
                   sex = "F",
                   study_year = "Second Winter",
                   week_integer = 14,
                   weight_std = seq(from = min(bat_active_data$weight_std),
                                           to = max(bat_active_data$weight_std),
                                       length = 100),
                   av_wind_speed_std = mean(bat_active_data$av_wind_speed_std),
                   daily_rain_std = mean(bat_active_data$daily_rain_std),
                   av_air_temp_std = mean(bat_active_data$av_air_temp_std),
                   av_rel_humidity_std = mean(bat_active_data$av_rel_humidity_std),
                   av_cloud_okta_std = mean(bat_active_data$av_cloud_okta_std),
                   meanMoonlightIntensity = mean(bat_active_data$meanMoonlightIntensity))

#' #' To compute population-level predictions for a given grouping
MyData$ring_no <- NA
MyData$time <- NA

#' Calculate the predicted  values.
Fit_mod <- predict(model_3,
                   newdata = MyData,
                   re.form = NULL,  #~0, Keep it to NULL..gives odd results with ~0
                   type = "link",
                   se = TRUE)

attr(Fit_mod$fit, "names") <- NULL
attr(Fit_mod$se.fit, "names") <- NULL

#' Extract the predicted values and standard errors, and calculate
#' standard errors (SE) for predicted values.
MyData$log_mu   <- Fit_mod$fit
MyData$mu <- exp(Fit_mod$fit)
MyData$SeUp <- exp(Fit_mod$fit + 1.96 * Fit_mod$se.fit)
MyData$SeLo <- exp(Fit_mod$fit - 1.96 * Fit_mod$se.fit)
MyData$se   <- Fit_mod$se.fit
head(MyData)

#back standardise
MyData$weight <- MyData$weight_std  * sd(bat_active_data$weight) +
  mean(bat_active_data$weight)


#plot the result
bat_plot <- ggplot(MyData, aes(x = weight, y = mu)) +
  geom_ribbon(aes(ymax = SeUp,
                  ymin = SeLo),
              alpha = 0.03,
              show.legend = FALSE) +
  geom_line(linewidth = 1) +
  geom_jitter(data = bat_active_data,
              aes(y = relative_emergence_time),
              alpha = 0.3,
              height = 0.01,
              width = 0.01) +
  labs(y = "Relative Emergence Time \n (minutes from sunset)",
       x = "Weight (g)") +
  theme_bw()
bat_plot 

ggsave(filename = paste0("figs/bat_daily_emergence_time_model_no_hibernation/emergence_time_model_trapping_weight.png"), plot = bat_plot, width = 10, height = 6, dpi = 300)



#eemeans for categorical variable ####
#age effect
m4_emm <- emmeans(model_3, ~ age, component = "cmean")

### pairwise comparisons
m4_contrast <- summary(contrast(m4_emm,
                                  method = "pairwise"),
                         infer = TRUE)

#Make into a dataframe
m4_emm_df <- data.frame(m4_emm)
m4_emm_contrast_df <- data.frame(m4_contrast)


# > m4_emm_df
#   age   emmean       SE  df asymp.LCL asymp.UCL
# 1   A 39.63865 2.837475 Inf  34.07730  45.19999
# 2   J 35.04174 3.831935 Inf  27.53128  42.55219
# 3  YA 57.45076 7.370149 Inf  43.00553  71.89599

# > m4_emm_contrast_df
#   contrast   estimate       SE  df  asymp.LCL asymp.UCL    z.ratio    p.value
# 1    A - J   4.596909 5.058640 Inf  -7.259028 16.452847  0.9087243 0.63477163
# 2   A - YA -17.812113 8.043531 Inf -36.663743  1.039516 -2.2144643 0.06874348
# 3   J - YA -22.409023 7.952519 Inf -41.047346 -3.770699 -2.8178522 0.01340975


#Change names
m4_emm_df$age <- factor(m4_emm_df$age , levels=c("J", "YA", "A"))
m4_emm_df$age <- plyr::revalue(m4_emm_df$age, c("J" = "Juvenile", "YA" = "Young Adult", "A" = "Adult"))

bat_active_data$age <- factor(bat_active_data$age , levels=c("J", "YA", "A"))
bat_active_data$age_rename <- plyr::revalue(bat_active_data$age, c("J" = "Juvenile", "YA" = "Young Adult", "A" = "Adult"))


#Make contrasting groups match variable names
m4_emm_contrast_df$group1 <- NA
m4_emm_contrast_df$group2 <- NA
m4_emm_contrast_df$contrast_pretty <- NA

m4_emm_contrast_df$group1 <- c("Adult", "Adult", "Juvenile")
m4_emm_contrast_df$group2 <- c("Juvenile", "Young Adult", "Young Adult")
m4_emm_contrast_df$contrast_pretty <- c("Adult - Juvenile", "Adult - Young Adult", "Juvenile - Young Adult")

m4_emm_contrast_df$p_rounded <- rstatix::p_round(m4_emm_contrast_df$p.value, digits = 2)
m4_emm_contrast_df$p_pretty <- rstatix::p_format(rstatix::p_round(m4_emm_contrast_df$p.value, digits = 2),
                                              accuracy = 1e-04, add.p = TRUE)

#plot the result
bat_plot <- ggplot(bat_active_data, aes(x = age_rename, y = relative_emergence_time)) +
  geom_jitter(alpha = 0.3,
              height = 0.1,
              width = 0.3) +
  geom_point(data = m4_emm_df,
             aes(y = emmean,
                 x = age), 
             size = 4,
             colour = "red") +
  geom_errorbar(data = m4_emm_df, 
                aes(y = emmean,
                    x = age,
                    ymin = asymp.LCL, 
                    ymax = asymp.UCL), 
                width = 0.05,
                colour = "red") +
  labs(y = "Relative Emergence Time \n (minutes from sunset)",
       x = "Age") +
  theme_pubr()
bat_plot


bat_effects <- ggplot(data = m4_emm_contrast_df, aes(x = estimate, y = contrast_pretty)) +
  geom_point(size = 4) +
  xlab("Effects (Emergence Time)") +
  geom_errorbar(aes(x = estimate,
                    xmin = asymp.LCL,
                    xmax = asymp.UCL),
                width = 0.1) +
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "blue") +
  annotate(geom = "text",
           x = m4_emm_contrast_df$estimate,
           y = 1:3 + 0.2,
           label = m4_emm_contrast_df$p_pretty,
           size = 3) +
  theme_pubr() +
  rremove("ylab")
bat_effects


# make int one plot
m4_bottom <- bat_plot +
  coord_flip() # rotate to horizontal

m4_top <- bat_effects + # assign to new plot object
  scale_x_continuous(position="top") # move x-axis to "top"

m4_plot <- cowplot::plot_grid(m4_top,
                     m4_bottom,
                     nrow = 2,
                     align = "v",
                     axis = "lr",
                     rel_heights = c(1, 1))
m4_plot

ggsave(filename = paste0("figs/bat_daily_emergence_time_model_no_hibernation/emergence_time_model_trapping_age.png"), plot = bat_plot, width = 10, height = 6, dpi = 300)

ggsave(filename = paste0("figs/bat_daily_emergence_time_model_no_hibernation/emergence_time_model_trapping_age_effect.png"), plot = m4_plot, width = 10, height = 6, dpi = 300)


#sex effect
m4_emm <- emmeans(model_3, ~ sex, component = "cmean")

### pairwise comparisons
m4_contrast <- summary(contrast(m4_emm,
                                  method = "pairwise"),
                         infer = TRUE)

#Make into a dataframe
m4_emm_df <- data.frame(m4_emm)
m4_emm_contrast_df <- data.frame(m4_contrast)


# > m4_emm_df
#   sex   emmean       SE  df asymp.LCL asymp.UCL
# 1   F 40.63369 3.398292 Inf  33.97316  47.29422
# 2   M 47.45373 4.397061 Inf  38.83565  56.07181

# > m4_emm_contrast_df
#   contrast  estimate       SE  df asymp.LCL asymp.UCL  z.ratio  p.value
# 1    F - M -6.820037 5.221321 Inf -17.05364  3.413565 -1.30619 0.191488


#Change names
m4_emm_df$sex <- factor(m4_emm_df$sex , levels=c("F", "M"))
m4_emm_df$sex <- plyr::revalue(m4_emm_df$sex, c("F" = "Female", "M" = "Male"))

bat_active_data$sex <- factor(bat_active_data$sex , levels=c("F", "M"))
bat_active_data$sex_rename <- plyr::revalue(bat_active_data$sex, c("F" = "Female", "M" = "Male"))

#Make contrasting groups match variable names
m4_emm_contrast_df$group1 <- NA
m4_emm_contrast_df$group2 <- NA
m4_emm_contrast_df$contrast_pretty <- NA

m4_emm_contrast_df$group1 <- c("Female")
m4_emm_contrast_df$group2 <- c("Male")
m4_emm_contrast_df$contrast_pretty <- c("Female - Male")

m4_emm_contrast_df$p_rounded <- rstatix::p_round(m4_emm_contrast_df$p.value, digits = 2)
m4_emm_contrast_df$p_pretty <- rstatix::p_format(rstatix::p_round(m4_emm_contrast_df$p.value, digits = 2),
                                              accuracy = 1e-04, add.p = TRUE)

#plot the result
bat_plot <- ggplot(bat_active_data, aes(x = sex_rename, y = relative_emergence_time)) +
  geom_jitter(alpha = 0.3,
              height = 0.1,
              width = 0.3) +
  geom_point(data = m4_emm_df,
             aes(y = emmean,
                 x = sex), 
             size = 4,
             colour = "red") +
  geom_errorbar(data = m4_emm_df, 
                aes(y = emmean,
                    x = sex,
                    ymin = asymp.LCL, 
                    ymax = asymp.UCL), 
                width = 0.05,
                colour = "red") +
  labs(y = "Relative Emergence Time \n (minutes from sunset)",
       x = "Sex") +
  theme_pubr()
bat_plot


bat_effects <- ggplot(data = m4_emm_contrast_df, aes(x = estimate, y = contrast_pretty)) +
  geom_point(size = 4) +
  xlab("Effects (Emergence Time)") +
  geom_errorbar(aes(x = estimate,
                    xmin = asymp.LCL,
                    xmax = asymp.UCL),
                width = 0.1) +
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "blue") +
  annotate(geom = "text",
           x = m4_emm_contrast_df$estimate,
           y = 1 + 0.2,
           label = m4_emm_contrast_df$p_pretty,
           size = 3) +
  theme_pubr() +
  rremove("ylab")
bat_effects


# make int one plot
m4_bottom <- bat_plot +
  coord_flip() # rotate to horizontal

m4_top <- bat_effects + # assign to new plot object
  scale_x_continuous(position="top") # move x-axis to "top"

m4_plot <- cowplot::plot_grid(m4_top,
                     m4_bottom,
                     nrow = 2,
                     align = "v",
                     axis = "lr",
                     rel_heights = c(1, 1))
m4_plot

ggsave(filename = paste0("figs/bat_daily_emergence_time_model_no_hibernation/emergence_time_model_trapping_sex.png"), plot = bat_plot, width = 10, height = 6, dpi = 300)

ggsave(filename = paste0("figs/bat_daily_emergence_time_model_no_hibernation/emergence_time_model_trapping_sex_effect.png"), plot = m4_plot, width = 10, height = 6, dpi = 300)


#year effect
m4_emm <- emmeans(model_3, ~ study_year, component = "cmean")

### pairwise comparisons
m4_contrast <- summary(contrast(m4_emm,
                                  method = "pairwise"),
                         infer = TRUE)

#Make into a dataframe
m4_emm_df <- data.frame(m4_emm)
m4_emm_contrast_df <- data.frame(m4_contrast)


# > m4_emm_df
#      study_year   emmean       SE  df asymp.LCL asymp.UCL
# 1  First Winter 42.73997 3.867344 Inf  35.16011  50.31982
# 2 Second Winter 45.34746 3.199767 Inf  39.07603  51.61889

# > m4_emm_contrast_df
#                       contrast estimate       SE  df asymp.LCL asymp.UCL    z.ratio   p.value
# 1 First Winter - Second Winter -2.60749 3.985829 Inf -10.41957  5.204591 -0.6541902 0.5129893


#Change names
m4_emm_df$study_year <- factor(m4_emm_df$study_year , levels=c("First Winter", "Second Winter"))

bat_active_data$study_year <- factor(bat_active_data$study_year , levels=c("First Winter", "Second Winter"))

#Make contrasting groups match variable names
m4_emm_contrast_df$group1 <- NA
m4_emm_contrast_df$group2 <- NA
m4_emm_contrast_df$contrast_pretty <- NA

m4_emm_contrast_df$group1 <-c("First Winter")
m4_emm_contrast_df$group2 <- c("Second Winter")
m4_emm_contrast_df$contrast_pretty <- c("First - Second")

m4_emm_contrast_df$p_rounded <- rstatix::p_round(m4_emm_contrast_df$p.value, digits = 2)
m4_emm_contrast_df$p_pretty <- rstatix::p_format(rstatix::p_round(m4_emm_contrast_df$p.value, digits = 2),
                                              accuracy = 1e-04, add.p = TRUE)

MyData <- crossing(age = "A",
                   sex = "F",
                   study_year = levels(bat_active_data$study_year),
                   week_integer = 14,
                   weight_std = mean(bat_active_data$weight_std),
                   av_wind_speed_std = mean(bat_active_data$av_wind_speed_std),
                   daily_rain_std = mean(bat_active_data$daily_rain_std),
                   av_air_temp_std = mean(bat_active_data$av_air_temp_std),
                   av_rel_humidity_std = mean(bat_active_data$av_rel_humidity_std),
                   av_cloud_okta_std = mean(bat_active_data$av_cloud_okta_std),
                   meanMoonlightIntensity = mean(bat_active_data$meanMoonlightIntensity))

#' dispersion parameter (σ)
Fit_mod <- predict(model_3,
                   newdata = MyData,
                   re.form = NA,
                   type = "disp")

attr(Fit_mod, "names") <- NULL

MyData$sigma   <- Fit_mod

# calculate variance ± SD
m4_emm_df <- m4_emm_df %>%
  left_join(MyData %>% select(study_year, sigma), by = "study_year") %>%
  mutate(
    variance = (emmean^2) * (sigma^2),
    SD_plus = emmean + 2*sqrt(variance),
    SD_minus = emmean - 2*sqrt(variance)
  )

m4_emm_df$study_year <- factor(m4_emm_df$study_year , levels=c("First Winter", "Second Winter"))


#plot the result
bat_plot <- ggplot(bat_active_data, aes(x = study_year, y = relative_emergence_time)) +
  geom_jitter(alpha = 0.3,
              height = 0.1,
              width = 0.3) +
  geom_point(data = m4_emm_df,
             aes(y = emmean,
                 x = study_year), 
             size = 4,
             colour = "red") +
  geom_errorbar(data = m4_emm_df, 
                aes(y = emmean,
                    x = study_year,
                    ymin = asymp.LCL, 
                    ymax = asymp.UCL), 
                width = 0.05,
                colour = "red") +
  geom_segment(data = m4_emm_df,
               aes(x = as.numeric(study_year) -  0.3,
                   xend = as.numeric(study_year) +  0.3,
                   y = SD_minus, yend = SD_minus),
             colour = "blue",
             linetype = "dashed",
             linewidth = 0.8,
             inherit.aes = FALSE) +
geom_segment(data = m4_emm_df,
             aes(x = as.numeric(study_year) -  0.3,
                 xend = as.numeric(study_year) +  0.3,
                 y = SD_plus, yend = SD_plus),
             colour = "blue",
             linetype = "dashed",
             linewidth = 0.8,
             inherit.aes = FALSE) +
  labs(y = "Relative Emergence Time \n (minutes from sunset)",
       x = "Study Year") +
  theme_pubr()
bat_plot


bat_effects <- ggplot(data = m4_emm_contrast_df, aes(x = estimate, y = contrast_pretty)) +
  geom_point(size = 4) +
  xlab("Effects (Emergence Time)") +
  geom_errorbar(aes(x = estimate,
                    xmin = asymp.LCL,
                    xmax = asymp.UCL),
                width = 0.1) +
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "blue") +
  annotate(geom = "text",
           x = m4_emm_contrast_df$estimate,
           y = 1 + 0.2,
           label = m4_emm_contrast_df$p_pretty,
           size = 3) +
  theme_pubr() +
  rremove("ylab")
bat_effects


# make int one plot
m4_bottom <- bat_plot +
  coord_flip() # rotate to horizontal

m4_top <- bat_effects + # assign to new plot object
  scale_x_continuous(position="top") # move x-axis to "top"

m4_plot <- cowplot::plot_grid(m4_top,
                     m4_bottom,
                     nrow = 2,
                     align = "v",
                     axis = "lr",
                     rel_heights = c(1, 1))
m4_plot

ggsave(filename = paste0("figs/bat_daily_emergence_time_model_no_hibernation/emergence_time_model_trapping_year.png"), plot = bat_plot, width = 10, height = 6, dpi = 300)

ggsave(filename = paste0("figs/bat_daily_emergence_time_model_no_hibernation/emergence_time_model_trapping_year_effect.png"), plot = m4_plot, width = 10, height = 6, dpi = 300)


#week effect
m4_emm <- emmeans(model_3, ~ week_integer, component = "cmean")

### pairwise comparisons
m4_contrast <- summary(contrast(m4_emm,
                                  method = "pairwise"),
                         infer = TRUE)

#Make into a dataframe
m4_emm_df <- data.frame(m4_emm)
m4_emm_contrast_df <- data.frame(m4_contrast)

MyData <- crossing(age = "A",
                   sex = "F",
                   study_year = "Second Winter",
                   week_integer = levels(bat_active_data$week_integer),
                   weight_std = mean(bat_active_data$weight_std),
                   av_wind_speed_std = mean(bat_active_data$av_wind_speed_std),
                   daily_rain_std = mean(bat_active_data$daily_rain_std),
                   av_air_temp_std = mean(bat_active_data$av_air_temp_std),
                   av_rel_humidity_std = mean(bat_active_data$av_rel_humidity_std),
                   av_cloud_okta_std = mean(bat_active_data$av_cloud_okta_std),
                   meanMoonlightIntensity = mean(bat_active_data$meanMoonlightIntensity))

#' dispersion parameter (σ)
Fit_mod <- predict(model_3,
                   newdata = MyData,
                   re.form = NA,
                   type = "disp")

attr(Fit_mod, "names") <- NULL

MyData$sigma   <- Fit_mod

# calculate variance ± SD
m4_emm_df <- m4_emm_df %>%
  left_join(MyData %>% select(week_integer, sigma), by = "week_integer") %>%
  mutate(
    variance = (emmean^2) * (sigma^2),
    SD_plus = emmean + 2*sqrt(variance),
    SD_minus = emmean - 2*sqrt(variance)
  )

#Change names
m4_emm_df$week_integer <- factor(m4_emm_df$week_integer , levels=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16"))
m4_emm_df$week_integer <- plyr::revalue(m4_emm_df$week_integer, c("1" = "38", "2" = "39", "3" = "40", "4" = "41", "5" = "42", "6" = "43", "7" = "44", "8" = "45", "9" = "46", "10" = "47", "11" = "48", "12" = "49", "13" = "50", "14" = "51", "15" = "52", "16" = "1"))

bat_active_data$week_integer <- factor(bat_active_data$week_integer , levels=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16"))
bat_active_data$week_integer_rename <- plyr::revalue(bat_active_data$week_integer, c("1" = "38", "2" = "39", "3" = "40", "4" = "41", "5" = "42", "6" = "43", "7" = "44", "8" = "45", "9" = "46", "10" = "47", "11" = "48", "12" = "49", "13" = "50", "14" = "51", "15" = "52", "16" = "1"))
                                                                                                                                                    

#Make contrasting groups match variable names
m4_emm_contrast_df$group1 <- NA
m4_emm_contrast_df$group2 <- NA
m4_emm_contrast_df$contrast_pretty <- NA

m4_emm_contrast_df$p_rounded <- rstatix::p_round(m4_emm_contrast_df$p.value, digits = 2)
m4_emm_contrast_df$p_pretty <- rstatix::p_format(rstatix::p_round(m4_emm_contrast_df$p.value, digits = 2),
                                              accuracy = 1e-04, add.p = TRUE)

#remove non-significant differences
m4_emm_contrast_df <- filter(m4_emm_contrast_df, p_rounded <= 0.05)

m4_emm_contrast_df$group1 <- c("44", "44", "44", "45", "46", "46", "48")
m4_emm_contrast_df$group2 <- c("45", "50", "52", "46", "50", "52", "50")
m4_emm_contrast_df$contrast_pretty <- c("44 - 45", "44 - 50", "44 - 52", "45 - 46", "46 - 50", "46 - 52", "48 - 50")


#plot the result
bat_plot <- ggplot(bat_active_data, aes(x = week_integer_rename, y = relative_emergence_time)) +
  geom_jitter(alpha = 0.3,
              height = 0.1,
              width = 0.3) +
  geom_point(data = m4_emm_df,
             aes(y = emmean,
                 x = week_integer), 
             size = 4,
             colour = "red") +
  geom_errorbar(data = m4_emm_df, 
                aes(y = emmean,
                    x = week_integer,
                    ymin = asymp.LCL, 
                    ymax = asymp.UCL), 
                width = 0.05,
                colour = "red") +
  geom_segment(data = m4_emm_df,
             aes(x = as.numeric(week_integer) -  0.3,
                 xend = as.numeric(week_integer) +  0.3,
                 y = SD_minus, yend = SD_minus),
             colour = "blue",
             linetype = "dashed",
             linewidth = 0.8,
             inherit.aes = FALSE) +
  geom_segment(data = m4_emm_df,
             aes(x = as.numeric(week_integer) -  0.3,
                 xend = as.numeric(week_integer) +  0.3,
                 y = SD_plus, yend = SD_plus),
             colour = "blue",
             linetype = "dashed",
             linewidth = 0.8,
             inherit.aes = FALSE) +
  ylim(c(0, (max(bat_active_data$relative_emergence_time)+15))) +
  labs(y = "Relative Emergence Time \n (minutes from sunset)",
       x = "Week in Year") +
  theme_pubr()
bat_plot


bat_effects <- ggplot(data = m4_emm_contrast_df, aes(x = estimate, y = contrast_pretty)) +
  geom_point(size = 4) +
  xlab("Effects (Emergence Time)") +
  geom_errorbar(aes(x = estimate,
                    xmin = asymp.LCL,
                    xmax = asymp.UCL),
                width = 0.1) +
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "blue") +
  annotate(geom = "text",
           x = m4_emm_contrast_df$estimate,
           y = 1:7 + 0.2,
           label = m4_emm_contrast_df$p_pretty,
           size = 3) +
  theme_pubr() +
  rremove("ylab")
bat_effects


# make int one plot
m4_bottom <- bat_plot +
  coord_flip() # rotate to horizontal

m4_top <- bat_effects + # assign to new plot object
  scale_x_continuous(position="top") # move x-axis to "top"

m4_plot <- cowplot::plot_grid(m4_top,
                     m4_bottom,
                     nrow = 2,
                     align = "v",
                     axis = "lr",
                     rel_heights = c(1, 1))
m4_plot

ggsave(filename = paste0("figs/bat_daily_emergence_time_model_no_hibernation/emergence_time_model_trapping_week.png"), plot = bat_plot, width = 10, height = 6, dpi = 300)

ggsave(filename = paste0("figs/bat_daily_emergence_time_model_no_hibernation/emergence_time_model_trapping_week_effeft.png"), plot = m4_plot, width = 10, height = 6, dpi = 300)
