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
bat_metrics <- mutate(bat_metrics,
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

#drop unused levels
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


#set priors for penalised regression
cprior <- data.frame(prior = rep("normal(0,3)", 2),
                     class = rep("fixef", 2),
                     coef = c("(Intercept)", ""))



# fitting model ####
print("Fitting Model")
#most complex model
model_1 <- glmmTMB(bat_emerged ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + weight_std + sex + age + survey_type + week_integer + study_year + sex*age + weight_std*av_air_temp_std + av_wind_speed_std*av_air_temp_std + (1 | ring_no),
                   data = bat_metrics,
                   family = binomial(link = "logit"),
                   priors = cprior) #add priors for penalised regression due to high 0's in some categories

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
plotResiduals(model_1_simulationOutput, form = bat_metrics$survey_type, xlab = "Rank-transformed Survey Type")


#make simpler mdoels without interactions
# sex*age
model_1.c <- glmmTMB(bat_emerged ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + weight_std + sex + age + survey_type + week_integer + study_year + weight_std*av_air_temp_std + av_wind_speed_std*av_air_temp_std + (1 | ring_no),
                   data = bat_metrics,
                   family = binomial(link = "logit"),
                   priors = cprior)

# weight_std*av_air_temp_std
model_1.d <- glmmTMB(bat_emerged ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + weight_std + sex + age + survey_type + week_integer + study_year + sex*age + av_wind_speed_std*av_air_temp_std + (1 | ring_no),
                   data = bat_metrics,
                   family = binomial(link = "logit"),
                   priors = cprior)

# av_wind_speed_std*av_air_temp_std
model_1.e <- glmmTMB(bat_emerged ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + weight_std + sex + age + survey_type + week_integer + study_year + sex*age + weight_std*av_air_temp_std + (1 | ring_no),
                   data = bat_metrics,
                   family = binomial(link = "logit"),
                   priors = cprior)


# test interactions
anova(model_1, model_1.c)
#           Df    AIC    BIC  logLik deviance Chisq Chi Df Pr(>Chisq)
# model_1.c 44 2042.7 2280.8 -977.36   1954.7
# model_1   46 2052.0 2300.9 -979.99   1960.0     0      2          1

anova(model_1, model_1.d)
#           Df    AIC    BIC  logLik deviance Chisq Chi Df Pr(>Chisq)
# model_1.d 45 2048 2291.4 -978.98     1958
# model_1   46 2052 2300.9 -979.99     1960     0      1          1

anova(model_1, model_1.e)
#           Df    AIC    BIC  logLik deviance Chisq Chi Df Pr(>Chisq)
# model_1.e 45 2046.5 2290.0 -978.23   1956.5
# model_1   46 2052.0 2300.9 -979.99   1960.0     0      1          1



#drop1 model selection ####
print("Dropping")
#model selection
drop1(model_1, test = "Chi")

#                                   Df    AIC      LRT Pr(>Chi)
# <none>                               2052.0
# meanMoonlightIntensity             1 2046.0  -4.0306   1.0000
# av_cloud_okta_std                  1 2050.3   0.3378   0.5611
# daily_rain_std                     1 2047.1  -2.9042   1.0000
# av_rel_humidity_std                1 2047.8  -2.1543   1.0000
# survey_type                        1 2049.4  -0.5362   1.0000
# week_integer                      28 1964.6 -31.3593   1.0000
# study_year                         1 2046.0  -3.9926   1.0000
# sex:age                            2 2042.7  -5.2596   1.0000
# av_air_temp_std:weight_std         1 2048.0  -2.0172   1.0000
# av_air_temp_std:av_wind_speed_std  1 2046.5  -3.5075   1.0000


# drop sex:age interaction
drop1(model_1.c, test = "Chi")

#                                   Df    AIC      LRT Pr(>Chi)
# <none>                               2042.7
# meanMoonlightIntensity             1 2036.7  -4.0279   1.0000
# av_cloud_okta_std                  1 2040.9   0.1781   0.6730
# daily_rain_std                     1 2037.8  -2.8909   1.0000
# av_rel_humidity_std                1 2038.7  -2.0545   1.0000
# sex                                1 2036.8  -3.9022   1.0000
# age                                2 2031.7  -7.0146   1.0000
# survey_type                        1 2042.3   1.5845   0.2081
# week_integer                      28 1957.3 -29.4416   1.0000
# study_year                         1 2036.9  -3.8585   1.0000
# av_air_temp_std:weight_std         1 2038.7  -2.0454   1.0000
# av_air_temp_std:av_wind_speed_std  1 2037.2  -3.5064   1.0000


# drop av_air_temp_std:av_wind_speed_std interaction
model_2 <- glmmTMB(bat_emerged ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + weight_std + sex + age + survey_type + week_integer + study_year + weight_std*av_air_temp_std + (1 | ring_no),
                   data = bat_metrics,
                   family = binomial(link = "logit"),
                   priors = cprior)

#model selection
drop1(model_2, test = "Chi")

#                            Df    AIC     LRT  Pr(>Chi)
# <none>                        2037.2
# av_wind_speed_std           1 2081.1  45.857 1.272e-11 ***
# meanMoonlightIntensity      1 2031.2  -3.976    1.0000
# av_cloud_okta_std           1 2035.2   0.011    0.9178
# daily_rain_std              1 2032.4  -2.766    1.0000
# av_rel_humidity_std         1 2033.5  -1.707    1.0000
# sex                         1 2031.3  -3.910    1.0000
# age                         2 2026.2  -7.002    1.0000
# survey_type                 1 2036.8   1.609    0.2046
# week_integer               28 1951.3 -29.917    1.0000
# study_year                  1 2031.3  -3.865    1.0000
# av_air_temp_std:weight_std  1 2033.1  -2.109    1.0000


# drop av_air_temp_std:weight_std interaction
model_3 <- glmmTMB(bat_emerged ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + weight_std + sex + age + survey_type + week_integer + study_year + (1 | ring_no),
                   data = bat_metrics,
                   family = binomial(link = "logit"),
                   priors = cprior)

#model selection
drop1(model_3, test = "Chi")

#                        Df    AIC     LRT  Pr(>Chi)
# <none>                    2033.1
# av_air_temp_std         1 2114.3  83.218 < 2.2e-16 ***
# av_wind_speed_std       1 2077.5  46.440 9.446e-12 ***
# meanMoonlightIntensity  1 2027.1  -3.974    1.0000
# av_cloud_okta_std       1 2031.3   0.204    0.6513
# daily_rain_std          1 2028.4  -2.681    1.0000
# av_rel_humidity_std     1 2029.5  -1.594    1.0000
# weight_std              1 2027.9  -3.217    1.0000
# sex                     1 2027.1  -3.967    1.0000
# age                     2 2022.2  -6.930    1.0000
# survey_type             1 2033.4   2.316    0.1280
# week_integer           28 1947.3 -29.785    1.0000
# study_year              1 2027.2  -3.903    1.0000



#check model fits ####
summary(model_3)
#  Family: binomial  ( logit )
# Formula:
# bat_emerged ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity +
#     av_cloud_okta_std + daily_rain_std + av_rel_humidity_std +
#     weight_std + sex + age + survey_type + week_integer + study_year +
#     (1 | ring_no)
# Data: bat_metrics

#       AIC       BIC    logLik -2*log(L)  df.resid
#    2033.1    2260.4    -974.6    1949.1      1612

# Random effects:

# Conditional model:
#  Groups  Name        Variance Std.Dev.
#  ring_no (Intercept) 0.5312   0.7289
# Number of obs: 1654, groups:  ring_no, 42

# Conditional model:
#                          Estimate Std. Error z value Pr(>|z|)
# (Intercept)              2.213595   1.355863   1.633   0.1026
# av_air_temp_std          1.192012   0.135877   8.773  < 2e-16 ***
# av_wind_speed_std       -4.901449   0.709619  -6.907 4.94e-12 ***
# meanMoonlightIntensity  -0.193640   0.784487  -0.247   0.8050
# av_cloud_okta_std       -0.640362   0.311784  -2.054   0.0400 *
# daily_rain_std           0.489457   0.422133   1.159   0.2463
# av_rel_humidity_std     -1.972320   1.262981  -1.562   0.1184
# weight_std              -0.203740   0.223652  -0.911   0.3623
# sexM                     0.090018   0.345383   0.261   0.7944
# ageJ                     0.015787   0.419566   0.038   0.9700
# ageYA                   -0.430267   0.440055  -0.978   0.3282
# survey_typetrapping      1.619131   0.645986   2.506   0.0122 *
# week_integer2           -1.763464   0.451262  -3.908 9.31e-05 ***
# week_integer3           -1.137530   0.480820  -2.366   0.0180 *
# week_integer4           -1.301389   0.598766  -2.173   0.0297 *
# week_integer5           -0.114222   0.800564  -0.143   0.8865
# week_integer6            0.139209   0.593321   0.235   0.8145
# week_integer7            0.631486   0.536527   1.177   0.2392
# week_integer8            0.504438   0.473380   1.066   0.2866
# week_integer9           -0.008588   0.444415  -0.019   0.9846
# week_integer10           0.261697   0.449823   0.582   0.5607
# week_integer11           0.018631   0.489490   0.038   0.9696
# week_integer12           0.378255   0.473592   0.799   0.4245
# week_integer13          -0.786964   0.468398  -1.680   0.0929 .
# week_integer14          -0.371994   0.473289  -0.786   0.4319
# week_integer15          -1.073741   0.510306  -2.104   0.0354 *
# week_integer16          -0.389071   0.550532  -0.707   0.4797
# week_integer17          -2.213005   1.742173  -1.270   0.2040
# week_integer18          -2.008405   1.789448  -1.122   0.2617
# week_integer19          -1.016012   1.125760  -0.903   0.3668
# week_integer20          -0.402847   0.944367  -0.427   0.6697
# week_integer21          -0.701528   1.114152  -0.630   0.5289
# week_integer22           0.272705   0.827504   0.330   0.7417
# week_integer23          -0.683760   1.111621  -0.615   0.5385
# week_integer24           0.868231   0.836260   1.038   0.2992
# week_integer25           0.537574   0.878635   0.612   0.5407
# week_integer26          -0.471509   1.161328  -0.406   0.6847
# week_integer27           0.022815   1.255066   0.018   0.9855
# week_integer28          -1.704449   1.912718  -0.891   0.3729
# week_integer29          -0.021384   1.259096  -0.017   0.9864
# study_yearSecond Winter -0.099552   0.273371  -0.364   0.7157
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Priors:
# fixef((Intercept)) ~ normal(0, 3)
# fixef ~ normal(0, 3)


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
plotResiduals(model_3_simulationOutput, form = bat_metrics$meanMoonlightIntensity, xlab = "Rank-transformed Standardised Moonlight Intensity")
plotResiduals(model_3_simulationOutput, form = bat_metrics$av_cloud_okta_std, xlab = "Rank-transformed Standardised Cloud Cover") 
plotResiduals(model_3_simulationOutput, form = bat_metrics$daily_rain_std, xlab = "Rank-transformed Standardised Rainfall") # slight deviation in 0.25 quartile
plotResiduals(model_3_simulationOutput, form = bat_metrics$av_rel_humidity_std, xlab = "Rank-transformed Standardised Humidity")
plotResiduals(model_3_simulationOutput, form = bat_metrics$weight_std, xlab = "Rank-transformed Standardised Weight")
plotResiduals(model_3_simulationOutput, form = bat_metrics$age, xlab = "Rank-transformed Age")
plotResiduals(model_3_simulationOutput, form = bat_metrics$sex, xlab = "Rank-transformed Sex")
plotResiduals(model_3_simulationOutput, form = bat_metrics$week_integer, xlab = "Rank-transformed Week")
plotResiduals(model_3_simulationOutput, form = bat_metrics$study_year, xlab = "Rank-transformed Study Year")
plotResiduals(model_3_simulationOutput, form = bat_metrics$survey_type, xlab = "Rank-transformed Survey Type")



#Test significance of main effects ####
model_3.b <- glmmTMB(bat_emerged ~ av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + weight_std + sex + age + survey_type + week_integer + study_year + (1 | ring_no),
                   data = bat_metrics,
                   family = binomial(link = "logit"),
                   priors = cprior)

model_3.c <- glmmTMB(bat_emerged ~ av_air_temp_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + weight_std + sex + age + survey_type + week_integer + study_year + (1 | ring_no),
                   data = bat_metrics,
                   family = binomial(link = "logit"),
                   priors = cprior)

model_3.d <- glmmTMB(bat_emerged ~ av_air_temp_std + av_wind_speed_std + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + weight_std + sex + age + survey_type + week_integer + study_year + (1 | ring_no),
                   data = bat_metrics,
                   family = binomial(link = "logit"),
                   priors = cprior)

model_3.e <- glmmTMB(bat_emerged ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + daily_rain_std + av_rel_humidity_std + weight_std + sex + age + survey_type + week_integer + study_year + (1 | ring_no),
                   data = bat_metrics,
                   family = binomial(link = "logit"),
                   priors = cprior)

model_3.f <- glmmTMB(bat_emerged ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + av_rel_humidity_std + weight_std + sex + age + survey_type + week_integer + study_year + (1 | ring_no),
                   data = bat_metrics,
                   family = binomial(link = "logit"),
                   priors = cprior)

model_3.g <- glmmTMB(bat_emerged ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + weight_std + sex + age + survey_type + week_integer + study_year + (1 | ring_no),
                   data = bat_metrics,
                   family = binomial(link = "logit"),
                   priors = cprior)

model_3.h <- glmmTMB(bat_emerged ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + sex + age + survey_type + week_integer + study_year + (1 | ring_no),
                   data = bat_metrics,
                   family = binomial(link = "logit"),
                   priors = cprior)

model_3.i <- glmmTMB(bat_emerged ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + weight_std + age + survey_type + week_integer + study_year + (1 | ring_no),
                   data = bat_metrics,
                   family = binomial(link = "logit"),
                   priors = cprior)

model_3.j <- glmmTMB(bat_emerged ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + weight_std + sex + survey_type + week_integer + study_year + (1 | ring_no),
                   data = bat_metrics,
                   family = binomial(link = "logit"),
                   priors = cprior)

model_3.k <- glmmTMB(bat_emerged ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + weight_std + sex + age + week_integer + study_year + (1 | ring_no),
                   data = bat_metrics,
                   family = binomial(link = "logit"),
                   priors = cprior)

model_3.l <- glmmTMB(bat_emerged ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + weight_std + sex + age + survey_type + study_year + (1 | ring_no),
                   data = bat_metrics,
                   family = binomial(link = "logit"),
                   priors = cprior)

model_3.m <- glmmTMB(bat_emerged ~ av_air_temp_std + av_wind_speed_std + meanMoonlightIntensity + av_cloud_okta_std + daily_rain_std + av_rel_humidity_std + weight_std + sex + age + survey_type + week_integer + (1 | ring_no),
                   data = bat_metrics,
                   family = binomial(link = "logit"),
                   priors = cprior)


anova(model_3, model_3.b) #air temp
#           Df    AIC    BIC   logLik deviance  Chisq Chi Df Pr(>Chisq)
# model_3.b 41 2114.3 2336.2 -1016.16   2032.3
# model_3   42 2033.1 2260.4  -974.55   1949.1 83.218      1  < 2.2e-16 ***
  
anova(model_3, model_3.c) #wind speed
#           Df    AIC    BIC  logLik deviance Chisq Chi Df Pr(>Chisq)
# model_3.c 41 2077.5 2299.4 -997.77   1995.5
# model_3   42 2033.1 2260.4 -974.55   1949.1 46.44      1  9.446e-12 ***
  
anova(model_3, model_3.d) #moonlight
#           Df    AIC    BIC  logLik deviance Chisq Chi Df Pr(>Chisq)
# model_3.d 41 2027.1 2249.0 -972.56   1945.1
# model_3   42 2033.1 2260.4 -974.55   1949.1     0      1          1

anova(model_3, model_3.e) #cloud
#           Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# model_3.e 41 2031.3 2253.2 -974.65   1949.3
# model_3   42 2033.1 2260.4 -974.55   1949.1 0.2042      1     0.6513

anova(model_3, model_3.f) #rain
#           Df    AIC    BIC  logLik deviance Chisq Chi Df Pr(>Chisq)
# model_3.f 41 2028.4 2250.3 -973.21   1946.4
# model_3   42 2033.1 2260.4 -974.55   1949.1     0      1          1

anova(model_3, model_3.g) #humidity
#           Df    AIC    BIC  logLik deviance Chisq Chi Df Pr(>Chisq)
# model_3.g 41 2029.5 2251.4 -973.75   1947.5
# model_3   42 2033.1 2260.4 -974.55   1949.1     0      1          1

anova(model_3, model_3.h) #weight
#           Df    AIC    BIC  logLik deviance Chisq Chi Df Pr(>Chisq)
# model_3.h 41 2027.9 2249.7 -972.94   1945.9
# model_3   42 2033.1 2260.4 -974.55   1949.1     0      1          1

anova(model_3, model_3.i) #sex
#           Df    AIC    BIC  logLik deviance Chisq Chi Df Pr(>Chisq)
# model_3.i 41 2027.1 2249.0 -972.57   1945.1
# model_3   42 2033.1 2260.4 -974.55   1949.1     0      1          1

anova(model_3, model_3.j) #age
#           Df    AIC    BIC  logLik deviance Chisq Chi Df Pr(>Chisq)
# model_3.j 40 2022.2 2238.6 -971.09   1942.2
# model_3   42 2033.1 2260.4 -974.55   1949.1     0      2          1

anova(model_3, model_3.k) #survey type
#           Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# model_3.k 41 2033.4 2255.3 -975.71   1951.4
# model_3   42 2033.1 2260.4 -974.55   1949.1 2.3164      1      0.128

anova(model_3, model_3.l) #week
#           Df    AIC    BIC  logLik deviance Chisq Chi Df Pr(>Chisq)
# model_3.l 14 1947.3 2023.1 -959.66   1919.3
# model_3   42 2033.1 2260.4 -974.55   1949.1     0     28          1

anova(model_3, model_3.m) #year
#           Df    AIC    BIC  logLik deviance Chisq Chi Df Pr(>Chisq)
# model_3.m 41 2027.2 2249.1 -972.60   1945.2
# model_3   42 2033.1 2260.4 -974.55   1949.1     0      1          1



# plotting ####
bat_plot <- ggplot(filter(bat_metrics, study_year == "First Winter"), aes(x = date, y = ring_no, colour = as.factor(bat_emerged))) +
  geom_point() +
  ggforce::facet_col(vars(survey_type),
                     scales = "free_y",
                     space="free") +
  labs(y = "Bat ID",
       x = "Date",
       colour = "Bat Emerged") +
  scale_color_manual(values = c("#56B4E9", "#D55E00")) +
  theme_bw() + theme(legend.position = "none")


ggsave(filename = paste0("figs/bat_daily_emergence_model_all_bats/emergence_model_all_hibernation_one_no_legend.png"), plot = bat_plot, width = 10, height = 6, dpi = 300)


bat_plot <- ggplot(filter(bat_metrics, study_year == "First Winter"), aes(x = date, y = ring_no, colour = as.factor(bat_emerged))) +
  geom_point() +
  ggforce::facet_col(vars(survey_type),
                     scales = "free_y",
                     space="free") +
  labs(y = "Bat ID",
       x = "Date",
       colour = "Bat Emerged") +
  scale_color_manual(values = c("#56B4E9", "#D55E00")) +
  theme_bw()


 ggsave(filename = paste0("figs/bat_daily_emergence_model_all_bats/emergence_model_all_hibernation_one.png"), plot = bat_plot, width = 10, height = 6, dpi = 300)


bat_plot <- ggplot(filter(bat_metrics, study_year == "Second Winter"), aes(x = date, y = ring_no, colour = as.factor(bat_emerged))) +
  geom_point() +
  ggforce::facet_col(vars(survey_type),
                     scales = "free_y",
                     space="free") +
  labs(y = "Bat ID",
       x = "Date",
       colour = "Bat Emerged") +
  scale_color_manual(values = c("#56B4E9", "#D55E00")) +
  theme_bw()


ggsave(filename = paste0("figs/bat_daily_emergence_model_all_bats/emergence_model_all_hibernation_two.png"), plot = bat_plot, width = 10, height = 6, dpi = 300)


#av_air_temp_std
MyData_temp <- crossing(age = "A",
                        sex = "F",
                        study_year = "Second Winter",
                        survey_type = "trapping",
                        week_integer = 14,
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
MyData_temp$time <- NA

#' Calculate the predicted  values.
Fit_temp <- predict(model_3,
                    newdata = MyData_temp,
                    re.form = NULL,  #~0, Keep it to NULL..gives odd results with ~0
                    type = "link",
                    se = TRUE)

#' Extract the predicted values and standard errors, and calculate
#' standard errors (SE) for predicted values.
MyData_temp$Pi   <- exp(Fit_temp$fit) / (1 +  exp(Fit_temp$fit))
MyData_temp$SeUp <- exp(Fit_temp$fit + 1.96 * Fit_temp$se.fit) / (1 + exp(Fit_temp$fit + 1.96 * Fit_temp$se.fit))
MyData_temp$SeLo <- exp(Fit_temp$fit - 1.96 * Fit_temp$se.fit) / (1 +  exp(Fit_temp$fit - 1.96 * Fit_temp$se.fit))
MyData_temp$se   <- Fit_temp$se.fit
head(MyData_temp)

#back standardise
MyData_temp$av_air_temp <- MyData_temp$av_air_temp_std * sd(bat_metrics$av_air_temp) +
  mean(bat_metrics$av_air_temp)


#plot the result
bat_plot <- ggplot(MyData_temp, aes(x = av_air_temp, y = Pi)) +
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
       x = "Average Air Temperature (DegC)") +
  scale_color_viridis(discrete=TRUE) +
  scale_fill_viridis(discrete=TRUE) +
  theme_bw()


ggsave(filename = paste0("figs/bat_daily_emergence_model_all_bats/emergence_model_all_temp.png"), plot = bat_plot, width = 10, height = 6, dpi = 300)


#wind
MyData_wind <- crossing(age = "A",
                        sex = "F",
                        study_year = "Second Winter",
                        survey_type = "trapping",
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
Fit_wind <- predict(model_3,
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
              alpha = 0.03,
              show.legend = FALSE) +
  geom_line(linewidth = 1) +
  geom_jitter(data = bat_metrics,
              aes(y = bat_emerged),
              alpha = 0.3,
              height = 0.01,
              width = 0.01) +
  labs(y = "Probability of Emergence",
       x = "Average Wind Speed (kn)") +
  scale_color_viridis(discrete=TRUE) +
  scale_fill_viridis(discrete=TRUE) +
  theme_bw()


ggsave(filename = paste0("figs/bat_daily_emergence_model_all_bats/emergence_model_all_wind.png"), plot = bat_plot, width = 10, height = 6, dpi = 300)

#survey effect
m4_emm <- emmeans(model_3, ~ survey_type, component = "cmean")

### pairwise comparisons
m4_contrast <- summary(contrast(m4_emm,
                                 method = "pairwise"),
                        infer = TRUE)

#Make into a dataframe
m4_emm_df <- data.frame(m4_emm)
m4_emm_contrast_df <- data.frame(m4_contrast)


# > m4_emm_df
#   survey_type    emmean         SE  df  asymp.LCL asymp.UCL
# 1 hibernation 0.1418823 0.05690891 Inf 0.03034289 0.2534217
# 2    trapping 0.4174284 0.06807317 Inf 0.28400742 0.5508493

# > m4_emm_contrast_df
#                 contrast   estimate         SE  df asymp.LCL   asymp.UCL  z.ratio     p.value
# 1 hibernation - trapping -0.2755461 0.09632826 Inf -0.464346 -0.08674615 -2.86049 0.004229863


#Change names
m4_emm_df$survey_type <- factor(m4_emm_df$survey_type , levels=c("hibernation", "trapping"))
m4_emm_df$survey_type <- plyr::revalue(m4_emm_df$survey_type, c("hibernation" = "Hibernation", "trapping" = "Trapping"))

bat_metrics$survey_type <- factor(bat_metrics$survey_type , levels=c("hibernation", "trapping"))
bat_metrics$survey_type <- plyr::revalue(bat_metrics$survey_type, c("hibernation" = "Hibernation", "trapping" = "Trapping"))


#Make contrasting groups match variable names
m4_emm_contrast_df$group1 <- NA
m4_emm_contrast_df$group2 <- NA
m4_emm_contrast_df$contrast_pretty <- NA

m4_emm_contrast_df$group1 <- c("hibernation")
m4_emm_contrast_df$group2 <- c("trapping")
m4_emm_contrast_df$contrast_pretty <- c("Hibernation - Trapping")

m4_emm_contrast_df$p_rounded <- rstatix::p_round(m4_emm_contrast_df$p.value, digits = 2)
m4_emm_contrast_df$p_pretty <- rstatix::p_format(rstatix::p_round(m4_emm_contrast_df$p.value, digits = 2),
                                                 accuracy = 1e-04, add.p = TRUE)

#plot the result
bat_plot <- ggplot(bat_metrics, aes(x = survey_type, y = bat_emerged)) +
  geom_jitter(alpha = 0.3,
              height = 0.1,
              width = 0.3) +
  geom_point(data = m4_emm_df,
             aes(y = emmean,
                 x = survey_type), 
             size = 4,
             colour = "red") +
  geom_errorbar(data = m4_emm_df, 
                aes(y = emmean,
                    x = survey_type,
                    ymin = asymp.LCL, 
                    ymax = asymp.UCL), 
                width = 0.05,
                colour = "red") +
  labs(y = "Probability of Emergence",
       x = "Survey Type") +
  ggpubr::theme_pubr()
bat_plot


bat_effects <- ggplot(data = m4_emm_contrast_df, aes(x = estimate, y = contrast_pretty)) +
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
           x = m4_emm_contrast_df$estimate,
           y = 1 + 0.2,
           label = m4_emm_contrast_df$p_pretty,
           size = 3) +
  ggpubr::theme_pubr() +
  ggpubr::rremove("ylab")
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

ggsave(filename = paste0("figs/bat_daily_emergence_model_all_bats/emergence_model_all_survey_effect.png"), plot = m4_plot, width = 10, height = 6, dpi = 300)

ggsave(filename = paste0("figs/bat_daily_emergence_model_all_bats/emergence_model_all_survey.png"), plot = bat_plot, width = 10, height = 6, dpi = 300)
