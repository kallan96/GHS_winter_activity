# ---------------------------
## load in packages that are required
library(moonlit)
library(tidyverse)

## load in custom functions
source("functions/UKBats_functions.R")



#Load and format data ####
#read csv of bat data
bat_metrics <-  read_csv("output/bat_data_filtered_class_join_red_night_metrics.csv")
bat_metrics <- distinct(bat_metrics, bat_name, bat_day, .keep_all = TRUE)

sun <-  read_csv("output/bat_tags.csv")
sun <- distinct(sun, Date, .keep_all = TRUE)

bat <-  read_csv("data/sun_data.csv")
bat <- distinct(bat, bat_name, .keep_all = TRUE)


# set data types
bat_metrics <- mutate(bat_metrics,
                      bat_name = as.factor(bat_name),
                      bat_day = as.Date(bat_day),
                      bat_day_next = as.Date(bat_day) + 1,
                      sunsetStart = as.POSIXct(sunsetStart),
                      sunriseEnd = as.POSIXct(sunriseEnd),
                      emergence_time = as.POSIXct(emergence_time),
                      reentry_time = as.POSIXct(reentry_time))


#order columns
bat_metrics <- select(bat_metrics, bat_name, bat_day, bat_day_next, sunsetStart, sunriseEnd, night_length, night_length_rounded, emergence_time, emergence_time_minutes, reentry_time, reentry_time_minutes, total_activity_minutes, total_activity_proportion, absolute_activity_minutes, absolute_activity_proportion)


#remove unused levels of factor
bat_metrics$bat_name <- droplevels(bat_metrics$bat_name)


#create column of season
bat_metrics <- mutate(bat_metrics,
                      study_year = ifelse(bat_day < '2023-06-01', 'First Winter', 'Second Winter'))


#arrange by bat names and date
bat_metrics <- arrange(bat_metrics,
                       bat_name, bat_day)


#filter bat info by bats in data
ind_bats <- levels(bat_metrics$bat_name)
bat <- filter(bat,
              bat_name %in% ind_bats)


#extract first and last records
start_FW <- pull(head(arrange(filter(bat, time < '2023-06-01'), first_time), 1), first_time)
end_FW <- pull(tail(arrange(filter(bat, time < '2023-06-01'), last_time), 1), last_time)
start_SW <- pull(head(arrange(filter(bat, time > '2023-06-01'), first_time), 1), first_time)
end_SW <- pull(tail(arrange(filter(bat, time > '2023-06-01'), last_time), 1), last_time)



#Create a custom color scale for locations ####
myColors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
              "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5",
              "#393b79", "#637939", "#8c6d31", "#843c39", "#7b4173", "#5254a3", "#8ca252", "#bd9e39", "#ad494a", "#a55194",
              "#6b6ecf", "#b5cf6b", "#e7ba52", "#d6616b", "#ce6dbd", "#9c9ede", "#cedb9c", "#e7cb94", "#f4a582", "#d5d5d5",
              "#c7e9c0", "#bcbddc", "#fdae6b", "#fdd0a2", "#e7969c", "#de9ed6", "#ff9896", "#fdae61", "#fee08b", "#9e0142",
              "#d53e4f", "#f46d43")
names(myColors) <- c("Bat Run", "Bovey Tracey", "Brookwood Mines", "Chudleigh", "Chudleigh Bat House", "Chudleigh Bowles", "Chudleigh Caves and Woods", "Chudleigh Cemetery", "Chudleigh Central", "Chudleigh Knighton", "Chudleigh Knighton Heath", "Clay Pits Way Bridge", "Clifford Business Park", "Cow Cave Wood", "Devon Express Way Layby", "Fosterville", "Gappah", "Gappah Cows", "Haytor and Smallacombe Iron Mines", "Highland Park", "Hillfort", "Huxbear Vineyard", "Kate Brook Farm", "Kate Brook Underpass", "Kiln Wood", "Kingsteington", "Lawell House", "Lawelldown Wood", "Little Lyndridge Farm", "Lower Upcott", "Olchard", "Olchard A", "Olchard B", "Olchard C", "Olchard D", "Olchard Lay-by", "Old Way", "Palace Meadow", "Penn Recca Tunnel", "Pipistrelle Close", "Rock Farm", "Sandford Orleigh House", "Sheepdip", "Speed Camera", "Stony Copse South", "Stover", "Stover School", "Teign Bridge", "Teigngrace Meadow", "The Granary Mill", "West Gappah", "Winstow Cottages")
colScale <- scale_colour_manual(name = "Receiver Locations",values = myColors)
fillScale <- scale_fill_manual(name = "Receiver Locations",values = myColors)



#create custome colour scale for bat_names ####
#scales::viridis_pal()(length(unique(beep$bat_name)))
myColors <- c("#440154FF", "#46075AFF", "#470D60FF", "#471365FF", "#48186AFF", "#481D6FFF", "#482374FF", "#482878FF", "#472D7BFF", "#46337EFF", "#453881FF", "#433C84FF", "#424186FF", "#404688FF", "#3E4A89FF", "#3C4F8AFF", "#3A538BFF", "#38588CFF", "#365C8DFF", "#34608DFF", "#32648EFF", "#31688EFF", "#2F6C8EFF", "#2D708EFF", "#2C738EFF", "#2A778EFF", "#297B8EFF", "#277F8EFF", "#26828EFF", "#24868EFF", "#238A8DFF", "#218E8DFF", "#20928CFF", "#1F968BFF", "#1F9A8AFF", "#1F9E89FF", "#1FA187FF", "#21A585FF", "#23A983FF", "#26AD81FF", "#2AB07FFF", "#2FB47CFF", "#35B779FF", "#3BBB75FF", "#42BE71FF", "#4AC16DFF", "#52C569FF", "#5AC864FF", "#64CB5FFF", "#6DCD59FF", "#76D054FF", "#80D34EFF", "#8AD547FF", "#94D840FF", "#9FDA3AFF", "#A9DC33FF", "#B4DE2CFF", "#BFDF25FF", "#C9E11FFF", "#D4E21AFF", "#DFE318FF", "#EAE51AFF", "#F4E61EFF", "#FDE725FF")

names(myColors) <- c("Abby", "Ana", "Andrea", "Bert", "Blob", "Bonnie", "Charlie", "Cindy", "Conchita", "Danny", "Derek", "Dorris", "Earl", "Edna", "Erik", "Felix", "Francine", "Fred",  "Gaston", "George", "Gerty",  "Heather",  "Henri",  "Ian",  "Iris", "Isla", "Jennet", "Jez",  "Joyce",  "Karen",  "Kay", "Kesha", "Larry",    "Laticia",  "Mary", "May",  "Moira", "Nancy", "Neil", "Norman", "Olga", "Olivia", "Pedro", "Polly", "Queenie", "Quincey", "Rex", "Ross", "Steve", "Strachy", "Tony", "Tracy", "Uri", "Usher", "Vera", "Vincent", "Wally", "Wilma", "Xanthi", "Xray", "Yasmin", "Yvette", "Zac", "Zeus" )

colbatScale <- scale_colour_manual(name = "Bat Names",values = myColors)
fillbatScale <- scale_fill_manual(name = "Bat Names",values = myColors)



# sorting weather data ####
#load in files
exeter_airport_2024 <- read_csv("data/weather_data/2024_exeter_airport_no_2_hourly_weather_data.csv")
exeter_aws_2024 <- read_csv("data/weather_data/2024_exeter_surfacenet_enclosure_aws_hourly_weather_data.csv")
exeter_airport_2023 <- read_csv("data/weather_data/2023_exeter_airport_no_2_hourly_weather_data.csv")
exeter_aws_2023 <- read_csv("data/weather_data/2023_exeter_surfacenet_enclosure_aws_hourly_weather_data.csv")
exeter_airport_2022 <- read_csv("data/weather_data/2022_exeter_airport_no_2_hourly_weather_data.csv")
exeter_aws_2022 <- read_csv("data/weather_data/2022_exeter_surfacenet_enclosure_aws_hourly_weather_data.csv")

#rain data
exeter_airport_rain_2022 <- read_csv("data/weather_data/2022_exeter_airport_no_2_hourly_rain.csv")
exeter_airport_rain_2023 <- read_csv("data/weather_data/2023_exeter_airport_no_2_hourly_rain.csv")


#join data frames
weather_data <- bind_rows(exeter_airport_2022, exeter_airport_2023, exeter_aws_2022, exeter_aws_2023)
rain_data_all <- bind_rows(exeter_airport_rain_2022, exeter_airport_rain_2023)


#make date column
weather_data <- mutate(weather_data, Date = as.Date(ob_time))


#filter rain data to keep only 1 obv count
rain_data <- filter(rain_data_all, ob_hour_count == 1)


#select needed columns
rain_data <- select(rain_data,
                    ob_end_time, prcp_amt)


#join matching columns
weather_data <- mutate(weather_data, 
                   time = ob_time,
                   Relative_Humidity = rltv_hum,
                   Hourly_Mean_Wind_Speed = wind_speed, 
                   Total_Cloud_Amount = cld_ttl_amt_id, 
                   Mean_Sea_Level_Pressure  = msl_pressure, 
                   Hourly_Air_Temperature = air_temperature)


#select columns 
weather_data <- select(weather_data,
                   Date, time, Relative_Humidity, Hourly_Mean_Wind_Speed, Total_Cloud_Amount, Mean_Sea_Level_Pressure, Hourly_Air_Temperature)


#change variable names of rain data
rain_data <- rename(rain_data,
                    Hourly_Rainfall = prcp_amt)


#join to full dataframe
weather_data <- left_join(weather_data, rain_data, by = c('time' = 'ob_end_time'))


#select columns 
weather_data <- select(weather_data,
                       Date, time, Relative_Humidity, Hourly_Mean_Wind_Speed, Total_Cloud_Amount, Mean_Sea_Level_Pressure, Hourly_Air_Temperature, Hourly_Rainfall)


#create date column
exeter_airport_2024 <- mutate(exeter_airport_2024, 
                              Date = as.Date(time))
exeter_aws_2024 <- mutate(exeter_aws_2024, 
                              Date = as.Date(time))

exeter_airport_2024 <- select(exeter_airport_2024,
                              time, Date, Relative_Humidity, Hourly_Mean_Wind_Speed, Total_Cloud_Amount, Mean_Sea_Level_Pressure, Hourly_Air_Temperature, Hourly_Rainfall)
exeter_aws_2024 <- select(exeter_aws_2024,
                          time, Date, Relative_Humidity, Total_Cloud_Amount, Mean_Sea_Level_Pressure, Hourly_Air_Temperature)


#join with 2024 data
weather_data <- bind_rows(weather_data, exeter_airport_2024, exeter_aws_2024)


#remove any columns where all data is na
weather_data <- weather_data[ , colSums(is.na(weather_data))!=nrow(weather_data)]


#filter to keep only data from study period
weather_FW <- filter(weather_data,
                     Date >= (as.Date(start_FW)-1) & Date <= as.Date(end_FW))
weather_SW <- filter(weather_data,
                     Date >= (as.Date(start_SW)-1) & Date <= as.Date(end_SW))
weather_data <- bind_rows(weather_FW, weather_SW)



# write tables to csv ####
write_csv(weather_data, "output/weather_data.csv")



# create moonlight variable ####
#summarise weather data dates to find first and last of each date
dates <- summarise(group_by(weather_data, Date),
                   first_time = min(time),
                   last_time = max(time))


#create dataframe of all times between these dates at 1 minute intervals
expanded_dates <- rename(select(unnest(mutate(rowwise(dates),
                                              time_seq = list(seq(from = first_time, 
                                                                  to = last_time, 
                                                                  by = "1 min"))),
                                       time_seq),
                                time_seq),
                         time = time_seq)


# Define location in chudleigh
latitude <- 50.590305 
longitude <- -3.612610 


#calculate moonlight data for each period and floor to each hour for grouping
moon_results <- calculateMoonlightIntensity(latitude, longitude, expanded_dates$time, 0.26)
moon_results <- mutate(moon_results, 
                       date_rounded = floor_date(date, "hour"))



# write tables to csv ####
write_csv(moon_results, "output/moon_data.csv")



#create hourly summary of moon data ####
#summarise per hour
moon_results_summary <- ungroup(summarise(group_by(moon_results, date_rounded),
                                          meanMoonlightIntensity = mean(moonlightModel),
                                          meanMoonPhase = mean(moonPhase)))


#create hourly summary of weather variables
weather_summary_hourly <- ungroup(summarise(group_by(weather_data, 
                                                     time, Date),
                                            av_rel_humidity = mean(Relative_Humidity, na.rm = TRUE),
                                            av_wind_speed = mean(Hourly_Mean_Wind_Speed, na.rm = TRUE),
                                            av_cloud_okta = mean(Total_Cloud_Amount, na.rm = TRUE),
                                            av_air_pressure = mean(Mean_Sea_Level_Pressure, na.rm = TRUE),
                                            av_air_temp = mean(Hourly_Air_Temperature, na.rm = TRUE),
                                            daily_rain = mean(Hourly_Rainfall, na.rm = TRUE)))


#find na_values for rainfall
na_rainfall <- filter(weather_summary_hourly, is.na(daily_rain) & time < "2024-01-01 00:00:00")


#keep accumulated readings
accumulated_rainfall <- filter(rain_data_all, ob_hour_count != 1)


#find times which match
match_times <- filter(accumulated_rainfall, ob_end_time %in% na_rainfall$time)


#arrange by time
match_times <- arrange(match_times, ob_end_time)


# fill missing rainfall values
weather_summary_hourly <- mutate(weather_summary_hourly,
                                 daily_rain = case_when(is.na(daily_rain) & time > match_times$ob_end_time[2] & time < match_times$ob_end_time[4] ~ (match_times$prcp_amt[3]/24),
                                                        is.na(daily_rain) & time < "2024-01-01 00:00:00" ~ 0,
                                                        TRUE ~ daily_rain))



#join moon data to weather data ####
#jopin dataframes
weather_summary_hourly <- left_join(weather_summary_hourly, moon_results_summary, by = c('time' = 'date_rounded'))


#create next day column
weather_summary_hourly <- mutate(weather_summary_hourly,
                                 Date_next = Date + 1)


#combine with sunrise/set data
weather_summary_hourly <- left_join(weather_summary_hourly, select(sun, Date, sunsetStart), by = c('Date' = 'Date'))
weather_summary_hourly <- left_join(weather_summary_hourly, select(sun, Date, sunriseEnd), by = c('Date_next' = 'Date'))


#round to nearest hour and set the Date
weather_summary_hourly <- mutate(weather_summary_hourly,
                                 sunriseEnd = ceiling_date(sunriseEnd, "hour"),
                                 sunsetStart = floor_date(sunsetStart, "hour"),
                                 Date = case_when(time <= sunriseEnd & time >= sunsetStart ~ Date,
                                                  TRUE ~ Date - 1))


#recreate next date column to ensuire correct date
weather_summary_hourly <- select(weather_summary_hourly,
                                 -Date_next, -sunsetStart, -sunriseEnd)

weather_summary_hourly <- mutate(weather_summary_hourly,
                                 Date_next = Date + 1)


#rejoin to sundata to ensure correct join
weather_summary_hourly <- left_join(weather_summary_hourly, select(sun, Date, sunsetStart), by = c('Date' = 'Date'))
weather_summary_hourly <- left_join(weather_summary_hourly, select(sun, Date, sunriseEnd), by = c('Date_next' = 'Date'))


#round again
weather_summary_hourly <- mutate(weather_summary_hourly,
                                 sunriseEnd = ceiling_date(sunriseEnd, "hour"),
                                 sunsetStart = floor_date(sunsetStart, "hour"))

#filter to keep nighttime data
weather_summary_hourly_night <- filter(arrange(weather_summary_hourly, time),
                                       time >= sunsetStart & time <= sunriseEnd)


#rejoin to sundata to ensure correct join
weather_summary_hourly <- left_join(select(weather_summary_hourly, -sunriseEnd, -sunsetStart), select(sun, Date, sunsetStart), by = c('Date' = 'Date'))
weather_summary_hourly <- left_join(weather_summary_hourly, select(sun, Date, sunriseEnd), by = c('Date_next' = 'Date'))

#rejoin to sundata to ensure correct join
weather_summary_hourly_night <- left_join(select(weather_summary_hourly_night, -sunriseEnd, -sunsetStart), select(sun, Date, sunsetStart), by = c('Date' = 'Date'))
weather_summary_hourly_night <- left_join(weather_summary_hourly_night, select(sun, Date, sunriseEnd), by = c('Date_next' = 'Date'))



# write tables to csv ####
write_csv(weather_summary_hourly, "output/weather_data_summary_hourly.csv")
write_csv(weather_summary_hourly_night, "output/weather_data_summary_hourly_night.csv")



#create daily summary of moon data ####
#create daily summary of weather variables
weather_summary_daily <- ungroup(summarise(group_by(weather_summary_hourly, 
                                                    Date),
                                           av_rel_humidity = mean(av_rel_humidity, na.rm = TRUE),
                                           av_wind_speed = mean(av_wind_speed, na.rm = TRUE),
                                           av_cloud_okta = mean(av_cloud_okta, na.rm = TRUE),
                                           av_air_pressure = mean(av_air_pressure, na.rm = TRUE),
                                           av_air_temp = mean(av_air_temp, na.rm = TRUE),
                                           daily_rain = sum(daily_rain, na.rm = TRUE),
                                           meanMoonlightIntensity = mean(meanMoonlightIntensity),
                                           meanMoonPhase = mean(meanMoonPhase)))

weather_summary_daily_night <- ungroup(summarise(group_by(weather_summary_hourly_night, 
                                                    Date),
                                                 av_rel_humidity = mean(av_rel_humidity, na.rm = TRUE),
                                                 av_wind_speed = mean(av_wind_speed, na.rm = TRUE),
                                                 av_cloud_okta = mean(av_cloud_okta, na.rm = TRUE),
                                                 av_air_pressure = mean(av_air_pressure, na.rm = TRUE),
                                                 av_air_temp = mean(av_air_temp, na.rm = TRUE),
                                                 daily_rain = sum(daily_rain, na.rm = TRUE),
                                                 meanMoonlightIntensity = mean(meanMoonlightIntensity),
                                                 meanMoonPhase = mean(meanMoonPhase)))


#create next day column
weather_summary_daily <- mutate(weather_summary_daily,
                                Date_next = Date + 1)

weather_summary_daily_night <- mutate(weather_summary_daily_night,
                                      Date_next = Date + 1)


#combine with sunrise/set data
weather_summary_daily <- left_join(weather_summary_daily, select(sun, Date, sunsetStart), by = c('Date' = 'Date'))
weather_summary_daily <- left_join(weather_summary_daily, select(sun, Date, sunriseEnd), by = c('Date_next' = 'Date'))

weather_summary_daily_night <- left_join(weather_summary_daily_night, select(sun, Date, sunsetStart), by = c('Date' = 'Date'))
weather_summary_daily_night <- left_join(weather_summary_daily_night, select(sun, Date, sunriseEnd), by = c('Date_next' = 'Date'))


rm(exeter_airport_2022, exeter_airport_2023, exeter_airport_2024, exeter_airport_rain_2022, exeter_airport_rain_2023, exeter_aws_2022, exeter_aws_2023, exeter_aws_2024, rain_data)



# write tables to csv ####
write_csv(weather_summary_daily, "output/weather_data_summary_daily.csv")
write_csv(weather_summary_daily_night, "output/weather_data_summary_daily_night.csv")



#join bat data to weather data ####
# Initialize an empty list to store results for each bat
bat_weather_today_list <- list()
bat_weather_previous_list <- list()

bat <- mutate(bat, 
              first_time_rounded = floor_date(first_time, "day"),
              last_time_rounded = ceiling_date(last_time, "day"))

# Loop over each unique bat
for (bat_id in ind_bats) {
  
  # Filter bat data for this specific bat
  bat_ind <- select(filter(bat_metrics, bat_name == bat_id),
                    -sunsetStart, -sunriseEnd)
  
  
  # Get the first and last observation date for the bat
  first_obs_date <- pull(filter(bat, bat_name == bat_id),
                         first_time_rounded)
  last_obs_date <- pull(filter(bat, bat_name == bat_id),
                        last_time_rounded)
                        
  
  # Create a complete sequence of dates from first_obs_date to last_obs_date
  complete_dates <- tibble(date = seq.Date(from = as.Date(first_obs_date), to = as.Date(last_obs_date), by = "day"))
  complete_dates <- mutate(complete_dates,
                           Date_prev = date -1)
  
  # Left join the complete date sequence with the weather data
  joined_data_today <- left_join(complete_dates, select(weather_summary_daily_night, -Date_next), by = c("date" = "Date"))
  joined_data_previous <- left_join(complete_dates, select(weather_summary_daily_night, -Date_next), by = c("Date_prev" = "Date"))
  
  # Left join bat data to weather data, adding `bat_id` and `bat_emerged`
  joined_data_today <- joined_data_today %>%
    left_join(bat_ind, by = c("date" = "bat_day")) %>%
    mutate(
      bat_emerged = if_else(!is.na(bat_name), 1, 0),  # 1 if bat data exists, 0 if not
      bat_name = bat_id
    )
  
  joined_data_previous <- joined_data_previous %>%
    left_join(bat_ind, by = c("date" = "bat_day")) %>%
    mutate(
      bat_emerged = if_else(!is.na(bat_name), 1, 0),  # 1 if bat data exists, 0 if not
      bat_name = bat_id
    )
  
  # Append result to list
  bat_weather_today_list[[bat_id]] <- joined_data_today
  bat_weather_previous_list[[bat_id]] <- joined_data_previous
  
}


# Combine all bat-specific dataframes into one final dataframe
bat_metrics_weather_today_daily <- bind_rows(bat_weather_today_list)
bat_metrics_weather_previous_daily <- bind_rows(bat_weather_previous_list)


#arrange by date
bat_metrics_weather_today_daily <- arrange(bat_metrics_weather_today_daily, date)
bat_metrics_weather_previous_daily <- arrange(bat_metrics_weather_previous_daily, date)


#recreate bat_day_next and season columns
bat_metrics_weather_today_daily <- mutate(bat_metrics_weather_today_daily,
                                          time = as.numeric(date),
                                          study_year = ifelse(date < '2023-06-01', 'First Winter', 'Second Winter'),
                                          week = isoweek(ymd(date)),
                                          month = as.numeric(format(date, "%m")),
                                          season = ifelse(month %in% c(3, 4, 5), 'spring',
                                                          ifelse(month %in% c(6, 7, 8), 'summer',
                                                                 ifelse( month %in% c(9, 10, 11), 'autumn', 'winter'))),
                                          week_integer = case_when(week < 38 ~ week + 15,
                                                                   TRUE ~ week -37),
                                          two_week = case_when(week_integer %% 2 == 0 ~ week_integer/2,
                                                               TRUE ~ (week_integer + 1)/2))

bat_metrics_weather_previous_daily <- mutate(bat_metrics_weather_previous_daily,
                                             time = as.numeric(date),
                                             study_year = ifelse(date < '2023-06-01', 'First Winter', 'Second Winter'),
                                             week = isoweek(ymd(date)),
                                             month = as.numeric(format(date, "%m")),
                                             season = ifelse(month %in% c(3, 4, 5), 'spring',
                                                             ifelse(month %in% c(6, 7, 8), 'summer',
                                                                    ifelse( month %in% c(9, 10, 11), 'autumn', 'winter'))),
                                             week_integer = case_when(week < 38 ~ week + 15,
                                                                      TRUE ~ week -37),
                                             two_week = case_when(week_integer %% 2 == 0 ~ week_integer/2,
                                                                  TRUE ~ (week_integer + 1)/2))


#replace NA's with 0's in specific column
bat_metrics_weather_today_daily <- replace_na(bat_metrics_weather_today_daily,
                                              list(emergence_time_minutes = 0, reentry_time_minutes = 0, total_activity_minutes = 0, total_activity_proportion = 0, absolute_activity_minutes = 0, absolute_activity_proportion = 0))

bat_metrics_weather_previous_daily <- replace_na(bat_metrics_weather_previous_daily,
                                                 list(emergence_time_minutes = 0, reentry_time_minutes = 0, total_activity_minutes = 0, total_activity_proportion = 0, absolute_activity_minutes = 0, absolute_activity_proportion = 0))


#select metrics from bat dataframe to be joined
bat <- select(bat, 
              bat_name, tag_id, location, survey_type, ring_no, forearm, sex, repro, age, weight)


#change age class YoY to Juvenile
bat$age[bat$age=="YoY"]<-"J"


#join into one dataframe
bat_metrics_weather_today_daily <- left_join(bat_metrics_weather_today_daily, bat)
bat_metrics_weather_previous_daily <- left_join(bat_metrics_weather_previous_daily, bat)


#reorder dataframes
bat_metrics_weather_today_daily <- select(bat_metrics_weather_today_daily, 
                                          time, study_year, season, month, week, week_integer, two_week, date, bat_name, tag_id, location, survey_type, ring_no, forearm, sex, repro, age, weight, bat_emerged, emergence_time, emergence_time_minutes, reentry_time, reentry_time_minutes, total_activity_minutes, total_activity_proportion, absolute_activity_minutes, absolute_activity_proportion, sunsetStart, sunriseEnd, night_length, night_length_rounded, av_rel_humidity, av_wind_speed, av_cloud_okta, av_air_pressure, av_air_temp, daily_rain, meanMoonlightIntensity, meanMoonPhase) 

bat_metrics_weather_previous_daily <- select(bat_metrics_weather_previous_daily, 
                                          time, study_year, season, month, week, week_integer, two_week, date, bat_name, tag_id, location, survey_type, ring_no, forearm, sex, repro, age, weight, bat_emerged, emergence_time, emergence_time_minutes, reentry_time, reentry_time_minutes, total_activity_minutes, total_activity_proportion, absolute_activity_minutes, absolute_activity_proportion, sunsetStart, sunriseEnd, night_length, night_length_rounded, av_rel_humidity, av_wind_speed, av_cloud_okta, av_air_pressure, av_air_temp, daily_rain, meanMoonlightIntensity, meanMoonPhase) 



# write tables to csv ####
write_csv(bat_metrics_weather_today_daily, "output/bat_data_fil_class_join_red_night_met_weather_today_daily.csv")
write_csv(bat_metrics_weather_previous_daily, "output/bat_data_fil_class_join_red_night_met_weather_previous_daily.csv")
