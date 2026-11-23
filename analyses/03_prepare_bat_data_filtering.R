# ---------------------------
## load in packages that are required
library(tidyverse)

## load in custom functions
source("functions/UKBats_functions.R")
source("functions/activity_predict_tRackIT.R")
source("functions/activity_vars_tRackIT.R")



#Load and format data ####
#read csv of bat data
bat_activity_class <-  read_csv("output/bat_data_filtered_classification.csv")
bat_data <-  read_csv("output/bat_data_filtered.csv")

# set data types and join together
bat_activity_class <- mutate(bat_activity_class,
                             time = as.POSIXct(timestamp),
                             location = as.factor(location),
                             radio_id = as.factor(receiver),
                             max_rssi = max,
                             activity_class = prediction,
                             bat_name = as.factor(bat_name),
                             week = isoweek(ymd(as.Date(timestamp))))

# keep only needed columns
bat_activity_class <- select(bat_activity_class,
                             radio_id, max_rssi, time, location, n_receivers, activity_class, bat_name, week)

# set column type and remove duplicated of specific columns
bat_data <- mutate(bat_data,
                   radio_id = as.factor(radio_id),
                   tag_id = as.factor(tag_id),
                   node_id = as.factor(node_id),
                   time = as.POSIXct(time),
                   station_id = as.factor(station_id),
                   method = as.factor(method),
                   location = as.factor(location),
                   bat_name = as.factor(bat_name),
                   ring_no = as.factor(ring_no),
                   capture_site = as.factor(capture_site),
                   Date = as.Date(time),
                   daytime = hms::as_hms(time),
                   week = isoweek(ymd(Date)))

# remove duplicate rows
bat_data <- distinct(bat_data,
                     radio_id, tag_rssi, time, location, bat_name, .keep_all = TRUE)


#remove bats from first year
bat_data <- filter(bat_data, 
                   !bat_name %in% c('Ana', 'Bonnie', 'Cindy', 'Danny', 'Earl', 'Francine', 'Gaston', 'Henri'))



#Set dataframes for classification ####
# join dataframes together to see what was classified
bat_activity_class <- left_join(bat_data, bat_activity_class, by = c("radio_id", "time", "location", "bat_name", "week"))

# replace na values with 'none'
bat_activity_class <- mutate(bat_activity_class,
                             activity_class = as.character(activity_class))

bat_activity_class <- replace_na(bat_activity_class, 
                                 list(activity_class = 'none'))

bat_activity_class <- mutate(bat_activity_class,
                             activity_class = as.factor(activity_class))

# make receiver type column
bat_activity_class <- mutate(bat_activity_class,
                             rec_type = case_when(is.na(node_id) ~ "Station",
                                                  TRUE ~ "Node"))



# Reduce dataset to keep only one activity class per timepoint ####
# Order factor
bat_activity_class <- mutate(bat_activity_class,
                             activity_class = factor(activity_class, levels = c("active", "passive", "none")))

# 1 classification, prioritising "active"
bat_activity_class <- bat_activity_class %>%
  group_by(bat_name, time) %>%
  arrange(bat_name, time, activity_class) %>%
  slice_head(n = 1) %>%
  ungroup()

#make factor again
bat_activity_class <- mutate(bat_activity_class,
                             activity_class = as.factor(activity_class),
                             bat_name = as.factor(bat_name))

#drop unused factors
bat_activity_class <- droplevels(bat_activity_class)
 
#remove bats that were not classified
bat_activity_class <- filter(bat_activity_class, !bat_name %in% c('Karen', 'Felix', 'Dorris', 'Edna', 'Kesha'))



# write tables to csv ####
write_csv(bat_activity_class, "output/bat_data_filtered_classification_joined.csv")
