## ---------------------------
## load in packages that are required
library(tidyverse)

## load in custom functions
source("functions/UKBats_functions.R")
source("functions/activity_predict_tRackIT.R")
source("functions/activity_vars_tRackIT.R")



#Load and format data ####
#read csv of bat data
beep <-  read_csv("data/bat_data.csv")

# set column type
beep <- mutate(beep,
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
               daytime = hms::as_hms(time))



# Filter records that are corrupt/wrong ####
# remove values of rssi from antenna which is corrupt
beep <- mutate(beep,
               tag_rssi = case_when(is.na(node_id) & station_id == 'V30B0154C8C9' & radio_id == 1 ~ NA,
                                    is.na(node_id) & station_id == 'V30B0154D690' & radio_id == 1 ~ NA,
                                    TRUE ~ tag_rssi))

# remove points that were before station firmware were fixed
beep <- mutate(beep,
               tag_rssi = case_when(is.na(node_id) & time > "2022-06-01 00:00:00" & time < "2023-03-30 00:00:00" ~ NA,
                                    TRUE ~ tag_rssi))


# Remove data from nodes with corrupted timestamps and small dataframe size
#create separate dataframe of kl_divergence data
small_data <- filter(beep, 
                     method == "kl_div")

#Remove KL Divergence calculated points smaller than 1000 rows
small_data <- filter(small_data, 
                     data_size >= 1000)

#filter to remove these observations
beep <- filter(beep, 
               method != "kl_div")

#combine dataframes together
beep <- bind_rows(beep, small_data)

rm(small_data)


# Remove data from bats with mistaken detections, e.g., known tag droppoffs or visual plotting data
#Henri
bat_data <- filter(beep, 
                   bat_name == "Henri")

bat_data <- filter(bat_data, 
                   time <= "2022-01-01 00:00:00")

beep <- filter(beep, 
               bat_name != "Henri")

beep <- bind_rows(beep, bat_data)


#Kay
bat_data <- filter(beep, bat_name == "Kay")

bat_data <- filter(bat_data, time <= "2023-10-03 06:00:00")

beep <- filter(beep, bat_name != "Kay")

beep <- bind_rows(beep, bat_data)


#Larry
bat_data <- filter(beep, bat_name == "Larry")

bat_data <- filter(bat_data, time <= "2022-10-11 00:00:00")

beep <- filter(beep, bat_name != "Larry")

beep <- bind_rows(beep, bat_data)


#May
bat_data <- filter(beep, bat_name == "May")

bat_data <- filter(bat_data, time <= "2023-10-02 17:00:00")

beep <- filter(beep, bat_name != "May")

beep <- bind_rows(beep, bat_data)


#Usher
bat_data <- filter(beep, bat_name == "Usher")

bat_data <- filter(bat_data, time <= "2023-01-05 07:00:00")

beep <- filter(beep, bat_name != "Usher")

beep <- bind_rows(beep, bat_data)


#Xanthi
bat_data <- filter(beep, bat_name == "Xanthi")

bat_data <- filter(bat_data, time <= "2023-03-14 22:00:00")

beep <- filter(beep, bat_name != "Xanthi")

beep <- bind_rows(beep, bat_data)

rm(bat_data)


# remove duplicate rows
beep <- distinct(beep,
                     radio_id, tag_id, tag_rssi, time, location, bat_name, ring_no, .keep_all = TRUE)



#Write filtered data to csv ####
write_csv(beep, "output/bat_data_filtered.csv")



# Remove bats from preliminary trials ####
#remove bats from first year
beep <- filter(beep, 
               !bat_name %in% c('Ana', 'Bonnie', 'Cindy', 'Danny', 'Earl', 'Francine', 'Gaston', 'Henri'))



# Classify activity state of data ####
# remove NA values for tag rssi
beep_no_na <- drop_na(beep, tag_rssi)

# remove levels of any variable no longer found in the dataframe
beep_no_na <- droplevels(beep_no_na)

# classify bat behaviour
bat_activity_classification <- bat_prediction_loop(beep_no_na) 



#Write classified data to csv ####
write_csv(bat_activity_classification, "output/bat_data_filtered_classification.csv")
