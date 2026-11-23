## ---------------------------
## load in packages that are required
library(tidyverse)

## load in custom functions
source("functions/UKBats_functions.R")



#Load and format data ####
#read csv of bat data
beep <-  read_csv("data/bat_data_filtered.csv")
bat_info <-  read_csv("data/bat_tags.csv")

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

#create column of rec type
beep <- mutate(beep,
               rec_type = case_when(is.na(node_id) ~ "Station",
                                    TRUE ~ "Node"))



# Remove data from bats when tag fall off ####
#Mary
mary_data <- filter(beep, bat_name == "Mary")

mary_data <- filter(mary_data, time <= "2022-10-26 12:00:00")

beep <- filter(beep, bat_name != "Mary")

beep <- bind_rows(beep, mary_data)


#Ross
ross_data <- filter(beep, bat_name == "Ross")

ross_data <- filter(ross_data, time <= "2023-10-09 00:00:00")

beep <- filter(beep, bat_name != "Ross")

beep <- bind_rows(beep, ross_data)


#Strachy
strachy_data <- filter(beep, bat_name == "Strachy")

strachy_data <- filter(strachy_data, time <= "2023-10-13 00:00:00")

beep <- filter(beep, bat_name != "Strachy")

beep <- bind_rows(beep, strachy_data)

rm(mary_data, ross_data, strachy_data)



# Calculate first and last detections ####
#remove unused levels of factor
beep$location <- droplevels(beep$location)

#remove any duplicate rows
beep <- distinct(beep, tag_id, time, location, .keep_all = TRUE)

# find the first and last detection of each bat
grouped_data <- beep %>%
  mutate(bat_name = as.character(bat_name)) %>%
  arrange(time) %>%
  group_by(bat_name) %>%
  summarise(first_time = min(time),
            last_time = max(time)) %>%
  data.frame()

bat_info <- bat_info %>%
  select(-first_time, -last_time)

# join together
bat_info <- left_join(bat_info, grouped_data)


ggplot(grouped_data, aes(x = first_time, y = bat_name, colour = bat_name)) +
  geom_point() +
  geom_point(aes(x = last_time)) +
  facet_wrap(~ year(first_time),
             scales = "free") +
  theme_bw() +
  theme(legend.position = "none")


#Write data to csv ####
write_csv(bat_info, "output/bat_tags.csv")
