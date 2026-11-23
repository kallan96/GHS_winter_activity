# ---------------------------
## load in packages that are required
library(glmmTMB)
library(tidyverse)

# set random seed
set.seed(180)



#Load and format data ####
#read csv of bat data
bat_metrics <-  read_csv("output/bat_data_fil_class_join_red_night_met_weather_today_daily.csv")
bat_metrics <- distinct(bat_metrics, bat_name, date, .keep_all = TRUE)


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




# Calculate number of detections per individual in models ####
#summaise raw counts
raw_counts <- bat_metrics %>% 
          group_by(study_year, week_integer, bat_emerged) %>%
          summarise(count = n(), .groups = "drop") %>%
          pivot_wider(names_from = bat_emerged,
                      values_from = count,
                      names_prefix = "count_",
                      values_fill = 0) %>%
          mutate(total = count_0 + count_1,
                 prop_active = count_1 / total)

#set week as factor
raw_counts$week_integer <- factor(raw_counts$week_integer , levels=c("1", "2", "3","4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29"))
raw_counts$week_integer <- plyr::revalue(raw_counts$week_integer, c("1" = "38", "2" = "39", "3" = "40", "4" = "41", "5" = "42", "6" = "43", "7" = "44", "8" = "45", "9" = "46", "10" = "47", "11" = "48", "12" = "49", "13" = "50", "14" = "51", "15" = "52", "16" = "1", "17" = "2", "18" = "3", "19" = "4", "20" = "5", "21" = "6", "22" = "7", "23" = "8", "24" = "9", "25" = "10", "26" = "11", "27" = "12", "28" = "13", "29" = "14"))

#fill values
raw_counts <- raw_counts %>%
  complete(study_year, week_integer,
           fill = list(count_0 = 0, count_1 = 0, total = 0, prop_active = 0))

#set weeks as months
raw_counts <- mutate(raw_counts,
                     month = case_when(week_integer %in% c("38", "39") ~ "Sep",
                                       week_integer %in% c("40", "41", "42", "43") ~ "Oct",
                                       week_integer %in% c("44", "45", "46", "47", "48") ~ "Nov",
                                       week_integer %in% c("49", "50", "51", "52") ~ "Dec",
                                       week_integer %in% c("1", "2", "3", "4") ~ "Jan",
                                       week_integer %in% c("5", "6", "7", "8") ~ "Feb",
                                       week_integer %in% c("9", "10", "11", "12", "13") ~ "Mar",
                                       TRUE ~ "Apr"))

#set as factor
raw_counts$month <- factor(raw_counts$month , levels=c("Sep", "Oct", "Nov","Dec", "Jan", "Feb", "Mar", "Apr"))                                    

# create pivot table of counts
pivot_long <- raw_counts %>%
                  select(study_year, month, week_integer, count_0, count_1) %>%
                    pivot_longer(cols=c('count_0', 'count_1'),
                                 names_to='emergence',
                                 values_to='count')

#set factor levels again
pivot_long$week_integer <- factor(pivot_long$week_integer , levels=c("38", "39", "40","41", "42", "43", "44", "45", "46", "47", "48", "49", "50", "51", "52", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14"))
pivot_long$month <- factor(pivot_long$month , levels=c("Sep", "Oct", "Nov","Dec", "Jan", "Feb", "Mar", "Apr"))                                    
pivot_long$emergence <- factor(pivot_long$emergence , levels=c("count_1", "count_0"))                                    
pivot_long$emergence <- plyr::revalue(pivot_long$emergence, c("count_1" = "Active", "count_0" = "Inactive"))


#plot results
bat_plot <- ggplot(pivot_long, aes(x = week_integer, y= count, fill = emergence)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Active" = "#D55E00", "Inactive" = "#56B4E9")) +
  facet_grid(cols = vars(month),
             rows = vars(study_year),
             scales = "free_x",
             space = "free_x",
             switch = "both",
             labeller = labeller(study_year = c("First Winter" = "September 2022 - April 2023",
                                                "Second Winter" = "September 2023 - April 2024"))) +
  labs(x = "Week of year", 
       y = "Number of Observations",
       fill = "Activity State") +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.spacing = unit(0, "lines"), 
        strip.background.x = element_blank(),
        strip.placement.x = "outside")


ggsave(filename = paste0("figs/bat_daily_emergence_model_all_bats/bat_days_count_stacked.png"), plot = bat_plot, width = 10, height = 6, dpi = 300)



# Plot of data distribution per individual bat ####
#remove bats that were not used in analyses
bat_red <- filter(bat_metrics, !bat_name %in% c('Helen', 'Lorna', 'Karen', 'Felix', 'Dorris', 'Edna', 'Kesha', 'Vera', 'Blob', 'Zeus', 'Abby', 'Gerty', 'Conchita', 'Jennet'))

# count number of detections per individual
raw_counts <- bat_red %>% 
  group_by(ring_no, bat_emerged) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = bat_emerged,
              values_from = count,
              names_prefix = "count_",
              values_fill = 0) %>%
  mutate(total = count_0 + count_1,
         prop_active = count_1 / total)

# create pivot table of counts
pivot_long <- raw_counts %>%
  select(ring_no, count_0, count_1) %>%
  pivot_longer(cols=c('count_0', 'count_1'),
               names_to='emergence',
               values_to='count')

# factor levels and set values
pivot_long$emergence <- factor(pivot_long$emergence , levels=c("count_1", "count_0"))                                    
pivot_long$emergence <- plyr::revalue(pivot_long$emergence, c("count_1" = "Active", "count_0" = "Inactive"))


# create plot of results
bat_plot <- ggplot(pivot_long, aes(x = ring_no, y= count, fill = emergence)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Active" = "#D55E00", "Inactive" = "#56B4E9")) +
  labs(x = "Bat ID", 
       y = "Number of Observations",
       fill = "Activity State") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, vjust = 0.75, hjust = 0.75))


ggsave(filename = paste0("figs/bat_daily_emergence_model_all_bats/bat_id_count_stacked.png"), plot = bat_plot, width = 10, height = 6, dpi = 300)


# create version split by year
raw_counts <- bat_red %>% 
  group_by(ring_no, study_year, bat_emerged) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = bat_emerged,
              values_from = count,
              names_prefix = "count_",
              values_fill = 0) %>%
  mutate(total = count_0 + count_1,
         prop_active = count_1 / total)

#make pivot table of counts
pivot_long <- raw_counts %>%
  select(study_year, ring_no, count_0, count_1) %>%
  pivot_longer(cols=c('count_0', 'count_1'),
               names_to='emergence',
               values_to='count')

# factor levels and set values
pivot_long$emergence <- factor(pivot_long$emergence , levels=c("count_1", "count_0"))                                    
pivot_long$emergence <- plyr::revalue(pivot_long$emergence, c("count_1" = "Active", "count_0" = "Inactive"))


# create plot of results
bat_plot <- ggplot(pivot_long, aes(x = ring_no, y= count, fill = emergence)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Active" = "#D55E00", "Inactive" = "#56B4E9")) +
  labs(x = "Bat ID", 
       y = "Count",
       fill = "Activity State") +
  facet_wrap(~study_year,
             scales = "free_x",
             space = "free_x",
             labeller = labeller(study_year = c("First Winter" = "September 2022 - April 2023",
                                                "Second Winter" = "September 2023 - April 2024"))) +
  labs(x = "Bat ID", 
       y = "Number of Observations",
       fill = "Activity State") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, vjust = 0.75, hjust = 0.75))


ggsave(filename = paste0("figs/bat_daily_emergence_model_all_bats/bat_id_count_year_stacked.png"), plot = bat_plot, width = 10, height = 6, dpi = 300)
