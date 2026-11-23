
#GPX file to tibble conversion ####
#A function to convert a GPX file to tibble format.GPX file must contain geometry column with corresponding latitude and longitude point.
#
#Returns database of waypoint locations in London timezone.
#
#Arguments
#route_gpx = A gpx file of transect

route_to_tibble <- function(route_gpx) {
  
  route.dataframe <- as.data.frame(route_gpx)
  
  #Separate geometry column into Latitude and Longitude
  route.dataframe <- separate(separate(separate(route.dataframe,
                                                col = geometry,
                                                into = c("Longitude", "Latitude"),
                                                sep = "\\,"),
                                       col = Longitude,
                                       into = c("c", "Longitude"),
                                       sep = "\\("),
                              col = Latitude,
                              into = c("Latitude"),
                              sep = "\\)")
  
  #Change Lat and Long from factor to numeric values
  route.dataframe <- mutate(route.dataframe,
                            Longitude = as.numeric(Longitude),
                            Latitude = as.numeric(Latitude))
  
  
  #Subset to meaningful data
  route.dataframe <- dplyr::select(route.dataframe, track_seg_point_id, ele, time, Longitude, Latitude)
  
  #Edit name of Long and Lat variable to correspond to route
  route.dataframe <- dplyr::rename(route.dataframe,
                                   R.Long = Longitude,
                                   R.Lat = Latitude)
  
  #Set timezone to match CTT datatime
  attr(route.dataframe$time, "tzone") <- "Europe/London"
  
  return(route.dataframe)
  
}



#Load in GPX files from folder ####
#A function to load all the GPX files from a folder and convert to individual tibble. The function also creates returns two separate databases, one of the start and end times of each transect, the other is a single database of all the gpx locations.
#
#Returns a list of the two dataframes.
#
#Arguments
#gpx_folder = pathway of the folder to be extracted from.

gpx_load <- function(gpx_folder) {
  
  #Get a list of all the file names
  gpx_filenames <- list.files(gpx_folder, pattern = "*.gpx")
  
  #Create blank data frames
  route_dataframe <- data.frame(matrix(ncol = 0, nrow = 0))
  route_timestamp <- data.frame(start_time = as.POSIXct(character()),
                                end_time = as.POSIXct(character()))
  
  #Create blank list
  individual_gpx <- list()
  
  #set transect number
  transect_number <- 1
  
  #Extract list into individual files
  for (filename in gpx_filenames) {
  
  #Assign filepath of individual filename
  filepath <- file.path(gpx_folder, filename)
  
  #Reading filename as GPX file
  gpx_file <- st_read(filepath, layer = "track_points")
  
  #Convert GPX file to tibble
  gpx_tibble <- route_to_tibble(gpx_file)
  
  #add transect number column
  gpx_tibble <- mutate(gpx_tibble,
                       transect_no = transect_number)
  
  #assign start and end time to dataframe
  new_time <- data.frame(start_time = as.POSIXct(format(head(gpx_tibble$time, 1)), tz = "Europe/London"),
                         end_time = as.POSIXct(format(tail(gpx_tibble$time, 1)), tz = "Europe/London"))
  
  #join into one dataframe
  route_dataframe <- bind_rows(route_dataframe, gpx_tibble)
  route_timestamp <- bind_rows(route_timestamp, new_time)
  
  #arrange by time
  route_dataframe <- arrange(route_dataframe, time)
  
  #Add individual GPX to list
  individual_gpx[[filename]] <- gpx_tibble
  
  #Increase transect number by 1
  transect_number <- transect_number + 1
  
  
  }
  
  #Combine into a list the route_dataframe and route_timestamp
  dataframes <- list(route = route_dataframe,
                     timestamp = route_timestamp)
  
  #Return output
  output <- list(summary_outputs = dataframes,
                 individual_ouputs = individual_gpx)
  
}



#Detection Waypoint Join Function ####
#A function to join CTT detection data to a corresponding waypoint by the time the record was taken. Both datasets must be in a tibble format and must include a variable "time" in POSIXct YYYY-MM-DD HH-MM-SS format.
#
#Returns database of positive detection and their corresponding waypoint position.
#
#Arguments
#detection_data = CTT detection database, must include time column
#route_data = transect route points, must include a latitude, longitude and time column
#type = the type of join you want to do with the data, see data.table for available types

detection_route_join <- function(detection_data, route_data, type) {
  
  #arranged data by time
  new_detection_data <- arrange(detection_data, time)
  new_route_data <- arrange(route_data, time)
  
  #change to datatable
  detection_data_dt <- data.table(new_detection_data, key = names(new_detection_data))
  route_data_dt <- data.table(new_route_data, key = names(new_route_data))
  
  #set column to match by (time)
  setkey(detection_data_dt, time)
  setkey(route_data_dt, time)
  
  #join datatables by nearest time
  all_trial_data <- route_data_dt[detection_data_dt, roll = type]
  
  #change to dataframe
  all_trial_data <- as.data.frame(all_trial_data)
  
  return(all_trial_data)
  
}



#Transect timeframe filter function ####
#A function to remove detection data from timepoints between transect periods of calibration data. Detection data from CTT sensor stations and nodes is the transect_data argument. timestamp corresponds to a dataframe of start and end times in as.POSIXct() format written as "YYYY-MM-DD HH:MM:SS" in whichever timezone you are in (see output of gpx_load[[1]][[2]] function). 
#
#Returns filtered database of positive detection, containing only points recorded during specified transect times.
#
#Arguments
#transect_data = CTT detection database, must include time column
#timestamp = timeframes of transect points - contains a "start" and "end" column (see gpx_load())

transect_data_filter <- function(transect_data, timestamp) {
  
  #Create blank data frames
  blank_dataframe <- data.frame(matrix(ncol = 0, nrow = 0))
  
  #filter points which lie between route times
  for (ntime in 1:nrow(timestamp)) {
    
    #Assign start and end values to object
    route_start <- timestamp[ntime,1]
    route_end <- timestamp[ntime,2]
    
    #filter data by start and end time
    route_data <- filter(transect_data,
                         time >= route_start & time <= route_end)
    
    #join into one dataframe
    blank_dataframe <- bind_rows(blank_dataframe, route_data)
    
    #arrange by time
    blank_dataframe <- arrange(blank_dataframe, time)
    
  }
  
  #Return output
  return(blank_dataframe)
  
}



#Function to calculate variables for activity classification ####
#
#This function loops through each receiver location in your dataframe and calculates the required explanatory variables for using in the Random Forest classifier. This function wraps those functions developed by Gottwald et al. 2023, available at https://github.com/Nature40/tRackIT. These functions have been modified to work on CTT data and are included in the repository. These functions require a trained random forest model to run. This function is wrapped by the 'bat_prediction_loop' function.
#
#Returns dataframe with explanatory variables for RF model.
#
#Arguments
#dataframe = CTT detection data carried over from 'bat_prediction_loop'

activity_loop <- function(dataframe) {
  
  # remove any unused levels in dataframe
  data <- droplevels(dataframe)
  
  # get unique locations
  data_locations <- unique(data$location)
  
  # create blank list to store results
  data_var_list <- list()
  
  # loop through each location to create model input
  for (loc in data_locations) {
    
    # filter dataframe by location
    df_filtered <- filter(data, location == loc)
    
    print(loc)
    print(nrow(df_filtered))
    
    # calculate variables
    df_variables <- activity_vars_tRackIT(df_filtered, tz = "UTC", rscale = 0)
    
    print(nrow(df_variables))
    
    # Append result to list
    data_var_list[[loc]] <- df_variables
    
  }
  
  # create dataframe of results
  data_var <- bind_rows(data_var_list)
  
  print(colnames(data_var))
  
  # check if there are any rows in the dataframe
  if(nrow(data_var) > 0) {
    
    # remove na_values
    data_var_drop <- drop_na(data_var,
                             timestamp, location, receiver, max, n_receivers, sumsq_max, sumsq_mean, std_hampel, max, sumsq_hampel, sumsq, skew)
    
    if(nrow(data_var_drop) > 0) {
      
      print('Making predictions.')
      #calculate predictions
      data_predict <- activity_predict_tRackIT(data = data_var_drop, get_data = TRUE)
      
      # return prediction results
      return(data_predict)
    }
  }
}
  


#Function to predict activity classification for bat dataframe ####
#
#This function loops through each bat ID in your dataframe and classifies the activity state of points using a Random Forest classifier. This function wraps those functions developed by Gottwald et al. 2023, available at https://github.com/Nature40/tRackIT. These functions have been modified to work on CTT data and are included in the repository. These functions require a trained random forest model to run. This function is a wrapper for the 'activity_loop' function.
#
#Returns dataframe with points classified as 'active' or 'inactive'.
#
#Arguments
#dataframe = CTT detection data of individuals bats

bat_prediction_loop <- function(dataframe) {
  
  # remove any unused levels in dataframe
  data <- droplevels(dataframe)
  
  # get unique bats
  unique_bats <- unique(data$bat_name)
  
  # create blank list to store results
  data_bat_list <- list()
  
  # loop through each bat to create model predictions
  for (bat in unique_bats) {
    
    print(bat)
    
    # filter dataframe by location
    df_filtered <- filter(data, bat_name == bat)
    
    # remove any unused levels in dataframe
    df_filtered <- droplevels(df_filtered)
    
    # calculate variables
    df_output <- activity_loop(df_filtered)
    
    # Skip if NULL or empty
    if (is.null(df_output) || nrow(df_output) == 0) {
      next
    }
    
    # add column with bat name
    df_output$bat_name <- bat
    
    # Append result to list
    data_bat_list[[bat]] <- df_output
    
  }
  
  # create dataframe of results
  data_bat <- bind_rows(data_bat_list)
  
  return(data_bat)
  
}
