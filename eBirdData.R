# 1: Initial prep

{# remotes::install_github("mstrimas/ebppackages")
library(auk)
library(lubridate)
library(sf)
library(gridExtra)
library(tidyverse)
library(raster)

memory.limit(200000)

# resolve namespace conflicts
select <- dplyr::select      

# function to convert time observation to hours since midnight
time_to_decimal <- function(x) {
  x <- hms(x, quiet = TRUE)
  hour(x) + minute(x) / 60 + second(x) / 3600
}

# Species list
list <- read_csv("~/ranges/species_list.csv")

# list of censored sp
censored <- c("Centrocercus minimus","Tympanuchus pallidicinctus",
              "Laterallus jamaicensis","Surnia ulula",
              "Strix occidentalis","Strix nebulosa","Falco rusticolus")

country="US"
}

##### Loop to create dataset for each country
for (month in 5){
  dir.create(paste0("~/data/processed/",country,"_",month))
  if(month %in% c(1,3,5,7,8,10,12)){end=31}
  if(month %in% c(4,6,9,11)){end=30}
  if(month==2){end=28}
  for (hhh in 0:67){
    if(hhh %in% 0:66){vals=((hhh*10)+1):((hhh*10)+10)}else{vals=671:677}
    # filter ebd 
    ebd <- auk_ebd(paste0("~/ebird/ebd_",country,"_relMay-2021.txt"),
                   file_sampling = "~/ebird/ebd_sampling_relMay-2021.txt")
    
    file <- "~/ebird/file.txt"
    checklists <- "~/ebird/checklists.txt"
    ebd_filters <- ebd %>% 
      auk_species(list[[1]][vals]) %>% 
      auk_country(c(country)) %>% 
      auk_date(date=c(paste0("*-",month,"-01"),
                      paste0("*-",month,"-",end))) %>% 
      auk_protocol(protocol = c("Stationary", "Traveling")) %>% 
      auk_distance(distance = c(0, 2)) %>% 
      auk_duration(duration = c(0, 180)) %>% 
      auk_complete()
    ebd_filters
    auk_filter(ebd_filters, file = file, 
               file_sampling = checklists,
               overwrite=T)
    
    # zero-fill
    ebird <- auk_zerofill(file, checklists, collapse = TRUE)
    
    # clean up variables
    ebird <- ebird %>% 
      mutate(
        # convert X to NA
        observation_count = if_else(observation_count == "X", 
                                    NA_character_, observation_count),
        observation_count = as.integer(observation_count),
        # effort_distance_km to 0 for non-traveling counts
        effort_distance_km = if_else(protocol_type != "Traveling", 
                                     0, effort_distance_km),
        # convert time to decimal hours since midnight
        time_observations_started = time_to_decimal(time_observations_started),
        # split date into year and day of year
        year = year(observation_date),
        day_of_year = yday(observation_date)
      )
    
    # additional filtering
    ebird <- ebird %>% 
      filter(
        # year
        year >= 2004,
        # 5 or fewer observers
        number_observers <= 5)
    
    # select final columns in dataset
    ebird <- ebird %>% 
      select(checklist_id, observer_id, locality_id,
             latitude, longitude,
             protocol_type, observation_date, year, 
             day_of_year, time_observations_started, 
             duration_minutes, effort_distance_km,
             number_observers, observation_count, scientific_name)

  # # remove inexperienced users
  # ebird_userlist <- ebird.merged %>%
  #   count(observer_id) %>%
  #   subset(n>=5)
  # ebird.merged <- ebird.merged %>%
  #   subset(observer_id %in% ebird_userlist$observer_id)
    
  # rearrange table
  ebird <- ebird %>% 
    pivot_wider(names_from = scientific_name, 
                values_from = observation_count)
  # save
  write_csv(ebird, paste0("~/data/processed/",country,"_",
                          month,"/temp",hhh,".csv"))

  }
}


# reload each piece and add zeroed sp with fixed list
# For any situations where the run was interrupted
month=10
file.list <- list.files(paste0("~/data/processed/US_",month,"/"))
for (lll in 1:length(file.list)){
 data <- read_csv(paste0("~/data/processed/US_",month,"/",file.list[lll]))
 if(lll==1){all <- data}else{
   all <- cbind(all, data[,14:ncol(data)])
 }
}
 zero.list <- list[[2]][which(!(list[[2]] %in% colnames(all)))]
 for (zzz in zero.list){
   all$new <- 0
   colnames(all)[ncol(all)] <- zzz
 }
 
 # standardize order
 all <- cbind(all[,1:13],all[,list[[2]]])
 # remove censored sp
 all <- all[,!(colnames(all) %in% censored)]
 # write
 write_csv(all, paste0("~/data/processed/US_",month,".csv"))
 
 
 
 
 
