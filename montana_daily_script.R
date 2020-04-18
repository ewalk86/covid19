# Pull daily updated State of Montana covid data, estimate R, save results/plots
# Ethan Walker
# 18 April 2020

library(tidyverse)
library(readxl)
library(naniar)
library(lubridate)
library(zoo)
library(eeptools)
library(knitr)
library(incidence)
library(EpiEstim)
library(earlyR)
library(projections)
library(distcrete)
library(epitrix)
library(jsonlite)
library(httr)
library(rlist)
library(taskscheduleR)


# Query and prep data
query_url <- "https://services.arcgis.com/qnjIrwR8z5Izc0ij/ArcGIS/rest/services/COVID_Cases_Production_View/FeatureServer/2/query?where=1%3D1&objectIds=&time=&resultType=none&outFields=OBJECTID%2C+Case_No%2C+Date_Reported_to_CDEpi%2C+County%2C+Age_Group%2C+Sex%2C+Hospitalization%2C+Outcome%2C+MT_case&returnIdsOnly=false&returnUniqueIdsOnly=false&returnCountOnly=false&returnDistinctValues=false&cacheHint=false&orderByFields=&groupByFieldsForStatistics=&outStatistics=&having=&resultOffset=&resultRecordCount=&sqlFormat=none&f=pjson&token="

initial_pull <- GET(query_url)

text_data <- content(initial_pull, as = "text")

parsed_data <- content(initial_pull, as = "parsed")

json_data <- fromJSON(text_data)

state_data <- as.data.frame(json_data$features$attributes) 

reg1 <- as.data.frame(c("Carter", "Custer", "Daniels", "Dawson", "Fallon", "Garfield", "McCone",
                        "Phillips", "Powder River", "Prairie", "Richland", "Roosevelt", "Rosebud",
                        "Sheridan", "Treasure", "Valley", "Wibaux")) %>% 
   rename(county = 1) %>% 
   mutate(region = 1)

reg2 <- as.data.frame(c("Blaine", "Cascade", "Chouteau", "Glacier", "Hill", "Liberty", "Pondera",
                        "Teton", "Toole")) %>% 
   rename(county = 1) %>% 
   mutate(region = 2)

reg3 <- as.data.frame(c("Big Horn", "Carbon", "Fergus", "Golden Valley", "Judith Basin",
                        "Musselshell", "Petroleum", "Stillwater", "Sweet Grass", "Wheatland",
                        "Yellowstone")) %>% 
   rename(county = 1) %>% 
   mutate(region = 3)

reg4 <- as.data.frame(c("Beaverhead", "Broadwater", "Deer Lodge", "Gallatin", "Granite", "Jefferson",
                        "Lewis and Clark", "Madison", "Meagher", "Park", "Powell", "Silver Bow")) %>% 
   rename(county = 1) %>% 
   mutate(region = 4)

reg5 <- as.data.frame(c("Flathead", "Lake", "Lincoln", "Mineral", "Missoula", 
                        "Ravalli", "Sanders")) %>% 
   rename(county = 1) %>% 
   mutate(region = 5)

counties_regions <- rbind(reg1, reg2, reg3, reg4, reg5)

state_data_clean <- state_data %>% 
   rename_all(tolower) %>% 
   mutate(date_reported_to_cdepi = date_reported_to_cdepi/1000,
          dates = as.POSIXct(date_reported_to_cdepi, origin = "1970-01-01")) %>% 
   separate(dates, c("dates", "trash"), sep = " ") %>% 
   mutate(dates = ymd(dates)) %>% 
   select(case_no, dates, county:mt_case) %>% 
   left_join(counties_regions, by = "county")

state_data_incidence <- state_data %>% 
   rename_all(tolower) %>% 
   mutate(date_reported_to_cdepi = date_reported_to_cdepi/1000,
          dates = as.POSIXct(date_reported_to_cdepi, origin = "1970-01-01")) %>% 
   separate(dates, c("dates", "trash"), sep = " ") %>% 
   mutate(dates = ymd(dates)) %>% 
   select(case_no, dates, county:mt_case) %>% 
   group_by(dates) %>% 
   mutate(I = n()) %>% 
   ungroup() %>% 
   distinct(dates, .keep_all = TRUE) %>% 
   select(dates, I)


# Run analysis and print results
serial_interval_mean <- 4
serial_interval_sd <- 4.75


## State results
state_incidence <- incidence(state_data_clean$dates)

state_results <- estimate_R(state_incidence, method="parametric_si", 
                            config = make_config(list(mean_si = serial_interval_mean, 
                                                      std_si = serial_interval_sd)))

state_plot <- estimate_R_plots(state_results, what = "R")

ggsave("C:/R/covid19/state_plot.png")

colomns <- c(1:4, 8)
state_r <- (state_results$R[,colomns]) %>% 
   mutate(region = "state")


## Reg1 results
reg1_data <- state_data_clean %>% 
   filter(region == 1)

reg1_incidence <- incidence(reg1_data$dates)

reg1_results <- estimate_R(reg1_incidence, method="parametric_si", 
                           config = make_config(list(mean_si = serial_interval_mean, 
                                                     std_si = serial_interval_sd)))

reg1_plot <- estimate_R_plots(reg1_results, what = "R")

ggsave("C:/R/covid19/reg1_plot.png")

colomns <- c(1:4, 8)
reg1_r <- (reg1_results$R[,colomns]) %>% 
   mutate(region = "1")


## Reg2 results
reg2_data <- state_data_clean %>% 
   filter(region == 2)

reg2_incidence <- incidence(reg2_data$dates)

reg2_results <- estimate_R(reg2_incidence, method="parametric_si", 
                           config = make_config(list(mean_si = serial_interval_mean, 
                                                     std_si = serial_interval_sd)))

reg2_plot <- estimate_R_plots(reg2_results, what = "R")

ggsave("C:/R/covid19/reg2_plot.png")

colomns <- c(1:4, 8)
reg2_r <- (reg2_results$R[,colomns]) %>% 
   mutate(region = "2")


## Reg3 results
reg3_data <- state_data_clean %>% 
   filter(region == 3)

reg3_incidence <- incidence(reg3_data$dates)

reg3_results <- estimate_R(reg3_incidence, method="parametric_si", 
                           config = make_config(list(mean_si = serial_interval_mean, 
                                                     std_si = serial_interval_sd)))

reg3_plot <- estimate_R_plots(reg3_results, what = "R")

ggsave("C:/R/covid19/reg3_plot.png")

colomns <- c(1:4, 8)
reg3_r <- (reg3_results$R[,colomns]) %>% 
   mutate(region = "3")


## Reg4 results
reg4_data <- state_data_clean %>% 
   filter(region == 4)

reg4_incidence <- incidence(reg4_data$dates)

reg4_results <- estimate_R(reg4_incidence, method="parametric_si", 
                           config = make_config(list(mean_si = serial_interval_mean, 
                                                     std_si = serial_interval_sd)))

reg4_plot <- estimate_R_plots(reg4_results, what = "R")

ggsave("C:/R/covid19/reg4_plot.png")

colomns <- c(1:4, 8)
reg4_r <- (reg4_results$R[,colomns]) %>% 
   mutate(region = "4")


## Reg5 results
reg5_data <- state_data_clean %>% 
   filter(region == 5)

reg5_incidence <- incidence(reg5_data$dates)

reg5_results <- estimate_R(reg5_incidence, method="parametric_si", 
                           config = make_config(list(mean_si = serial_interval_mean, 
                                                     std_si = serial_interval_sd)))

reg5_plot <- estimate_R_plots(reg5_results, what = "R")

ggsave("C:/R/covid19/reg5_plot.png")

colomns <- c(1:4, 8)
reg5_r <- (reg5_results$R[,colomns]) %>% 
   mutate(region = "5")



# Bind files and save
all_regions_r <- rbind(state_r, reg1_r, reg2_r, reg3_r, reg4_r, reg5_r)

write_csv(all_regions_r, "C:/R/covid19/all_regions_r.csv")


