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

query_url <- "https://services.arcgis.com/qnjIrwR8z5Izc0ij/ArcGIS/rest/services/COVID_Cases_Production_View/FeatureServer/2/query?where=1%3D1&objectIds=&time=&resultType=none&outFields=OBJECTID%2C+Case_No%2C+Date_Reported_to_CDEpi%2C+County%2C+Age_Group%2C+Sex%2C+Hospitalization%2C+Outcome%2C+MT_case&returnIdsOnly=false&returnUniqueIdsOnly=false&returnCountOnly=false&returnDistinctValues=false&cacheHint=false&orderByFields=&groupByFieldsForStatistics=&outStatistics=&having=&resultOffset=&resultRecordCount=&sqlFormat=none&f=pjson&token="

initial_pull <- GET(query_url)

text_data <- content(initial_pull, as = "text")

parsed_data <- content(initial_pull, as = "parsed")

json_data <- fromJSON(text_data)

state_data <- as.data.frame(json_data$features$attributes) 

state_data_clean <- state_data %>% 
   rename_all(tolower) %>% 
   mutate(date_reported_to_cdepi = date_reported_to_cdepi/1000,
          dates = as.POSIXct(date_reported_to_cdepi, origin = "1970-01-01")) %>% 
   separate(dates, c("dates", "trash"), sep = " ") %>% 
   mutate(dates = ymd(dates)) %>% 
   select(case_no, dates, county:mt_case) 

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



serial_interval_mean <- 4
serial_interval_sd <- 4.75

state_incidence <- incidence(state_data_clean$dates)

estimate_results <- estimate_R(state_incidence, method="parametric_si", 
                               config = make_config(list(mean_si = serial_interval_mean, 
                                                         std_si = serial_interval_sd)))

estimate_plots <- estimate_R_plots(estimate_results, what = "R")
#estimate_plots

colomns <- c(1:4, 8)
save_r <- tail(estimate_results$R[,colomns], 10)

write_csv(save_r, "C:/R/covid19/save_r.csv")
