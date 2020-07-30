# Pull daily updated State of Montana covid data
# Print results for counties, hospitalizations, and case data
# Ethan Walker
# 1 July 2020

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
library(Rsftp)
#library(sftp)


############## Query and prep data
query_url <- "https://services.arcgis.com/qnjIrwR8z5Izc0ij/ArcGIS/rest/services/COVID_Cases_Production_View/FeatureServer/2/query?where=1%3D1&objectIds=&time=&resultType=none&outFields=OBJECTID%2C+Case_No%2C+Date_Reported_to_CDEpi%2C+County%2C+Age_Group%2C+Sex%2C+Hospitalization%2C+Outcome%2C+MT_case&returnIdsOnly=false&returnUniqueIdsOnly=false&returnCountOnly=false&returnDistinctValues=false&cacheHint=false&orderByFields=&groupByFieldsForStatistics=&outStatistics=&having=&resultOffset=&resultRecordCount=&sqlFormat=none&f=pjson&token="

initial_pull <- GET(query_url)

text_data <- content(initial_pull, as = "text")

parsed_data <- content(initial_pull, as = "parsed")

json_data <- fromJSON(text_data)

state_data <- as.data.frame(json_data$features$attributes) 

reg1 <- as.data.frame(c("Carter", "Custer", "Daniels", "Dawson", "Fallon", "Garfield", "McCone",
                        "Phillips", "Powder River", "Prairie", "Richland", "Roosevelt", "Rosebud",
                        "Sheridan", "Treasure", "Valley", "Wibaux")) %>% 
   rename_at(1, ~"county") %>% 
   mutate(region = 1)

reg2 <- as.data.frame(c("Blaine", "Cascade", "Chouteau", "Glacier", "Hill", "Liberty", "Pondera",
                        "Teton", "Toole")) %>% 
   rename_at(1, ~"county") %>% 
   mutate(region = 2)

reg3 <- as.data.frame(c("Big Horn", "Carbon", "Fergus", "Golden Valley", "Judith Basin",
                        "Musselshell", "Petroleum", "Stillwater", "Sweet Grass", "Wheatland",
                        "Yellowstone")) %>% 
   rename_at(1, ~"county") %>% 
   mutate(region = 3)

reg4 <- as.data.frame(c("Beaverhead", "Broadwater", "Deer Lodge", "Gallatin", "Granite", "Jefferson",
                        "Lewis and Clark", "Madison", "Meagher", "Park", "Powell", "Silver Bow")) %>% 
   rename_at(1, ~"county") %>% 
   mutate(region = 4)

reg5 <- as.data.frame(c("Flathead", "Lake", "Lincoln", "Mineral", "Missoula", 
                        "Ravalli", "Sanders")) %>% 
   rename_at(1, ~"county") %>% 
   mutate(region = 5)

counties_regions <- rbind(reg1, reg2, reg3, reg4, reg5)

state_data_clean <- state_data %>% 
   rename_all(tolower) %>% 
   mutate(date_reported_to_cdepi = date_reported_to_cdepi/1000,
          dates = as.POSIXct(date_reported_to_cdepi, origin = "1970-01-01")) %>% 
   separate(dates, c("dates", "trash"), sep = " ") %>% 
   mutate(dates = ymd(dates)) %>% 
   select(case_no, dates, county:mt_case) %>% 
   left_join(counties_regions, by = "county") %>% 
   mutate(case = 1) %>% 
   mutate(age_group2 = if_else(age_group == "80-89" | age_group == "90-99",
                               "80+", age_group),
          age_group_new = factor(age_group2, 
                                 levels = c("0-9", "10-19", "20-29", 
                                            "30-39", "40-49", "50-59", 
                                            "60-69", "70-79", "80+"),
                                 labels = c("0 to 9, 0.12", "10 to 19, 0.13", "20 to 29, 0.13", 
                                            "30 to 39, 0.13", "40 to 49, 0.11", "50 to 59, 0.13", 
                                            "60 to 69, 0.14", "70 to 79, 0.08", "80+, 0.04"))) %>% 
   separate(age_group_new, c("age_group_new", "age_group_new_percent"), sep = ",") %>% 
   mutate(age_group_new_percent = as.numeric(age_group_new_percent),
          state_pop = as.numeric(1062000)) %>% 
   select(-age_group2, -age_group_new, -age_group_new_percent, -state_pop) %>% 
   mutate(hospitalization = factor(hospitalization,
                                   levels = c("Y", "N", "P", "U"),
                                   labels = c("Hosp: Yes", "Hosp: No", 
                                              "Hosp: Past", "Hosp: Unknown"))) %>% 
   select(-case_no) %>% 
   rownames_to_column(var = "case_no")

output_path <- c("C:/Users/ethan.walker/Box/Ethan Walker UM/R/covid19/")

write_rds(state_data_clean, paste0(output_path, "Output/state_data_clean.rds"))
write_csv(state_data_clean, paste0(output_path, "Output/state_data_clean.csv"), na = " ")


########## Function to get county-specific outcomes
county_function <- function(county_name, data = state_data_clean){
   
   county_data <- data %>% 
      filter(mt_case != "X") %>% 
      filter(county == county_name) %>% 
      select(-case_no) %>% 
      rownames_to_column(var = "case_no")
   
   county_incidence <- incidence(county_data$dates, 
                                 first_date = "2020-03-13", 
                                 last_date = Sys.Date(),
                                 standard = FALSE)
   
   county_dates <- as.data.frame(county_incidence$dates) %>% 
      rename_at(1, ~"dates")
   
   county_cases <- as.data.frame(county_incidence$counts) %>% 
      rename_at(1, ~"Daily total cases")
   
   county_data_new <<- county_dates %>% 
      left_join(county_data, by = "dates") %>% 
      mutate(case_new = case) %>% 
      group_by(dates) %>% 
      mutate(hospitalization = as.character(hospitalization)) %>% 
      pivot_wider(names_from = "hospitalization", values_from = "case") %>% 
      select(-"NA") %>% 
      group_by(dates) %>% 
      pivot_wider(names_from = "outcome", values_from = "case_new") %>% 
      select(-"NA") %>% 
      group_by(dates) %>% 
      mutate_at(vars(one_of("Hosp: No")), sum, na.rm = TRUE) %>% 
      mutate_at(vars(one_of("Hosp: Yes")), sum, na.rm = TRUE) %>% 
      mutate_at(vars(one_of("Hosp: Unkown")), sum, na.rm = TRUE) %>% 
      mutate_at(vars(one_of("Hosp: Past")), sum, na.rm = TRUE) %>% 
      mutate_at(vars(one_of("Active")), sum, na.rm = TRUE) %>% 
      mutate_at(vars(one_of("Deceased")), sum, na.rm = TRUE) %>% 
      mutate_at(vars(one_of("Recovered")), sum, na.rm = TRUE) %>% 
      distinct(dates, .keep_all = TRUE) %>% 
      ungroup() %>% 
      fill(county, .direction = c("downup")) %>% 
      fill(region, .direction = c("downup")) %>% 
      select(-case_no, -age_group, -sex, -mt_case)
   
   county_data_final <- cbind(county_data_new, county_cases) 
   
}

missoula_county <- county_function("Missoula")
gallatin_county <- county_function("Gallatin")
yellowstone_county <- county_function("Yellowstone")
richland_county <- county_function("Richland")
roosevelt_county <- county_function("Roosevelt")
cascade_county <- county_function("Cascade")
glacier_county <- county_function("Glacier")
hill_county <- county_function("Hill")
liberty_county <- county_function("Liberty")
pondera_county <- county_function("Pondera")
toole_county <- county_function("Toole")
bighorn_county <- county_function("Big Horn")
carbon_county <- county_function("Carbon")
fergus_county <- county_function("Fergus")
gvalley_county <- county_function("Golden Valley")
musselshell_county <- county_function("Musselshell")
stillwater_county <- county_function("Stillwater")
wheatland_county <- county_function("Wheatland")
beaverhead_county <- county_function("Beaverhead")
broadwater_county <- county_function("Broadwater")
deerlodge_county <- county_function("Deer Lodge")
jefferson_county <- county_function("Jefferson")
lewisandclark_county <- county_function("Lewis and Clark")
madison_county <- county_function("Madison")
meagher_county <- county_function("Meagher")
park_county <- county_function("Park")
silverbow_county <- county_function("Silver Bow")
flathead_county <- county_function("Flathead")
lake_county <- county_function("Lake")
lincoln_county <- county_function("Lincoln")
ravalli_county <- county_function("Ravalli")
rosebud_county <- county_function("Rosebud")
custer_county <- county_function("Custer")
dawson_county <- county_function("Dawson")
valley_county <- county_function("Valley")
treasure_county <- county_function("Treasure")
granite_county <- county_function("Granite")
fallon_county <- county_function("Fallon")
garield_county <- county_function("Garfield")
sheridan_county <- county_function("Sheridan")
teton_county <- county_function("Teton")
daniels_county <- county_function("Daniels")
mccone_county <- county_function("McCone")
priver_county <- county_function("Powder River")
wibaux_county <- county_function("Wibaux")
blaine_county <- county_function("Blaine")
chouteau_county <- county_function("Chouteau")
jbasin_county <- county_function("Judith Basin")
sweetgrass_county <- county_function("Sweet Grass")
sanders_county <- county_function("Sanders")
powell_county <- county_function("Powell")




county_data_combined <- plyr::rbind.fill(missoula_county, gallatin_county,
                                         yellowstone_county, richland_county,
                                         roosevelt_county, cascade_county,
                                         glacier_county, hill_county,
                                         liberty_county, pondera_county,
                                         toole_county, bighorn_county,
                                         carbon_county, fergus_county,
                                         gvalley_county, musselshell_county,
                                         stillwater_county, wheatland_county,
                                         beaverhead_county, broadwater_county,
                                         deerlodge_county, jefferson_county,
                                         lewisandclark_county, madison_county,
                                         meagher_county, park_county,
                                         silverbow_county, flathead_county,
                                         lake_county, lincoln_county,
                                         ravalli_county, rosebud_county,
                                         custer_county, dawson_county,
                                         valley_county, treasure_county,
                                         granite_county, fallon_county,
                                         garield_county, sheridan_county,
                                         teton_county, daniels_county,
                                         mccone_county, priver_county, 
                                         wibaux_county, blaine_county,
                                         chouteau_county, jbasin_county,
                                         sweetgrass_county, sanders_county,
                                         powell_county) %>% 
   group_by(county) %>% 
   mutate("Cumulative total cases" = cumsum(`Daily total cases`),
          "Cumulative active cases" = cumsum(Active),
          "Cumulative recovered cases" = cumsum(Recovered),
          "Cumulative deseased cases" = cumsum(Deceased))

write_csv(county_data_combined, "C:/R/covid19/state_daily_results/county_data_combined.csv")


#################### Run and save state hospitalization data
state_hosp_data <- state_data %>% 
   rename_all(tolower) %>% 
   mutate(date_reported_to_cdepi = date_reported_to_cdepi/1000,
          dates = as.POSIXct(date_reported_to_cdepi, origin = "1970-01-01")) %>% 
   separate(dates, c("dates", "trash"), sep = " ") %>% 
   mutate(dates = ymd(dates)) %>% 
   select(case_no, dates, county:mt_case) %>% 
   left_join(counties_regions, by = "county") %>% 
   mutate(case = 1) %>% 
   mutate(age_group2 = if_else(age_group == "80-89" | age_group == "90-99",
                               "80+", age_group),
          age_group_new = factor(age_group2, 
                                 levels = c("0-9", "10-19", "20-29", 
                                            "30-39", "40-49", "50-59", 
                                            "60-69", "70-79", "80+"),
                                 labels = c("0 to 9", "10 to 19", "20 to 29", 
                                            "30 to 39", "40 to 49", "50 to 59", 
                                            "60 to 69", "70 to 79", "80+"))) %>% 
   select(-age_group2, -age_group) %>% 
   mutate(hospitalization = factor(hospitalization,
                                   levels = c("Y", "N", "P", "U"),
                                   labels = c("Hosp: Yes", "Hosp: No", 
                                              "Hosp: Past", "Hosp: Unknown")))

state_hosp <- state_hosp_data %>% 
   filter(mt_case != "X") %>% 
   group_by(age_group_new) %>% 
   mutate(age_group_cases = sum(case)) %>% 
   mutate(hosp = if_else(hospitalization == "Hosp: Yes" | hospitalization == "Hosp: Past", 1, 0),
          hosp_yes = sum(hosp, na.rm = TRUE),
          hosp_no = age_group_cases - hosp_yes,
          hosp_percent = round(hosp_yes/age_group_cases*100, digits = 2)) %>% 
   ungroup() %>% 
   distinct(age_group_new, .keep_all = TRUE) %>% 
   select(age_group_new, age_group_cases, hosp_yes, hosp_no, hosp_percent) %>% 
   rename(age_group = age_group_new) %>% 
   arrange(age_group)

write_csv(state_hosp, "C:/R/covid19/state_daily_results/state_hosp.csv", na = " ")


#################### Run and save daily case, hosp, outcome, test data
query_url <- "https://services.arcgis.com/qnjIrwR8z5Izc0ij/ArcGIS/rest/services/COVID_Cases_Production_View/FeatureServer/1/query?where=1%3D1&objectIds=&time=&resultType=none&outFields=OBJECTID%2C+Total_Tests_Completed%2C+New_Tests_Completed%2C+Test_Date%2C+ScriptRunDate&returnIdsOnly=false&returnUniqueIdsOnly=false&returnCountOnly=false&returnDistinctValues=false&cacheHint=false&orderByFields=&groupByFieldsForStatistics=&outStatistics=&having=&resultOffset=&resultRecordCount=&sqlFormat=none&f=pjson&token="

initial_pull <- GET(query_url)

text_data <- content(initial_pull, as = "text")

parsed_data <- content(initial_pull, as = "parsed")

json_data <- fromJSON(text_data)

state_test_data <- as.data.frame(json_data$features$attributes) 

state_test_data_clean <- state_test_data %>% 
   rename_all(tolower) %>% 
   mutate(dates = test_date/1000,
          dates = as.POSIXct(dates, origin = "1970-01-01")) %>% 
   separate(dates, c("dates", "trash"), sep = " ") %>% 
   mutate(dates = ymd(dates)) %>% 
   select(dates, total_tests_completed, new_tests_completed) %>% 
   arrange(dates, desc(total_tests_completed)) %>% 
   distinct(dates, .keep_all = TRUE)


hosp_data_daily <- state_data_clean %>% 
   filter(mt_case != "X") %>% 
   select(dates, hospitalization, case) %>% 
   mutate(hospitalization = factor(hospitalization,
                                   labels = c("hosp_active",
                                              "hosp_no",
                                              "hosp_past",
                                              "hosp_unknown"))) %>% 
   group_by(dates, hospitalization) %>% 
   mutate(sum_hosp = sum(case)) %>% 
   group_by(dates) %>% 
   select(-case) %>% 
   distinct(dates, hospitalization, .keep_all = TRUE) %>% 
   pivot_wider(names_from = "hospitalization", values_from = "sum_hosp") %>% 
   mutate(hosp_active = if_else(is.na(hosp_active), 0, hosp_active),
          hosp_no = if_else(is.na(hosp_no), 0, hosp_no),
          hosp_past = if_else(is.na(hosp_past), 0, hosp_past),
          hosp_unknown = if_else(is.na(hosp_unknown), 0, hosp_unknown)) %>% 
   group_by(dates) %>% 
   mutate(daily_hosp = hosp_active + hosp_past) %>% 
   ungroup() %>% 
   mutate(total_hosp = cumsum(daily_hosp)) %>% 
   select(dates, daily_hosp, total_hosp)


outcome_data_daily <- state_data_clean %>% 
   filter(mt_case != "X") %>% 
   select(dates, outcome, case) %>% 
   mutate(outcome = factor(outcome,
                           labels = c("active",
                                      "deceased",
                                      "recovered"))) %>% 
   group_by(dates, outcome) %>% 
   mutate(total_outcome = sum(case)) %>% 
   group_by(dates) %>% 
   select(-case) %>% 
   distinct(dates, outcome, .keep_all = TRUE) %>% 
   pivot_wider(names_from = "outcome", values_from = "total_outcome") %>% 
   mutate(active = if_else(is.na(active), 0, active),
          deceased = if_else(is.na(deceased), 0, deceased),
          recovered = if_else(is.na(recovered), 0, recovered)) %>% 
   group_by(dates) %>% 
   mutate(daily_deceased = deceased,
          daily_active = active,
          daily_recovered = recovered) %>% 
   ungroup() %>% 
   mutate(total_deceased = cumsum(daily_deceased),
          total_active = cumsum(daily_active),
          total_recovered = cumsum(daily_recovered)) %>% 
   select(dates, daily_deceased, total_deceased, daily_active,
          total_active, daily_recovered, total_recovered)


case_data_daily <- state_data_clean %>% 
   filter(mt_case != "X") %>% 
   select(dates, case) %>% 
   group_by(dates) %>% 
   mutate(daily_cases = sum(case)) %>% 
   ungroup() %>% 
   distinct(dates, .keep_all = TRUE) %>% 
   mutate(total_cases = cumsum(daily_cases)) %>% 
   select(-case)


case_hosp_test_outcome <- case_data_daily %>% 
   left_join(hosp_data_daily, by = "dates") %>% 
   left_join(outcome_data_daily, by = "dates") %>% 
   full_join(state_test_data_clean, by = "dates") %>% 
   arrange(dates) %>% 
   mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>% 
   mutate(total_cases = cumsum(daily_cases),
          total_hosp = cumsum(daily_hosp),
          total_deceased = cumsum(daily_deceased),
          total_active = cumsum(daily_active),
          total_recovered = cumsum(daily_recovered),
          total_tests_completed = if_else(total_tests_completed == 0, 
                                          lag(total_tests_completed), total_tests_completed),
          total_tests_completed = if_else(total_tests_completed == 0, 
                                          lag(total_tests_completed), total_tests_completed),
          total_tests_completed = if_else(is.na(total_tests_completed), 
                                          0, total_tests_completed)) %>% 
   rename(total_tests = total_tests_completed,
          daily_tests = new_tests_completed) 

write_csv(case_hosp_test_outcome, "C:/R/covid19/state_daily_results/mt_case_outcome_hosp_data.csv", na = " ")


# Send files to the sftp server
# host = 'elbastion.dbs.umt.edu'
# port = 22
# username = 'celftp'
# password  = 'celftp'
# remotepath = '/celFtpFiles/covid19/Rt/incoming/'

sftpUpload("elbastion.dbs.umt.edu", "celftp", "celftp",
           "/celFtpFiles/covid19/Rt/incoming/state_hosp.csv",
           "C:/R/covid19/state_daily_results/state_hosp.csv")

sftpUpload("elbastion.dbs.umt.edu", "celftp", "celftp",
           "/celFtpFiles/covid19/Rt/incoming/county_data_combined.csv",
           "C:/R/covid19/state_daily_results/county_data_combined.csv")

sftpUpload("elbastion.dbs.umt.edu", "celftp", "celftp",
           "/celFtpFiles/covid19/Rt/incoming/mt_case_outcome_hosp_data.csv",
           "C:/R/covid19/state_daily_results/mt_case_outcome_hosp_data.csv")
