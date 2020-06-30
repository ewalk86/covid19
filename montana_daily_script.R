# Pull daily updated State of Montana covid data, estimate R, save results/plots
# Ethan Walker
# 16 June 2020

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
                                              "Hosp: Past", "Hosp: Unknown")))

state_data_wide <- state_data_clean %>% 
   mutate(age_group2 = if_else(age_group == "80-89" | age_group == "90-99",
                               "80+", age_group),
          age_group_new = factor(age_group2, 
                                 levels = c("0-9", "10-19", "20-29", 
                                            "30-39", "40-49", "50-59", 
                                            "60-69", "70-79", "80+"),
                                 labels = c("0 to 9", "10 to 19", "20 to 29", 
                                            "30 to 39", "40 to 49", "50 to 59", 
                                            "60 to 69", "70 to 79", "80+")),
          sex = factor(sex, labels = c("Female", "Male"))) %>% 
   select(-age_group, -age_group2) %>% 
   mutate(case = 1) %>% 
   pivot_wider(names_from = "age_group_new", values_from = "case") %>% 
   mutate(case = 1) %>% 
   pivot_wider(names_from = "sex", values_from = "case") %>% 
   mutate(case = 1) %>% 
   pivot_wider(names_from = "hospitalization", values_from = "case") %>% 
   mutate(case = 1) %>% 
   pivot_wider(names_from = "county", values_from = "case")

state_wide_date <- state_data_wide %>% 
   mutate(region = "state") %>% 
   group_by(dates) %>% 
   mutate_at(c(6:51), sum, na.rm = TRUE) %>% 
   mutate(case = 1,
          daily_cases = sum(case)) %>% 
   arrange(dates) %>% 
   distinct(dates, .keep_all = TRUE) %>% 
   ungroup() %>% 
   mutate(cumulative_cases = cumsum(daily_cases)) %>% 
   select(region, dates, daily_cases, cumulative_cases, `70 to 79`:Fergus)

reg1_wide_date <- state_data_wide %>% 
   filter(region == 1) %>% 
   group_by(dates) %>% 
   mutate_at(c(6:51), sum, na.rm = TRUE) %>% 
   mutate(case = 1,
          daily_cases = sum(case)) %>% 
   arrange(dates) %>% 
   distinct(dates, .keep_all = TRUE) %>% 
   ungroup() %>% 
   mutate(cumulative_cases = cumsum(daily_cases)) %>% 
   select(region, dates, daily_cases, cumulative_cases, `70 to 79`:Fergus)

reg2_wide_date <- state_data_wide %>% 
   filter(region == 2) %>% 
   group_by(dates) %>% 
   mutate_at(c(6:51), sum, na.rm = TRUE) %>% 
   mutate(case = 1,
          daily_cases = sum(case)) %>% 
   arrange(dates) %>% 
   distinct(dates, .keep_all = TRUE) %>% 
   ungroup() %>% 
   mutate(cumulative_cases = cumsum(daily_cases)) %>% 
   select(region, dates, daily_cases, cumulative_cases, `70 to 79`:Fergus)

reg3_wide_date <- state_data_wide %>% 
   filter(region == 3) %>% 
   group_by(dates) %>% 
   mutate_at(c(6:51), sum, na.rm = TRUE) %>% 
   mutate(case = 1,
          daily_cases = sum(case)) %>% 
   arrange(dates) %>% 
   distinct(dates, .keep_all = TRUE) %>% 
   ungroup() %>% 
   mutate(cumulative_cases = cumsum(daily_cases)) %>% 
   select(region, dates, daily_cases, cumulative_cases, `70 to 79`:Fergus)

reg4_wide_date <- state_data_wide %>% 
   filter(region == 4) %>% 
   group_by(dates) %>% 
   mutate_at(c(6:51), sum, na.rm = TRUE) %>% 
   mutate(case = 1,
          daily_cases = sum(case)) %>% 
   arrange(dates) %>% 
   distinct(dates, .keep_all = TRUE) %>% 
   ungroup() %>% 
   mutate(cumulative_cases = cumsum(daily_cases)) %>% 
   select(region, dates, daily_cases, cumulative_cases, `70 to 79`:Fergus)

reg5_wide_date <- state_data_wide %>% 
   filter(region == 5) %>% 
   group_by(dates) %>% 
   mutate_at(c(6:51), sum, na.rm = TRUE) %>% 
   mutate(case = 1,
          daily_cases = sum(case)) %>% 
   arrange(dates) %>% 
   distinct(dates, .keep_all = TRUE) %>% 
   ungroup() %>% 
   mutate(cumulative_cases = cumsum(daily_cases)) %>% 
   select(region, dates, daily_cases, cumulative_cases, `70 to 79`:Fergus)

all_data_wide <- rbind(state_wide_date, reg1_wide_date, reg2_wide_date,
                       reg3_wide_date, reg4_wide_date, reg5_wide_date)

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


########## Function to get county-specific outcomes
county_function <- function(county_name, data = state_data_clean){
   
   county_data <- data %>% 
      filter(county == county_name)
   
   county_incidence <- incidence(county_data$dates, 
                                 first_date = "2020-03-13", 
                                 last_date = Sys.Date(),
                                 standard = FALSE)
   county_dates <- as.data.frame(county_incidence$dates) %>% 
      rename_at(1, ~"dates")
   county_cases <- county_incidence$counts
   
   county_data_new <<- county_dates %>% 
      left_join(county_data, by = "dates") %>% 
      mutate(case_new = case) %>% 
      group_by(dates) %>% 
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
                                         granite_county) 

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
   filter(mt_case != "N" & mt_case != "X") %>% 
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
   filter(mt_case == "Y") %>% 
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
   filter(mt_case == "Y") %>% 
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
   filter(mt_case == "Y") %>% 
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
                                          0, total_tests_completed)) 

write_csv(case_hosp_test_outcome, "C:/R/covid19/state_daily_results/mt_case_outcome_hosp_data.csv", na = " ")



#################### Run analysis and print results
serial_interval_mean <- 4
serial_interval_sd <- 4.75


## State results
current_date <- format(Sys.Date() - 1, "%Y-%m-%d")
state_incidence <- incidence(state_data_clean$dates, last_date = current_date)

state_n <- as.data.frame(state_incidence$dates)
time_var <- nrow(state_n)
time_start <- seq(2, time_var-13)
time_end <- time_start + 13

state_results <- estimate_R(state_incidence, method="parametric_si", 
                            config = make_config(list(mean_si = serial_interval_mean, 
                                                      std_si = serial_interval_sd,
                                                      t_start =  time_start,
                                                      t_end = time_end)))

#state_plot <- estimate_R_plots(state_results, what = "R")

colomns <- c(1:4, 8)
state_r <- (state_results$R[,colomns]) %>% 
   mutate(region = "state") %>% 
   rename(mean_r = `Mean(R)`,
          sd_r = `Std(R)`,
          median_r = `Median(R)`)

state_dates <- as.data.frame(state_results$dates)
state_i <- as.data.frame(state_results$I)
state_cil <- as.data.frame(state_results$R$`Quantile.0.025(R)`)
state_cih <- as.data.frame(state_results$R$`Quantile.0.975(R)`)
state_dates_new <- cbind(state_dates, state_i) %>% 
   rename(dates = 1,
          incidence = 2) %>% 
   mutate(dates = ymd(dates))
state_dates_new <- state_dates_new[-(1:14), 1:2]
state_dates_new <- cbind(state_dates_new, state_cil, state_cih) %>% 
   rename(cl_low = 3,
          cl_high = 4)
state_r_clean <- cbind(state_r, state_dates_new)


total_cases <- sum(state_data_clean$case)
date_today <- format(Sys.Date(), "%d %b %Y")
date_10 <- format(Sys.Date() - 10, "%d %b %Y")


state_r_plot <- state_r_clean %>% 
   ggplot() +
   geom_line(aes(dates, mean_r), size = 1.5, color = "black") +
   geom_line(aes(dates, cl_low), size = 1.5, color = "grey") +
   geom_line(aes(dates, cl_high), size = 1.5, color = "grey") +
   labs(title = "COVID-19 Rolling 14-day R-values, State of Montana, 2020",
        subtitle = paste0("Data current as of ", date_today),
        color = "") +
   ylab("R-value") +
   xlab("") +
   geom_hline(yintercept = 1, color = "black", size = 1.2) +
   scale_x_date(date_breaks = "2 days", date_labels = "%d-%b") +
   scale_y_continuous(breaks = seq(0, 10, 0.5), labels = seq(0, 10, 0.5)) +
   theme_minimal() +
   theme(strip.text = element_text(size = 16, colour = "black"),
         title = element_text(size = 12, colour = "black"),
         panel.grid = element_blank(),
         panel.grid.major.y = element_line(colour = "grey"),
         axis.text.x = element_text(size = 12, colour = "black", 
                                    angle = 90, vjust = 0.4),
         axis.text.y = element_text(size = 12, colour = "black"),
         legend.text = element_text(size = 12, colour = "black"),
         axis.title.y = element_text(size = 12, colour = "black",
                                     margin = unit(c(0, 5, 0, 0), "mm")),
         axis.title.x = element_text(size = 12, colour = "black",
                                     margin = unit(c(5, 0, 0, 0), "mm")),
         axis.line.x = element_blank(), 
         axis.line.y = element_blank(), 
         axis.ticks = element_blank()) +
   scale_color_manual(values = c("black")) 
#state_r_plot
ggsave("C:/R/covid19/state_daily_results/state_r_plot.png", width = 10, height = 8)

state_inc_plot <- state_data_clean %>% 
   ggplot() +
   geom_col(aes(dates, case), 
            fill = "steelblue") +
   labs(title = paste0("COVID-19 Cases in State of Montana, 2020", 
                       " (Total cases = ", total_cases, ")"),
        subtitle = paste0("Data current as of ", date_today)) +
   ylab("Number of Cases") +
   xlab("") +
   scale_x_date(date_breaks = "2 day", date_labels = "%d-%b") +
   theme_minimal() +
   theme(strip.text = element_text(size = 16, colour = "black"),
         title = element_text(size = 12, colour = "black"),
         panel.grid = element_blank(),
         panel.grid.major.y = element_line(colour = "grey"),
         axis.text.x = element_text(size = 12, colour = "black", 
                                    angle = 90, vjust = 0.4),
         axis.text.y = element_text(size = 12, colour = "black"),
         axis.title.y = element_text(size = 12, colour = "black",
                                     margin = unit(c(0, 5, 0, 0), "mm")),
         axis.title.x = element_text(size = 12, colour = "black",
                                     margin = unit(c(5, 0, 0, 0), "mm")),
         axis.line.x = element_blank(), 
         axis.line.y = element_blank(), 
         axis.ticks = element_blank()) 
#state_inc_plot

ggsave("C:/R/covid19/state_daily_results/state_inc_plot.png", width = 10, height = 8)



## Reg1 results
reg1_data <- state_data_clean %>% 
   filter(region == 1)

reg1_incidence <- incidence(reg1_data$dates, last_date = current_date)

reg1_n <- as.data.frame(reg1_incidence$dates)
time_var <- nrow(reg1_n)
time_start <- seq(2, time_var-13)
time_end <- time_start + 13

reg1_results <- estimate_R(reg1_incidence, method="parametric_si", 
                           config = make_config(list(mean_si = serial_interval_mean, 
                                                     std_si = serial_interval_sd,
                                                     t_start =  time_start,
                                                     t_end = time_end)))

#reg1_plot <- estimate_R_plots(reg1_results, what = "R")


colomns <- c(1:4, 8)
reg1_r <- (reg1_results$R[,colomns]) %>% 
   mutate(region = "1") %>% 
   rename(mean_r = `Mean(R)`,
          sd_r = `Std(R)`,
          median_r = `Median(R)`)

reg1_dates <- as.data.frame(reg1_results$dates)
reg1_i <- as.data.frame(reg1_results$I)
reg1_cil <- as.data.frame(reg1_results$R$`Quantile.0.025(R)`)
reg1_cih <- as.data.frame(reg1_results$R$`Quantile.0.975(R)`)
reg1_dates_new <- cbind(reg1_dates, reg1_i) %>% 
   rename(dates = 1,
          incidence = 2) %>% 
   mutate(dates = ymd(dates))
reg1_dates_new <- reg1_dates_new[-(1:14), 1:2]
reg1_dates_new <- cbind(reg1_dates_new, reg1_cil, reg1_cih) %>% 
   rename(cl_low = 3,
          cl_high = 4)
reg1_r_clean <- cbind(reg1_r, reg1_dates_new)

reg1_data_clean <- state_data_clean %>% 
   filter(region == 1)
total_cases <- sum(reg1_data_clean$case)


reg1_r_plot <- reg1_r_clean %>% 
   ggplot() +
   geom_line(aes(dates, mean_r), size = 1.5, color = "black") +
   geom_line(aes(dates, cl_low), size = 1.5, color = "grey") +
   geom_line(aes(dates, cl_high), size = 1.5, color = "grey") +
   labs(title = "COVID-19 Rolling 14-day R-values, Montana Region 1, 2020",
        subtitle = paste0("Data current as of ", date_today),
        color = "") +
   ylab("R-value") +
   xlab("") +
   geom_hline(yintercept = 1, color = "black", size = 1.2) +
   scale_x_date(date_breaks = "2 days", date_labels = "%d-%b") +
   #scale_y_continuous(breaks = seq(0, 20, 2), labels = seq(0, 20, 2)) +
   theme_minimal() +
   theme(strip.text = element_text(size = 16, colour = "black"),
         title = element_text(size = 12, colour = "black"),
         panel.grid = element_blank(),
         panel.grid.major.y = element_line(colour = "grey"),
         axis.text.x = element_text(size = 12, colour = "black", 
                                    angle = 90, vjust = 0.4),
         axis.text.y = element_text(size = 12, colour = "black"),
         legend.text = element_text(size = 12, colour = "black"),
         axis.title.y = element_text(size = 12, colour = "black",
                                     margin = unit(c(0, 5, 0, 0), "mm")),
         axis.title.x = element_text(size = 12, colour = "black",
                                     margin = unit(c(5, 0, 0, 0), "mm")),
         axis.line.x = element_blank(), 
         axis.line.y = element_blank(), 
         axis.ticks = element_blank()) +
   scale_color_manual(values = c("black")) 
#reg1_r_plot
ggsave("C:/R/covid19/state_daily_results/reg1_r_plot.png", width = 10, height = 8)

reg1_inc_plot <- reg1_data_clean %>% 
   ggplot() +
   geom_col(aes(dates, case), 
            fill = "steelblue") +
   labs(title = paste0("COVID-19 Cases in Montana Region 1, 2020", 
                       " (Total cases = ", total_cases, ")"),
        subtitle = paste0("Data current as of ", date_today)) +
   ylab("Number of Cases") +
   xlab("") +
   scale_x_date(date_breaks = "2 day", date_labels = "%d-%b") +
   theme_minimal() +
   theme(strip.text = element_text(size = 16, colour = "black"),
         title = element_text(size = 12, colour = "black"),
         panel.grid = element_blank(),
         panel.grid.major.y = element_line(colour = "grey"),
         axis.text.x = element_text(size = 12, colour = "black", 
                                    angle = 90, vjust = 0.4),
         axis.text.y = element_text(size = 12, colour = "black"),
         axis.title.y = element_text(size = 12, colour = "black",
                                     margin = unit(c(0, 5, 0, 0), "mm")),
         axis.title.x = element_text(size = 12, colour = "black",
                                     margin = unit(c(5, 0, 0, 0), "mm")),
         axis.line.x = element_blank(), 
         axis.line.y = element_blank(), 
         axis.ticks = element_blank()) 
#reg1_inc_plot

ggsave("C:/R/covid19/state_daily_results/reg1_inc_plot.png", width = 10, height = 8)



## Reg2 results
reg2_data <- state_data_clean %>% 
   filter(region == 2)

reg2_incidence <- incidence(reg2_data$dates, last_date = current_date)

reg2_n <- as.data.frame(reg2_incidence$dates)
time_var <- nrow(reg2_n)
time_start <- seq(2, time_var-13)
time_end <- time_start + 13

reg2_results <- estimate_R(reg2_incidence, method="parametric_si", 
                           config = make_config(list(mean_si = serial_interval_mean, 
                                                     std_si = serial_interval_sd,
                                                     t_start =  time_start,
                                                     t_end = time_end)))

#reg2_plot <- estimate_R_plots(reg2_results, what = "R")


colomns <- c(1:4, 8)
reg2_r <- (reg2_results$R[,colomns]) %>% 
   mutate(region = "2") %>% 
   rename(mean_r = `Mean(R)`,
          sd_r = `Std(R)`,
          median_r = `Median(R)`)

reg2_dates <- as.data.frame(reg2_results$dates)
reg2_i <- as.data.frame(reg2_results$I)
reg2_cil <- as.data.frame(reg2_results$R$`Quantile.0.025(R)`)
reg2_cih <- as.data.frame(reg2_results$R$`Quantile.0.975(R)`)
reg2_dates_new <- cbind(reg2_dates, reg2_i) %>% 
   rename(dates = 1,
          incidence = 2) %>% 
   mutate(dates = ymd(dates))
reg2_dates_new <- reg2_dates_new[-(1:14), 1:2]
reg2_dates_new <- cbind(reg2_dates_new, reg2_cil, reg2_cih) %>% 
   rename(cl_low = 3,
          cl_high = 4)
reg2_r_clean <- cbind(reg2_r, reg2_dates_new)


reg2_data_clean <- state_data_clean %>% 
   filter(region == 2)
total_cases <- sum(reg2_data_clean$case)


reg2_r_plot <- reg2_r_clean %>% 
   ggplot() +
   geom_line(aes(dates, mean_r), size = 1.5, color = "black") +
   geom_line(aes(dates, cl_low), size = 1.5, color = "grey") +
   geom_line(aes(dates, cl_high), size = 1.5, color = "grey") +
   labs(title = "COVID-19 Rolling 14-day R-values, Montana Region 2, 2020",
        subtitle = paste0("Data current as of ", date_today),
        color = "") +
   ylab("R-value") +
   xlab("") +
   geom_hline(yintercept = 1, color = "black", size = 1.2) +
   scale_x_date(date_breaks = "2 days", date_labels = "%d-%b") +
   #scale_y_continuous(breaks = seq(0, 15, 0.5), labels = seq(0, 15, 0.5)) +
   theme_minimal() +
   theme(strip.text = element_text(size = 16, colour = "black"),
         title = element_text(size = 12, colour = "black"),
         panel.grid = element_blank(),
         panel.grid.major.y = element_line(colour = "grey"),
         axis.text.x = element_text(size = 12, colour = "black", 
                                    angle = 90, vjust = 0.4),
         axis.text.y = element_text(size = 12, colour = "black"),
         legend.text = element_text(size = 12, colour = "black"),
         axis.title.y = element_text(size = 12, colour = "black",
                                     margin = unit(c(0, 5, 0, 0), "mm")),
         axis.title.x = element_text(size = 12, colour = "black",
                                     margin = unit(c(5, 0, 0, 0), "mm")),
         axis.line.x = element_blank(), 
         axis.line.y = element_blank(), 
         axis.ticks = element_blank()) +
   scale_color_manual(values = c("black")) 
#reg2_r_plot
ggsave("C:/R/covid19/state_daily_results/reg2_r_plot.png", width = 10, height = 8)

reg2_inc_plot <- reg2_data_clean %>% 
   ggplot() +
   geom_col(aes(dates, case), 
            fill = "steelblue") +
   labs(title = paste0("COVID-19 Cases in Montana Region 2, 2020", 
                       " (Total cases = ", total_cases, ")"),
        subtitle = paste0("Data current as of ", date_today)) +
   ylab("Number of Cases") +
   xlab("") +
   scale_x_date(date_breaks = "2 day", date_labels = "%d-%b") +
   theme_minimal() +
   theme(strip.text = element_text(size = 16, colour = "black"),
         title = element_text(size = 12, colour = "black"),
         panel.grid = element_blank(),
         panel.grid.major.y = element_line(colour = "grey"),
         axis.text.x = element_text(size = 12, colour = "black", 
                                    angle = 90, vjust = 0.4),
         axis.text.y = element_text(size = 12, colour = "black"),
         axis.title.y = element_text(size = 12, colour = "black",
                                     margin = unit(c(0, 5, 0, 0), "mm")),
         axis.title.x = element_text(size = 12, colour = "black",
                                     margin = unit(c(5, 0, 0, 0), "mm")),
         axis.line.x = element_blank(), 
         axis.line.y = element_blank(), 
         axis.ticks = element_blank()) 
#reg2_inc_plot

ggsave("C:/R/covid19/state_daily_results/reg2_inc_plot.png", width = 10, height = 8)


## Reg3 results
reg3_data <- state_data_clean %>% 
   filter(region == 3)

reg3_incidence <- incidence(reg3_data$dates, last_date = current_date)

reg3_n <- as.data.frame(reg3_incidence$dates)
time_var <- nrow(reg3_n)
time_start <- seq(2, time_var-13)
time_end <- time_start + 13

reg3_results <- estimate_R(reg3_incidence, method="parametric_si", 
                           config = make_config(list(mean_si = serial_interval_mean, 
                                                     std_si = serial_interval_sd,
                                                     t_start =  time_start,
                                                     t_end = time_end)))

#reg3_plot <- estimate_R_plots(reg3_results, what = "R")

colomns <- c(1:4, 8)
reg3_r <- (reg3_results$R[,colomns]) %>% 
   mutate(region = "3") %>% 
   rename(mean_r = `Mean(R)`,
          sd_r = `Std(R)`,
          median_r = `Median(R)`)

reg3_dates <- as.data.frame(reg3_results$dates)
reg3_i <- as.data.frame(reg3_results$I)
reg3_cil <- as.data.frame(reg3_results$R$`Quantile.0.025(R)`)
reg3_cih <- as.data.frame(reg3_results$R$`Quantile.0.975(R)`)
reg3_dates_new <- cbind(reg3_dates, reg3_i) %>% 
   rename(dates = 1,
          incidence = 2) %>% 
   mutate(dates = ymd(dates))
reg3_dates_new <- reg3_dates_new[-(1:14), 1:2]
reg3_dates_new <- cbind(reg3_dates_new, reg3_cil, reg3_cih) %>% 
   rename(cl_low = 3,
          cl_high = 4)
reg3_r_clean <- cbind(reg3_r, reg3_dates_new)


reg3_data_clean <- state_data_clean %>% 
   filter(region == 3)
total_cases <- sum(reg3_data_clean$case)


reg3_r_plot <- reg3_r_clean %>% 
   ggplot() +
   geom_line(aes(dates, mean_r), size = 1.5, color = "black") +
   geom_line(aes(dates, cl_low), size = 1.5, color = "grey") +
   geom_line(aes(dates, cl_high), size = 1.5, color = "grey") +
   labs(title = "COVID-19 Rolling 14-day R-values, Montana Region 3, 2020",
        subtitle = paste0("Data current as of ", date_today),
        color = "") +
   ylab("R-value") +
   xlab("") +
   geom_hline(yintercept = 1, color = "black", size = 1.2) +
   scale_x_date(date_breaks = "2 days", date_labels = "%d-%b") +
   #scale_y_continuous(breaks = seq(0, 15, 0.5), labels = seq(0, 15, 0.5)) +
   theme_minimal() +
   theme(strip.text = element_text(size = 16, colour = "black"),
         title = element_text(size = 12, colour = "black"),
         panel.grid = element_blank(),
         panel.grid.major.y = element_line(colour = "grey"),
         axis.text.x = element_text(size = 12, colour = "black", 
                                    angle = 90, vjust = 0.4),
         axis.text.y = element_text(size = 12, colour = "black"),
         legend.text = element_text(size = 12, colour = "black"),
         axis.title.y = element_text(size = 12, colour = "black",
                                     margin = unit(c(0, 5, 0, 0), "mm")),
         axis.title.x = element_text(size = 12, colour = "black",
                                     margin = unit(c(5, 0, 0, 0), "mm")),
         axis.line.x = element_blank(), 
         axis.line.y = element_blank(), 
         axis.ticks = element_blank()) +
   scale_color_manual(values = c("black")) 
#reg3_r_plot
ggsave("C:/R/covid19/state_daily_results/reg3_r_plot.png", width = 10, height = 8)

reg3_inc_plot <- reg3_data_clean %>% 
   ggplot() +
   geom_col(aes(dates, case), 
            fill = "steelblue") +
   labs(title = paste0("COVID-19 Cases in Montana Region 3, 2020", 
                       " (Total cases = ", total_cases, ")"),
        subtitle = paste0("Data current as of ", date_today)) +
   ylab("Number of Cases") +
   xlab("") +
   scale_x_date(date_breaks = "2 day", date_labels = "%d-%b") +
   theme_minimal() +
   theme(strip.text = element_text(size = 16, colour = "black"),
         title = element_text(size = 12, colour = "black"),
         panel.grid = element_blank(),
         panel.grid.major.y = element_line(colour = "grey"),
         axis.text.x = element_text(size = 12, colour = "black", 
                                    angle = 90, vjust = 0.4),
         axis.text.y = element_text(size = 12, colour = "black"),
         axis.title.y = element_text(size = 12, colour = "black",
                                     margin = unit(c(0, 5, 0, 0), "mm")),
         axis.title.x = element_text(size = 12, colour = "black",
                                     margin = unit(c(5, 0, 0, 0), "mm")),
         axis.line.x = element_blank(), 
         axis.line.y = element_blank(), 
         axis.ticks = element_blank()) 
#reg3_inc_plot

ggsave("C:/R/covid19/state_daily_results/reg3_inc_plot.png", width = 10, height = 8)


## Reg4 results
reg4_data <- state_data_clean %>% 
   filter(region == 4)

reg4_incidence <- incidence(reg4_data$dates, last_date = current_date)

reg4_n <- as.data.frame(reg4_incidence$dates)
time_var <- nrow(reg4_n)
time_start <- seq(2, time_var-13)
time_end <- time_start + 13

reg4_results <- estimate_R(reg4_incidence, method="parametric_si", 
                           config = make_config(list(mean_si = serial_interval_mean, 
                                                     std_si = serial_interval_sd,
                                                     t_start =  time_start,
                                                     t_end = time_end)))

#reg4_plot <- estimate_R_plots(reg4_results, what = "R")

colomns <- c(1:4, 8)
reg4_r <- (reg4_results$R[,colomns]) %>% 
   mutate(region = "4") %>% 
   rename(mean_r = `Mean(R)`,
          sd_r = `Std(R)`,
          median_r = `Median(R)`)

reg4_dates <- as.data.frame(reg4_results$dates)
reg4_i <- as.data.frame(reg4_results$I)
reg4_cil <- as.data.frame(reg4_results$R$`Quantile.0.025(R)`)
reg4_cih <- as.data.frame(reg4_results$R$`Quantile.0.975(R)`)
reg4_dates_new <- cbind(reg4_dates, reg4_i) %>% 
   rename(dates = 1,
          incidence = 2) %>% 
   mutate(dates = ymd(dates))
reg4_dates_new <- reg4_dates_new[-(1:14), 1:2]
reg4_dates_new <- cbind(reg4_dates_new, reg4_cil, reg4_cih) %>% 
   rename(cl_low = 3,
          cl_high = 4)
reg4_r_clean <- cbind(reg4_r, reg4_dates_new)


reg4_data_clean <- state_data_clean %>% 
   filter(region == 4)
total_cases <- sum(reg4_data_clean$case)


reg4_r_plot <- reg4_r_clean %>% 
   ggplot() +
   geom_line(aes(dates, mean_r), size = 1.5, color = "black") +
   geom_line(aes(dates, cl_low), size = 1.5, color = "grey") +
   geom_line(aes(dates, cl_high), size = 1.5, color = "grey") +
   labs(title = "COVID-19 Rolling 14-day R-values, Montana Region 4, 2020",
        subtitle = paste0("Data current as of ", date_today),
        color = "") +
   ylab("R-value") +
   xlab("") +
   geom_hline(yintercept = 1, color = "black", size = 1.2) +
   scale_x_date(date_breaks = "2 days", date_labels = "%d-%b") +
   #scale_y_continuous(breaks = seq(0, 15, 0.5), labels = seq(0, 15, 0.5)) +
   theme_minimal() +
   theme(strip.text = element_text(size = 16, colour = "black"),
         title = element_text(size = 12, colour = "black"),
         panel.grid = element_blank(),
         panel.grid.major.y = element_line(colour = "grey"),
         axis.text.x = element_text(size = 12, colour = "black", 
                                    angle = 90, vjust = 0.4),
         axis.text.y = element_text(size = 12, colour = "black"),
         legend.text = element_text(size = 12, colour = "black"),
         axis.title.y = element_text(size = 12, colour = "black",
                                     margin = unit(c(0, 5, 0, 0), "mm")),
         axis.title.x = element_text(size = 12, colour = "black",
                                     margin = unit(c(5, 0, 0, 0), "mm")),
         axis.line.x = element_blank(), 
         axis.line.y = element_blank(), 
         axis.ticks = element_blank()) +
   scale_color_manual(values = c("black")) 
#reg4_r_plot
ggsave("C:/R/covid19/state_daily_results/reg4_r_plot.png", width = 10, height = 8)

reg4_inc_plot <- reg4_data_clean %>% 
   ggplot() +
   geom_col(aes(dates, case), 
            fill = "steelblue") +
   labs(title = paste0("COVID-19 Cases in Montana Region 4, 2020", 
                       " (Total cases = ", total_cases, ")"),
        subtitle = paste0("Data current as of ", date_today)) +
   ylab("Number of Cases") +
   xlab("") +
   scale_x_date(date_breaks = "2 day", date_labels = "%d-%b") +
   theme_minimal() +
   theme(strip.text = element_text(size = 16, colour = "black"),
         title = element_text(size = 12, colour = "black"),
         panel.grid = element_blank(),
         panel.grid.major.y = element_line(colour = "grey"),
         axis.text.x = element_text(size = 12, colour = "black", 
                                    angle = 90, vjust = 0.4),
         axis.text.y = element_text(size = 12, colour = "black"),
         axis.title.y = element_text(size = 12, colour = "black",
                                     margin = unit(c(0, 5, 0, 0), "mm")),
         axis.title.x = element_text(size = 12, colour = "black",
                                     margin = unit(c(5, 0, 0, 0), "mm")),
         axis.line.x = element_blank(), 
         axis.line.y = element_blank(), 
         axis.ticks = element_blank()) 
#reg4_inc_plot

ggsave("C:/R/covid19/state_daily_results/reg4_inc_plot.png", width = 10, height = 8)


## Reg5 results
reg5_data <- state_data_clean %>% 
   filter(region == 5)

reg5_incidence <- incidence(reg5_data$dates, last_date = current_date)

reg5_n <- as.data.frame(reg5_incidence$dates)
time_var <- nrow(reg5_n)
time_start <- seq(2, time_var-13)
time_end <- time_start + 13

reg5_results <- estimate_R(reg5_incidence, method="parametric_si", 
                           config = make_config(list(mean_si = serial_interval_mean, 
                                                     std_si = serial_interval_sd,
                                                     t_start =  time_start,
                                                     t_end = time_end)))

#reg5_plot <- estimate_R_plots(reg5_results, what = "R")

colomns <- c(1:4, 8)
reg5_r <- (reg5_results$R[,colomns]) %>% 
   mutate(region = "5") %>% 
   rename(mean_r = `Mean(R)`,
          sd_r = `Std(R)`,
          median_r = `Median(R)`)

reg5_dates <- as.data.frame(reg5_results$dates)
reg5_i <- as.data.frame(reg5_results$I)
reg5_cil <- as.data.frame(reg5_results$R$`Quantile.0.025(R)`)
reg5_cih <- as.data.frame(reg5_results$R$`Quantile.0.975(R)`)
reg5_dates_new <- cbind(reg5_dates, reg5_i) %>% 
   rename(dates = 1,
          incidence = 2) %>% 
   mutate(dates = ymd(dates))
reg5_dates_new <- reg5_dates_new[-(1:14), 1:2]
reg5_dates_new <- cbind(reg5_dates_new, reg5_cil, reg5_cih) %>% 
   rename(cl_low = 3,
          cl_high = 4)
reg5_r_clean <- cbind(reg5_r, reg5_dates_new)


reg5_data_clean <- state_data_clean %>% 
   filter(region == 5)
total_cases <- sum(reg5_data_clean$case)


reg5_r_plot <- reg5_r_clean %>% 
   ggplot() +
   geom_line(aes(dates, mean_r), size = 1.5, color = "black") +
   geom_line(aes(dates, cl_low), size = 1.5, color = "grey") +
   geom_line(aes(dates, cl_high), size = 1.5, color = "grey") +
   labs(title = "COVID-19 Rolling 14-day R-values, Montana Region 5, 2020",
        subtitle = paste0("Data current as of ", date_today),
        color = "") +
   ylab("R-value") +
   xlab("") +
   geom_hline(yintercept = 1, color = "black", size = 1.2) +
   scale_x_date(date_breaks = "2 days", date_labels = "%d-%b") +
   #scale_y_continuous(breaks = seq(0, 20, 1), labels = seq(0, 20, 1)) +
   theme_minimal() +
   theme(strip.text = element_text(size = 16, colour = "black"),
         title = element_text(size = 12, colour = "black"),
         panel.grid = element_blank(),
         panel.grid.major.y = element_line(colour = "grey"),
         axis.text.x = element_text(size = 12, colour = "black", 
                                    angle = 90, vjust = 0.4),
         axis.text.y = element_text(size = 12, colour = "black"),
         legend.text = element_text(size = 12, colour = "black"),
         axis.title.y = element_text(size = 12, colour = "black",
                                     margin = unit(c(0, 5, 0, 0), "mm")),
         axis.title.x = element_text(size = 12, colour = "black",
                                     margin = unit(c(5, 0, 0, 0), "mm")),
         axis.line.x = element_blank(), 
         axis.line.y = element_blank(), 
         axis.ticks = element_blank()) +
   scale_color_manual(values = c("black")) 
#reg5_r_plot

ggsave("C:/R/covid19/state_daily_results/reg5_r_plot.png", width = 10, height = 8)

reg5_inc_plot <- reg5_data_clean %>% 
   ggplot() +
   geom_col(aes(dates, case), 
            fill = "steelblue") +
   labs(title = paste0("COVID-19 Cases in Montana Region 5, 2020", 
                       " (Total cases = ", total_cases, ")"),
        subtitle = paste0("Data current as of ", date_today)) +
   ylab("Number of Cases") +
   xlab("") +
   scale_x_date(date_breaks = "2 day", date_labels = "%d-%b") +
   theme_minimal() +
   theme(strip.text = element_text(size = 16, colour = "black"),
         title = element_text(size = 12, colour = "black"),
         panel.grid = element_blank(),
         panel.grid.major.y = element_line(colour = "grey"),
         axis.text.x = element_text(size = 12, colour = "black", 
                                    angle = 90, vjust = 0.4),
         axis.text.y = element_text(size = 12, colour = "black"),
         axis.title.y = element_text(size = 12, colour = "black",
                                     margin = unit(c(0, 5, 0, 0), "mm")),
         axis.title.x = element_text(size = 12, colour = "black",
                                     margin = unit(c(5, 0, 0, 0), "mm")),
         axis.line.x = element_blank(), 
         axis.line.y = element_blank(), 
         axis.ticks = element_blank()) 
#reg5_inc_plot

ggsave("C:/R/covid19/state_daily_results/reg5_inc_plot.png", width = 10, height = 8)



# Bind files and save
all_regions_r <- rbind(state_r_clean, reg1_r_clean, reg2_r_clean, reg3_r_clean, 
                       reg4_r_clean, reg5_r_clean) %>% 
   left_join(all_data_wide, by = c("region", "dates"))

write_csv(all_regions_r, "C:/R/covid19/state_daily_results/all_regions_r.csv", na = " ")



# Send files to the sftp server
# host = 'elbastion.dbs.umt.edu'
# port = 22
# username = 'celftp'
# password  = 'celftp'
# remotepath = '/celFtpFiles/covid19/Rt/incoming/'

sftpUpload("elbastion.dbs.umt.edu", "celftp", "celftp",
           "/celFtpFiles/covid19/Rt/incoming/all_regions_r.csv",
           "C:/R/covid19/state_daily_results/all_regions_r.csv")

sftpUpload("elbastion.dbs.umt.edu", "celftp", "celftp",
           "/celFtpFiles/covid19/Rt/incoming/state_hosp.csv",
           "C:/R/covid19/state_daily_results/state_hosp.csv")

sftpUpload("elbastion.dbs.umt.edu", "celftp", "celftp",
           "/celFtpFiles/covid19/Rt/incoming/county_data_combined.csv",
           "C:/R/covid19/state_daily_results/county_data_combined.csv")

sftpUpload("elbastion.dbs.umt.edu", "celftp", "celftp",
           "/celFtpFiles/covid19/Rt/incoming/mt_case_outcome_hosp_data.csv",
           "C:/R/covid19/state_daily_results/mt_case_outcome_hosp_data.csv")

sftpUpload("elbastion.dbs.umt.edu", "celftp", "celftp",
           "/celFtpFiles/covid19/Rt/incoming/state_r_plot.png",
           "C:/R/covid19/state_daily_results/state_r_plot.png")

sftpUpload("elbastion.dbs.umt.edu", "celftp", "celftp",
           "/celFtpFiles/covid19/Rt/incoming/reg1_r_plot.png",
           "C:/R/covid19/state_daily_results/reg1_r_plot.png")

sftpUpload("elbastion.dbs.umt.edu", "celftp", "celftp",
           "/celFtpFiles/covid19/Rt/incoming/reg2_r_plot.png",
           "C:/R/covid19/state_daily_results/reg2_r_plot.png")

sftpUpload("elbastion.dbs.umt.edu", "celftp", "celftp",
           "/celFtpFiles/covid19/Rt/incoming/reg3_r_plot.png",
           "C:/R/covid19/state_daily_results/reg3_r_plot.png")

sftpUpload("elbastion.dbs.umt.edu", "celftp", "celftp",
           "/celFtpFiles/covid19/Rt/incoming/reg4_r_plot.png",
           "C:/R/covid19/state_daily_results/reg4_r_plot.png")

sftpUpload("elbastion.dbs.umt.edu", "celftp", "celftp",
           "/celFtpFiles/covid19/Rt/incoming/reg5_r_plot.png",
           "C:/R/covid19/state_daily_results/reg5_r_plot.png")

sftpUpload("elbastion.dbs.umt.edu", "celftp", "celftp",
           "/celFtpFiles/covid19/Rt/incoming/state_inc_plot.png",
           "C:/R/covid19/state_daily_results/state_inc_plot.png")

sftpUpload("elbastion.dbs.umt.edu", "celftp", "celftp",
           "/celFtpFiles/covid19/Rt/incoming/reg1_inc_plot.png",
           "C:/R/covid19/state_daily_results/reg1_inc_plot.png")

sftpUpload("elbastion.dbs.umt.edu", "celftp", "celftp",
           "/celFtpFiles/covid19/Rt/incoming/reg2_inc_plot.png",
           "C:/R/covid19/state_daily_results/reg2_inc_plot.png")

sftpUpload("elbastion.dbs.umt.edu", "celftp", "celftp",
           "/celFtpFiles/covid19/Rt/incoming/reg3_inc_plot.png",
           "C:/R/covid19/state_daily_results/reg3_inc_plot.png")

sftpUpload("elbastion.dbs.umt.edu", "celftp", "celftp",
           "/celFtpFiles/covid19/Rt/incoming/reg4_inc_plot.png",
           "C:/R/covid19/state_daily_results/reg4_inc_plot.png")

sftpUpload("elbastion.dbs.umt.edu", "celftp", "celftp",
           "/celFtpFiles/covid19/Rt/incoming/reg5_inc_plot.png",
           "C:/R/covid19/state_daily_results/reg5_inc_plot.png")



# Test to see if I can pull file back from server

#sftpDownload("elbastion.dbs.umt.edu", "celftp", "celftp",
#           "/celFtpFiles/covid19/Rt/incoming/reg5_plot.png",
#           "C:/R/covid19/reg5_plot.png")



