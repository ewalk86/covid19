# Estimate R, save results/plots
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
library(googledrive)
library(googlesheets4)



file_path <- c("C:/Users/ethan.walker/Box/Ethan Walker UM/R/covid19/")

##### Data prep

# Make dataframes for 5 health regions
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



# Load/format case data
mt_case_data <- read_xlsx(paste0(file_path, "Input/SI_Local_v_Import Data_11.12.2020.xlsx"),
                                 sheet = 2) %>% 
   rename_all(tolower) %>% 
   select(-case_no) %>% 
   rownames_to_column(var = "case_no") %>% 
   mutate(local = if_else(local_import == 0, 1, 0),
          imported = local_import,
          date_reported = ymd(date_reported),
          dates = as.Date(symptom_onset_date),
          dates = ymd(dates),
          case = 1) %>% 
   rename(hospitalization = "ever_hospitalized") %>% 
   left_join(counties_regions, by = "county") %>% 
   mutate(age_group_new = if_else(age_group == "80-89" | age_group == "90-99" | age_group == "100" | age_group == "100-110",
                               "80+", age_group),
          age_group_new = factor(age_group_new, 
                                 levels = c("0-9", "10-19", "20-29", 
                                            "30-39", "40-49", "50-59", 
                                            "60-69", "70-79", "80+"),
                                 labels = c("0 to 9", "10 to 19", "20 to 29", 
                                            "30 to 39", "40 to 49", "50 to 59", 
                                            "60 to 69", "70 to 79", "80+"))) %>% 
   mutate(hospitalization = factor(hospitalization,
                                   levels = c("Y", "N", "P", "U"),
                                   labels = c("Hosp: Yes", "Hosp: No", 
                                              "Hosp: Past", "Hosp: Unknown")),
          sex = fct_collapse(sex, "M" = c("m", "M"), "F" = c("f", "F"))) %>% 
   filter(!is.na(dates)) %>% 
   arrange(dates) %>% 
   ungroup() %>% 
   select(case_no, mt_case, county:hospitalization, local_import:age_group_new)



# Change case data to wide format
state_data_wide <- mt_case_data %>% 
   mutate(sex = factor(sex, labels = c("Female", "Male"))) %>% 
   select(-age_group) %>% 
   mutate(case = 1) %>% 
   pivot_wider(names_from = "age_group_new", values_from = "case") %>% 
   mutate(case = 1) %>% 
   select(-"NA") %>% 
   pivot_wider(names_from = "sex", values_from = "case") %>% 
   mutate(case = 1) %>% 
   select(-"NA") %>% 
   pivot_wider(names_from = "hospitalization", values_from = "case") %>% 
   mutate(case = 1) %>% 
   select(-"NA") %>% 
   pivot_wider(names_from = "county", values_from = "case") 


# Filter and format wide data for state and 5 health regions
state_wide_date <- state_data_wide %>% 
   select(region, dates, 8:78) %>% 
   mutate(region = "state") %>% 
   group_by(region, dates) %>% 
   mutate_all(sum, na.rm = TRUE) %>% 
   mutate(case = 1,
          daily_cases = sum(case)) %>% 
   arrange(dates) %>% 
   distinct(dates, .keep_all = TRUE) %>% 
   ungroup() %>% 
   mutate(cumulative_cases = cumsum(daily_cases)) %>% 
   select(region, dates, daily_cases, cumulative_cases, 3:75)

reg1_wide_date <- state_data_wide %>% 
   filter(region == 1) %>% 
   select(region, dates, 8:78) %>% 
   group_by(region, dates) %>% 
   mutate_all(sum, na.rm = TRUE) %>% 
   mutate(case = 1,
          daily_cases = sum(case)) %>% 
   arrange(dates) %>% 
   distinct(dates, .keep_all = TRUE) %>% 
   ungroup() %>% 
   mutate(cumulative_cases = cumsum(daily_cases)) %>% 
   select(region, dates, daily_cases, cumulative_cases, 3:75)

reg2_wide_date <- state_data_wide %>% 
   filter(region == 2) %>% 
   select(region, dates, 8:78) %>% 
   group_by(region, dates) %>% 
   mutate_all(sum, na.rm = TRUE) %>% 
   mutate(case = 1,
          daily_cases = sum(case)) %>% 
   arrange(dates) %>% 
   distinct(dates, .keep_all = TRUE) %>% 
   ungroup() %>% 
   mutate(cumulative_cases = cumsum(daily_cases)) %>% 
   select(region, dates, daily_cases, cumulative_cases, 3:75)

reg3_wide_date <- state_data_wide %>% 
   filter(region == 3) %>% 
   select(region, dates, 8:78) %>% 
   group_by(region, dates) %>% 
   mutate_all(sum, na.rm = TRUE) %>% 
   mutate(case = 1,
          daily_cases = sum(case)) %>% 
   arrange(dates) %>% 
   distinct(dates, .keep_all = TRUE) %>% 
   ungroup() %>% 
   mutate(cumulative_cases = cumsum(daily_cases)) %>% 
   select(region, dates, daily_cases, cumulative_cases, 3:75)

reg4_wide_date <- state_data_wide %>% 
   filter(region == 4) %>% 
   select(region, dates, 8:78) %>% 
   group_by(region, dates) %>% 
   mutate_all(sum, na.rm = TRUE) %>% 
   mutate(case = 1,
          daily_cases = sum(case)) %>% 
   arrange(dates) %>% 
   distinct(dates, .keep_all = TRUE) %>% 
   ungroup() %>% 
   mutate(cumulative_cases = cumsum(daily_cases)) %>% 
   select(region, dates, daily_cases, cumulative_cases, 3:75)

reg5_wide_date <- state_data_wide %>% 
   filter(region == 5) %>% 
   select(region, dates, 8:78) %>% 
   group_by(region, dates) %>% 
   mutate_all(sum, na.rm = TRUE) %>% 
   mutate(case = 1,
          daily_cases = sum(case)) %>% 
   arrange(dates) %>% 
   distinct(dates, .keep_all = TRUE) %>% 
   ungroup() %>% 
   mutate(cumulative_cases = cumsum(daily_cases)) %>% 
   select(region, dates, daily_cases, cumulative_cases, 3:75)

all_data_wide <- rbind(state_wide_date, reg1_wide_date, reg2_wide_date,
                       reg3_wide_date, reg4_wide_date, reg5_wide_date)




#################### Run analysis and print results ########################


##### State results #####

# Format analysis data
mt_li_analysis_data <- mt_case_data %>% 
   group_by(dates) %>% 
   mutate(local = sum(local),
          imported = sum(imported)) %>% 
   select(dates, local, imported) %>% 
   distinct(dates, .keep_all = TRUE) %>% 
   arrange(dates)


mt_data <- mt_li_analysis_data %>% 
   mutate(I = local + imported) %>% 
   select(dates, I) 

latest_date <- format(Sys.Date() - 14, "%Y-%m-%d")

mt_incidence_data <- incidence(mt_data$dates, last_date = latest_date)

mt_li_inc_data <- as.data.frame(mt_incidence_data$dates) %>% 
   mutate(dates = ymd(mt_incidence_data$dates)) %>% 
   select(dates) %>% 
   left_join(mt_li_analysis_data, by = "dates") %>% 
   mutate(local = if_else(is.na(local), 0, local),
          imported = if_else(is.na(imported), 0, imported))


# Set dates for 14-day rolling averages
state_n <- as.data.frame(mt_li_inc_data$dates)
time_var <- nrow(state_n)
time_start <- seq(2, time_var-13)
time_end <- time_start + 13


# Serial Interval derived from State of Montana paired case data
serial_interval_mean <- 5.29
serial_interval_sd <- 4.45

mt_r_results <- estimate_R(mt_li_inc_data, method="parametric_si", 
                           config = make_config(list(mean_si = serial_interval_mean, 
                                                     std_si = serial_interval_sd,
                                                     t_start =  time_start,
                                                     t_end = time_end)))

# Format and save analysis results
colomns <- c(1:4, 8)
state_r <- (mt_r_results$R[,colomns]) %>% 
   mutate(region = "state") %>% 
   rename(mean_r = `Mean(R)`,
          sd_r = `Std(R)`,
          median_r = `Median(R)`)

state_dates <- as.data.frame(mt_r_results$dates)
state_i <- as.data.frame(mt_r_results$I)
state_cil <- as.data.frame(mt_r_results$R$`Quantile.0.025(R)`)
state_cih <- as.data.frame(mt_r_results$R$`Quantile.0.975(R)`)
state_dates_new <- cbind(state_dates, state_i) %>% 
   rename(dates = 1,
          incidence = 2) %>% 
   mutate(dates = ymd(dates))
state_dates_new <- state_dates_new[-(1:14), 1:2]
state_dates_new <- cbind(state_dates_new, state_cil, state_cih) %>% 
   rename(cl_low = 3,
          cl_high = 4)
state_r <- cbind(state_r, state_dates_new)



##### Function for running R results for 5 health regions #####

region_analysis_function <- function(filter_var, label_var, data = mt_case_data) {
   
   mt_analysis_data <- mt_case_data %>% 
      filter(region == filter_var) %>% 
      group_by(dates) %>% 
      mutate(local = sum(local),
             imported = sum(imported)) %>% 
      select(dates, local, imported) %>% 
      distinct(dates, .keep_all = TRUE) %>% 
      arrange(dates)
   
   
   mt_data <- mt_analysis_data %>% 
      mutate(I = local + imported) %>% 
      select(dates, I) 
   
   mt_incidence_data <- incidence(mt_data$dates, last_date = latest_date)
   
   mt_li_inc_data <- as.data.frame(mt_incidence_data$dates) %>% 
      mutate(dates = ymd(mt_incidence_data$dates)) %>% 
      select(dates) %>% 
      left_join(mt_analysis_data, by = "dates") %>% 
      mutate(local = if_else(is.na(local), 0, local),
             imported = if_else(is.na(imported), 0, imported)) 
   
   
   state_n <- as.data.frame(mt_li_inc_data$dates)
   time_var <- nrow(state_n)
   time_start <- seq(2, time_var-13)
   time_end <- time_start + 13
   
   
   # Serial Interval derived from State of Montana case data
   serial_interval_mean <- 5.29
   serial_interval_sd <- 4.45
   
   mt_r_results <- estimate_R(mt_li_inc_data, method="parametric_si", 
                              config = make_config(list(mean_si = serial_interval_mean, 
                                                        std_si = serial_interval_sd,
                                                        t_start =  time_start,
                                                        t_end = time_end)))
   
   
   colomns <- c(1:4, 8)
   state_r <- (mt_r_results$R[,colomns]) %>% 
      mutate(region = label_var) %>% 
      rename(mean_r = `Mean(R)`,
             sd_r = `Std(R)`,
             median_r = `Median(R)`)
   
   state_dates <- as.data.frame(mt_r_results$dates)
   state_i <- as.data.frame(mt_r_results$I)
   state_cil <- as.data.frame(mt_r_results$R$`Quantile.0.025(R)`)
   state_cih <- as.data.frame(mt_r_results$R$`Quantile.0.975(R)`)
   state_dates_new <- cbind(state_dates, state_i) %>% 
      rename(dates = 1,
             incidence = 2) %>% 
      mutate(dates = ymd(dates))
   state_dates_new <- state_dates_new[-(1:14), 1:2]
   state_dates_new <- cbind(state_dates_new, state_cil, state_cih) %>% 
      rename(cl_low = 3,
             cl_high = 4)
   r_clean <<- cbind(state_r, state_dates_new)
   
}

reg1_r <- region_analysis_function("1", "1")
reg2_r <- region_analysis_function("2", "2")
reg3_r <- region_analysis_function("3", "3")
reg4_r <- region_analysis_function("4", "4")
reg5_r <- region_analysis_function("5", "5")



##### Function for running R results for counties #####

county_analysis_function <- function(filter_var, label_var, data = mt_case_data) {
   
   mt_analysis_data <- mt_case_data %>% 
      filter(county == filter_var) %>% 
      group_by(dates) %>% 
      mutate(local = sum(local),
             imported = sum(imported)) %>% 
      select(dates, local, imported) %>% 
      distinct(dates, .keep_all = TRUE) %>% 
      arrange(dates)
   
   
   mt_data <- mt_analysis_data %>% 
      mutate(I = local + imported) %>% 
      select(dates, I) 
   
   mt_incidence_data <- incidence(mt_data$dates, last_date = latest_date)
   
   mt_li_inc_data <- as.data.frame(mt_incidence_data$dates) %>% 
      mutate(dates = ymd(mt_incidence_data$dates)) %>% 
      select(dates) %>% 
      left_join(mt_analysis_data, by = "dates") %>% 
      mutate(local = if_else(is.na(local), 0, local),
             imported = if_else(is.na(imported), 0, imported)) 
   
   
   state_n <- as.data.frame(mt_li_inc_data$dates)
   time_var <- nrow(state_n)
   time_start <- seq(2, time_var-13)
   time_end <- time_start + 13
   
   
   # Serial Interval derived from State of Montana case data
   serial_interval_mean <- 5.29
   serial_interval_sd <- 4.45
   
   mt_r_results <- estimate_R(mt_li_inc_data, method="parametric_si", 
                              config = make_config(list(mean_si = serial_interval_mean, 
                                                        std_si = serial_interval_sd,
                                                        t_start =  time_start,
                                                        t_end = time_end)))
   
   
   colomns <- c(1:4, 8)
   state_r <- (mt_r_results$R[,colomns]) %>% 
      mutate(region = label_var) %>% 
      rename(mean_r = `Mean(R)`,
             sd_r = `Std(R)`,
             median_r = `Median(R)`)
   
   state_dates <- as.data.frame(mt_r_results$dates)
   state_i <- as.data.frame(mt_r_results$I)
   state_cil <- as.data.frame(mt_r_results$R$`Quantile.0.025(R)`)
   state_cih <- as.data.frame(mt_r_results$R$`Quantile.0.975(R)`)
   state_dates_new <- cbind(state_dates, state_i) %>% 
      rename(dates = 1,
             incidence = 2) %>% 
      mutate(dates = ymd(dates))
   state_dates_new <- state_dates_new[-(1:14), 1:2]
   state_dates_new <- cbind(state_dates_new, state_cil, state_cih) %>% 
      rename(cl_low = 3,
             cl_high = 4)
   r_clean <<- cbind(state_r, state_dates_new)
   
}

missoula_r <- county_analysis_function("Missoula", "Missoula County")
gallatin_r <- county_analysis_function("Gallatin", "Gallatin County")
yellowstone_r <- county_analysis_function("Yellowstone", "Yellowstone County")
bighorn_r <- county_analysis_function("Big Horn", "Big Horn County")
lake_r <- county_analysis_function("Lake", "Lake County")
lewisandclark_r <- county_analysis_function("Lewis and Clark", "Lewis and Clark County")
flathead_r <- county_analysis_function("Flathead", "Flathead County")
cascade_r <- county_analysis_function("Cascade", "Cascade County")
silverbow_r <- county_analysis_function("Silver Bow", "Silver Bow County")
rosebud_r <- county_analysis_function("Rosebud", "Rosebud County")
glacier_r <- county_analysis_function("Glacier", "Glacier County")
roosevelt_r <- county_analysis_function("Roosevelt", "Roosevelt County")




# Bind files and save
all_regions_r <- rbind(state_r, reg1_r, reg2_r, reg3_r, reg4_r, reg5_r,
                       missoula_r, gallatin_r, yellowstone_r, bighorn_r,
                       lake_r, lewisandclark_r, flathead_r, cascade_r,
                       silverbow_r, rosebud_r, glacier_r, roosevelt_r) %>% 
   left_join(all_data_wide, by = c("region", "dates")) %>% 
   mutate(daily_cases = incidence) %>% 
   mutate(mean_r = round(mean_r, digits = 2),
          median_r = round(median_r, digits = 2),
          sd_r = round(sd_r, digits = 2),
          cl_low = round(cl_low, digits = 2),
          cl_high = round(cl_high, digits = 2))

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



# Test to see if I can pull file back from server

#sftpDownload("elbastion.dbs.umt.edu", "celftp", "celftp",
#           "/celFtpFiles/covid19/Rt/incoming/reg5_plot.png",
#           "C:/R/covid19/reg5_plot.png")



############ Push files to Google ##################3

options(gargle_oauth_email = "ethanwalker86@gmail.com")
drive_auth(email = "ethanwalker86@gmail.com")
gs4_auth(token = drive_token())


r_filter_dates <- all_regions_r %>% 
   filter(dates > Sys.Date() - 60)


# Write State/region R data to google
mt_data <- r_filter_dates %>% 
   filter(region == "state") %>% 
   select(dates, cl_low, mean_r, cl_high) %>% 
   rename("Dates" = dates,
          "Lower_Confidence_Limit" = cl_low,
          "Mean_R" = mean_r,
          "Upper_Confidence_Limit" = cl_high)

sheet_write(mt_data, 
            ss = "https://docs.google.com/spreadsheets/d/1q-QUjewIBdPu_LNR8E53EPQMtw_nEaiOojYkJ3_cSYc/edit#gid=441849247",
            sheet = 1)

reg1_data <- r_filter_dates %>% 
   filter(region == "1") %>% 
   select(dates, cl_low, mean_r, cl_high) %>% 
   rename("Dates" = dates,
          "Lower_Confidence_Limit" = cl_low,
          "Mean_R" = mean_r,
          "Upper_Confidence_Limit" = cl_high)

sheet_write(reg1_data, 
            ss = "https://docs.google.com/spreadsheets/d/1q-QUjewIBdPu_LNR8E53EPQMtw_nEaiOojYkJ3_cSYc/edit#gid=441849247",
            sheet = 2)

reg2_data <- r_filter_dates %>% 
   filter(region == "2") %>% 
   select(dates, cl_low, mean_r, cl_high) %>% 
   rename("Dates" = dates,
          "Lower_Confidence_Limit" = cl_low,
          "Mean_R" = mean_r,
          "Upper_Confidence_Limit" = cl_high)

sheet_write(reg2_data, 
            ss = "https://docs.google.com/spreadsheets/d/1q-QUjewIBdPu_LNR8E53EPQMtw_nEaiOojYkJ3_cSYc/edit#gid=441849247",
            sheet = 3)

reg3_data <- r_filter_dates %>% 
   filter(region == "3") %>% 
   select(dates, cl_low, mean_r, cl_high) %>% 
   rename("Dates" = dates,
          "Lower_Confidence_Limit" = cl_low,
          "Mean_R" = mean_r,
          "Upper_Confidence_Limit" = cl_high)

sheet_write(reg3_data, 
            ss = "https://docs.google.com/spreadsheets/d/1q-QUjewIBdPu_LNR8E53EPQMtw_nEaiOojYkJ3_cSYc/edit#gid=441849247",
            sheet = 4)

reg4_data <- r_filter_dates %>% 
   filter(region == "4") %>% 
   select(dates, cl_low, mean_r, cl_high) %>% 
   rename("Dates" = dates,
          "Lower_Confidence_Limit" = cl_low,
          "Mean_R" = mean_r,
          "Upper_Confidence_Limit" = cl_high)

sheet_write(reg4_data, 
            ss = "https://docs.google.com/spreadsheets/d/1q-QUjewIBdPu_LNR8E53EPQMtw_nEaiOojYkJ3_cSYc/edit#gid=441849247",
            sheet = 5)

reg5_data <- r_filter_dates %>% 
   filter(region == "5") %>% 
   select(dates, cl_low, mean_r, cl_high) %>% 
   rename("Dates" = dates,
          "Lower_Confidence_Limit" = cl_low,
          "Mean_R" = mean_r,
          "Upper_Confidence_Limit" = cl_high)

sheet_write(reg5_data, 
            ss = "https://docs.google.com/spreadsheets/d/1q-QUjewIBdPu_LNR8E53EPQMtw_nEaiOojYkJ3_cSYc/edit#gid=441849247",
            sheet = 6)



# Write State/region case data to google
mt_data <- all_regions_r %>% 
   filter(region == "state") %>% 
   select(dates, incidence) %>% 
   rename("Dates" = dates,
          "Cases" = incidence)

sheet_write(mt_data, 
            ss = "https://docs.google.com/spreadsheets/d/11adsCi8okCW26lyfkZSXItlBCMGuvnysSxM6lVB8pQY/edit#gid=0",
            sheet = 1)

reg1_data <- all_regions_r %>% 
   filter(region == "1") %>% 
   select(dates, incidence) %>% 
   rename("Dates" = dates,
          "Cases" = incidence)

sheet_write(reg1_data, 
            ss = "https://docs.google.com/spreadsheets/d/11adsCi8okCW26lyfkZSXItlBCMGuvnysSxM6lVB8pQY/edit#gid=0",
            sheet = 2)

reg2_data <- all_regions_r %>% 
   filter(region == "2") %>% 
   select(dates, incidence) %>% 
   rename("Dates" = dates,
          "Cases" = incidence)

sheet_write(reg2_data, 
            ss = "https://docs.google.com/spreadsheets/d/11adsCi8okCW26lyfkZSXItlBCMGuvnysSxM6lVB8pQY/edit#gid=0",
            sheet = 3)

reg3_data <- all_regions_r %>% 
   filter(region == "3") %>% 
   select(dates, incidence) %>% 
   rename("Dates" = dates,
          "Cases" = incidence)

sheet_write(reg3_data, 
            ss = "https://docs.google.com/spreadsheets/d/11adsCi8okCW26lyfkZSXItlBCMGuvnysSxM6lVB8pQY/edit#gid=0",
            sheet = 4)

reg4_data <- all_regions_r %>% 
   filter(region == "4") %>% 
   select(dates, incidence) %>% 
   rename("Dates" = dates,
          "Cases" = incidence)

sheet_write(reg4_data, 
            ss = "https://docs.google.com/spreadsheets/d/11adsCi8okCW26lyfkZSXItlBCMGuvnysSxM6lVB8pQY/edit#gid=0",
            sheet = 5)

reg5_data <- all_regions_r %>% 
   filter(region == "5") %>% 
   select(dates, incidence) %>% 
   rename("Dates" = dates,
          "Cases" = incidence)

sheet_write(reg5_data, 
            ss = "https://docs.google.com/spreadsheets/d/11adsCi8okCW26lyfkZSXItlBCMGuvnysSxM6lVB8pQY/edit#gid=0",
            sheet = 6)



# Write county R data to google
miss_data <- r_filter_dates %>% 
   filter(region == "Missoula County") %>% 
   select(dates, cl_low, mean_r, cl_high) %>% 
   rename("Dates" = dates,
          "Lower_Confidence_Limit" = cl_low,
          "Mean_R" = mean_r,
          "Upper_Confidence_Limit" = cl_high)

sheet_write(miss_data, 
            ss = "https://docs.google.com/spreadsheets/d/1L1pBAp0e5RvU7x4IfnfyeJ5uaQs0OTuYXYmgJB0r_oE/edit#gid=0",
            sheet = 1)

gall_data <- r_filter_dates %>% 
   filter(region == "Gallatin County") %>% 
   select(dates, cl_low, mean_r, cl_high) %>% 
   rename("Dates" = dates,
          "Lower_Confidence_Limit" = cl_low,
          "Mean_R" = mean_r,
          "Upper_Confidence_Limit" = cl_high)

sheet_write(gall_data, 
            ss = "https://docs.google.com/spreadsheets/d/1L1pBAp0e5RvU7x4IfnfyeJ5uaQs0OTuYXYmgJB0r_oE/edit#gid=0",
            sheet = 2)

ystn_data <- r_filter_dates %>% 
   filter(region == "Yellowstone County") %>% 
   select(dates, cl_low, mean_r, cl_high) %>% 
   rename("Dates" = dates,
          "Lower_Confidence_Limit" = cl_low,
          "Mean_R" = mean_r,
          "Upper_Confidence_Limit" = cl_high)

sheet_write(ystn_data, 
            ss = "https://docs.google.com/spreadsheets/d/1L1pBAp0e5RvU7x4IfnfyeJ5uaQs0OTuYXYmgJB0r_oE/edit#gid=0",
            sheet = 3)

bh_data <- r_filter_dates %>% 
   filter(region == "Big Horn County") %>% 
   select(dates, cl_low, mean_r, cl_high) %>% 
   rename("Dates" = dates,
          "Lower_Confidence_Limit" = cl_low,
          "Mean_R" = mean_r,
          "Upper_Confidence_Limit" = cl_high)

sheet_write(bh_data, 
            ss = "https://docs.google.com/spreadsheets/d/1L1pBAp0e5RvU7x4IfnfyeJ5uaQs0OTuYXYmgJB0r_oE/edit#gid=0",
            sheet = 4)

lake_data <- r_filter_dates %>% 
   filter(region == "Lake County") %>% 
   select(dates, cl_low, mean_r, cl_high) %>% 
   rename("Dates" = dates,
          "Lower_Confidence_Limit" = cl_low,
          "Mean_R" = mean_r,
          "Upper_Confidence_Limit" = cl_high)

sheet_write(lake_data, 
            ss = "https://docs.google.com/spreadsheets/d/1L1pBAp0e5RvU7x4IfnfyeJ5uaQs0OTuYXYmgJB0r_oE/edit#gid=0",
            sheet = 5)

lac_data <- r_filter_dates %>% 
   filter(region == "Lewis and Clark County") %>% 
   select(dates, cl_low, mean_r, cl_high) %>% 
   rename("Dates" = dates,
          "Lower_Confidence_Limit" = cl_low,
          "Mean_R" = mean_r,
          "Upper_Confidence_Limit" = cl_high)

sheet_write(lac_data, 
            ss = "https://docs.google.com/spreadsheets/d/1L1pBAp0e5RvU7x4IfnfyeJ5uaQs0OTuYXYmgJB0r_oE/edit#gid=0",
            sheet = 6)

flat_data <- r_filter_dates %>% 
   filter(region == "Flathead County") %>% 
   select(dates, cl_low, mean_r, cl_high) %>% 
   rename("Dates" = dates,
          "Lower_Confidence_Limit" = cl_low,
          "Mean_R" = mean_r,
          "Upper_Confidence_Limit" = cl_high)

sheet_write(flat_data, 
            ss = "https://docs.google.com/spreadsheets/d/1L1pBAp0e5RvU7x4IfnfyeJ5uaQs0OTuYXYmgJB0r_oE/edit#gid=0",
            sheet = 7)

casc_data <- r_filter_dates %>% 
   filter(region == "Cascade County") %>% 
   select(dates, cl_low, mean_r, cl_high) %>% 
   rename("Dates" = dates,
          "Lower_Confidence_Limit" = cl_low,
          "Mean_R" = mean_r,
          "Upper_Confidence_Limit" = cl_high)

sheet_write(casc_data, 
            ss = "https://docs.google.com/spreadsheets/d/1L1pBAp0e5RvU7x4IfnfyeJ5uaQs0OTuYXYmgJB0r_oE/edit#gid=0",
            sheet = 8)

sb_data <- r_filter_dates %>% 
   filter(region == "Silver Bow County") %>% 
   select(dates, cl_low, mean_r, cl_high) %>% 
   rename("Dates" = dates,
          "Lower_Confidence_Limit" = cl_low,
          "Mean_R" = mean_r,
          "Upper_Confidence_Limit" = cl_high)

sheet_write(sb_data, 
            ss = "https://docs.google.com/spreadsheets/d/1L1pBAp0e5RvU7x4IfnfyeJ5uaQs0OTuYXYmgJB0r_oE/edit#gid=0",
            sheet = 9)

rosebud_data <- r_filter_dates %>% 
   filter(region == "Rosebud County") %>% 
   select(dates, cl_low, mean_r, cl_high) %>% 
   rename("Dates" = dates,
          "Lower_Confidence_Limit" = cl_low,
          "Mean_R" = mean_r,
          "Upper_Confidence_Limit" = cl_high)

sheet_write(rosebud_data, 
            ss = "https://docs.google.com/spreadsheets/d/1L1pBAp0e5RvU7x4IfnfyeJ5uaQs0OTuYXYmgJB0r_oE/edit#gid=0",
            sheet = 10)

glacier_data <- r_filter_dates %>% 
   filter(region == "Glacier County") %>% 
   select(dates, cl_low, mean_r, cl_high) %>% 
   rename("Dates" = dates,
          "Lower_Confidence_Limit" = cl_low,
          "Mean_R" = mean_r,
          "Upper_Confidence_Limit" = cl_high)

sheet_write(glacier_data, 
            ss = "https://docs.google.com/spreadsheets/d/1L1pBAp0e5RvU7x4IfnfyeJ5uaQs0OTuYXYmgJB0r_oE/edit#gid=0",
            sheet = 11)

roosevelt_data <- r_filter_dates %>% 
   filter(region == "Roosevelt County") %>% 
   select(dates, cl_low, mean_r, cl_high) %>% 
   rename("Dates" = dates,
          "Lower_Confidence_Limit" = cl_low,
          "Mean_R" = mean_r,
          "Upper_Confidence_Limit" = cl_high)

sheet_write(roosevelt_data, 
            ss = "https://docs.google.com/spreadsheets/d/1L1pBAp0e5RvU7x4IfnfyeJ5uaQs0OTuYXYmgJB0r_oE/edit#gid=0",
            sheet = 12)



# Write county case data to google
miss_data <- all_regions_r %>% 
   filter(region == "Missoula County") %>% 
   select(dates, incidence) %>% 
   rename("Dates" = dates,
          "Cases" = incidence)

sheet_write(miss_data, 
            ss = "https://docs.google.com/spreadsheets/d/1X0vxDLVxQT_XrPgNMtyRwT1kTVQSN2PzQPYDWi2oxBo/edit#gid=0",
            sheet = 1)

gall_data <- all_regions_r %>% 
   filter(region == "Gallatin County") %>% 
   select(dates, incidence) %>% 
   rename("Dates" = dates,
          "Cases" = incidence)

sheet_write(gall_data, 
            ss = "https://docs.google.com/spreadsheets/d/1X0vxDLVxQT_XrPgNMtyRwT1kTVQSN2PzQPYDWi2oxBo/edit#gid=0",
            sheet = 2)

ystn_data <- all_regions_r %>% 
   filter(region == "Yellowstone County") %>% 
   select(dates, incidence) %>% 
   rename("Dates" = dates,
          "Cases" = incidence)

sheet_write(ystn_data, 
            ss = "https://docs.google.com/spreadsheets/d/1X0vxDLVxQT_XrPgNMtyRwT1kTVQSN2PzQPYDWi2oxBo/edit#gid=0",
            sheet = 3)

bh_data <- all_regions_r %>% 
   filter(region == "Big Horn County") %>% 
   select(dates, incidence) %>% 
   rename("Dates" = dates,
          "Cases" = incidence)

sheet_write(bh_data, 
            ss = "https://docs.google.com/spreadsheets/d/1X0vxDLVxQT_XrPgNMtyRwT1kTVQSN2PzQPYDWi2oxBo/edit#gid=0",
            sheet = 4)

lake_data <- all_regions_r %>% 
   filter(region == "Lake County") %>% 
   select(dates, incidence) %>% 
   rename("Dates" = dates,
          "Cases" = incidence)

sheet_write(lake_data, 
            ss = "https://docs.google.com/spreadsheets/d/1X0vxDLVxQT_XrPgNMtyRwT1kTVQSN2PzQPYDWi2oxBo/edit#gid=0",
            sheet = 5)

lac_data <- all_regions_r %>% 
   filter(region == "Lewis and Clark County") %>% 
   select(dates, incidence) %>% 
   rename("Dates" = dates,
          "Cases" = incidence)

sheet_write(lac_data, 
            ss = "https://docs.google.com/spreadsheets/d/1X0vxDLVxQT_XrPgNMtyRwT1kTVQSN2PzQPYDWi2oxBo/edit#gid=0",
            sheet = 6)

flat_data <- all_regions_r %>% 
   filter(region == "Flathead County") %>% 
   select(dates, incidence) %>% 
   rename("Dates" = dates,
          "Cases" = incidence)

sheet_write(flat_data, 
            ss = "https://docs.google.com/spreadsheets/d/1X0vxDLVxQT_XrPgNMtyRwT1kTVQSN2PzQPYDWi2oxBo/edit#gid=0",
            sheet = 7)

casc_data <- all_regions_r %>% 
   filter(region == "Cascade County") %>% 
   select(dates, incidence) %>% 
   rename("Dates" = dates,
          "Cases" = incidence)

sheet_write(casc_data, 
            ss = "https://docs.google.com/spreadsheets/d/1X0vxDLVxQT_XrPgNMtyRwT1kTVQSN2PzQPYDWi2oxBo/edit#gid=0",
            sheet = 8)

sb_data <- all_regions_r %>% 
   filter(region == "Silver Bow County") %>% 
   select(dates, incidence) %>% 
   rename("Dates" = dates,
          "Cases" = incidence)

sheet_write(sb_data, 
            ss = "https://docs.google.com/spreadsheets/d/1X0vxDLVxQT_XrPgNMtyRwT1kTVQSN2PzQPYDWi2oxBo/edit#gid=0",
            sheet = 9)

rosebud_data <- all_regions_r %>% 
   filter(region == "Rosebud County") %>% 
   select(dates, incidence) %>% 
   rename("Dates" = dates,
          "Cases" = incidence)

sheet_write(rosebud_data, 
            ss = "https://docs.google.com/spreadsheets/d/1X0vxDLVxQT_XrPgNMtyRwT1kTVQSN2PzQPYDWi2oxBo/edit#gid=0",
            sheet = 10)

glacier_data <- all_regions_r %>% 
   filter(region == "Glacier County") %>% 
   select(dates, incidence) %>% 
   rename("Dates" = dates,
          "Cases" = incidence)

sheet_write(glacier_data, 
            ss = "https://docs.google.com/spreadsheets/d/1X0vxDLVxQT_XrPgNMtyRwT1kTVQSN2PzQPYDWi2oxBo/edit#gid=0",
            sheet = 11)

roosevelt_data <- all_regions_r %>% 
   filter(region == "Roosevelt County") %>% 
   select(dates, incidence) %>% 
   rename("Dates" = dates,
          "Cases" = incidence)

sheet_write(roosevelt_data, 
            ss = "https://docs.google.com/spreadsheets/d/1X0vxDLVxQT_XrPgNMtyRwT1kTVQSN2PzQPYDWi2oxBo/edit#gid=0",
            sheet = 12)
