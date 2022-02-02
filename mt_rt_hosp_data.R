# Estimate R, save results/plots
# Ethan Walker
# 7 Dec 2020

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


mt_county_fips <- read_csv(paste0(file_path, "Input/mt_county_fips.csv")) %>% 
   rename(county_fips = fips) %>% 
   select(-county_id) %>% 
   mutate(county_fips = county_fips + 30000) %>% 
   full_join(counties_regions, by = "county")


# Load/format case data -- change file name in next line
mt_case_data <- read_xlsx(paste0(file_path, "Input/uom_covid_01312022.xlsx"),
                          sheet = 1, skip = 1,
                          col_names = c("midis_add_datetime", 
                                        "inv_start_date", 
                                        "state_num", "county_fips", 
                                        "sex", "spatial_age", "age", "age_unit", "hospitalization", 
                                        "hosp_dur", "hosp_dur_unit", "deceased",
                                        "onset_date", "spec_coll_date", "diagnosis_date",
                                        "reinfection",
                                        "breakthrough", "variant_class",
                                        "variant_type", "zipcode", "inv_rpt_date"),
                          col_types = c("date", 
                                        "date", "numeric", "numeric",
                                        "text", "text", "text", "text", "text",  
                                        "numeric", "text", "text", 
                                        "date", "date", "date",
                                        "numeric",
                                        "text", "text",
                                        "text", "text", "date")) %>% 
   rownames_to_column(var = "case_no") %>% 
   mutate(case = 1) %>% 
   left_join(mt_county_fips, by = "county_fips") %>% 
   filter(!is.na(county)) %>% 
   mutate(age_group = fct_collapse(age,
                                   "0 to 9" = "0-9",
                                   "10 to 19" = "10-19", 
                                   "20 to 29" = "20-29", 
                                   "30 to 39" = "30-39", 
                                   "40 to 49" = "40-49", 
                                   "50 to 59" = "50-59", 
                                   "60 to 69" = "60-69", 
                                   "70 to 79" = "70-79", 
                                   "80+" = c("80-89", "90-99", "100-109"))) %>% 
   mutate(sex = factor(sex,
                       levels = c("F", "M", "U"),
                       labels = c("Female", "Male", "Unknown")),
          deceased = if_else(is.na(deceased), "N", deceased),
          deceased = fct_collapse(deceased,
                                  "Deceased: Yes" = "Y",
                                  "Deceased: No" = "N",
                                  "Deceased: Unknown" = "UNK"),
          hospitalization = if_else(is.na(hospitalization), "N", hospitalization),
          hospitalization = fct_collapse(hospitalization,
                                  "Hosp: Yes" = "Y",
                                  "Hosp: No" = "N",
                                  "Hosp: Unknown" = "UNK")) %>% 
   filter(diagnosis_date > "2020-03-10" | is.na(diagnosis_date)) %>% 
   filter(spec_coll_date > "2020-03-10" | is.na(spec_coll_date)) %>% 
   filter(onset_date > "2020-02-13" | is.na(onset_date)) %>% 
   mutate(onset_inv_diff = as.duration(interval(onset_date, inv_start_date)),
          onset_inv_diff = as.numeric(onset_inv_diff)/86400,
          mean_onset_inv_diff = round(mean(onset_inv_diff, na.rm = TRUE)) * 86400,
          onset_coll_diff = as.duration(interval(onset_date, spec_coll_date)),
          onset_coll_diff = as.numeric(onset_coll_diff)/86400,
          mean_onset_coll_diff = round(mean(onset_coll_diff, na.rm = TRUE)) * 86400,
          onset_diag_diff = as.duration(interval(onset_date, diagnosis_date)),
          onset_diag_diff = as.numeric(onset_diag_diff)/86400,
          mean_onset_diag_diff = round(mean(onset_diag_diff, na.rm = TRUE)) * 86400,
          onset_date_2 = if_else(is.na(onset_date) & !is.na(spec_coll_date), 
                                 spec_coll_date - mean_onset_coll_diff, onset_date),
          onset_date_3 = if_else(is.na(onset_date_2) & !is.na(diagnosis_date), 
                                 diagnosis_date - mean_onset_diag_diff, onset_date_2),
          onset_date_4 = if_else(is.na(onset_date_3) & !is.na(inv_start_date), 
                                 inv_start_date - mean_onset_inv_diff, onset_date_3)) %>% 
   separate(onset_date_4, into = c("dates", "time"), sep = " ") %>% 
   mutate(dates = ymd(dates),
          region = as.factor(region),
          zipcode = as.factor(zipcode)) %>% 
   select(-time) %>% 
   filter(!is.na(dates)) %>% 
   arrange(dates) %>% 
   ungroup() 

# summary(mt_case_data)

#################### Run analysis and print results ########################


##### State results #####

# Format analysis data
mt_analysis_data <- mt_case_data %>% 
   group_by(dates) %>% 
   mutate(I = sum(case)) %>% 
   select(dates, I) %>% 
   arrange(dates)

latest_date <- format(Sys.Date() - 14, "%Y-%m-%d")

mt_incidence_data <- incidence(mt_analysis_data$dates, last_date = latest_date)

# Set dates for 14-day rolling averages
state_n <- as.data.frame(mt_incidence_data$dates)
time_var <- nrow(state_n)
time_start <- seq(2, time_var-13)
time_end <- time_start + 13


# Serial Interval derived from State of Montana paired case data
serial_interval_mean <- 4.42
serial_interval_sd <- 3.51

mt_r_results <- estimate_R(mt_incidence_data, method="parametric_si", 
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
      mutate(I = sum(case)) %>% 
      select(dates, I) %>% 
      arrange(dates)
   
   mt_incidence_data <- incidence(mt_analysis_data$dates, last_date = latest_date)
   
   state_n <- as.data.frame(mt_incidence_data$dates)
   time_var <- nrow(state_n)
   time_start <- seq(2, time_var-13)
   time_end <- time_start + 13
   
   
   # Serial Interval derived from State of Montana case data
   serial_interval_mean <- 4.42
   serial_interval_sd <- 3.51
   
   mt_r_results <- estimate_R(mt_incidence_data, method="parametric_si", 
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
      mutate(I = sum(case)) %>% 
      select(dates, I) %>% 
      arrange(dates)
   
   mt_incidence_data <- incidence(mt_analysis_data$dates, last_date = latest_date)
   
   state_n <- as.data.frame(mt_incidence_data$dates)
   time_var <- nrow(state_n)
   time_start <- seq(2, time_var-13)
   time_end <- time_start + 13
   
   
   # Serial Interval derived from State of Montana case data
   serial_interval_mean <- 4.42
   serial_interval_sd <- 3.51
   
   mt_r_results <- estimate_R(mt_incidence_data, method="parametric_si", 
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
ravalli_r <- county_analysis_function("Ravalli", "Ravalli County")


# Bind files and save
all_regions_r <- rbind(state_r, reg1_r, reg2_r, reg3_r, reg4_r, reg5_r,
                       missoula_r, gallatin_r, yellowstone_r, bighorn_r,
                       lake_r, lewisandclark_r, flathead_r, cascade_r,
                       silverbow_r, rosebud_r, glacier_r, roosevelt_r,
                       ravalli_r) %>% 
   mutate(daily_cases = incidence) %>% 
   mutate(mean_r = round(mean_r, digits = 2),
          median_r = round(median_r, digits = 2),
          sd_r = round(sd_r, digits = 2),
          cl_low = round(cl_low, digits = 2),
          cl_high = round(cl_high, digits = 2))

write_csv(all_regions_r, "C:/R/covid19/state_daily_results/all_regions_r_new.csv", na = " ")



# Send files to the sftp server
# host = 'elbastion.dbs.umt.edu'
# port = 22
# username = 'celftp'
# password  = 'celftp'
# remotepath = '/celFtpFiles/covid19/Rt/incoming/'

#sftpUpload("elbastion.dbs.umt.edu", "celftp", "celftp",
#           "/celFtpFiles/covid19/Rt/incoming/all_regions_r_new.csv",
#           "C:/R/covid19/state_daily_results/all_regions_r_new.csv")

sftpUpload("celftp.nephelai.net", "celftp", "celftp@umt",
           "/home/celftp/data/all_regions_r_new.csv",
           "C:/R/covid19/state_daily_results/all_regions_r_new.csv")



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


options(gargle_oauth_email = "ethanwalker86@gmail.com")
drive_auth(email = "ethanwalker86@gmail.com")
gs4_auth(token = drive_token())


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

ravalli_data <- r_filter_dates %>% 
   filter(region == "Ravalli County") %>% 
   select(dates, cl_low, mean_r, cl_high) %>% 
   rename("Dates" = dates,
          "Lower_Confidence_Limit" = cl_low,
          "Mean_R" = mean_r,
          "Upper_Confidence_Limit" = cl_high)

sheet_write(ravalli_data, 
            ss = "https://docs.google.com/spreadsheets/d/1L1pBAp0e5RvU7x4IfnfyeJ5uaQs0OTuYXYmgJB0r_oE/edit#gid=0",
            sheet = 13)


options(gargle_oauth_email = "ethanwalker86@gmail.com")
drive_auth(email = "ethanwalker86@gmail.com")
gs4_auth(token = drive_token())


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

ravalli_data <- all_regions_r %>% 
   filter(region == "Ravalli County") %>% 
   select(dates, incidence) %>% 
   rename("Dates" = dates,
          "Cases" = incidence)

sheet_write(ravalli_data, 
            ss = "https://docs.google.com/spreadsheets/d/1X0vxDLVxQT_XrPgNMtyRwT1kTVQSN2PzQPYDWi2oxBo/edit#gid=0",
            sheet = 13)



#################### Run and save state hospitalization data

state_hosp <- mt_case_data %>% 
   filter(!is.na(age_group)) %>% 
   filter(dates < Sys.Date()-30) %>% 
   group_by(age_group) %>% 
   mutate(age_group_cases = sum(case)) %>% 
   mutate(hosp = if_else(hospitalization == "Hosp: Yes", 1, 0),
          death = if_else(deceased == "Deceased: Yes", 1, 0),
          hosp_yes = sum(hosp, na.rm = TRUE),
          hosp_no = age_group_cases - hosp_yes,
          hosp_percent = round(hosp_yes/age_group_cases*100, digits = 2),
          deaths = sum(death, na.rm = TRUE),
          deaths_percent = round(deaths/age_group_cases*100, digits = 2)) %>% 
   ungroup() %>% 
   distinct(age_group, .keep_all = TRUE) %>% 
   select(age_group, age_group_cases, hosp_yes, hosp_no, hosp_percent) %>% 
   arrange(age_group)

write_csv(state_hosp, "C:/R/covid19/state_daily_results/state_hosp.csv", na = " ")


state_deaths <- mt_case_data %>% 
   filter(!is.na(age_group)) %>% 
   filter(dates < Sys.Date()-14) %>% 
   group_by(age_group) %>% 
   mutate(age_group_cases = sum(case)) %>% 
   mutate(hosp = if_else(hospitalization == "Hosp: Yes", 1, 0),
          death = if_else(deceased == "Deceased: Yes", 1, 0),
          hosp_yes = sum(hosp, na.rm = TRUE),
          hosp_no = age_group_cases - hosp_yes,
          hosp_percent = round(hosp_yes/age_group_cases*100, digits = 2),
          deaths = sum(death, na.rm = TRUE),
          deaths_percent = round(deaths/age_group_cases*100, digits = 2)) %>% 
   ungroup() %>% 
   distinct(age_group, .keep_all = TRUE) %>% 
   select(age_group, age_group_cases, deaths, deaths_percent) %>% 
   arrange(age_group)

write_csv(state_deaths, "C:/R/covid19/state_daily_results/state_deaths.csv", na = " ")


state_hosp_month <- mt_case_data %>% 
   filter(!is.na(age_group)) %>% 
   filter(dates < Sys.Date()-30) %>% 
   mutate(onset_month = lubridate::month(dates, label = TRUE)) %>% 
   group_by(age_group, onset_month) %>% 
   mutate(age_group_cases = sum(case)) %>% 
   mutate(hosp = if_else(hospitalization == "Hosp: Yes", 1, 0),
          death = if_else(deceased == "Deceased: Yes", 1, 0),
          hosp_yes = sum(hosp, na.rm = TRUE),
          hosp_no = age_group_cases - hosp_yes,
          hosp_percent = round(hosp_yes/age_group_cases*100, digits = 2),
          deaths = sum(death, na.rm = TRUE),
          deaths_percent = round(deaths/age_group_cases*100, digits = 2)) %>% 
   ungroup() %>% 
   distinct(onset_month, age_group, .keep_all = TRUE) %>% 
   select(onset_month, age_group, age_group_cases, hosp_yes, hosp_no, hosp_percent) %>% 
   arrange(onset_month, age_group)

write_csv(state_hosp_month, "C:/R/covid19/state_daily_results/state_hosp_month.csv", na = " ")


state_deaths_month <- mt_case_data %>% 
   filter(!is.na(age_group)) %>% 
   filter(dates < Sys.Date()-14) %>% 
   mutate(onset_month = lubridate::month(dates, label = TRUE)) %>% 
   group_by(age_group, onset_month) %>% 
   mutate(age_group_cases = sum(case)) %>% 
   mutate(hosp = if_else(hospitalization == "Hosp: Yes", 1, 0),
          death = if_else(deceased == "Deceased: Yes", 1, 0),
          hosp_yes = sum(hosp, na.rm = TRUE),
          hosp_no = age_group_cases - hosp_yes,
          hosp_percent = round(hosp_yes/age_group_cases*100, digits = 2),
          deaths = sum(death, na.rm = TRUE),
          deaths_percent = round(deaths/age_group_cases*100, digits = 2)) %>% 
   ungroup() %>% 
   distinct(onset_month, age_group, .keep_all = TRUE) %>% 
   select(onset_month, age_group, age_group_cases, deaths, deaths_percent) %>% 
   arrange(onset_month, age_group)

write_csv(state_deaths_month, "C:/R/covid19/state_daily_results/state_deaths_month.csv", na = " ")


reg_hosp <- mt_case_data %>% 
   filter(!is.na(age_group)) %>% 
   filter(dates < Sys.Date()-30) %>% 
   group_by(age_group, region) %>% 
   mutate(age_group_cases = sum(case)) %>% 
   mutate(hosp = if_else(hospitalization == "Hosp: Yes", 1, 0),
          death = if_else(deceased == "Deceased: Yes", 1, 0),
          hosp_yes = sum(hosp, na.rm = TRUE),
          hosp_no = age_group_cases - hosp_yes,
          hosp_percent = round(hosp_yes/age_group_cases*100, digits = 2),
          deaths = sum(death, na.rm = TRUE),
          deaths_percent = round(deaths/age_group_cases*100, digits = 2)) %>% 
   ungroup() %>% 
   distinct(age_group, region, .keep_all = TRUE) %>% 
   select(region, age_group, age_group_cases, hosp_yes, hosp_no, hosp_percent) %>% 
   arrange(region, age_group)

write_csv(reg_hosp, "C:/R/covid19/state_daily_results/reg_hosp.csv", na = " ")


reg_deaths <- mt_case_data %>% 
   filter(!is.na(age_group)) %>% 
   filter(dates < Sys.Date()-14) %>% 
   group_by(age_group, region) %>% 
   mutate(age_group_cases = sum(case)) %>% 
   mutate(hosp = if_else(hospitalization == "Hosp: Yes", 1, 0),
          death = if_else(deceased == "Deceased: Yes", 1, 0),
          hosp_yes = sum(hosp, na.rm = TRUE),
          hosp_no = age_group_cases - hosp_yes,
          hosp_percent = round(hosp_yes/age_group_cases*100, digits = 2),
          deaths = sum(death, na.rm = TRUE),
          deaths_percent = round(deaths/age_group_cases*100, digits = 2)) %>% 
   ungroup() %>% 
   distinct(age_group, region, .keep_all = TRUE) %>% 
   select(region, age_group, age_group_cases, deaths, deaths_percent) %>% 
   arrange(region, age_group)

write_csv(reg_deaths, "C:/R/covid19/state_daily_results/reg_deaths.csv", na = " ")



hosp_data_initial <- mt_case_data %>% 
   mutate(age_group_new = factor(age_group, 
                                 levels = c("0 to 9", "10 to 19", "20 to 29", 
                                            "30 to 39", "40 to 49", "50 to 59", 
                                            "60 to 69", "70 to 79", "80+"),
                                 labels = c("0 to 9, 0.12", "10 to 19, 0.12", "20 to 29, 0.13", 
                                            "30 to 39, 0.13", "40 to 49, 0.11", "50 to 59, 0.12", 
                                            "60 to 69, 0.14", "70 to 79, 0.08", "80+, 0.04"))) %>% 
   separate(age_group_new, c("age_group_new", "age_group_new_percent"), sep = ",") %>% 
   mutate(age_group_new_percent = as.numeric(age_group_new_percent),
          state_pop = as.numeric(1068778)) 


age_rates <- hosp_data_initial %>% 
   select(case, age_group_new, age_group_new_percent, state_pop) %>% 
   mutate(total_cases = n()) %>% 
   group_by(age_group_new) %>% 
   mutate(group_cases = n(),
          group_prop = group_cases/total_cases*100,
          group_pop = age_group_new_percent*state_pop,
          group_rate = group_cases/group_pop*100000) %>% 
   ungroup() %>% 
   distinct(age_group_new, .keep_all = TRUE) %>% 
   mutate(region = "Montana") %>% 
   rename(age_group = age_group_new,
          age_group_cases = group_cases) %>% 
   mutate(age_group_rate = round(group_rate, digits = 0)) %>% 
   select(age_group, age_group_cases, age_group_rate) %>% 
   filter(!is.na(age_group)) %>% 
   arrange(age_group)

hosp_data <- hosp_data_initial %>% 
   select(case, age_group_new, age_group_new_percent, state_pop, hospitalization) %>% 
   mutate(total_cases = n(),
          hospitalization_status = if_else(hospitalization == "Hosp: Yes", 
                                           "Yes", "No")) %>% 
   group_by(age_group_new, hospitalization_status) %>% 
   mutate(group_cases = n(),
          group_prop = group_cases/total_cases*100,
          group_pop = age_group_new_percent*state_pop,
          group_rate = group_cases/group_pop*100000) %>% 
   ungroup() %>% 
   distinct(age_group_new, hospitalization_status, .keep_all = TRUE) 

hosp_cases <- hosp_data %>% 
   arrange(age_group_new, hospitalization_status) %>% 
   pivot_wider(names_from = hospitalization_status, values_from = group_cases) %>% 
   rename(age_group = age_group_new,
          Hospitalized = Yes,
          Not_Hospitalized = No) %>% 
   select(age_group, Hospitalized, Not_Hospitalized) %>% 
   filter(!is.na(age_group)) 

hosp_rates <- hosp_data %>% 
   arrange(age_group_new, hospitalization_status) %>% 
   mutate(age_group_rate = round(group_rate, digits = 0)) %>% 
   pivot_wider(names_from = hospitalization_status, values_from = age_group_rate) %>% 
   rename(age_group = age_group_new,
          Hospitalized = Yes,
          Not_Hospitalized = No) %>% 
   select(age_group, Hospitalized, Not_Hospitalized) %>% 
   filter(!is.na(age_group)) 



# Send files to the sftp server
# host = 'elbastion.dbs.umt.edu'
# port = 22
# username = 'celftp'
# password  = 'celftp'
# remotepath = '/celFtpFiles/covid19/Rt/incoming/'

#sftpUpload("elbastion.dbs.umt.edu", "celftp", "celftp",
#           "/celFtpFiles/covid19/Rt/incoming/state_hosp.csv",
#           "C:/R/covid19/state_daily_results/state_hosp.csv")

sftpUpload("celftp.nephelai.net", "celftp", "celftp@umt",
           "/home/celftp/data/state_hosp.csv",
           "C:/R/covid19/state_daily_results/state_hosp.csv")

#sftpUpload("elbastion.dbs.umt.edu", "celftp", "celftp",
#           "/celFtpFiles/covid19/Rt/incoming/state_deaths.csv",
#           "C:/R/covid19/state_daily_results/state_deaths.csv")

sftpUpload("celftp.nephelai.net", "celftp", "celftp@umt",
           "/home/celftp/data/state_deaths.csv",
           "C:/R/covid19/state_daily_results/state_deaths.csv")

#sftpUpload("elbastion.dbs.umt.edu", "celftp", "celftp",
#           "/celFtpFiles/covid19/Rt/incoming/reg_hosp.csv",
#           "C:/R/covid19/state_daily_results/reg_hosp.csv")

sftpUpload("celftp.nephelai.net", "celftp", "celftp@umt",
           "/home/celftp/data/reg_hosp.csv",
           "C:/R/covid19/state_daily_results/reg_hosp.csv")

#sftpUpload("elbastion.dbs.umt.edu", "celftp", "celftp",
#           "/celFtpFiles/covid19/Rt/incoming/reg_death.csv",
#           "C:/R/covid19/state_daily_results/reg_death.csv")

sftpUpload("celftp.nephelai.net", "celftp", "celftp@umt",
           "/home/celftp/data/reg_death.csv",
           "C:/R/covid19/state_daily_results/reg_death.csv")

#sftpUpload("elbastion.dbs.umt.edu", "celftp", "celftp",
#           "/celFtpFiles/covid19/Rt/incoming/state_hosp_month.csv",
#           "C:/R/covid19/state_daily_results/state_hosp_month.csv")

sftpUpload("celftp.nephelai.net", "celftp", "celftp@umt",
           "/home/celftp/data/state_hosp_month.csv",
           "C:/R/covid19/state_daily_results/state_hosp_month.csv")

#sftpUpload("elbastion.dbs.umt.edu", "celftp", "celftp",
#           "/celFtpFiles/covid19/Rt/incoming/state_deaths_month.csv",
#           "C:/R/covid19/state_daily_results/state_deaths_month.csv")

sftpUpload("celftp.nephelai.net", "celftp", "celftp@umt",
           "/home/celftp/data/state_deaths_month.csv",
           "C:/R/covid19/state_daily_results/state_deaths_month.csv")


########## Push hosp data by age group to Google sheets

options(gargle_oauth_email = "ethanwalker86@gmail.com")
drive_auth(email = "ethanwalker86@gmail.com")
gs4_auth(token = drive_token())


# Write age case/rate data to google
mt_data <- age_rates %>% 
   select(age_group, age_group_cases) %>% 
   rename("Age_Group" = age_group,
          "Cases" = age_group_cases)

sheet_write(mt_data, 
            ss = "https://docs.google.com/spreadsheets/d/1H5e3OPlxzlCacDAD_foj72EqZEzyPZ66FUyh-fxokag/edit#gid=0",
            sheet = 1)

mt_data <- age_rates %>% 
   select(age_group, age_group_rate) %>% 
   rename("Age_Group" = age_group,
          "Rate per 100,000 population" = age_group_rate)

sheet_write(mt_data, 
            ss = "https://docs.google.com/spreadsheets/d/1H5e3OPlxzlCacDAD_foj72EqZEzyPZ66FUyh-fxokag/edit#gid=0",
            sheet = 2)

mt_data <- hosp_cases %>% 
   select(age_group, Hospitalized, Not_Hospitalized) %>% 
   rename("Age_Group" = age_group)

sheet_write(mt_data, 
            ss = "https://docs.google.com/spreadsheets/d/1H5e3OPlxzlCacDAD_foj72EqZEzyPZ66FUyh-fxokag/edit#gid=0",
            sheet = 3)

mt_data <- hosp_rates %>% 
   select(age_group, Hospitalized, Not_Hospitalized) %>% 
   rename("Age_Group" = age_group)

sheet_write(mt_data, 
            ss = "https://docs.google.com/spreadsheets/d/1H5e3OPlxzlCacDAD_foj72EqZEzyPZ66FUyh-fxokag/edit#gid=0",
            sheet = 4)


# County level data for spatial analysis
spatial_data <- mt_case_data %>% 
   select(dates, county, spatial_age, sex, hospitalization, deceased, case) %>% 
   mutate(spatial_age = factor(spatial_age,
                               levels = c("0-4", "5-24", "25-64", "65+"),
                               labels = c("0 to 4", "5 to 24", "25 to 64", "65+")),
          hospitalization_status = if_else(hospitalization == "Hosp: Yes", 
                                           1, 0),
          deceased_status = if_else(deceased == "Deceased: Yes", 
                                    1, 0)) %>% 
   group_by(county, spatial_age) %>% 
   mutate(latest_date = max(dates)) %>% 
   summarize(county_age_case = sum(case),
             county_age_hosp = sum(hospitalization_status),
             county_age_deaths = sum(deceased_status),
             latest_onset_date = latest_date) %>% 
   distinct(county, spatial_age, .keep_all = T) %>% 
   ungroup() %>% 
   arrange(county, spatial_age) 

write_csv(spatial_data, "C:/R/covid19/state_daily_results/spatial_data.csv", na = " ")


sftpUpload("celftp.nephelai.net", "celftp", "celftp@umt",
           "/home/celftp/data/spatial_data.csv",
           "C:/R/covid19/state_daily_results/spatial_data.csv")

