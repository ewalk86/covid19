# Estimate R, save results/plots
# Ethan Walker
# 24 Nov 2020

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


# Load/format case data
mt_case_data <- read_xlsx(paste0(file_path, "Input/uom_covid_11232020_2.xlsx"),
                          sheet = 1, skip = 1,
                          col_names = c("jurisdiction", "spec_coll_date", "inv_start_date", 
                                        "diagnosis_date", "county_fips", "state_num", 
                                        "age", "age_unit", "sex", "hosp_adm_date", 
                                        "hospitalization", "onset_date"),
                          col_types = c("text", "date", "date",
                                        "date", "numeric", "numeric", 
                                        "numeric", "text", "text", "date", 
                                        "text", "date")) %>% 
   rownames_to_column(var = "case_no") %>% 
   mutate(jurisdiction = str_to_title(jurisdiction),
          case = 1) %>% 
   left_join(mt_county_fips, by = "county_fips") %>% 
   filter(jurisdiction != "Out Of State") %>% 
   filter(!is.na(county)) %>% 
   mutate(age = as.numeric(age),
          age_group = cut(age, breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 130),
                          labels = c("0 to 9", "10 to 19", "20 to 29", 
                                     "30 to 39", "40 to 49", "50 to 59", 
                                     "60 to 69", "70 to 79", "80+"), right = FALSE)) %>% 
   mutate(hospitalization = factor(hospitalization,
                                   levels = c("Y", "N", "UNK"),
                                   labels = c("Hosp: Yes", "Hosp: No", 
                                              "Hosp: Unknown")),
          sex = factor(sex,
                       levels = c("F", "M", "U"),
                       labels = c("Female", "Male", "Unknown"))) %>% 
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
                                 inv_start_date - mean_onset_inv_diff, onset_date_3),
          onset_date_5 = if_else(is.na(onset_date) & !is.na(inv_start_date), 
                                 inv_start_date - mean_onset_inv_diff, onset_date)) %>% 
   separate(onset_date_4, into = c("dates", "time"), sep = " ") %>% 
   separate(onset_date_5, into = c("dates_2", "time"), sep = " ")

dates_onset_coll_diag_inv <- mt_case_data %>% 
   mutate(dates = ymd(dates)) %>% 
   select(-time) %>% 
   filter(!is.na(dates)) %>% 
   arrange(dates) %>% 
   ungroup() 

dates_onset_inv <- mt_case_data %>% 
   mutate(dates = ymd(dates_2)) %>% 
   select(-time) %>% 
   filter(!is.na(dates)) %>% 
   arrange(dates) %>% 
   ungroup() 

dates_state_pull <- read_rds(paste0(file_path, "Output/state_data_clean.rds"))


#################### Run analysis and print results ########################


##### State results - dates_onset_coll_diag_inv #####

# Format analysis data
mt_analysis_data <- dates_onset_coll_diag_inv %>% 
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
serial_interval_mean <- 5.29
serial_interval_sd <- 4.45

mt_r_results <- estimate_R(mt_incidence_data, method="parametric_si", 
                           config = make_config(list(mean_si = serial_interval_mean, 
                                                     std_si = serial_interval_sd,
                                                     t_start =  time_start,
                                                     t_end = time_end)))

# Format and save analysis results
colomns <- c(1:4, 8)
state_r <- (mt_r_results$R[,colomns]) %>% 
   mutate(method = "Onset_Coll_Diag_Inv") %>% 
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
r_dates_onset_coll_diag_inv <- cbind(state_r, state_dates_new)


##### State results - dates_onset_inv #####

# Format analysis data
mt_analysis_data <- dates_onset_inv %>% 
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
serial_interval_mean <- 5.29
serial_interval_sd <- 4.45

mt_r_results <- estimate_R(mt_incidence_data, method="parametric_si", 
                           config = make_config(list(mean_si = serial_interval_mean, 
                                                     std_si = serial_interval_sd,
                                                     t_start =  time_start,
                                                     t_end = time_end)))

# Format and save analysis results
colomns <- c(1:4, 8)
state_r <- (mt_r_results$R[,colomns]) %>% 
   mutate(method = "Onset_Investigation") %>% 
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
r_dates_onset_inv <- cbind(state_r, state_dates_new)



##### State results - dates_state_pull #####

# Format analysis data
mt_analysis_data <- dates_state_pull %>% 
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
serial_interval_mean <- 5.29
serial_interval_sd <- 4.45

mt_r_results <- estimate_R(mt_incidence_data, method="parametric_si", 
                           config = make_config(list(mean_si = serial_interval_mean, 
                                                     std_si = serial_interval_sd,
                                                     t_start =  time_start,
                                                     t_end = time_end)))

# Format and save analysis results
colomns <- c(1:4, 8)
state_r <- (mt_r_results$R[,colomns]) %>% 
   mutate(method = "State_Pull") %>% 
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
r_dates_state_pull <- cbind(state_r, state_dates_new)




# Bind files and save
date_comparison_r <- rbind(r_dates_onset_coll_diag_inv, r_dates_state_pull) %>% 
   mutate(daily_cases = incidence) %>% 
   mutate(mean_r = round(mean_r, digits = 2),
          median_r = round(median_r, digits = 2),
          sd_r = round(sd_r, digits = 2),
          cl_low = round(cl_low, digits = 2),
          cl_high = round(cl_high, digits = 2))

write_csv(date_comparison_r, "C:/R/covid19/state_daily_results/date_comparison_r.csv", na = " ")


# Plot results

date_comparison_r_plot <- date_comparison_r %>% 
   filter(dates > "2020-07-31") %>% 
   ggplot() +
   geom_line(aes(dates, mean_r, color = method), size = 1.5) +
   #geom_line(aes(dates, cl_low, color = method), size = 0.5) +
   #geom_line(aes(dates, cl_high, color = method), size = 0.5) +
   labs(title = "COVID-19 Rolling 14-day R-number, Montana",
        color = "") +
   ylab("R-number") +
   xlab("") +
   geom_hline(yintercept = 1, color = "red", size = 1.2) +
   scale_x_date(date_breaks = "3 days", date_labels = "%d-%b") +
   scale_y_continuous(breaks = seq(0, 5, 0.1), labels = seq(0, 5, 0.1)) +
   theme_minimal() +
   theme(strip.text = element_text(size = 16, colour = "black"),
         title = element_text(size = 18, colour = "black"),
         panel.grid = element_blank(),
         panel.grid.major.y = element_line(colour = "grey"),
         axis.text.x = element_text(size = 16, colour = "black", 
                                    angle = 90, vjust = 0.4),
         axis.text.y = element_text(size = 16, colour = "black"),
         legend.text = element_text(size = 16, colour = "black"),
         axis.title.y = element_text(size = 16, colour = "black",
                                     margin = unit(c(0, 5, 0, 0), "mm")),
         axis.title.x = element_text(size = 16, colour = "black",
                                     margin = unit(c(5, 0, 0, 0), "mm")),
         axis.line.x = element_blank(), 
         axis.line.y = element_blank(), 
         axis.ticks = element_blank())  
date_comparison_r_plot

ggsave("C:/R/covid19/date_comparison_r_plot.png", width = 12, height = 6)


############# Compare error in the estimates of onset date

mt_case_data <- read_xlsx(paste0(file_path, "Input/uom_covid_11232020_2.xlsx"),
                          sheet = 1, skip = 1,
                          col_names = c("jurisdiction", "spec_coll_date", "inv_start_date", 
                                        "diagnosis_date", "county_fips", "state_num", 
                                        "age", "age_unit", "sex", "hosp_adm_date", 
                                        "hospitalization", "onset_date"),
                          col_types = c("text", "date", "date",
                                        "date", "numeric", "numeric", 
                                        "numeric", "text", "text", "date", 
                                        "text", "date")) %>% 
   rownames_to_column(var = "case_no") %>% 
   mutate(jurisdiction = str_to_title(jurisdiction),
          case = 1) %>% 
   left_join(mt_county_fips, by = "county_fips") %>% 
   filter(jurisdiction != "Out Of State") %>% 
   filter(!is.na(county)) %>% 
   mutate(age = as.numeric(age),
          age_group = cut(age, breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 130),
                          labels = c("0 to 9", "10 to 19", "20 to 29", 
                                     "30 to 39", "40 to 49", "50 to 59", 
                                     "60 to 69", "70 to 79", "80+"), right = FALSE)) %>% 
   mutate(hospitalization = factor(hospitalization,
                                   levels = c("Y", "N", "UNK"),
                                   labels = c("Hosp: Yes", "Hosp: No", 
                                              "Hosp: Unknown")),
          sex = factor(sex,
                       levels = c("F", "M", "U"),
                       labels = c("Female", "Male", "Unknown"))) %>% 
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
          onset_date_coll = spec_coll_date - mean_onset_coll_diff,
          onset_coll_diff2 = as.duration(interval(onset_date, onset_date_coll)),
          onset_coll_diff2 = as.numeric(onset_coll_diff2)/86400,
          mean_onset_coll_diff2 = mean(onset_coll_diff2, na.rm = TRUE),
          sd_onset_coll_diff2 = sd(onset_coll_diff2, na.rm = TRUE),
          onset_date_diag = diagnosis_date - mean_onset_diag_diff,
          onset_diag_diff2 = as.duration(interval(onset_date, onset_date_diag)),
          onset_diag_diff2 = as.numeric(onset_diag_diff2)/86400,
          mean_onset_diag_diff2 = mean(onset_diag_diff2, na.rm = TRUE),
          sd_onset_diag_diff2 = sd(onset_diag_diff2, na.rm = TRUE),
          onset_date_inv = inv_start_date - mean_onset_inv_diff,
          onset_inv_diff2 = as.duration(interval(onset_date, onset_date_inv)),
          onset_inv_diff2 = as.numeric(onset_inv_diff2)/86400,
          mean_onset_inv_diff2 = mean(onset_inv_diff2, na.rm = TRUE),
          sd_onset_inv_diff2 = sd(onset_inv_diff2, na.rm = TRUE)) 
