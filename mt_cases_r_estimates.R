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


# Load serial interval and case data
mt_si_data <- read_xlsx(paste0(file_path, "Input/SI_Local_v_Import Data_07.01.2020.xlsx"),
                          sheet = 1) %>% 
   select(-Pair_No) %>% 
   rename(EL = EL_Infector_Lower,
          ER = ER_Infector_Upper,
          SL = SL_Infected_Lower,
          SR = SR_Infected_Upper) %>% 
   mutate(SL2 = if_else(SL<1, SL + (1 - SL), SL),
          SR2 = if_else(SL<1, SR + (1 - SL), SR),
          type = 0) %>% 
   select(EL, ER, SL2, SR2, type) %>% 
   rename(SL = SL2,
          SR = SR2) %>% 
   mutate_all(as.integer)


# Load/format case data
mt_case_data <- read_xlsx(paste0(file_path, "Input/SI_Local_v_Import Data_07.01.2020.xlsx"),
                                 sheet = 2) %>% 
   rename_all(tolower) %>% 
   select(-case_no) %>% 
   rownames_to_column(var = "case_no") %>% 
   mutate(local = if_else(local_import == 0, 1, 0),
          imported = local_import,
          date_reported = ymd(date_reported),
          dates = ymd(symptom_onset_date),
          case = 1) %>% 
   rename(hospitalization = "ever_hospitalized") %>% 
   left_join(counties_regions, by = "county") %>% 
   mutate(age_group_new = if_else(age_group == "80-89" | age_group == "90-99",
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
                                              "Hosp: Past", "Hosp: Unknown"))) %>% 
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
   pivot_wider(names_from = "county", values_from = "case") %>% 
   select(-"NA")


# Filter and format wide data for state and 5 health regions
state_wide_date <- state_data_wide %>% 
   select(region, dates, 8:59) %>% 
   mutate(region = "state") %>% 
   group_by(region, dates) %>% 
   mutate_all(sum, na.rm = TRUE) %>% 
   mutate(case = 1,
          daily_cases = sum(case)) %>% 
   arrange(dates) %>% 
   distinct(dates, .keep_all = TRUE) %>% 
   ungroup() %>% 
   mutate(cumulative_cases = cumsum(daily_cases)) %>% 
   select(region, dates, daily_cases, cumulative_cases, 3:53)

reg1_wide_date <- state_data_wide %>% 
   filter(region == 1) %>% 
   select(region, dates, 8:59) %>% 
   group_by(region, dates) %>% 
   mutate_all(sum, na.rm = TRUE) %>% 
   mutate(case = 1,
          daily_cases = sum(case)) %>% 
   arrange(dates) %>% 
   distinct(dates, .keep_all = TRUE) %>% 
   ungroup() %>% 
   mutate(cumulative_cases = cumsum(daily_cases)) %>% 
   select(region, dates, daily_cases, cumulative_cases, 3:53)

reg2_wide_date <- state_data_wide %>% 
   filter(region == 2) %>% 
   select(region, dates, 8:59) %>% 
   group_by(region, dates) %>% 
   mutate_all(sum, na.rm = TRUE) %>% 
   mutate(case = 1,
          daily_cases = sum(case)) %>% 
   arrange(dates) %>% 
   distinct(dates, .keep_all = TRUE) %>% 
   ungroup() %>% 
   mutate(cumulative_cases = cumsum(daily_cases)) %>% 
   select(region, dates, daily_cases, cumulative_cases, 3:53)

reg3_wide_date <- state_data_wide %>% 
   filter(region == 3) %>% 
   select(region, dates, 8:59) %>% 
   group_by(region, dates) %>% 
   mutate_all(sum, na.rm = TRUE) %>% 
   mutate(case = 1,
          daily_cases = sum(case)) %>% 
   arrange(dates) %>% 
   distinct(dates, .keep_all = TRUE) %>% 
   ungroup() %>% 
   mutate(cumulative_cases = cumsum(daily_cases)) %>% 
   select(region, dates, daily_cases, cumulative_cases, 3:53)

reg4_wide_date <- state_data_wide %>% 
   filter(region == 4) %>% 
   select(region, dates, 8:59) %>% 
   group_by(region, dates) %>% 
   mutate_all(sum, na.rm = TRUE) %>% 
   mutate(case = 1,
          daily_cases = sum(case)) %>% 
   arrange(dates) %>% 
   distinct(dates, .keep_all = TRUE) %>% 
   ungroup() %>% 
   mutate(cumulative_cases = cumsum(daily_cases)) %>% 
   select(region, dates, daily_cases, cumulative_cases, 3:53)

reg5_wide_date <- state_data_wide %>% 
   filter(region == 5) %>% 
   select(region, dates, 8:59) %>% 
   group_by(region, dates) %>% 
   mutate_all(sum, na.rm = TRUE) %>% 
   mutate(case = 1,
          daily_cases = sum(case)) %>% 
   arrange(dates) %>% 
   distinct(dates, .keep_all = TRUE) %>% 
   ungroup() %>% 
   mutate(cumulative_cases = cumsum(daily_cases)) %>% 
   select(region, dates, daily_cases, cumulative_cases, 3:53)

all_data_wide <- rbind(state_wide_date, reg1_wide_date, reg2_wide_date,
                       reg3_wide_date, reg4_wide_date, reg5_wide_date)




#################### Run analysis and print results
## State results
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


state_n <- as.data.frame(mt_li_inc_data$dates)
time_var <- nrow(state_n)
time_start <- seq(2, time_var-13)
time_end <- time_start + 13


MCMC_seed <- 1
overall_seed <- 2
mcmc_control <- make_mcmc_control(seed = MCMC_seed, burnin = 1000, thin = 10)
dist <- "G"  # fitting a Gamma distribution for the SI
si_config <- make_config(list(si_parametric_distr = dist, 
                              mcmc_control = mcmc_control, seed = overall_seed, n1 = 500, n2 = 50,
                              t_start = time_start, t_end = time_end))

mt_r_results <- estimate_R(mt_li_inc_data, method = "si_from_data", 
                       si_data = mt_si_data, config = si_config)


# calculate mean serial interval using SI data
mean_si <- mt_r_results$SI.Moments %>% 
   summarize(mean_si = mean(Mean), sd_si = mean(Std)) 



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
state_r_clean <- cbind(state_r, state_dates_new)

total_cases <- sum(mt_case_data$case)
date_today <- format(Sys.Date(), "%d %b %Y")


state_r_plot <- state_r_clean %>% 
   ggplot() +
   geom_line(aes(dates, mean_r), size = 1.5, color = "black") +
   geom_line(aes(dates, cl_low), size = 1.5, color = "grey") +
   geom_line(aes(dates, cl_high), size = 1.5, color = "grey") +
   labs(title = "COVID-19 Rolling 14-day R-values, State of Montana, 2020",
        color = "") +
   ylab("R-value") +
   xlab("") +
   geom_hline(yintercept = 1, color = "black", size = 1.2) +
   scale_x_date(date_breaks = "3 days", date_labels = "%d-%b") +
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

state_inc_plot <- mt_case_data %>% 
   ggplot() +
   geom_col(aes(dates, case), 
            fill = "steelblue") +
   labs(title = paste0("COVID-19 Cases in State of Montana, 2020", 
                       " (Total cases = ", total_cases, ")"),
        subtitle = paste0("Data current as of ", date_today)) +
   ylab("Number of Cases") +
   xlab("") +
   scale_x_date(date_breaks = "3 day", date_labels = "%d-%b") +
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



## Region 1 results
# Format analysis data
mt_li_analysis_data <- mt_case_data %>% 
   filter(region == 1) %>% 
   group_by(dates) %>% 
   mutate(local = sum(local),
          imported = sum(imported)) %>% 
   select(dates, local, imported) %>% 
   distinct(dates, .keep_all = TRUE) %>% 
   arrange(dates)


reg1_data <- mt_li_analysis_data 

mt_data <- reg1_data %>% 
   mutate(I = local + imported) %>% 
   select(dates, I) 

mt_incidence_data <- incidence(mt_data$dates, last_date = latest_date)

mt_li_inc_data <- as.data.frame(mt_incidence_data$dates) %>% 
   mutate(dates = ymd(mt_incidence_data$dates)) %>% 
   select(dates) %>% 
   left_join(reg1_data, by = "dates") %>% 
   mutate(local = if_else(is.na(local), 0, local),
          imported = if_else(is.na(imported), 0, imported)) 


state_n <- as.data.frame(mt_li_inc_data$dates)
time_var <- nrow(state_n)
time_start <- seq(2, time_var-13)
time_end <- time_start + 13


MCMC_seed <- 1
overall_seed <- 2
mcmc_control <- make_mcmc_control(seed = MCMC_seed, burnin = 1000, thin = 10)
dist <- "G"  # fitting a Gamma distribution for the SI
si_config <- make_config(list(si_parametric_distr = dist, 
                              mcmc_control = mcmc_control, seed = overall_seed, n1 = 500, n2 = 50,
                              t_start = time_start, t_end = time_end))

mt_r_results <- estimate_R(mt_li_inc_data, method = "si_from_data", 
                           si_data = mt_si_data, config = si_config)



colomns <- c(1:4, 8)
state_r <- (mt_r_results$R[,colomns]) %>% 
   mutate(region = "1") %>% 
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
reg1_r_clean <- cbind(state_r, state_dates_new)

reg1_data_clean <- mt_case_data %>% 
   filter(region == 1)
total_cases <- sum(reg1_data_clean$case)


state_r_plot <- reg1_r_clean %>% 
   ggplot() +
   geom_line(aes(dates, mean_r), size = 1.5, color = "black") +
   geom_line(aes(dates, cl_low), size = 1.5, color = "grey") +
   geom_line(aes(dates, cl_high), size = 1.5, color = "grey") +
   labs(title = "COVID-19 Rolling 14-day R-values, Montana Region 1, 2020",
        color = "") +
   ylab("R-value") +
   xlab("") +
   geom_hline(yintercept = 1, color = "black", size = 1.2) +
   scale_x_date(date_breaks = "3 days", date_labels = "%d-%b") +
   #scale_y_continuous(breaks = seq(0, 10, 0.5), labels = seq(0, 10, 0.5)) +
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
ggsave("C:/R/covid19/state_daily_results/reg1_r_plot.png", width = 10, height = 8)

state_inc_plot <- reg1_data_clean %>% 
   ggplot() +
   geom_col(aes(dates, case), 
            fill = "steelblue") +
   labs(title = paste0("COVID-19 Cases in Montana Region 1, 2020", 
                       " (Total cases = ", total_cases, ")"),
        subtitle = paste0("Data current as of ", date_today)) +
   ylab("Number of Cases") +
   xlab("") +
   scale_x_date(date_breaks = "3 day", date_labels = "%d-%b") +
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

ggsave("C:/R/covid19/state_daily_results/reg1_inc_plot.png", width = 10, height = 8)



## Region 2 results
# Format analysis data
mt_li_analysis_data <- mt_case_data %>% 
   filter(region == 2) %>% 
   group_by(dates) %>% 
   mutate(local = sum(local),
          imported = sum(imported)) %>% 
   select(dates, local, imported) %>% 
   distinct(dates, .keep_all = TRUE) %>% 
   arrange(dates)


reg2_data <- mt_li_analysis_data 

mt_data <- reg2_data %>% 
   mutate(I = local + imported) %>% 
   select(dates, I) 

mt_incidence_data <- incidence(mt_data$dates, last_date = latest_date)

mt_li_inc_data <- as.data.frame(mt_incidence_data$dates) %>% 
   mutate(dates = ymd(mt_incidence_data$dates)) %>% 
   select(dates) %>% 
   left_join(reg2_data, by = "dates") %>% 
   mutate(local = if_else(is.na(local), 0, local),
          imported = if_else(is.na(imported), 0, imported)) 


state_n <- as.data.frame(mt_li_inc_data$dates)
time_var <- nrow(state_n)
time_start <- seq(2, time_var-13)
time_end <- time_start + 13


MCMC_seed <- 1
overall_seed <- 2
mcmc_control <- make_mcmc_control(seed = MCMC_seed, burnin = 1000, thin = 10)
dist <- "G"  # fitting a Gamma distribution for the SI
si_config <- make_config(list(si_parametric_distr = dist, 
                              mcmc_control = mcmc_control, seed = overall_seed, n1 = 500, n2 = 50,
                              t_start = time_start, t_end = time_end))

mt_r_results <- estimate_R(mt_li_inc_data, method = "si_from_data", 
                           si_data = mt_si_data, config = si_config)



colomns <- c(1:4, 8)
state_r <- (mt_r_results$R[,colomns]) %>% 
   mutate(region = "2") %>% 
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
reg2_r_clean <- cbind(state_r, state_dates_new)

reg2_data_clean <- mt_case_data %>% 
   filter(region == 2)
total_cases <- sum(reg2_data_clean$case)


state_r_plot <- reg2_r_clean %>% 
   ggplot() +
   geom_line(aes(dates, mean_r), size = 1.5, color = "black") +
   geom_line(aes(dates, cl_low), size = 1.5, color = "grey") +
   geom_line(aes(dates, cl_high), size = 1.5, color = "grey") +
   labs(title = "COVID-19 Rolling 14-day R-values, Montana Region 2, 2020",
        color = "") +
   ylab("R-value") +
   xlab("") +
   geom_hline(yintercept = 1, color = "black", size = 1.2) +
   scale_x_date(date_breaks = "3 days", date_labels = "%d-%b") +
   #scale_y_continuous(breaks = seq(0, 10, 0.5), labels = seq(0, 10, 0.5)) +
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
ggsave("C:/R/covid19/state_daily_results/reg2_r_plot.png", width = 10, height = 8)

state_inc_plot <- reg2_data_clean %>% 
   ggplot() +
   geom_col(aes(dates, case), 
            fill = "steelblue") +
   labs(title = paste0("COVID-19 Cases in Montana Region 2, 2020", 
                       " (Total cases = ", total_cases, ")"),
        subtitle = paste0("Data current as of ", date_today)) +
   ylab("Number of Cases") +
   xlab("") +
   scale_x_date(date_breaks = "3 day", date_labels = "%d-%b") +
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

ggsave("C:/R/covid19/state_daily_results/reg2_inc_plot.png", width = 10, height = 8)




## Region 3 results
# Format analysis data
mt_li_analysis_data <- mt_case_data %>% 
   filter(region == 3) %>% 
   group_by(dates) %>% 
   mutate(local = sum(local),
          imported = sum(imported)) %>% 
   select(dates, local, imported) %>% 
   distinct(dates, .keep_all = TRUE) %>% 
   arrange(dates)


reg3_data <- mt_li_analysis_data 

mt_data <- reg3_data %>% 
   mutate(I = local + imported) %>% 
   select(dates, I) 

mt_incidence_data <- incidence(mt_data$dates, last_date = latest_date)

mt_li_inc_data <- as.data.frame(mt_incidence_data$dates) %>% 
   mutate(dates = ymd(mt_incidence_data$dates)) %>% 
   select(dates) %>% 
   left_join(reg3_data, by = "dates") %>% 
   mutate(local = if_else(is.na(local), 0, local),
          imported = if_else(is.na(imported), 0, imported)) 


state_n <- as.data.frame(mt_li_inc_data$dates)
time_var <- nrow(state_n)
time_start <- seq(2, time_var-13)
time_end <- time_start + 13


MCMC_seed <- 1
overall_seed <- 2
mcmc_control <- make_mcmc_control(seed = MCMC_seed, burnin = 1000, thin = 10)
dist <- "G"  # fitting a Gamma distribution for the SI
si_config <- make_config(list(si_parametric_distr = dist, 
                              mcmc_control = mcmc_control, seed = overall_seed, n1 = 500, n2 = 50,
                              t_start = time_start, t_end = time_end))

mt_r_results <- estimate_R(mt_li_inc_data, method = "si_from_data", 
                           si_data = mt_si_data, config = si_config)



colomns <- c(1:4, 8)
state_r <- (mt_r_results$R[,colomns]) %>% 
   mutate(region = "3") %>% 
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
reg3_r_clean <- cbind(state_r, state_dates_new)

reg3_data_clean <- mt_case_data %>% 
   filter(region == 3)
total_cases <- sum(reg3_data_clean$case)


state_r_plot <- reg3_r_clean %>% 
   ggplot() +
   geom_line(aes(dates, mean_r), size = 1.5, color = "black") +
   geom_line(aes(dates, cl_low), size = 1.5, color = "grey") +
   geom_line(aes(dates, cl_high), size = 1.5, color = "grey") +
   labs(title = "COVID-19 Rolling 14-day R-values, Montana Region 3, 2020",
        color = "") +
   ylab("R-value") +
   xlab("") +
   geom_hline(yintercept = 1, color = "black", size = 1.2) +
   scale_x_date(date_breaks = "3 days", date_labels = "%d-%b") +
   #scale_y_continuous(breaks = seq(0, 10, 0.5), labels = seq(0, 10, 0.5)) +
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
ggsave("C:/R/covid19/state_daily_results/reg3_r_plot.png", width = 10, height = 8)

state_inc_plot <- reg3_data_clean %>% 
   ggplot() +
   geom_col(aes(dates, case), 
            fill = "steelblue") +
   labs(title = paste0("COVID-19 Cases in Montana Region 3, 2020", 
                       " (Total cases = ", total_cases, ")"),
        subtitle = paste0("Data current as of ", date_today)) +
   ylab("Number of Cases") +
   xlab("") +
   scale_x_date(date_breaks = "3 day", date_labels = "%d-%b") +
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

ggsave("C:/R/covid19/state_daily_results/reg3_inc_plot.png", width = 10, height = 8)




## Region 4 results
# Format analysis data
mt_li_analysis_data <- mt_case_data %>% 
   filter(region == 4) %>% 
   group_by(dates) %>% 
   mutate(local = sum(local),
          imported = sum(imported)) %>% 
   select(dates, local, imported) %>% 
   distinct(dates, .keep_all = TRUE) %>% 
   arrange(dates)


reg4_data <- mt_li_analysis_data 

mt_data <- reg4_data %>% 
   mutate(I = local + imported) %>% 
   select(dates, I) 

mt_incidence_data <- incidence(mt_data$dates, last_date = latest_date)

mt_li_inc_data <- as.data.frame(mt_incidence_data$dates) %>% 
   mutate(dates = ymd(mt_incidence_data$dates)) %>% 
   select(dates) %>% 
   left_join(reg4_data, by = "dates") %>% 
   mutate(local = if_else(is.na(local), 0, local),
          imported = if_else(is.na(imported), 0, imported)) 


state_n <- as.data.frame(mt_li_inc_data$dates)
time_var <- nrow(state_n)
time_start <- seq(2, time_var-13)
time_end <- time_start + 13


MCMC_seed <- 1
overall_seed <- 2
mcmc_control <- make_mcmc_control(seed = MCMC_seed, burnin = 1000, thin = 10)
dist <- "G"  # fitting a Gamma distribution for the SI
si_config <- make_config(list(si_parametric_distr = dist, 
                              mcmc_control = mcmc_control, seed = overall_seed, n1 = 500, n2 = 50,
                              t_start = time_start, t_end = time_end))

mt_r_results <- estimate_R(mt_li_inc_data, method = "si_from_data", 
                           si_data = mt_si_data, config = si_config)



colomns <- c(1:4, 8)
state_r <- (mt_r_results$R[,colomns]) %>% 
   mutate(region = "4") %>% 
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
reg4_r_clean <- cbind(state_r, state_dates_new)

reg4_data_clean <- mt_case_data %>% 
   filter(region == 4)
total_cases <- sum(reg4_data_clean$case)


state_r_plot <- reg4_r_clean %>% 
   ggplot() +
   geom_line(aes(dates, mean_r), size = 1.5, color = "black") +
   geom_line(aes(dates, cl_low), size = 1.5, color = "grey") +
   geom_line(aes(dates, cl_high), size = 1.5, color = "grey") +
   labs(title = "COVID-19 Rolling 14-day R-values, Montana Region 4, 2020",
        color = "") +
   ylab("R-value") +
   xlab("") +
   geom_hline(yintercept = 1, color = "black", size = 1.2) +
   scale_x_date(date_breaks = "3 days", date_labels = "%d-%b") +
   #scale_y_continuous(breaks = seq(0, 10, 0.5), labels = seq(0, 10, 0.5)) +
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
ggsave("C:/R/covid19/state_daily_results/reg4_r_plot.png", width = 10, height = 8)

state_inc_plot <- reg4_data_clean %>% 
   ggplot() +
   geom_col(aes(dates, case), 
            fill = "steelblue") +
   labs(title = paste0("COVID-19 Cases in Montana Region 4, 2020", 
                       " (Total cases = ", total_cases, ")"),
        subtitle = paste0("Data current as of ", date_today)) +
   ylab("Number of Cases") +
   xlab("") +
   scale_x_date(date_breaks = "3 day", date_labels = "%d-%b") +
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

ggsave("C:/R/covid19/state_daily_results/reg4_inc_plot.png", width = 10, height = 8)




## Region 5 results
# Format analysis data
mt_li_analysis_data <- mt_case_data %>% 
   filter(region == 5) %>% 
   group_by(dates) %>% 
   mutate(local = sum(local),
          imported = sum(imported)) %>% 
   select(dates, local, imported) %>% 
   distinct(dates, .keep_all = TRUE) %>% 
   arrange(dates)


reg5_data <- mt_li_analysis_data 

mt_data <- reg5_data %>% 
   mutate(I = local + imported) %>% 
   select(dates, I) 

mt_incidence_data <- incidence(mt_data$dates, last_date = latest_date)

mt_li_inc_data <- as.data.frame(mt_incidence_data$dates) %>% 
   mutate(dates = ymd(mt_incidence_data$dates)) %>% 
   select(dates) %>% 
   left_join(reg5_data, by = "dates") %>% 
   mutate(local = if_else(is.na(local), 0, local),
          imported = if_else(is.na(imported), 0, imported)) 


state_n <- as.data.frame(mt_li_inc_data$dates)
time_var <- nrow(state_n)
time_start <- seq(2, time_var-13)
time_end <- time_start + 13


MCMC_seed <- 1
overall_seed <- 2
mcmc_control <- make_mcmc_control(seed = MCMC_seed, burnin = 1000, thin = 10)
dist <- "G"  # fitting a Gamma distribution for the SI
si_config <- make_config(list(si_parametric_distr = dist, 
                              mcmc_control = mcmc_control, seed = overall_seed, n1 = 500, n2 = 50,
                              t_start = time_start, t_end = time_end))

mt_r_results <- estimate_R(mt_li_inc_data, method = "si_from_data", 
                           si_data = mt_si_data, config = si_config)



colomns <- c(1:4, 8)
state_r <- (mt_r_results$R[,colomns]) %>% 
   mutate(region = "5") %>% 
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
reg5_r_clean <- cbind(state_r, state_dates_new)

reg5_data_clean <- mt_case_data %>% 
   filter(region == 5)
total_cases <- sum(reg5_data_clean$case)


state_r_plot <- reg5_r_clean %>% 
   ggplot() +
   geom_line(aes(dates, mean_r), size = 1.5, color = "black") +
   geom_line(aes(dates, cl_low), size = 1.5, color = "grey") +
   geom_line(aes(dates, cl_high), size = 1.5, color = "grey") +
   labs(title = "COVID-19 Rolling 14-day R-values, Montana Region 5, 2020",
        color = "") +
   ylab("R-value") +
   xlab("") +
   geom_hline(yintercept = 1, color = "black", size = 1.2) +
   scale_x_date(date_breaks = "3 days", date_labels = "%d-%b") +
   #scale_y_continuous(breaks = seq(0, 10, 0.5), labels = seq(0, 10, 0.5)) +
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
ggsave("C:/R/covid19/state_daily_results/reg5_r_plot.png", width = 10, height = 8)

state_inc_plot <- reg5_data_clean %>% 
   ggplot() +
   geom_col(aes(dates, case), 
            fill = "steelblue") +
   labs(title = paste0("COVID-19 Cases in Montana Region 5, 2020", 
                       " (Total cases = ", total_cases, ")"),
        subtitle = paste0("Data current as of ", date_today)) +
   ylab("Number of Cases") +
   xlab("") +
   scale_x_date(date_breaks = "3 day", date_labels = "%d-%b") +
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

