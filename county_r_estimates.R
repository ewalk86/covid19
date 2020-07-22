# Estimate R, save results/plots for specific counties
# Ethan Walker
# 16 July 2020

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
mt_si_data <- read_xlsx(paste0(file_path, "Input/SI_Local_v_Import Data_07.22.2020.xlsx"),
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
mt_case_data <- read_xlsx(paste0(file_path, "Input/SI_Local_v_Import Data_07.22.2020.xlsx"),
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


#################### Run analysis and print results
## Missoula county results
# Format analysis data
missoula_li_analysis_data <- mt_case_data %>% 
   filter(county == "Missoula") %>% 
   group_by(dates) %>% 
   mutate(local = sum(local),
          imported = sum(imported)) %>% 
   select(dates, local, imported) %>% 
   distinct(dates, .keep_all = TRUE) %>% 
   arrange(dates)


missoula_data <- missoula_li_analysis_data 

inc_data <- missoula_data %>% 
   mutate(I = local + imported) %>% 
   select(dates, I) 

latest_date <- format(Sys.Date() - 14, "%Y-%m-%d")

missoula_inc_data <- incidence(inc_data$dates, last_date = latest_date)

missoula_li_inc_data <- as.data.frame(missoula_inc_data$dates) %>% 
   mutate(dates = ymd(missoula_inc_data$dates)) %>% 
   select(dates) %>% 
   left_join(missoula_data, by = "dates") %>% 
   mutate(local = if_else(is.na(local), 0, local),
          imported = if_else(is.na(imported), 0, imported)) 


missoula_n <- as.data.frame(missoula_li_inc_data$dates)
time_var <- nrow(missoula_n)
time_start <- seq(2, time_var-13)
time_end <- time_start + 13


MCMC_seed <- 1
overall_seed <- 2
mcmc_control <- make_mcmc_control(seed = MCMC_seed, burnin = 1000, thin = 10)
dist <- "G"  # fitting a Gamma distribution for the SI
si_config <- make_config(list(si_parametric_distr = dist, 
                              mcmc_control = mcmc_control, seed = overall_seed, n1 = 500, n2 = 50,
                              t_start = time_start, t_end = time_end))

missoula_r_results <- estimate_R(missoula_li_inc_data, method = "si_from_data", 
                           si_data = mt_si_data, config = si_config)

# calculate mean serial interval using SI data
mean_si <- missoula_r_results$SI.Moments %>% 
   summarize(mean_si = mean(Mean), sd_si = mean(Std)) 



colomns <- c(1:4, 8)
missoula_r <- (missoula_r_results$R[,colomns]) %>% 
   mutate(region = "Missoula") %>% 
   rename(mean_r = `Mean(R)`,
          sd_r = `Std(R)`,
          median_r = `Median(R)`)

missoula_dates <- as.data.frame(missoula_r_results$dates)
missoula_i <- as.data.frame(missoula_r_results$I)
missoula_cil <- as.data.frame(missoula_r_results$R$`Quantile.0.025(R)`)
missoula_cih <- as.data.frame(missoula_r_results$R$`Quantile.0.975(R)`)
missoula_dates_new <- cbind(missoula_dates, missoula_i) %>% 
   rename(dates = 1,
          incidence = 2) %>% 
   mutate(dates = ymd(dates))
missoula_dates_new <- missoula_dates_new[-(1:14), 1:2]
missoula_dates_new <- cbind(missoula_dates_new, missoula_cil, missoula_cih) %>% 
   rename(cl_low = 3,
          cl_high = 4)
missoula_r_clean <- cbind(missoula_r, missoula_dates_new)

missoula_data_clean <- mt_case_data %>% 
   filter(county == "Missoula")
total_cases <- sum(missoula_data_clean$case)


missoula_r_plot <- missoula_r_clean %>% 
   filter(dates > Sys.Date() - 35) %>% 
   ggplot() +
   geom_line(aes(dates, mean_r), size = 1.5, color = "black") +
   geom_line(aes(dates, cl_low), size = 1.5, color = "grey") +
   geom_line(aes(dates, cl_high), size = 1.5, color = "grey") +
   labs(title = "COVID-19 Rolling 14-day R-values, Missoula County, 2020",
        color = "") +
   ylab("R-value") +
   xlab("") +
   geom_hline(yintercept = 1, color = "red", size = 1.2) +
   scale_x_date(date_breaks = "1 days", date_labels = "%d-%b") +
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
#missoula_r_plot
ggsave("C:/R/covid19/missoula_r_plot.png", width = 10, height = 8)

missoula_inc_plot <- missoula_data_clean %>% 
   ggplot() +
   geom_col(aes(dates, case), 
            fill = "steelblue") +
   labs(title = paste0("COVID-19 Cases in Missoula County, 2020", 
                       " (Total cases = ", total_cases, ")")) +
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
#missoula_inc_plot

ggsave("C:/R/covid19/missoula_inc_plot.png", width = 10, height = 8)



## Gallatin county results
# Format analysis data
gallatin_li_analysis_data <- mt_case_data %>% 
   filter(county == "Gallatin") %>% 
   group_by(dates) %>% 
   mutate(local = sum(local),
          imported = sum(imported)) %>% 
   select(dates, local, imported) %>% 
   distinct(dates, .keep_all = TRUE) %>% 
   arrange(dates)


gallatin_data <- gallatin_li_analysis_data 

inc_data <- gallatin_data %>% 
   mutate(I = local + imported) %>% 
   select(dates, I) 

latest_date <- format(Sys.Date() - 14, "%Y-%m-%d")

gallatin_inc_data <- incidence(inc_data$dates, last_date = latest_date)

gallatin_li_inc_data <- as.data.frame(gallatin_inc_data$dates) %>% 
   mutate(dates = ymd(gallatin_inc_data$dates)) %>% 
   select(dates) %>% 
   left_join(gallatin_data, by = "dates") %>% 
   mutate(local = if_else(is.na(local), 0, local),
          imported = if_else(is.na(imported), 0, imported)) 


gallatin_n <- as.data.frame(gallatin_li_inc_data$dates)
time_var <- nrow(gallatin_n)
time_start <- seq(2, time_var-13)
time_end <- time_start + 13


MCMC_seed <- 1
overall_seed <- 2
mcmc_control <- make_mcmc_control(seed = MCMC_seed, burnin = 1000, thin = 10)
dist <- "G"  # fitting a Gamma distribution for the SI
si_config <- make_config(list(si_parametric_distr = dist, 
                              mcmc_control = mcmc_control, seed = overall_seed, n1 = 500, n2 = 50,
                              t_start = time_start, t_end = time_end))

gallatin_r_results <- estimate_R(gallatin_li_inc_data, method = "si_from_data", 
                                 si_data = mt_si_data, config = si_config)



colomns <- c(1:4, 8)
gallatin_r <- (gallatin_r_results$R[,colomns]) %>% 
   mutate(region = "Gallatin") %>% 
   rename(mean_r = `Mean(R)`,
          sd_r = `Std(R)`,
          median_r = `Median(R)`)

gallatin_dates <- as.data.frame(gallatin_r_results$dates)
gallatin_i <- as.data.frame(gallatin_r_results$I)
gallatin_cil <- as.data.frame(gallatin_r_results$R$`Quantile.0.025(R)`)
gallatin_cih <- as.data.frame(gallatin_r_results$R$`Quantile.0.975(R)`)
gallatin_dates_new <- cbind(gallatin_dates, gallatin_i) %>% 
   rename(dates = 1,
          incidence = 2) %>% 
   mutate(dates = ymd(dates))
gallatin_dates_new <- gallatin_dates_new[-(1:14), 1:2]
gallatin_dates_new <- cbind(gallatin_dates_new, gallatin_cil, gallatin_cih) %>% 
   rename(cl_low = 3,
          cl_high = 4)
gallatin_r_clean <- cbind(gallatin_r, gallatin_dates_new)

gallatin_data_clean <- mt_case_data %>% 
   filter(county == "Gallatin")
total_cases <- sum(gallatin_data_clean$case)


gallatin_r_plot <- gallatin_r_clean %>% 
   filter(dates > Sys.Date() - 35) %>% 
   ggplot() +
   geom_line(aes(dates, mean_r), size = 1.5, color = "black") +
   geom_line(aes(dates, cl_low), size = 1.5, color = "grey") +
   geom_line(aes(dates, cl_high), size = 1.5, color = "grey") +
   labs(title = "COVID-19 Rolling 14-day R-values, Gallatin County, 2020",
        color = "") +
   ylab("R-value") +
   xlab("") +
   geom_hline(yintercept = 1, color = "red", size = 1.2) +
   scale_x_date(date_breaks = "1 days", date_labels = "%d-%b") +
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
#gallatin_r_plot
ggsave("C:/R/covid19/gallatin_r_plot.png", width = 10, height = 8)

gallatin_inc_plot <- gallatin_data_clean %>% 
   ggplot() +
   geom_col(aes(dates, case), 
            fill = "steelblue") +
   labs(title = paste0("COVID-19 Cases in Gallatin County, 2020", 
                       " (Total cases = ", total_cases, ")")) +
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
#gallatin_inc_plot

ggsave("C:/R/covid19/gallatin_inc_plot.png", width = 10, height = 8)




## Yellowstone county results
# Format analysis data
yellowstone_li_analysis_data <- mt_case_data %>% 
   filter(county == "Yellowstone") %>% 
   group_by(dates) %>% 
   mutate(local = sum(local),
          imported = sum(imported)) %>% 
   select(dates, local, imported) %>% 
   distinct(dates, .keep_all = TRUE) %>% 
   arrange(dates)


yellowstone_data <- yellowstone_li_analysis_data 

inc_data <- yellowstone_data %>% 
   mutate(I = local + imported) %>% 
   select(dates, I) 

latest_date <- format(Sys.Date() - 14, "%Y-%m-%d")

yellowstone_inc_data <- incidence(inc_data$dates, last_date = latest_date)

yellowstone_li_inc_data <- as.data.frame(yellowstone_inc_data$dates) %>% 
   mutate(dates = ymd(yellowstone_inc_data$dates)) %>% 
   select(dates) %>% 
   left_join(yellowstone_data, by = "dates") %>% 
   mutate(local = if_else(is.na(local), 0, local),
          imported = if_else(is.na(imported), 0, imported)) 


yellowstone_n <- as.data.frame(yellowstone_li_inc_data$dates)
time_var <- nrow(yellowstone_n)
time_start <- seq(2, time_var-13)
time_end <- time_start + 13


MCMC_seed <- 1
overall_seed <- 2
mcmc_control <- make_mcmc_control(seed = MCMC_seed, burnin = 1000, thin = 10)
dist <- "G"  # fitting a Gamma distribution for the SI
si_config <- make_config(list(si_parametric_distr = dist, 
                              mcmc_control = mcmc_control, seed = overall_seed, n1 = 500, n2 = 50,
                              t_start = time_start, t_end = time_end))

yellowstone_r_results <- estimate_R(yellowstone_li_inc_data, method = "si_from_data", 
                                 si_data = mt_si_data, config = si_config)



colomns <- c(1:4, 8)
yellowstone_r <- (yellowstone_r_results$R[,colomns]) %>% 
   mutate(region = "Yellowstone") %>% 
   rename(mean_r = `Mean(R)`,
          sd_r = `Std(R)`,
          median_r = `Median(R)`)

yellowstone_dates <- as.data.frame(yellowstone_r_results$dates)
yellowstone_i <- as.data.frame(yellowstone_r_results$I)
yellowstone_cil <- as.data.frame(yellowstone_r_results$R$`Quantile.0.025(R)`)
yellowstone_cih <- as.data.frame(yellowstone_r_results$R$`Quantile.0.975(R)`)
yellowstone_dates_new <- cbind(yellowstone_dates, yellowstone_i) %>% 
   rename(dates = 1,
          incidence = 2) %>% 
   mutate(dates = ymd(dates))
yellowstone_dates_new <- yellowstone_dates_new[-(1:14), 1:2]
yellowstone_dates_new <- cbind(yellowstone_dates_new, yellowstone_cil, yellowstone_cih) %>% 
   rename(cl_low = 3,
          cl_high = 4)
yellowstone_r_clean <- cbind(yellowstone_r, yellowstone_dates_new)

yellowstone_data_clean <- mt_case_data %>% 
   filter(county == "Yellowstone")
total_cases <- sum(yellowstone_data_clean$case)


yellowstone_r_plot <- yellowstone_r_clean %>% 
   filter(dates > Sys.Date() - 35) %>% 
   ggplot() +
   geom_line(aes(dates, mean_r), size = 1.5, color = "black") +
   geom_line(aes(dates, cl_low), size = 1.5, color = "grey") +
   geom_line(aes(dates, cl_high), size = 1.5, color = "grey") +
   labs(title = "COVID-19 Rolling 14-day R-values, Yellowstone County, 2020",
        color = "") +
   ylab("R-value") +
   xlab("") +
   geom_hline(yintercept = 1, color = "red", size = 1.2) +
   scale_x_date(date_breaks = "1 days", date_labels = "%d-%b") +
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
#yellowstone_r_plot
ggsave("C:/R/covid19/yellowstone_r_plot.png", width = 10, height = 8)

yellowstone_inc_plot <- yellowstone_data_clean %>% 
   ggplot() +
   geom_col(aes(dates, case), 
            fill = "steelblue") +
   labs(title = paste0("COVID-19 Cases in Yellowstone County, 2020", 
                       " (Total cases = ", total_cases, ")")) +
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
#yellowstone_inc_plot

ggsave("C:/R/covid19/yellowstone_inc_plot.png", width = 10, height = 8)




## Bighorn county results
# Format analysis data
bighorn_li_analysis_data <- mt_case_data %>% 
   filter(county == "Big Horn") %>% 
   group_by(dates) %>% 
   mutate(local = sum(local),
          imported = sum(imported)) %>% 
   select(dates, local, imported) %>% 
   distinct(dates, .keep_all = TRUE) %>% 
   arrange(dates)


bighorn_data <- bighorn_li_analysis_data 

inc_data <- bighorn_data %>% 
   mutate(I = local + imported) %>% 
   select(dates, I) 

latest_date <- format(Sys.Date() - 14, "%Y-%m-%d")

bighorn_inc_data <- incidence(inc_data$dates, last_date = latest_date)

bighorn_li_inc_data <- as.data.frame(bighorn_inc_data$dates) %>% 
   mutate(dates = ymd(bighorn_inc_data$dates)) %>% 
   select(dates) %>% 
   left_join(bighorn_data, by = "dates") %>% 
   mutate(local = if_else(is.na(local), 0, local),
          imported = if_else(is.na(imported), 0, imported)) 


bighorn_n <- as.data.frame(bighorn_li_inc_data$dates)
time_var <- nrow(bighorn_n)
time_start <- seq(2, time_var-13)
time_end <- time_start + 13


MCMC_seed <- 1
overall_seed <- 2
mcmc_control <- make_mcmc_control(seed = MCMC_seed, burnin = 1000, thin = 10)
dist <- "G"  # fitting a Gamma distribution for the SI
si_config <- make_config(list(si_parametric_distr = dist, 
                              mcmc_control = mcmc_control, seed = overall_seed, n1 = 500, n2 = 50,
                              t_start = time_start, t_end = time_end))

bighorn_r_results <- estimate_R(bighorn_li_inc_data, method = "si_from_data", 
                                 si_data = mt_si_data, config = si_config)



colomns <- c(1:4, 8)
bighorn_r <- (bighorn_r_results$R[,colomns]) %>% 
   mutate(region = "Big Horn") %>% 
   rename(mean_r = `Mean(R)`,
          sd_r = `Std(R)`,
          median_r = `Median(R)`)

bighorn_dates <- as.data.frame(bighorn_r_results$dates)
bighorn_i <- as.data.frame(bighorn_r_results$I)
bighorn_cil <- as.data.frame(bighorn_r_results$R$`Quantile.0.025(R)`)
bighorn_cih <- as.data.frame(bighorn_r_results$R$`Quantile.0.975(R)`)
bighorn_dates_new <- cbind(bighorn_dates, bighorn_i) %>% 
   rename(dates = 1,
          incidence = 2) %>% 
   mutate(dates = ymd(dates))
bighorn_dates_new <- bighorn_dates_new[-(1:14), 1:2]
bighorn_dates_new <- cbind(bighorn_dates_new, bighorn_cil, bighorn_cih) %>% 
   rename(cl_low = 3,
          cl_high = 4)
bighorn_r_clean <- cbind(bighorn_r, bighorn_dates_new)

bighorn_data_clean <- mt_case_data %>% 
   filter(county == "Big Horn")
total_cases <- sum(bighorn_data_clean$case)


bighorn_r_plot <- bighorn_r_clean %>% 
   filter(dates > Sys.Date() - 35) %>% 
   ggplot() +
   geom_line(aes(dates, mean_r), size = 1.5, color = "black") +
   geom_line(aes(dates, cl_low), size = 1.5, color = "grey") +
   geom_line(aes(dates, cl_high), size = 1.5, color = "grey") +
   labs(title = "COVID-19 Rolling 14-day R-values, Big Horn County, 2020",
        color = "") +
   ylab("R-value") +
   xlab("") +
   geom_hline(yintercept = 1, color = "red", size = 1.2) +
   scale_x_date(date_breaks = "1 days", date_labels = "%d-%b") +
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
#bighorn_r_plot
ggsave("C:/R/covid19/bighorn_r_plot.png", width = 10, height = 8)

bighorn_inc_plot <- bighorn_data_clean %>% 
   ggplot() +
   geom_col(aes(dates, case), 
            fill = "steelblue") +
   labs(title = paste0("COVID-19 Cases in Big Horn County, 2020", 
                       " (Total cases = ", total_cases, ")")) +
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
#bighorn_inc_plot

ggsave("C:/R/covid19/bighorn_inc_plot.png", width = 10, height = 8)




## Lake county results
# Format analysis data
lake_li_analysis_data <- mt_case_data %>% 
   filter(county == "Lake") %>% 
   group_by(dates) %>% 
   mutate(local = sum(local),
          imported = sum(imported)) %>% 
   select(dates, local, imported) %>% 
   distinct(dates, .keep_all = TRUE) %>% 
   arrange(dates)


lake_data <- lake_li_analysis_data 

inc_data <- lake_data %>% 
   mutate(I = local + imported) %>% 
   select(dates, I) 

latest_date <- format(Sys.Date() - 14, "%Y-%m-%d")

lake_inc_data <- incidence(inc_data$dates, last_date = latest_date)

lake_li_inc_data <- as.data.frame(lake_inc_data$dates) %>% 
   mutate(dates = ymd(lake_inc_data$dates)) %>% 
   select(dates) %>% 
   left_join(lake_data, by = "dates") %>% 
   mutate(local = if_else(is.na(local), 0, local),
          imported = if_else(is.na(imported), 0, imported)) 


lake_n <- as.data.frame(lake_li_inc_data$dates)
time_var <- nrow(lake_n)
time_start <- seq(2, time_var-13)
time_end <- time_start + 13


MCMC_seed <- 1
overall_seed <- 2
mcmc_control <- make_mcmc_control(seed = MCMC_seed, burnin = 1000, thin = 10)
dist <- "G"  # fitting a Gamma distribution for the SI
si_config <- make_config(list(si_parametric_distr = dist, 
                              mcmc_control = mcmc_control, seed = overall_seed, n1 = 500, n2 = 50,
                              t_start = time_start, t_end = time_end))

lake_r_results <- estimate_R(lake_li_inc_data, method = "si_from_data", 
                                si_data = mt_si_data, config = si_config)



colomns <- c(1:4, 8)
lake_r <- (lake_r_results$R[,colomns]) %>% 
   mutate(region = "Lake") %>% 
   rename(mean_r = `Mean(R)`,
          sd_r = `Std(R)`,
          median_r = `Median(R)`)

lake_dates <- as.data.frame(lake_r_results$dates)
lake_i <- as.data.frame(lake_r_results$I)
lake_cil <- as.data.frame(lake_r_results$R$`Quantile.0.025(R)`)
lake_cih <- as.data.frame(lake_r_results$R$`Quantile.0.975(R)`)
lake_dates_new <- cbind(lake_dates, lake_i) %>% 
   rename(dates = 1,
          incidence = 2) %>% 
   mutate(dates = ymd(dates))
lake_dates_new <- lake_dates_new[-(1:14), 1:2]
lake_dates_new <- cbind(lake_dates_new, lake_cil, lake_cih) %>% 
   rename(cl_low = 3,
          cl_high = 4)
lake_r_clean <- cbind(lake_r, lake_dates_new)

lake_data_clean <- mt_case_data %>% 
   filter(county == "Lake")
total_cases <- sum(lake_data_clean$case)


lake_r_plot <- lake_r_clean %>% 
   filter(dates > Sys.Date() - 35) %>% 
   ggplot() +
   geom_line(aes(dates, mean_r), size = 1.5, color = "black") +
   geom_line(aes(dates, cl_low), size = 1.5, color = "grey") +
   geom_line(aes(dates, cl_high), size = 1.5, color = "grey") +
   labs(title = "COVID-19 Rolling 14-day R-values, Lake County, 2020",
        color = "") +
   ylab("R-value") +
   xlab("") +
   geom_hline(yintercept = 1, color = "red", size = 1.2) +
   scale_x_date(date_breaks = "1 days", date_labels = "%d-%b") +
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
#lake_r_plot
ggsave("C:/R/covid19/lake_r_plot.png", width = 10, height = 8)

lake_inc_plot <- lake_data_clean %>% 
   ggplot() +
   geom_col(aes(dates, case), 
            fill = "steelblue") +
   labs(title = paste0("COVID-19 Cases in Lake County, 2020", 
                       " (Total cases = ", total_cases, ")")) +
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
#lake_inc_plot

ggsave("C:/R/covid19/lake_inc_plot.png", width = 10, height = 8)




## Lewis and Clark county results
# Format analysis data
lewisandclark_li_analysis_data <- mt_case_data %>% 
   filter(county == "Lewis and Clark") %>% 
   group_by(dates) %>% 
   mutate(local = sum(local),
          imported = sum(imported)) %>% 
   select(dates, local, imported) %>% 
   distinct(dates, .keep_all = TRUE) %>% 
   arrange(dates)


lewisandclark_data <- lewisandclark_li_analysis_data 

inc_data <- lewisandclark_data %>% 
   mutate(I = local + imported) %>% 
   select(dates, I) 

latest_date <- format(Sys.Date() - 14, "%Y-%m-%d")

lewisandclark_inc_data <- incidence(inc_data$dates, last_date = latest_date)

lewisandclark_li_inc_data <- as.data.frame(lewisandclark_inc_data$dates) %>% 
   mutate(dates = ymd(lewisandclark_inc_data$dates)) %>% 
   select(dates) %>% 
   left_join(lewisandclark_data, by = "dates") %>% 
   mutate(local = if_else(is.na(local), 0, local),
          imported = if_else(is.na(imported), 0, imported)) 


lewisandclark_n <- as.data.frame(lewisandclark_li_inc_data$dates)
time_var <- nrow(lewisandclark_n)
time_start <- seq(2, time_var-13)
time_end <- time_start + 13


MCMC_seed <- 1
overall_seed <- 2
mcmc_control <- make_mcmc_control(seed = MCMC_seed, burnin = 1000, thin = 10)
dist <- "G"  # fitting a Gamma distribution for the SI
si_config <- make_config(list(si_parametric_distr = dist, 
                              mcmc_control = mcmc_control, seed = overall_seed, n1 = 500, n2 = 50,
                              t_start = time_start, t_end = time_end))

lewisandclark_r_results <- estimate_R(lewisandclark_li_inc_data, method = "si_from_data", 
                                si_data = mt_si_data, config = si_config)



colomns <- c(1:4, 8)
lewisandclark_r <- (lewisandclark_r_results$R[,colomns]) %>% 
   mutate(region = "Lewis and Clark") %>% 
   rename(mean_r = `Mean(R)`,
          sd_r = `Std(R)`,
          median_r = `Median(R)`)

lewisandclark_dates <- as.data.frame(lewisandclark_r_results$dates)
lewisandclark_i <- as.data.frame(lewisandclark_r_results$I)
lewisandclark_cil <- as.data.frame(lewisandclark_r_results$R$`Quantile.0.025(R)`)
lewisandclark_cih <- as.data.frame(lewisandclark_r_results$R$`Quantile.0.975(R)`)
lewisandclark_dates_new <- cbind(lewisandclark_dates, lewisandclark_i) %>% 
   rename(dates = 1,
          incidence = 2) %>% 
   mutate(dates = ymd(dates))
lewisandclark_dates_new <- lewisandclark_dates_new[-(1:14), 1:2]
lewisandclark_dates_new <- cbind(lewisandclark_dates_new, lewisandclark_cil, lewisandclark_cih) %>% 
   rename(cl_low = 3,
          cl_high = 4)
lewisandclark_r_clean <- cbind(lewisandclark_r, lewisandclark_dates_new)

lewisandclark_data_clean <- mt_case_data %>% 
   filter(county == "Lewis and Clark")
total_cases <- sum(lewisandclark_data_clean$case)


lewisandclark_r_plot <- lewisandclark_r_clean %>% 
   filter(dates > Sys.Date() - 35) %>% 
   ggplot() +
   geom_line(aes(dates, mean_r), size = 1.5, color = "black") +
   geom_line(aes(dates, cl_low), size = 1.5, color = "grey") +
   geom_line(aes(dates, cl_high), size = 1.5, color = "grey") +
   labs(title = "COVID-19 Rolling 14-day R-values, Lewis and Clark County, 2020",
        color = "") +
   ylab("R-value") +
   xlab("") +
   geom_hline(yintercept = 1, color = "red", size = 1.2) +
   scale_x_date(date_breaks = "1 days", date_labels = "%d-%b") +
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
#lewisandclark_r_plot
ggsave("C:/R/covid19/lewisandclark_r_plot.png", width = 10, height = 8)

lewisandclark_inc_plot <- lewisandclark_data_clean %>% 
   ggplot() +
   geom_col(aes(dates, case), 
            fill = "steelblue") +
   labs(title = paste0("COVID-19 Cases in Lewis and Clark County, 2020", 
                       " (Total cases = ", total_cases, ")")) +
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
#lewisandclark_inc_plot

ggsave("C:/R/covid19/lewisandclark_inc_plot.png", width = 10, height = 8)





# Bind files and save
all_counties_r <- rbind(missoula_r_clean, gallatin_r_clean, yellowstone_r_clean, 
                        bighorn_r_clean, lake_r_clean, lewisandclark_r_clean) %>% 
   mutate(daily_cases = incidence) 

write_csv(all_counties_r, "C:/R/covid19/all_counties_r.csv", na = " ")

# Send files to the sftp server
# host = 'elbastion.dbs.umt.edu'
# port = 22
# username = 'celftp'
# password  = 'celftp'
# remotepath = '/celFtpFiles/covid19/Rt/incoming/'

sftpUpload("elbastion.dbs.umt.edu", "celftp", "celftp",
           "/celFtpFiles/covid19/Rt/incoming/all_counties_r.csv",
           "C:/R/covid19/all_counties_r.csv")


# Combine county R file and region/state R file; push to server

all_counties_r <- read_csv("C:/R/covid19/all_counties_r.csv") 

all_regions_r <- read_csv("C:/R/covid19/state_daily_results/all_regions_r.csv")

all_regions_cos_r <- plyr::rbind.fill(all_counties_r, all_regions_r)


write_csv(all_regions_cos_r, "C:/R/covid19/all_regions_cos_r.csv", na = " ")


sftpUpload("elbastion.dbs.umt.edu", "celftp", "celftp",
           "/celFtpFiles/covid19/Rt/incoming/all_regions_r.csv",
           "C:/R/covid19/all_regions_cos_r.csv")
