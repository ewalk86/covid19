# Pull daily updated State of Montana covid data, estimate R, save results/plots
# Ethan Walker
# 21 April 2020

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
   left_join(counties_regions, by = "county") %>% 
   mutate(case = 1)

state_data_wide <- state_data_clean %>% 
   mutate(age_group = factor(age_group,
                             labels = c("<10 yrs", "10 to 19 yrs",
                                        "20 to 29 yrs", "30 to 39 yrs",
                                        "40 to 49 yrs", "50 to 59 yrs",
                                        "60 to 69 yrs", "70 to 79 yrs",
                                        "80 to 89 yrs", "90 to 99 yrs")),
          sex = factor(sex, labels = c("Female", "Male")),
          hospitalization = factor(hospitalization,
                                   labels = c("Hosp: No", "Hosp: Possible",
                                              "Hosp: Unknown", "Hosp: Yes"))) %>% 
   mutate(case = 1) %>% 
   pivot_wider(names_from = "age_group", values_from = "case") %>% 
   mutate(case = 1) %>% 
   pivot_wider(names_from = "sex", values_from = "case") %>% 
   mutate(case = 1) %>% 
   pivot_wider(names_from = "hospitalization", values_from = "case") %>% 
   mutate(case = 1) %>% 
   pivot_wider(names_from = "county", values_from = "case")

state_wide_date <- state_data_wide %>% 
   mutate(region = "state") %>% 
   group_by(dates) %>% 
   mutate_at(c(6:50), sum, na.rm = TRUE) %>% 
   mutate(case = 1,
          daily_cases = sum(case)) %>% 
   arrange(dates) %>% 
   distinct(dates, .keep_all = TRUE) %>% 
   ungroup() %>% 
   mutate(cumulative_cases = cumsum(daily_cases)) %>% 
   select(region, dates, daily_cases, cumulative_cases, `70 to 79 yrs`:Pondera)

reg1_wide_date <- state_data_wide %>% 
   filter(region == 1) %>% 
   group_by(dates) %>% 
   mutate_at(c(6:50), sum, na.rm = TRUE) %>% 
   mutate(case = 1,
          daily_cases = sum(case)) %>% 
   arrange(dates) %>% 
   distinct(dates, .keep_all = TRUE) %>% 
   ungroup() %>% 
   mutate(cumulative_cases = cumsum(daily_cases)) %>% 
   select(region, dates, daily_cases, cumulative_cases, `70 to 79 yrs`:Pondera)

reg2_wide_date <- state_data_wide %>% 
   filter(region == 2) %>% 
   group_by(dates) %>% 
   mutate_at(c(6:50), sum, na.rm = TRUE) %>% 
   mutate(case = 1,
          daily_cases = sum(case)) %>% 
   arrange(dates) %>% 
   distinct(dates, .keep_all = TRUE) %>% 
   ungroup() %>% 
   mutate(cumulative_cases = cumsum(daily_cases)) %>% 
   select(region, dates, daily_cases, cumulative_cases, `70 to 79 yrs`:Pondera)

reg3_wide_date <- state_data_wide %>% 
   filter(region == 3) %>% 
   group_by(dates) %>% 
   mutate_at(c(6:50), sum, na.rm = TRUE) %>% 
   mutate(case = 1,
          daily_cases = sum(case)) %>% 
   arrange(dates) %>% 
   distinct(dates, .keep_all = TRUE) %>% 
   ungroup() %>% 
   mutate(cumulative_cases = cumsum(daily_cases)) %>% 
   select(region, dates, daily_cases, cumulative_cases, `70 to 79 yrs`:Pondera)

reg4_wide_date <- state_data_wide %>% 
   filter(region == 4) %>% 
   group_by(dates) %>% 
   mutate_at(c(6:50), sum, na.rm = TRUE) %>% 
   mutate(case = 1,
          daily_cases = sum(case)) %>% 
   arrange(dates) %>% 
   distinct(dates, .keep_all = TRUE) %>% 
   ungroup() %>% 
   mutate(cumulative_cases = cumsum(daily_cases)) %>% 
   select(region, dates, daily_cases, cumulative_cases, `70 to 79 yrs`:Pondera)

reg5_wide_date <- state_data_wide %>% 
   filter(region == 5) %>% 
   group_by(dates) %>% 
   mutate_at(c(6:50), sum, na.rm = TRUE) %>% 
   mutate(case = 1,
          daily_cases = sum(case)) %>% 
   arrange(dates) %>% 
   distinct(dates, .keep_all = TRUE) %>% 
   ungroup() %>% 
   mutate(cumulative_cases = cumsum(daily_cases)) %>% 
   select(region, dates, daily_cases, cumulative_cases, `70 to 79 yrs`:Pondera)

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


# Run analysis and print results
serial_interval_mean <- 4
serial_interval_sd <- 4.75


## State results
state_incidence <- incidence(state_data_clean$dates)

state_results <- estimate_R(state_incidence, method="parametric_si", 
                            config = make_config(list(mean_si = serial_interval_mean, 
                                                      std_si = serial_interval_sd)))

#state_plot <- estimate_R_plots(state_results, what = "R")

colomns <- c(1:4, 8)
state_r <- (state_results$R[,colomns]) %>% 
   mutate(region = "state") %>% 
   rename(mean_r = `Mean(R)`,
          sd_r = `Std(R)`,
          median_r = `Median(R)`)

state_dates <- as.data.frame(state_results$dates)
state_i <- as.data.frame(state_results$I)
state_dates_new <- cbind(state_dates, state_i) %>% 
   rename(dates = 1,
          incidence = 2) %>% 
   mutate(dates = ymd(dates))
state_dates_new <- state_dates_new[-(1:7), 1:2]
state_r_clean <- cbind(state_r, state_dates_new) 


total_cases <- sum(state_data_clean$case)
date_today <- format(Sys.Date(), "%d %b %Y")
date_10 <- format(Sys.Date() - 10, "%d %b %Y")


state_r_plot <- state_r_clean %>% 
   ggplot() +
   geom_line(aes(dates, mean_r), size = 1.5, color = "black") +
   labs(title = "COVID-19 Rolling 7-day R-values, State of Montana, 2020",
        subtitle = paste0("Data current as of ", date_today),
        color = "") +
   ylab("R-value") +
   xlab("") +
   geom_hline(yintercept = 1, color = "black", size = 1.2) +
   scale_x_date(date_breaks = "2 days", date_labels = "%d-%b") +
   scale_y_continuous(breaks = seq(0, 5, 0.25), labels = seq(0, 5, 0.25)) +
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
ggsave("C:/R/covid19/state_daily_results/state_r_plot.png")

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

ggsave("C:/R/covid19/state_daily_results/state_inc_plot.png")



## Reg1 results
reg1_data <- state_data_clean %>% 
   filter(region == 1)

reg1_incidence <- incidence(reg1_data$dates)

reg1_results <- estimate_R(reg1_incidence, method="parametric_si", 
                           config = make_config(list(mean_si = serial_interval_mean, 
                                                     std_si = serial_interval_sd)))

#reg1_plot <- estimate_R_plots(reg1_results, what = "R")


colomns <- c(1:4, 8)
reg1_r <- (reg1_results$R[,colomns]) %>% 
   mutate(region = "1") %>% 
   rename(mean_r = `Mean(R)`,
          sd_r = `Std(R)`,
          median_r = `Median(R)`)

reg1_dates <- as.data.frame(reg1_results$dates)
reg1_i <- as.data.frame(reg1_results$I)
reg1_dates_new <- cbind(reg1_dates, reg1_i) %>% 
   rename(dates = 1,
          incidence = 2) %>% 
   mutate(dates = ymd(dates))
reg1_dates_new <- reg1_dates_new[-(1:7), 1:2]
reg1_r_clean <- cbind(reg1_r, reg1_dates_new)

reg1_data_clean <- state_data_clean %>% 
   filter(region == 1)
total_cases <- sum(reg1_data_clean$case)


reg1_r_plot <- reg1_r_clean %>% 
   ggplot() +
   geom_line(aes(dates, mean_r), size = 1.5, color = "black") +
   labs(title = "COVID-19 Rolling 7-day R-values, Montana Region 1, 2020",
        subtitle = paste0("Data current as of ", date_today),
        color = "") +
   ylab("R-value") +
   xlab("") +
   geom_hline(yintercept = 1, color = "black", size = 1.2) +
   scale_x_date(date_breaks = "2 days", date_labels = "%d-%b") +
   scale_y_continuous(breaks = seq(0, 20, 2), labels = seq(0, 20, 2)) +
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
ggsave("C:/R/covid19/state_daily_results/reg1_r_plot.png")

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

ggsave("C:/R/covid19/state_daily_results/reg1_inc_plot.png")



## Reg2 results
reg2_data <- state_data_clean %>% 
   filter(region == 2)

reg2_incidence <- incidence(reg2_data$dates)

reg2_results <- estimate_R(reg2_incidence, method="parametric_si", 
                           config = make_config(list(mean_si = serial_interval_mean, 
                                                     std_si = serial_interval_sd)))

#reg2_plot <- estimate_R_plots(reg2_results, what = "R")


colomns <- c(1:4, 8)
reg2_r <- (reg2_results$R[,colomns]) %>% 
   mutate(region = "2") %>% 
   rename(mean_r = `Mean(R)`,
          sd_r = `Std(R)`,
          median_r = `Median(R)`)

reg2_dates <- as.data.frame(reg2_results$dates)
reg2_i <- as.data.frame(reg2_results$I)
reg2_dates_new <- cbind(reg2_dates, reg2_i) %>% 
   rename(dates = 1,
          incidence = 2) %>% 
   mutate(dates = ymd(dates))
reg2_dates_new <- reg2_dates_new[-(1:7), 1:2]
reg2_r_clean <- cbind(reg2_r, reg2_dates_new)


reg2_data_clean <- state_data_clean %>% 
   filter(region == 2)
total_cases <- sum(reg2_data_clean$case)


reg2_r_plot <- reg2_r_clean %>% 
   ggplot() +
   geom_line(aes(dates, mean_r), size = 1.5, color = "black") +
   labs(title = "COVID-19 Rolling 7-day R-values, Montana Region 2, 2020",
        subtitle = paste0("Data current as of ", date_today),
        color = "") +
   ylab("R-value") +
   xlab("") +
   geom_hline(yintercept = 1, color = "black", size = 1.2) +
   scale_x_date(date_breaks = "2 days", date_labels = "%d-%b") +
   scale_y_continuous(breaks = seq(0, 10, 0.25), labels = seq(0, 10, 0.25)) +
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
ggsave("C:/R/covid19/state_daily_results/reg2_r_plot.png")

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

ggsave("C:/R/covid19/state_daily_results/reg2_inc_plot.png")


## Reg3 results
reg3_data <- state_data_clean %>% 
   filter(region == 3)

reg3_incidence <- incidence(reg3_data$dates)

reg3_results <- estimate_R(reg3_incidence, method="parametric_si", 
                           config = make_config(list(mean_si = serial_interval_mean, 
                                                     std_si = serial_interval_sd)))

#reg3_plot <- estimate_R_plots(reg3_results, what = "R")

colomns <- c(1:4, 8)
reg3_r <- (reg3_results$R[,colomns]) %>% 
   mutate(region = "3") %>% 
   rename(mean_r = `Mean(R)`,
          sd_r = `Std(R)`,
          median_r = `Median(R)`)

reg3_dates <- as.data.frame(reg3_results$dates)
reg3_i <- as.data.frame(reg3_results$I)
reg3_dates_new <- cbind(reg3_dates, reg3_i) %>% 
   rename(dates = 1,
          incidence = 2) %>% 
   mutate(dates = ymd(dates))
reg3_dates_new <- reg3_dates_new[-(1:7), 1:2]
reg3_r_clean <- cbind(reg3_r, reg3_dates_new)


reg3_data_clean <- state_data_clean %>% 
   filter(region == 3)
total_cases <- sum(reg3_data_clean$case)


reg3_r_plot <- reg3_r_clean %>% 
   ggplot() +
   geom_line(aes(dates, mean_r), size = 1.5, color = "black") +
   labs(title = "COVID-19 Rolling 7-day R-values, Montana Region 3, 2020",
        subtitle = paste0("Data current as of ", date_today),
        color = "") +
   ylab("R-value") +
   xlab("") +
   geom_hline(yintercept = 1, color = "black", size = 1.2) +
   scale_x_date(date_breaks = "2 days", date_labels = "%d-%b") +
   scale_y_continuous(breaks = seq(0, 10, 0.25), labels = seq(0, 10, 0.25)) +
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
ggsave("C:/R/covid19/state_daily_results/reg3_r_plot.png")

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

ggsave("C:/R/covid19/state_daily_results/reg3_inc_plot.png")


## Reg4 results
reg4_data <- state_data_clean %>% 
   filter(region == 4)

reg4_incidence <- incidence(reg4_data$dates)

reg4_results <- estimate_R(reg4_incidence, method="parametric_si", 
                           config = make_config(list(mean_si = serial_interval_mean, 
                                                     std_si = serial_interval_sd)))

#reg4_plot <- estimate_R_plots(reg4_results, what = "R")

colomns <- c(1:4, 8)
reg4_r <- (reg4_results$R[,colomns]) %>% 
   mutate(region = "4") %>% 
   rename(mean_r = `Mean(R)`,
          sd_r = `Std(R)`,
          median_r = `Median(R)`)

reg4_dates <- as.data.frame(reg4_results$dates)
reg4_i <- as.data.frame(reg4_results$I)
reg4_dates_new <- cbind(reg4_dates, reg4_i) %>% 
   rename(dates = 1,
          incidence = 2) %>% 
   mutate(dates = ymd(dates))
reg4_dates_new <- reg4_dates_new[-(1:7), 1:2]
reg4_r_clean <- cbind(reg4_r, reg4_dates_new)


reg4_data_clean <- state_data_clean %>% 
   filter(region == 4)
total_cases <- sum(reg4_data_clean$case)


reg4_r_plot <- reg4_r_clean %>% 
   ggplot() +
   geom_line(aes(dates, mean_r), size = 1.5, color = "black") +
   labs(title = "COVID-19 Rolling 7-day R-values, Montana Region 4, 2020",
        subtitle = paste0("Data current as of ", date_today),
        color = "") +
   ylab("R-value") +
   xlab("") +
   geom_hline(yintercept = 1, color = "black", size = 1.2) +
   scale_x_date(date_breaks = "2 days", date_labels = "%d-%b") +
   scale_y_continuous(breaks = seq(0, 10, 0.25), labels = seq(0, 10, 0.25)) +
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
ggsave("C:/R/covid19/state_daily_results/reg4_r_plot.png")

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

ggsave("C:/R/covid19/state_daily_results/reg4_inc_plot.png")


## Reg5 results
reg5_data <- state_data_clean %>% 
   filter(region == 5)

reg5_incidence <- incidence(reg5_data$dates)

reg5_results <- estimate_R(reg5_incidence, method="parametric_si", 
                           config = make_config(list(mean_si = serial_interval_mean, 
                                                     std_si = serial_interval_sd)))

#reg5_plot <- estimate_R_plots(reg5_results, what = "R")

colomns <- c(1:4, 8)
reg5_r <- (reg5_results$R[,colomns]) %>% 
   mutate(region = "5") %>% 
   rename(mean_r = `Mean(R)`,
          sd_r = `Std(R)`,
          median_r = `Median(R)`)

reg5_dates <- as.data.frame(reg5_results$dates)
reg5_i <- as.data.frame(reg5_results$I)
reg5_dates_new <- cbind(reg5_dates, reg5_i) %>% 
   rename(dates = 1,
          incidence = 2) %>% 
   mutate(dates = ymd(dates))
reg5_dates_new <- reg5_dates_new[-(1:7), 1:2]
reg5_r_clean <- cbind(reg5_r, reg5_dates_new)


reg5_data_clean <- state_data_clean %>% 
   filter(region == 5)
total_cases <- sum(reg5_data_clean$case)


reg5_r_plot <- reg5_r_clean %>% 
   ggplot() +
   geom_line(aes(dates, mean_r), size = 1.5, color = "black") +
   labs(title = "COVID-19 Rolling 7-day R-values, Montana Region 5, 2020",
        subtitle = paste0("Data current as of ", date_today),
        color = "") +
   ylab("R-value") +
   xlab("") +
   geom_hline(yintercept = 1, color = "black", size = 1.2) +
   scale_x_date(date_breaks = "2 days", date_labels = "%d-%b") +
   scale_y_continuous(breaks = seq(0, 10, 0.25), labels = seq(0, 10, 0.25)) +
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

ggsave("C:/R/covid19/state_daily_results/reg5_r_plot.png")

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

ggsave("C:/R/covid19/state_daily_results/reg5_inc_plot.png")



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



