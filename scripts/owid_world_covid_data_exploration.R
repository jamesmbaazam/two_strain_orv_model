#package
#' packages
library(openxlsx)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(scales)
library(lubridate)
library(conflicted)
library(janitor)

#' resolve package function conflicts
conflict_prefer('read.xlsx', 'openxlsx')
conflict_prefer('filter', 'dplyr')
conflict_prefer('mutate', 'dplyr')
conflict_prefer('summarise', 'dplyr')
conflict_prefer('week', 'lubridate')


#data
covid_world_data <- openxlsx::read.xlsx('./data/owid-covid-data.xlsx')


#' Subset of African countries to analyse. #The choices represent African countries 
#' with high enough population sizes that have vaccinated less than <1%, 1-10%, 
#' and 50-100% respectively
covid_africa_data <- covid_world_data %>% 
    filter(location == 'Chad' | location == 'South Africa' | location == 'Morocco' ) %>% 
    remove_empty(which = 'cols') %>% 
    as_tibble() %>% 
    mutate(month = month(date, label = TRUE, abbr = TRUE))


stringency_index_and_vax_per_hundred_plot <- ggplot(data = covid_africa_data) + 
    geom_line(aes(x = date, 
                   y = stringency_index,
                   group = location,
                   color = location
                   ),
              size = 2
               ) +
    geom_area(data = covid_africa_data, 
              aes(x = date, 
                  y = total_vaccinations_per_hundred,
                  group = location,
                  fill = location
                  ),
              size = 2,
              alpha = 0.35
              ) +
    scale_x_discrete(breaks = seq(min(ymd(covid_africa_data$date)), 
                                  max(ymd(covid_africa_data$date)), 
                                  by = '2 weeks'
                                  ),
                       labels = seq(min(ymd(covid_africa_data$date)), 
                                    max(ymd(covid_africa_data$date)), 
                                    by = '2 weeks'
                                    )
                       ) +
    labs(title = 'Stringency index and population vaccinated per hundred', 
         y = 'Per hundred'
         ) +
    theme_minimal()

plot(stringency_index_and_vax_per_hundred_plot)
