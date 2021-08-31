library(deSolve)
library(ggplot2)

source('./scripts/Kissler_two_strain_model/Kissler_two_strain_model.R')
start_date <- as.Date('2020-01-01')
end_date <- as.Date('2030-01-01')
plot_dates <- seq.Date(start_date, end_date, by = 'year')
seed_I <- 50 / 58000000
step_size <- 0.1
MAXTIME <- 52 * 10
tmp <- readxl::read_excel('./scripts/sample_R_code/Kissler_two_strain_model/parameters.xlsx')
pp <- as.list(tmp$value)
names(pp) <- tmp$par_name
init <- readRDS('./scripts/Kissler_two_strain_model/init.RDS')

## Scenarios with moderate cross immunity from seasonal coronaviruses

### With 40 week duration of immunity to SARS-CoV-2:

ts_cross <- data.frame(lsoda(
  y = init,                  # Initial conditions for population
  times = seq(0, MAXTIME, step_size), # Timepoints for evaluation
  func = corona_model,              # Function to evaluate
  parms = pp,             # Vector of parameters
))  
ts_cross$date <- as.Date(ts_cross$time*7, origin = start_date)
ts_cross$I1 <- ts_cross$I1S2 + ts_cross$I1E2 + ts_cross$I1I2 + ts_cross$I1R2
ts_cross$I2 <- ts_cross$S1I2 + ts_cross$E1I2 + ts_cross$I1I2 + ts_cross$R1I2

(ggplot(ts_cross)
  + geom_line(aes(x = date, y = I1*1000), color = '#E03127')
  + geom_line(aes(x = date, y = I2*1000), color = '#F7941E')
  + xlab('')
  + scale_x_date(breaks = plot_dates, labels = format(plot_dates, format = '%Y'))
  + ylab('Prevalence per 1,000 people')
  + ylim(0,100)
  + theme_bw(base_size = 12)
)

### With 2 year duration of immunity to SARS-CoV-2:

pp_lti <- pp
pp_lti$sigma_1 <- 1 / 104

ts_cross <- data.frame(lsoda(
  y = init,                  # Initial conditions for population
  times = seq(0, MAXTIME, step_size), # Timepoints for evaluation
  func = corona_model,              # Function to evaluate
  parms = pp_lti,             # Vector of parameters
))  
ts_cross$date <- as.Date(ts_cross$time*7, origin = start_date)
ts_cross$I1 <- ts_cross$I1S2 + ts_cross$I1E2 + ts_cross$I1I2 + ts_cross$I1R2
ts_cross$I2 <- ts_cross$S1I2 + ts_cross$E1I2 + ts_cross$I1I2 + ts_cross$R1I2

(ggplot(ts_cross)
  + geom_line(aes(x = date, y = I1*1000), color = '#E03127')
  + geom_line(aes(x = date, y = I2*1000), color = '#F7941E')
  + xlab('')
  + scale_x_date(breaks = plot_dates, labels = format(plot_dates, format = '%Y'))
  + ylab('Prevalence per 1,000 people')
  + ylim(0,100)
  + theme_bw(base_size = 12)
)
