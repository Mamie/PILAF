# forecast summary statistics
library(dplyr)
library(ggplot2)
library(PILAF)

H3N2_tree <- ape::read.nexus(system.file("extdata/H3N2/MCC.tree",
                                         package = "PILAF"))
plot(ape::ladderize(H3N2_tree),show.tip.label=FALSE)
H3N2_last_time <- 2019.0739726027398
H3N2_root_time <- H3N2_last_time - max(phylodyn::summarize_phylo(H3N2_tree)$coal_times)
ar2_formula <- Y ~ -1 + beta0 +
  f(time, model="ar", order = 2, hyper = hyper, constr = FALSE) +
  f(time2, w, copy="time", fixed=FALSE, param=c(0, beta1_prec))


seasons <- seq(2014, 2018)
week_start <- c(40:52, 1:15)
pred <- 4

# Prediction tasks
tasks <- expand.grid(year = seasons, wk_start = week_start) %>%
  arrange(year, wk_start) %>%
  filter(!(year == 2014 & wk_start < 40)) %>%
  rbind(data.frame(year = 2019, wk_start = 1:15)) %>%
  mutate(pred = pred)

# forecast for each task
forecast_res <- purrr::pmap(tasks,
                            ~forecast_starting(H3N2_tree,
                                               last_time = H3N2_last_time,
                                               week_start = ..2,
                                               year_start = ..1,
                                               label = "H3N2_Ne",
                                               formula = ar2_formula,
                                               pred = ..3,
                                               include_trunc = F))

tasks$start <- purrr::pmap(tasks, ~ifelse(..2 < 21, ..1 - 1, ..1))

# obtain forecast results
H3N2_forecasts <- purrr::pmap(list(res = forecast_res, name = 1:length(forecast_res),
                                   year = tasks$start),
                              ~mutate(..1$df[1:pred,], forecast = ..2)) %>%
  dplyr::bind_rows()
H3N2_forecasts$week <- purrr::map_dbl(H3N2_forecasts$Time, ~ lubridate::week(lubridate::round_date(lubridate::date_decimal(.x))))

H3N2_forecasts <- H3N2_forecasts %>%
  arrange(forecast, week) %>%
  group_by(forecast) %>%
  mutate(i = 1:pred)

# Use whole data AR(2) as reference
H3N2_Ne_forecast <- BNPR_forecast(H3N2_tree, last_time = H3N2_last_time, formula = ar2_formula, pred = 0, pref = T)
H3N2_truth <- with(H3N2_Ne_forecast,
                     data.frame(time = H3N2_last_time - result$summary.random$time$ID,
                                week = weeks,
                                H3N2_Ne = effpopmean,
                                lwr = effpop025,
                                upr = effpop975))
H3N2_truth$year <- floor(H3N2_truth$time)

# Compute summary statistics
summary_stats <- H3N2_forecasts %>%
  mutate(year = floor(Time)) %>%
  select(year, Median, Upper, Lower, forecast, week, i) %>%
  left_join(H3N2_truth %>% select(year, week, H3N2_Ne), by = c("year", "week")) %>%
  tidyr::drop_na() %>%
  group_by(i) %>%
  summarize(n = n(), # number of forecasts
            MRE = mean_relative_error(Median, H3N2_Ne),
            COV = mean_coverage(Lower, Upper, H3N2_Ne),
            MW = mean_width(Lower, Upper, H3N2_Ne)) %>%
  rename(week = i)

plot_nice_table <- function(table, digits = 2) {
  kableExtra::kable(table, digits = digits) %>% kableExtra::kable_styling()
}
plot_nice_table(summary_stats)

# repeat for H1N1

H1N1_tree <- ape::read.nexus(system.file("extdata/H1N1/MCC.tree",
                                         package = "PILAF"))
H1N1_last_time <- 2019.0931506849315

forecasts <- purrr::pmap(tasks,
                         ~try(forecast_starting(H1N1_tree,
                                            last_time = H1N1_last_time,
                                            week_start = ..2,
                                            year_start = ..1,
                                            label = "H1N1_Ne",
                                            formula = ar2_formula,
                                            pred = ..3,
                                            include_trunc = F)))

# find failed forecasts
failed <- which(purrr::map_chr(forecasts, class) != "list")
H1N1_forecasts <- purrr::pmap(list(res = forecasts[-failed], name = (1:length(forecasts))[-failed],
                                   year = tasks$start[-failed]),
                              ~mutate(..1$df[1:pred,], forecast = ..2)) %>%
  dplyr::bind_rows()
H1N1_forecasts$week <- purrr::map_dbl(H1N1_forecasts$Time, ~lubridate::week(lubridate::round_date(lubridate::date_decimal(.x))))

H1N1_forecasts <- H1N1_forecasts %>%
  arrange(forecast, week) %>%
  group_by(forecast) %>%
  mutate(i = 1:pred)

H1N1_Ne_forecast <- BNPR_forecast(H1N1_tree, last_time = H1N1_last_time, formula = ar2_formula, pred = 0, pref = T)
H1N1_truth <- with(H1N1_Ne_forecast,
                   data.frame(time = H1N1_last_time - result$summary.random$time$ID,
                              week = weeks,
                              H1N1_Ne = effpopmean,
                              lwr = effpop025,
                              upr = effpop975))
H1N1_truth$year <- floor(H1N1_truth$time)

# Compute summary statistics
summary_stats <- H1N1_forecasts %>%
  mutate(year = floor(Time)) %>%
  select(year, Median, Upper, Lower, forecast, week, i) %>%
  left_join(H1N1_truth %>% select(year, week, H1N1_Ne), by = c("year", "week")) %>%
  tidyr::drop_na() %>%
  group_by(i) %>%
  summarize(n = n(), # number of forecasts
            MRE = mean_relative_error(Median, H1N1_Ne),
            COV = mean_coverage(Lower, Upper, H1N1_Ne),
            MW = mean_width(Lower, Upper, H1N1_Ne)) %>%
  rename(week = i)

plot_nice_table(summary_stats)
