#' Create a list of tasks to run for each season
#'
#' @param season_list A data frame of the start and end time of train and predict
#' @param start_time Starting time of the data
#' @param start_week Starting week number
#' @param end_week Ending week number
#' @param min_pred Minimum number of weeks to predict
#' @return A data frame of tasks to run
#' @export
create_season_task_list <- function(season_list, starting_time, start_week, end_week, min_pred) {
  predict_start <- season_list$start_time - (start_week - season_list$start_weeks)
  predict_end <- season_list$end_time + (season_list$end_weeks - end_week)
  n_windows <- predict_start - predict_end + 1 - min_pred
  task_list <- data.frame()
  task_list <- cbind(train_start = rep(starting_time, n_windows))
  task_list <- data.frame(task_list)
  task_list$predict_start <- seq(predict_start, predict_end + min_pred, by = -1)
  task_list$train_end <- task_list$predict_start - 1
  task_list$predict_end <- rep(predict_end, n_windows)
  task_list$season <- season_list$season
  task_list$task <- seq(nrow(task_list))
  return(task_list)
}

#' Generate task list
#' @param data A data frame of data available
#' @param start_week Starting week for prediction of each season
#' @param end_week Ending week for prediction of each season
#' @param min_pred Minimum number of week to predict
#' @return A data frame of tasks to run
#' @export
generate_task_list <- function(data, start_week = 40, end_week = 22, min_pred = 4) {
  starting_time <- max(data$time)
  data_avail_table <- create_data_avail_table(data) %>%
    filter(start_weeks < start_week & end_weeks > end_week & start_time - end_time > 0)

  all_task_list <- data.frame()
  for (i in seq(nrow(data_avail_table))) {
    all_task_list <- rbind(all_task_list, create_season_task_list(data_avail_table[i, ], starting_time, start_week, end_week, min_pred))
  }
  return(all_task_list)
}

#' Create data availability table
#' @param data A PILAF object
#' @return A data frame of available data for each season
#' @export
create_data_avail_table <- function(data) {
  time_table <- data %>%
    arrange(year, week) %>%
    mutate(season = ifelse(week >= 30, paste0(year, '-', year + 1), paste0(year - 1, '-', year))) %>%
    group_by(season)
  start_weeks <- time_table %>% top_n(1, time) %>% select(season, time, week) %>% rename(start_weeks = week, start_time = time)
  end_weeks <- time_table %>% top_n(-1, time) %>% select(season, time, week) %>% rename(end_weeks = week, end_time = time)
  time_table <- left_join(start_weeks, end_weeks, by = 'season')
  return(time_table)
}

#' Leave one season out evaluation
#' @param data A PILAF object
#' @param task_list A data frame of tasks to run
#' @param method Method used for forecasting (count or joint)
#' @param formula Model formula
#' @return A list of prediction mean and 95% BCI
#' @export
LOSO_forecast <- function(data, task_list, method = 'count', formula = NULL) {
  eval <- list()
  n <- nrow(task_list)
  p <- dplyr::progress_estimated(n)
  for (i in seq(n)) {
    # browser()
    task_entry <- task_list[i, ]
    train_time <- seq(task_entry$train_start, task_entry$train_end, by = -1)
    test_time <- seq(task_entry$predict_start, task_entry$predict_end, by = -1)
    train_data <- data[data$time %in% train_time, ]
    test_data <- data[data$time %in% test_time, ]
    try({
      fit <- forecast(train_data, test_data$time, test_data$week, return_model = F, verbose = F, method = method, formula = formula)
      eval[[i]] <- fit$forecast})
    p$pause(0.1)$tick()$print()
  }
  return(eval)
}

