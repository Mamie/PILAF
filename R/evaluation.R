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
  models <- list()
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
      fit <- forecast(train_data, test_data$time, test_data$week, return_model = T, verbose = F, method = method, formula = formula)
      eval[[i]] <- fit$forecast
      models[[i]] <- fit$models[[1]]})
    p$pause(0.1)$tick()$print()
  }
  return(list(forecast = eval, models = models))
}


#' return the peak week for each model  by simulation
compute_peak_week <- function(task_list, models, task_season_range, n = 100) {
  peak_week <- data.frame(time = integer(), peak_week = integer(), quant0.025 = double(), quant0.975 = double(), median = double(), mean = double())
  peak_ILI <- data.frame(time = integer(), peak_ILI = integer(), quant0.025 = double(), quant0.975 = double(), median = double(), mean = double())
  p <- dplyr::progress_estimated(length(models))
  for (i in seq_along(models)) {
    try({
      posterior_samples <- INLA::inla.posterior.sample(n, models[[i]])
      task_entry <- task_list[i, ]
      seaons_range <- task_season_range[task_season_range$season == task_entry$season, ]
      peak_week_samples <- data.frame(peak_week = integer(), peak_ILI = double())

      for (j in seq_along(posterior_samples)) {
        pred <- find_peak_week(posterior_samples[[j]], task_entry, season_range)
        peak_week_samples <- rbind(peak_week_samples, pred)
      }

      peak_week <- rbind(peak_week, cbind(time = i, get_summary_stats(peak_week_samples$peak_week)))
      peak_ILI <- rbind(peak_ILI, cbind(time = i, get_summary_stats(peak_week_samples$peak_ILI)))
    })
    p$pause(0.1)$tick()$print()
  }
  return(list(peak_week = peak_week, peak_ILI = peak_ILI))
}

get_summary_stats <- function(samples) {
  quants <- unname(quantile(samples, c(0.025, 0.975)))
  quant0.025 <- quants[1]
  quant0.975 <- quants[2]
  mean <- mean(samples)
  median <- median(samples)
  return(data.frame(quant0.025 = quant0.025, quant0.975 = quant0.975, mean = mean, median = median, stringsAsFactors = F))
}

find_peak_week <- function(posterior_sample, task_entry, season_range) {
  actual_times <- seq(task_entry$predict_end, task_entry$train_start)
  latent_sample <- posterior_sample$latent
  predictor_idx <- which(grepl('Predictor', rownames(latent_sample), fixed = F))[seq_along(actual_times)]
  latent_sample_predictor <- exp(latent_sample[predictor_idx])
  time_latent_map <- hashmap::hashmap(actual_times, latent_sample_predictor)

  season_range_seq <- seq(season_range$time_end, season_range$time_start)
  ILI <- sapply(season_range_seq, function(x) time_latent_map[[x]])
  return(data.frame(peak_week = season_range_seq[which.max(ILI)[1]], peak_ILI = max(ILI)))
}