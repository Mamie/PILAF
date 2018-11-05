create_season_task_list <- function(season_list, start_time, start_week, end_week, min_pred) {
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

generate_task_list <- function(data, start_week = 40, end_week = 22, min_pred = 4) {
  starting_time <- max(data$time)
  data_avail_table <- create_data_avail_table(data) %>%
    filter(start_weeks < start_week & end_weeks > end_week & start_time - end_time > 0)

  all_task_list <- data.frame()
  for (i in seq(nrow(data_avail_table))) {
    all_task_list <- rbind(all_task_list, create_season_task_list(data_avail_table[i, ], data, start_week, end_week, min_pred))
  }

  return(all_task_list)
}

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
