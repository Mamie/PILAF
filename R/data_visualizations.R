paste_dash <- function(...) paste(..., sep = "-")

#' Visualize the time series for each Flu Season
#'
#' Plot the time series for each flu season (week 40 - week 21).
#'
#' @param data A data frame contain three columns (year, week, name_of_time_series)
#' @param time_series_name The column name for the time series
#' @param subset_seasons Seasons to keep
#' @return A ggplot object visualizing from week 40 to week 21
#' @export
#' @import ggplot2
visualize_flu_season <- function(data, time_series_name, subset_seasons = NULL) {
  flu_season <- data %>%
    filter(week >= 40 | week <= 21)
  colnames(flu_season)[colnames(flu_season) == time_series_name] <- "time_series"
  flu_season$season <- purrr::pmap_chr(flu_season[, c("year", "week")],
                                        ~ ifelse(..2 <= 21, paste_dash(..1 - 1, ..1), paste_dash(..1, ..1 + 1)))
  flu_season$week <- factor(flu_season$week, levels = c(40:53, 1:21))

  if (!is.null(subset_seasons))
    flu_season <- dplyr::filter(flu_season, season %in% subset_seasons)

  ggplot(data = flu_season, aes(x = week, y = time_series, group = season, color = season)) +
    geom_line() +
    geom_text(data = group_by(flu_season, season) %>%
                top_n(1, time_series), aes(label = season), position = "dodge") +
    scale_y_log10() +
    scale_x_discrete(breaks = c(40, 45, 50, 1, 5, 10, 15, 20)) +
    theme_classic() +
    theme(legend.position = "none") +
    ylab(time_series_name)
}
