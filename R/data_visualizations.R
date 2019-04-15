paste_dash <- function(...) paste(..., sep = "-")

#' Visualize the time series for each Flu Season
#'
#' Plot the time series for each flu season (week 40 - week 21).
#'
#' @param data A data frame contain three columns (year, week, the time series);
#' the column names for year and week must be "year" and "week"
#' @param datalist A list of data frame to plot multiple time series (alternative to data)
#' @param time_series_names The column name for the time series (same in length as data + datalist)
#' @param subset_seasons Seasons to keep
#' @param ncol Number of columns
#' @param
#' @return A ggplot object visualizing from week 40 to week 21
#' @export
#' @import ggplot2
visualize_flu_season <- function(..., time_series_names, datalist = NULL,
                                 subset_seasons = NULL, ncol = 1, no_y_axis = FALSE) {
  datalist <- c(list(...), datalist)
  stopifnot(length(time_series_names) == length(datalist))
  data <- purrr::map2(datalist, time_series_names, ~ mutate(.x, name = .y) %>%
                         rename_(.dots = list("time_series" = .y)) %>%
                         select(year, week, time_series, name)) %>%
    dplyr::bind_rows() %>%
    mutate(name = factor(name, levels = time_series_names))

  flu_season <- data %>%
    filter(week >= 40 | week <= 21)

  flu_season$season <- purrr::pmap_chr(flu_season[, c("year", "week")],
                                        ~ ifelse(..2 <= 21, paste_dash(..1 - 1, ..1), paste_dash(..1, ..1 + 1)))
  flu_season$week <- factor(flu_season$week, levels = c(40:53, 1:21))

  if (!is.null(subset_seasons))
    flu_season <- dplyr::filter(flu_season, season %in% subset_seasons)

  p <- ggplot(data = flu_season, aes(x = week, y = time_series, group = season, color = season)) +
    geom_line() +
    geom_text(data = group_by(flu_season, season, name) %>%
                top_n(1, time_series), aes(label = season)) +
    scale_y_log10() +
    scale_x_discrete(breaks = c(40, 45, 50, 1, 5, 10, 15, 20)) +
    theme_classic() +
    theme(legend.position = "none", strip.background = element_blank())
  if (no_y_axis) p <- p + theme(axis.title.y = element_blank())
  if (length(datalist) > 1) p <- p + facet_wrap(~name, ncol = ncol, scale = "free_y")
  return(p)
}
