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
#' @param no_y_axis Whether to emove y axis label
#' @param season_label Whether to remove season label
#' @param title The title of the plot
#' @param scales Parameter for facet_wrap ("free_x", "free_y", "fixed")
#' @param log_y Whether to log y axis
#' @param error_band_names A data frame that record the variable names for lower and upper bound
#' @param correlate_ts A data frame containing three columns (time, week, time series) for annotation; must be in order
#' @param label_x x coordinate of correlation annotation
#' @param label_y y coordinate of correlation annotation
#' @param add_forecasts A data frame containing the forecast values (output df from forecast_starting)
#' @return A ggplot object visualizing from week 40 to week 21
#' @export
#' @import ggplot2
visualize_flu_season <- function(..., time_series_names, datalist = NULL,
                                 subset_seasons = NULL, ncol = 1, no_y_axis = TRUE,
                                 season_label = T, title = NULL, scales = "fixed",
                                 log_y = T, error_band_names = NULL, correlate_ts = NULL,
                                 label_x = NULL, label_y = NULL, add_forecasts = NULL) {
  datalist <- c(list(...), datalist)
  stopifnot(length(time_series_names) == length(datalist))
  if (is.null(error_band_names)) {
    data <- purrr::map2(datalist, time_series_names, ~ mutate(.x, name = .y) %>%
                          rename_(.dots = list("time_series" = .y)) %>%
                          select(year, week, time_series, name))
  } else {
    stopifnot(ncol(error_band_names) == 2 & nrow(error_band_names) == length(datalist))
    data <- purrr::pmap(list(data = datalist, ts = time_series_names,
                             lwr = error_band_names[,1], upr = error_band_names[,2]),
                        ~ mutate(..1, name = ..2) %>%
                          rename_(.dots = list("time_series" = ..2,
                                               "quant025" = ..3,
                                               "quant975" = ..4)) %>%
                          select(year, week, time_series, name, quant025, quant975)
                        )

  }

  data <- dplyr::bind_rows(data) %>%
    mutate(name = factor(name, levels = time_series_names))
  if(!is.null(correlate_ts)) {
    colnames(correlate_ts) <- c("year", "week", "ts")
    data <- left_join(data, correlate_ts, by = c("year", "week"))
  }

  flu_season <- data %>%
    filter(week >= 40 | week <= 21)

  flu_season$season <- purrr::pmap_chr(flu_season[, c("year", "week")],
                                        ~ ifelse(..2 <= 21, paste_dash(..1 - 1, ..1), paste_dash(..1, ..1 + 1)))
  #flu_season$week <- factor(flu_season$week, levels = c(40:53, 1:21))


  if (!is.null(add_forecasts)) {
    input <- data.frame(year = floor(add_forecasts$Time),
                        week = add_forecasts$week,
                        time_series = add_forecasts$Mean,
                        name = add_forecasts$Label,
                        quant025 = add_forecasts$Lower,
                        quant975 = add_forecasts$Upper,
                        forecast = add_forecasts$forecast)
    input$season <- purrr::pmap_chr(input[, c("year", "week")],
                                    ~ ifelse(..2 <= 21, paste_dash(..1 - 1, ..1), paste_dash(..1, ..1 + 1)))
    #View(input)
    flu_season$forecast <- "0"

    flu_season <- rbind(flu_season, input)
    flu_season$is_forecast <- as.character(flu_season$forecast) != "0"

  }

  flu_season$week <- factor(flu_season$week, levels = c(40:53, 1:21))
  if (!is.null(subset_seasons))
    flu_season <- dplyr::filter(flu_season, season %in% subset_seasons)

  if (!is.null(add_forecasts)) {
    p <- ggplot(data = flu_season, aes(x = week, y = time_series, group = forecast, color = is_forecast)) +
      scale_color_manual(values = c("black", "turquoise"))
  }
  else p <- ggplot(data = flu_season, aes(x = week, y = time_series, group = season))
  if (!is.null(error_band_names)) p <- p + geom_ribbon(aes(ymin = quant025, ymax = quant975), fill = "lightgray", color = "white", alpha = 0.7)
  p <- p + geom_line(size = 0.7)

  if (!is.null(correlate_ts)) {
    annot <- flu_season %>%
      select(name, season, year, week, ts, time_series) %>%
      group_by(name, season) %>%
      summarize(rho = cor.test(ts, time_series)$estimate,
             p = cor.test(ts, time_series)$p.value,
             max_y = max(time_series))

    if (is.null(label_x)) label_x <- 0
    if (is.null(label_y)) label_y <- max(annot$max_y)
    p <- p +
      geom_text(data = annot, aes(x = label_x, y = label_y, label = round(rho, 2)),
                colour="black", inherit.aes=FALSE, parse=FALSE)
  }


  if (log_y) p <- p + scale_y_log10()

  if (length(unique(flu_season$name)) > 1)
    p <- p + facet_wrap(~name + season, ncol = ncol, scales = scales)
  else
    p <- p + facet_wrap(~season, ncol = ncol, scales = scales)

  if (season_label) p <- p + geom_text(data = group_by(flu_season, season, name) %>%
                                top_n(1, time_series), aes(label = season)) +
  if (no_y_axis) p <- p + theme(axis.title.y = element_blank())
  if (!is.null(title)) p <- p + ggtitle(title)

  p <- p +
    scale_x_discrete(breaks = c(40, 45, 50, 1, 5, 10, 15, 20)) +
    theme_classic() +
    theme(legend.position = "none", strip.background = element_blank())
  return(p)
}


