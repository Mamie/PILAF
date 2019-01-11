#' Constructor to Forecast object
#'
#' Constructor method to create an object from Forecast class
#'
#' @param time A numeric vector of time.
#' @param mean A numeric vector of forecast mean.
#' @param quant0.025 A numeric vector of 0.025 quantile.
#' @param quant0.975 A numeric vector 0.975 quantile.
#' @param iter A numeric vector to indicate group of trajectories.
#' @return A Forecast object.
#' @export
Forecast <- function(time=integer(), week=integer(), mean=double(), `quant0.025`=double(),
                     `quant0.975`=double(), iter=integer()) {
  forecast <- data.frame(time = time, mean = mean, `quant0.025` = `quant0.025`,
                         `quant0.975` = `quant0.975`, iter=iter)
  class(forecast) <- c("Forecast", class(forecast))
  return(forecast)
}

#' Check if an Object is Forecast
#'
#' Check if an object is of Forecast class
#'
#' @param x An object.
#' @return A logical scalar indicating whether the object belongs to Forecast class.
#' @export
is.Forecast = function(x) {
  inherits(x, "Forecast")
}


#' Plot the Forecast Results
#'
#' Plot the forecasting results
#'
#' @param forecast A Forecast object.
#' @param x A PILAF object.
#' @param truth An optional PILAF object that shows the truth.
#' @return A ggplot object for the forecasting results.
#' @export
plot.Forecast = function(x, pilaf, truth=NULL, ...) {
  x$type='ILI'
  pilaf = with(pilaf, data.frame(time=time, coalescent=coal, sampling=samp,
                                 ILI=ILI, iter=iter))
  pilaf = tidyr::gather(pilaf, type, counts, -c(time, iter))
  p = ggplot(data=pilaf) +
    geom_line(aes(x=time, y=counts, group=iter, color=iter), alpha=0.8, size=0.1) +
    geom_line(data=x, aes(x=time, y=mean, group=iter, color=iter), alpha=0.8, size=0.1) +
    geom_ribbon(data=x, aes(x=time, ymin=quant0.025, ymax=quant0.975, group=iter), ...) +
    geom_vline(aes(xintercept=max(x$time), linetype='dotted'), alpha=0.1) +
    facet_wrap(~type, scales='free_y', ncol=1) +
    xlab('time to present') +
    theme_classic() +
    theme(axis.ticks.x=element_blank(), legend.position='none',
          strip.background=element_blank())
  if(!is.null(truth)) {
    truth$type = 'ILI'
    p = p +
      geom_line(data=truth, aes(x=time, y=ILI, group=iter, color=iter),
                linetype='dotdash', alpha=0.3, size=0.3)
  }
  return(p)
}

#' Compute absolute error and width between forecast and truth.
computeAEW = function(x, truth) {
  x %>%
    dplyr::left_join(truth, by = c("time", "iter")) %>%
    group_by(time) %>%
    dplyr::summarize(`mean relative error` = mean(abs(mean - ILI)/ILI),
                  `mean relative width` = mean((`quant0.975` - `quant0.025`)/ILI),
                  `coverage` = mean((`quant0.025` <= ILI & ILI <= `quant0.975`))) %>%
    dplyr::select(time, `mean relative error`, `mean relative width`, `coverage`)
}

#' Generic function for comparison
#' @export
compare <- function(x, ...) {
  UseMethod('compare', x)
}

#' Compare performance of two models by AE and W
#'
#' @param x, y Two Forecast object.
#' @return A ggplot object comparing the AE and W between the two methods.
#' @export
compare.Forecast = function(x, y, truth, method = c(1, 2)) {
  x = computeAEW(x, truth)
  x$method = method[1]
  y = computeAEW(y, truth)
  y$method = method[2]
  rbind(x, y) %>%
    tidyr::gather(type, value, -c(time, method)) %>%
    ggplot(data = ., aes(x = time, y = value, color = method)) +
    geom_point() +
    geom_line() +
    facet_wrap(~type, ncol = 1, scale = 'free_y') +
    scale_colour_manual(values=c("#0083C3", "#EB975F")) +
    theme_classic() +
    theme(strip.background = element_blank(), axis.title.y = element_blank()) +
    scale_x_reverse(breaks = c(-2, -4, -6, -8, -10)) +
    xlab('time')
}

#' Summarize the Forecast
#'
#' Summarize the forecast by  mean absolute error and mean relative
#' width.
#' @param x A forecast object.
#' @param truth A data frame containing time, true ILI.
#' @return A data frame containing the absolute error of forecast.
#' @export
#' @import dplyr
summary.Forecast <- function(x, truth) {
  mean_narm <- partial(mean, na.rm = T)
  merged <- computeAEW(x, truth) %>%
    dplyr::group_by(time) %>%
    dplyr::summarize(MRE = mean_narm(`mean relative error`),
                     MRW = mean_narm(`mean relative width`),
                     CV = mean_narm(`coverage`))
  return(merged)
}

#' Compute SE for given variable.
summarySE = function(data, measure_var, group_var) {
  data %>%
    dplyr::group_by_(.dots = group_var) %>%
    dplyr::mutate(n = dplyr::n()) %>%
    dplyr::mutate_(.dots = setNames(paste0('mean(`', measure_var, '`)'), 'mean')) %>%
    dplyr::mutate_(.dots = setNames(paste0('sd(`', measure_var, '`)'), 'sd')) %>%
    dplyr::select(-dplyr::matches(measure_var)) %>%
    dplyr::distinct() %>%
    dplyr::mutate(SE = sd/sqrt(n))
}

#' Generic function for evaluation
#' @export
evaluate = function(x, ...) {
  UseMethod('evaluate', x)
}

#' Evaluate the Forecast
#'
#' Evaluate the forecast by  absolute error and BCI
#' width. Return a plot of AE and width.
#' @param x A forecast object.
#' @param truth A data frame containing time, true ILI.
#' @return A ggplot object
#' @export
#' @import dplyr
evaluate.Forecast = function(x, truth) {
  merged = computeAEW(x, truth)
  merged_AE = merged %>%
    dplyr::select(time, `mean relative error`)
  merged_W = merged %>%
    dplyr::select(time, `mean relative width`)
  merged_CV <- merged %>%
    dplyr::select(time, `coverage`)
  AE = summarySE(merged_AE, 'mean relative error', 'time')
  AE$type = 'mean relative error'
  W = summarySE(merged_W, 'mean relative width', 'time')
  W$type = 'mean relative width'
  CV = summarySE(merged_CV, 'coverage', 'time')
  CV$type = 'coverage'

  rbind(AE, W, CV) %>%
    dplyr::ungroup() %>%
    ggplot(data=., aes(x = time, y = mean)) +
    geom_line(size=0.3, alpha=0.8) +
    geom_errorbar(aes(ymin = mean - SE, ymax = mean + SE), size = 0.3, alpha = 0.8, width = 0.3) +
    facet_wrap(~type, scale='free_y') +
    theme_classic() +
    theme(strip.background = element_blank(), axis.title.y = element_blank()) +
    xlab('time') +
    scale_x_reverse()
}





