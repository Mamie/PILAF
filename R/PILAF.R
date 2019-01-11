#' Constructor to PILAF object
#'
#' Constructor method to create an object from PILAF class
#'
#' @param time A numeric vector of time.
#' @param week A numeric vector of week.
#' @param coal A numeric vector of coalescent event counts.
#' @param samp A numeric vector of sampling event counts.
#' @param ILI A numeric vector ILI counts.
#' @param coal.E A numeric vector of offset to coal.
#' @param samp.E A numeric vector of offset to samp.
#' @param ILI.E A numeric vector of offset to ILI.
#' @param iter A numeric vector to indicate group of trajectories.
#' @return A PILAF object.
#' @export
PILAF <- function(time, week, coal, samp, ILI, coal.E, samp.E, ILI.E, iter) {
  pilaf <- data.frame(time = time,
                     week = week,
                     coal = coal,
                     samp = samp,
                     ILI = ILI,
                     coal.E = coal.E,
                     samp.E = samp.E,
                     ILI.E = ILI.E,
                     iter = iter)
  pilaf = dplyr::arrange(pilaf, iter, time)
  class(pilaf) = c("PILAF", class(pilaf))
  return(pilaf)
}

#' Check if an Object is PILAF
#'
#' Check if an object is of PILAF class
#'
#' @param x An object.
#' @return A logical scalar indicating whether the object belongs to PILAF class.
#' @export
is.PILAF = function(x) {
  inherits(x, "PILAF")
}

#' Plot PILAF trajectory
#'
#' Plot the count trajectory of PILAF class.
#'
#' @param x A PILAF object
#' @export
#' @import ggplot2
plot.PILAF = function(x, ...) {
  x = with(x,
          dplyr::data_frame(time = time, `coalescent event counts`=coal, `sampling event counts`=samp, `ILI counts`=ILI, iter=iter))
  x = tidyr::gather(x, type, counts, -c(time, iter))
  ggplot(data=x) +
    geom_line(aes(x=time, y=counts, group=iter, color=iter), ...) +
    facet_wrap(~type, scales='free', ncol=1) +
    theme_classic() +
    xlab('Time to present (weeks)') +
    theme(axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position='none',
          strip.background=element_blank())
}

#' Generic function for forecast
#' @export
forecast = function(x, ...) {
  UseMethod('forecast', x)
}

#' Forecast a PILAF Object.
#'
#' Forecast a PILAF object using one of the three methods: count only, joint
#' modeling of count and coalescent events, joint modeling of count, coalescent
#' and sampling envets.
#'
#' @param x A PILAF object.
#' @param time.forecast A negative numeric vector of timepoints to forecast.
#' @param week.forecast A numeric vector of week corresponding to timepoints to forecast.
#' @param formula An optional character scalar of formula (debugging).
#' @export
forecast.PILAF = function(x, time.forecast, week.forecast, method='count', formula='NULL',
                          return_model = F, verbose = F) {
  if (!method %in% c('count', 'joint', 'ps')) {
    warnings('Forecast method not recognized. Use default="count"')
    method = 'count'
  }
  iter.ids = unique(x$iter)
  forecast.all = Forecast()
  n = length(iter.ids)
  if (verbose) p = dplyr::progress_estimated(n)
  models = list()
  for (iter.id in iter.ids) {
    #if (verbose) print(iter.id)
    x.iter = dplyr::filter(x, iter==iter.id)
    x.forecast = rbind(data.frame(time=time.forecast,
                                  week=week.forecast,
                                  coal=NA, samp=NA, ILI=NA,
                                  coal.E=1, samp.E=1, ILI.E=1,
                                  iter=x.iter$iter[1]), x.iter)
    # time is with respect to present, future has smaller value of time
    x.forecast = x.forecast[order(x.forecast[,'iter'], x.forecast[,'time']),]
    try({
      invisible({
        forecast.inla = eval(parse(text=paste0('forecast.PILAF.', method, '(x.forecast, formula=',
                                               formula, ')')))
      })
      if (return_model) models[[iter.id]] = forecast.inla
      forecast = forecast.inla$summary.fitted.values[1:length(time.forecast),]
      forecast = Forecast(time = time.forecast,
                          week = week.forecast,
                          mean = forecast$mean,
                          quant0.025 = forecast$`0.025quant`,
                          quant0.975 = forecast$`0.975quant`,
                          iter = iter.id)
      forecast.all = rbind(forecast.all, forecast)
    })
    if (verbose) p$pause(0.1)$tick()$print()
  }
  return(list(forecast = forecast.all, models = models))
}

#' Forecast using Only Count
#'
#' @param x A PILAF object with NA for values to forecast.
#' @return A PILAF object with forecasted ILI counts and 95 % BCI.
#' @export
forecast.PILAF.count = function(x, formula=NULL) {
  link <- rep(1, nrow(x))
  if(is.null(formula)) {
    formula = as.formula('ILI ~ -1 + f(time, model="ar", order=2)')
  }
  forecast = INLA::inla(formula, family = "poisson", data = x,
                        control.predictor = list(compute = T, link = link),
                        E = x$ILI.E, control.compute = list(config = T))
  return(forecast)
}

#' Forecast using Count and Coalescent events
#'
#' @param x A PILAF object with NA for values to forecast.
#' @param formula A INLA formula.
#' @return A PILAF object with forecasted ILI counts and 95 % BCI.
#' @export
forecast.PILAF.joint = function(x, formula=NULL, verbose = F) {
  n <- nrow(x)
  link <- c(rep(1, n), rep(2, n))
  E <- c(x$ILI.E, x$coal.E)
  time <- rep(x$time, 2)
  week <- rep(x$week, 2)
  beta0 <- c(rep(0, n), rep(1, n))
  w0 <- c(rep(1, n), rep(0, n))
  w <- c(rep(0, n), rep(-1, n))
  Y <- matrix(NA, nrow=2*n, ncol=2)
  Y[1:n, 1] <- x$ILI
  Y[(1:n)+n, 1] <- x$coal
  r <- rep(1:2, each = n)
  X <- list(time=time, time2=time, week=week, week2=week, beta0=beta0, w0=w0, w=w, Y=Y)
  if(is.null(formula)) {
    formula <- as.formula(paste0("Y ~ -1 + beta0 + f(time, w0, model='ar', order=2) +",
                                "f(time2, w, copy='time', fixed=F)"))
  }
  forecast <- INLA::inla(formula, family = c("poisson", "poisson"), data = X,
                        control.predictor = list(compute = T, link = link),
                        E=X$E, verbose = verbose, control.compute = list(config = T))
  return(forecast)
}
