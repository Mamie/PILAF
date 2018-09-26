#' Constructor to PILAF object
#'
#' Constructor method to create an object from PILAF class
#'
#' @param time A numeric vector of time.
#' @param coal A numeric vector of coalescent event counts.
#' @param samp A numeric vector of sampling event counts.
#' @param ILI A numeric vector ILI counts.
#' @param coal.E A numeric vector of offset to coal.
#' @param samp.E A numeric vector of offset to samp.
#' @param ILI.E A numeric vector of offset to ILI.
#' @param iter A numeric vector to indicate group of trajectories.
#' @return A PILAF object.
#' @export
PILAF = function(time, coal, samp, ILI, coal.E, samp.E, ILI.E, iter) {
  pilaf = data.frame(time=time,
                     coal=coal,
                     samp=samp,
                     ILI=ILI,
                     coal.E=coal.E,
                     samp.E=samp.E,
                     ILI.E=ILI.E,
                     iter=iter)
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
plot.PILAF = function(x) {
  x = with(x,
          data.frame(time = time, coalescent=coal, sampling=samp,ILI=ILI,iter=iter))
  x = tidyr::gather(x, type, counts, -c(time, iter))
  ggplot(data=x) +
    geom_line(aes(x=time, y=counts, group=iter, color=iter, alpha=0.8), size=0.1) +
    facet_wrap(~type, scales='free', ncol=1) +
    theme_classic() +
    xlab('time to present') +
    theme(axis.ticks.x=element_blank(),
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
#' @param formula An optional character scalar of formula (debugging).
#' @export
forecast.PILAF = function(x, time.forecast, method='count', formula='NULL') {
  if (!method %in% c('count', 'joint', 'ps')) {
    warnings('Forecast method not recognized. Use default="count"')
    method = 'count'
  }
  iter.ids = unique(x$iter)
  forecast.all = c()
  n = length(iter.ids)
  p = dplyr::progress_estimated(n)
  for (iter.id in iter.ids) {
    x.iter = dplyr::filter(x, iter==iter.id)
    x.forecast = rbind(data.frame(time=time.forecast,
                                     coal=NA, samp=NA, ILI=NA,
                                     coal.E=1, samp.E=1, ILI.E=1,
                                     iter=x.iter$iter[1]), x.iter)
    x.forecast = x.forecast[order(x.forecast[,'iter'], x.forecast[,'time']),]
    invisible({
      forecast.inla = eval(parse(text=paste0('forecast.PILAF.', method, '(x.forecast, formula=',
                                           formula, ')')))
    })
    forecast = forecast.inla$summary.fitted.values[1:length(time.forecast),]
    forecast = cbind(time=time.forecast, forecast)
    forecast$iter = iter.id
    forecast.all = rbind(forecast.all, forecast)
    p$pause(0.1)$tick()$print()
  }

  class(forecast.all) = c("forecast", class(forecast.all))
  return(forecast.all)
}

#' Forecast using Only Count
#'
#' @param x A PILAF object with NA for values to forecast.
#' @return A PILAF object with forecasted ILI counts and 95 % BCI.
#' @export
forecast.PILAF.count = function(x, formula=NULL) {
  link = rep(1, nrow(x))
  if(is.null(formula)) {
    formula = as.formula('ILI ~ -1 + f(time, model="ar", order=3)')
  }
  forecast = INLA::inla(formula,
                     family="poisson", data=x,
                     control.predictor=list(compute=T, link=link),
                     E=x$ILI.E)
  return(forecast)
}

#' Forecast using Count and Coalescent events
#'
#' @param x A PILAF object with NA for values to forecast.
#' @param formula A INLA formula.
#' @return A PILAF object with forecasted ILI counts and 95 % BCI.
#' @export
forecast.PILAF.joint = function(x, formula=NULL) {
  n = nrow(x)
  link = c(rep(1, n), rep(2, n))
  E = c(x$ILI.E, x$coal.E)
  time = rep(x$time, 2)
  beta0 = c(rep(0, n), rep(1, n))
  w0 = c(rep(1, n), rep(0, n))
  w = c(rep(0, n), rep(-1, n))
  Y = matrix(NA, nrow=2*n, ncol=2)
  Y[1:n, 1] = x$ILI
  Y[(1:n)+n, 1] = x$coal
  X = list(time=time, time2=time, beta0=beta0, w0=w0, w=w, Y=Y)
  if(is.null(formula)) {
    formula = as.formula(paste0("Y ~ -1 + beta0 + f(time, w0, model='ar', order=3) +",
                                "f(time2, w, copy='time', fixed=F)"))
  }
  forecast = INLA::inla(formula,
                        family=c("poisson", "poisson"), data=X,
                        control.predictor=list(compute=T, link=link),
                        E=X$E)
  return(forecast)
}

#' Plot the Forecast Results
#'
#' Plot the forecasting results
#'
#' @param forecast A forecast object.
#' @param x A PILAF object.
#' @param truth An optional PILAF object that shows the truth.
#' @return A ggplot object for the forecasting results.
#' @export
plot.forecast = function(x, pilaf, truth=NULL) {
  x = with(x, data.frame(time=time, mean=mean, quant0.025=`0.025quant`,
                         quant0.975=`0.975quant`, iter=iter))
  x$type='ILI'
  pilaf = with(pilaf,
           data.frame(time=time, coalescent=coal, sampling=samp, ILI=ILI,
                      iter=iter))
  pilaf = tidyr::gather(pilaf, type, counts, -c(time, iter))
  p = ggplot(data=pilaf) +
    geom_line(aes(x=time, y=counts, group=iter, color=iter), alpha=0.8, size=0.1) +
    geom_line(data=x, aes(x=time, y=mean, group=iter, color=iter), alpha=0.8, size=0.1) +
    geom_ribbon(data=x, aes(x=time, ymin=quant0.025, ymax=quant0.975, group=iter), alpha=0.01) +
    geom_vline(aes(xintercept=max(x$time), linetype='dotted'), alpha=0.1) +
    facet_wrap(~type, scales='free_y', ncol=1) +
    xlab('time to present') +
    theme_classic() +
    theme(axis.ticks.x=element_blank(),
          legend.position='none',
          strip.background=element_blank())
  if(!is.null(truth)) {
    truth$type = 'ILI'
    p = p +
      geom_line(data=truth, aes(x=time, y=ILI, group=iter, color=iter),
                                linetype='dotdash', alpha=0.3, size=0.3)
  }
  return(p)
}


#' Generic function for evaluation
#' @export
evaluate = function(x, ...) {
  UseMethod('evaluate', x)
}

#' Evaluate the Forecast
#'
#' Evaluate the forecast by pointwise mean absolute error and mean relative
#' width.
#' @param x A forecast object.
#' @param truth A data frame containing time, true ILI.
#' @return A data frame containing the absolute error of forecast.
#' @export
#' @import dplyr
evaluate.forecast = function(x, truth) {
  merged = dplyr::left_join(x, truth, by=c('time', 'iter')) %>%
    dplyr::group_by(time) %>%
    dplyr::summarize(MAE = mean(abs(mean - ILI), na.rm=T),
                     MW = mean(`0.975quant` - `0.025quant`, na.rm=T))
  return(merged)
}
