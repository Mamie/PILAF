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
                     ILI.E=ILI.E,
                     iter=iter)
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
#' @import magrittr
plot.PILAF = function(x) {
  x = with(x,
          data.frame(time = time, coalescent=coal, sampling=samp,ILI=ILI,iter=iter))
  x = tidyr::gather(x, type, counts, -c(time, iter))
  ggplot(data=x) +
    geom_line(aes(x=time, y=counts, group=iter, color=iter, alpha=0.4), size=0.1) +
    facet_wrap(~type, scales='free', ncol=1) +
    theme_classic() +
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
#' @param time.forecast A numeric vector of timepoints to forecast.
#' @param formula An optional character scalar of formula (debugging).
#' @export
forecast.PILAF = function(x, time.forecast, method='count', formula='NULL') {
  if (!method %in% c('count', 'joint', 'ps')) {
    warnings('Forecast method not recognized. Use default="count"')
    method = 'count'
  }
  x.forecast = rbind(x, data.frame(time=time.forecast,
                                   coal=NA, samp=NA, ILI=NA,
                                   coal.E=1, samp.E=1, ILI.E=1))
  forecast = eval(paste0('forecast.PILAF.', method, '(x.forecast, formula=)', formula))
  return(forecast)
}

forecast.PILAF.count = function(x) {
}
