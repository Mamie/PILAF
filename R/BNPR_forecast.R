#' Forecast using Bayesian Nonparametric phylodynamics
#'
#' @param data An ape phylo object
#' @param last_time The time of last sampling point
#' @param formula INLA formula to model Ne (y) with respect to time and week
#' @inheritParams phylodyn::BNPR
#' @export
BNPR_forecast <- function (data, last_time, formula, lengthout = 100, pref = FALSE, prec_alpha = 0.01,
          prec_beta = 0.01, beta1_prec = 0.001, fns = NULL, log_fns = TRUE,
          simplify = TRUE, derivative = FALSE, forward = TRUE, pred = 4)
{
  # browser()
  if (class(data) == "phylo") {
    phy <- phylodyn::summarize_phylo(data)
  }
  else if (all(c("coal_times", "samp_times", "n_sampled") %in%
               names(data))) {
    phy <- with(data, list(samp_times = samp_times, coal_times = coal_times,
                           n_sampled = n_sampled))
  }
  result <- infer_coal_samp_pred(samp_times = phy$samp_times, coal_times = phy$coal_times,
                                 last_time = last_time, n_sampled = phy$n_sampled,
                                 formula = formula, fns = fns, lengthout = lengthout,
                            prec_alpha = prec_alpha, prec_beta = prec_beta, beta1_prec = beta1_prec,
                            use_samp = pref, log_fns = log_fns, simplify = simplify,
                            derivative = derivative, pred = pred)

  result$samp_times <- phy$samp_times
  result$n_sampled <- phy$n_sampled
  result$coal_times <- phy$coal_times

  weeks <- result$grid_week$train$week
  if (pred > 0) weeks <- c(weeks, result$grid_week$test$week)
  result$weeks <- weeks

  result$effpop <- result$result$summary.random$time$`0.5quant`
  result$effpopmean <- result$result$summary.random$time$mean
  result$effpop975 <- result$result$summary.random$time$`0.025quant`
  result$effpop025 <- result$result$summary.random$time$`0.975quant`


  if (grepl("week", as.character(formula)[3])) {
    effpop_map <- with(result$result$summary.random$week,
                       hashmap::hashmap(ID, `0.5quant`))
    effpopmean_map <- with(result$result$summary.random$week,
                           hashmap::hashmap(ID, mean))
    effpop975_map <- with(result$result$summary.random$week,
                          hashmap::hashmap(ID, `0.975quant`))
    effpop025_map <- with(result$result$summary.random$week,
                          hashmap::hashmap(ID, `0.025quant`))
    effpop_map <- with(result$result$summary.random$week,
                       hashmap::hashmap(ID, `0.5quant`))

    result$effpop <- result$effpop + purrr::map_dbl(weeks, ~effpop_map[[.x]])
    result$effpopmean <- result$effpopmean + purrr::map_dbl(weeks, ~effpopmean_map[[.x]])
    result$effpop975 <- result$effpop975 + purrr::map_dbl(weeks, ~effpop025_map[[.x]])
    result$effpop025 <- result$effpop025 + purrr::map_dbl(weeks, ~effpop975_map[[.x]])
  }

  if (grepl("seasonal", as.character(formula)[3])) {
    result$effpop <- result$effpop + result$result$summary.random$seasonal$`0.5quant`
    result$effpopmean <- result$effpopmean + result$result$summary.random$seasonal$mean
    result$effpop975 <- result$effpop975 + result$result$summary.random$seasonal$`0.025quant`
    result$effpop025 <- result$effpop025 + result$result$summary.random$seasonal$`0.975quant`
  }

  result$effpop <- exp(-result$effpop)
  result$effpopmean <- exp(-result$effpopmean)
  result$effpop975 <- exp(-result$effpop975)
  result$effpop025 <- exp(-result$effpop025)


  # result$summary <- with(result,
  #                        data.frame(time = ID, mean = effpopmean, sd = sd * effpopmean,
  #                                   quant0.025 = exp(-`0.975quant`), quant0.5 = exp(-`0.5quant`),
  #                                   quant0.975 = exp(-`0.025quant`)))

  if (derivative) {
    if (forward)
      ind <- c(1:(lengthout - 1), (lengthout - 1))
    else ind <- c(1, 1:(lengthout - 1))
    result$derivative <- with(result$result$summary.lincomb,
                              data.frame(time = result$x, mean = -mean[ind], sd = sd[ind],
                                         quant0.025 = -`0.975quant`[ind], quant0.5 = -`0.5quant`[ind],
                                         quant0.975 = -`0.025quant`[ind]))
  }
  if (pref) {
    result$beta0 <- result$result$summary.fixed["beta0",
                                                "0.5quant"]
    result$beta0summ <- result$result$summary.fixed["beta0",
                                                    ]
    rownames(result$beta0summ) <- "Beta0"
    result$beta1 <- result$result$summary.hyperpar[2, "0.5quant"]
    result$beta1summ <- result$result$summary.hyperpar[2,
                                                       ]
    rownames(result$beta1summ) <- "Beta1"
  }
  return(result)
}

lastmon <- function(x) 7 * floor(as.numeric(x-1+4)/7) + as.Date(1-4, origin="1970-01-01")

time2grid <- function(time) {
  time <- sort(time)
  diff <- abs(time[2] - time[1]) / 2
  return(c(time - diff, time[length(time)] + diff))
}

#' Create grid of time and week
#'
#' @param samp_times Sampling times; the last sampling time is set to 0
#' @param coal_times Coalescent times; the same reference as samp_times
#' @param last_time The last time point
#' @param pred Integer; the number of week to forecast ahead
#' @return A list containing the training and testing time and week
#' @export
create_grid <- function(samp_times, coal_times, last_time, pred){
  # browser()
  grid_week <- list()
  samp_times_real <- last_time - samp_times
  coal_times_real <- last_time - coal_times
  min_time <- min(coal_times_real)
  #max_time <- max(samp_times_real)
  min_date <- as.character(lubridate::round_date(lubridate::date_decimal(min_time)))
  max_date <- as.character(lubridate::round_date(lubridate::date_decimal(last_time)))
  min_mon_date <- lastmon(as.Date(min_date))
  max_mon_date <- lastmon(as.Date(max_date))

  train_time <- seq(min_mon_date, max_mon_date, by = 7)
  train_week <- rev(as.numeric(strftime(train_time, format = "%V")))
  train_grid_0 <- rev(last_time - lubridate::decimal_date(train_time))
  if (pred == 0) {
    test_time_0 <- NULL
    test_week <- NULL
  } else {
    test_time <- seq(max_mon_date + 7, max_mon_date + 7 * pred, by = 7)
    test_week <- rev(as.numeric(strftime(test_time, format = "%V")))
    test_time_0 <- rev(last_time - lubridate::decimal_date(test_time))
  }

  grid_week$train$grid <- train_grid_0
  grid_week$train$week <- train_week
  grid_week$test$time <- test_time_0
  grid_week$test$week <- test_week

  return(grid_week)
}


infer_coal_samp_pred <- function (samp_times, coal_times, last_time, n_sampled = NULL, formula = formula, fns = NULL,
          lengthout = 100, prec_alpha = 0.01, prec_beta = 0.01, beta1_prec = 0.001,
          use_samp = FALSE, log_fns = TRUE, simplify = FALSE, events_only = FALSE,
          derivative = FALSE, pred = 4)
{
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop("INLA needed for this function to work. Use install.packages(\"INLA\", repos=\"https://www.math.ntnu.no/inla/R/stable\").",
         call. = FALSE)
  }
  if (min(coal_times) < min(samp_times))
    stop("First coalescent time occurs before first sampling time")
  if (max(samp_times) > max(coal_times))
    stop("Last sampling time occurs after last coalescent time")
  #browser()

  grid_week <- create_grid(samp_times, coal_times, last_time, pred = pred)
  grid <- grid_week$train$grid
  week <- grid_week$train$week
  time_pred <- grid_week$test$time
  week_pred <- grid_week$test$week

  if (is.null(n_sampled))
    n_sampled <- rep(1, length(samp_times))
  coal_data <- phylodyn:::coal_stats(grid = time2grid(grid), samp_times = samp_times,
                          n_sampled = n_sampled, coal_times = coal_times)
  coal_data <- with(coal_data, phylodyn:::condense_stats(time = time,
                                                event = event, E = E))
  coal_data$week <- week
  coal_data$seasonal <- coal_data$time
  hyper <- list(prec = list(param = c(prec_alpha, prec_beta)))
  if (!use_samp) {
    data <- with(coal_data,
                 data.frame(y = event, time = time, E_log = E_log,
                            week = week, seasonal = seasonal))

    if (pred != 0) {
      data <- rbind(data, data.frame(y = rep(NA, pred),
                                     time = time_pred,
                                     E_log = rep(NA, pred),
                                     week = week_pred,
                                     seasonal = time_pred))
    }
    family <- "poisson"
  }
  else if (use_samp) {
    if (events_only)
      samp_data <- samp_stats(grid = time2grid(grid), samp_times = samp_times)
    else samp_data <- samp_stats(grid = time2grid(grid), samp_times = samp_times, n_sampled = n_sampled)
    data <- joint_stats(coal_data = coal_data, samp_data = samp_data, pred = pred, time_pred = time_pred, week_pred = week_pred)
    #View(data)
    family <- c("poisson", "poisson")
  }
  else stop("Invalid use_samp value, should be boolean.")
  if (derivative) {
    Imat <- diag(lengthout)
    A <- utils::head(Imat, -1) - utils::tail(Imat, -1)
    field <- grid[-1] - diff(grid)/2
    A <- diag(1/diff(field)) %*% A
    A[A == 0] <- NA
    lc_many <- INLA::inla.make.lincombs(time = A)
  }
  else {
    lc_many <- NULL
  }

  mod <- INLA::inla(formula, family = family, data = data,
                    lincomb = lc_many, offset = data$E_log,
                    control.predictor = list(compute = TRUE),
                    control.inla = list(lincomb.derived.only = FALSE))
  return(list(result = mod, data = data, grid = grid, x = coal_data$time, grid_week = grid_week))
}

samp_stats <- function (grid, samp_times, n_sampled = NULL, trim_end = FALSE) {
  lengthout <- length(grid) - 1
  field <- grid[-1] - diff(grid)/2
  E <- diff(grid)
  bins <- cut(x = samp_times, breaks = grid, include.lowest = TRUE)
  if (is.null(n_sampled))
    count <- as.vector(table(bins))
  else {
    tab <- stats::aggregate(n_sampled ~ bins, FUN = sum,
                            labels = FALSE)
    count <- rep(0, lengthout)
    count[as.numeric(tab$bins)] <- tab$n_sampled
  }
  count[utils::head(grid, -1) >= max(samp_times)] <- NA
  result <- data.frame(time = field, count = count, E = E,
                       E_log = log(E))
  if (trim_end)
    result <- result[stats::complete.cases(result), ]
  return(result)
}

joint_stats <- function (coal_data, samp_data, pred = 0,
                         time_pred = NULL, week_pred = NULL) {
  n1 <- length(coal_data$time)
  n2 <- length(samp_data$time)
  beta0 <- c(rep(0, n1 + pred), rep(1, n2))
  E_log <- c(rep(NA, pred), coal_data$E_log, samp_data$E_log)
  Y <- matrix(c(rep(NA, pred), coal_data$event, rep(NA, n2),
                rep(NA, n1 + pred), samp_data$count),
              nrow = n1 + n2 + pred, byrow = FALSE)
  w <- c(rep(1, n1 + pred), rep(-1, n2))
  time <- c(time_pred, coal_data$time, rep(NA, n2))
  time2 <- c(rep(NA, n1 + pred), samp_data$time)
  week <- c(week_pred, coal_data$week, rep(NA, n2))
  week2 <- c(rep(NA, n1 + pred), coal_data$week)
  return(list(Y = Y, beta0 = beta0, time = time, time2 = time2,
              week = week, week2 = week2, w = w, E_log = E_log))
}

find_nearest_later_time <- function(ref, time, n = 2, thresh = 0.01) {
  dist <- time - ref
  ref <- ref[dist >= 0]
  dist <- dist[dist >= 0]
  res <- ref[order(dist)[1:n]]
  if(sum(abs(res - time) > thresh) > 0) {
    stop(paste0("The nearest timepoints are > ", thresh, " from ", time))
  }
  return(res)
}

#' Truncates a phylogentic tree
#'
#' @param tree a phylo object
#' @param truncation_time The time of last sampling point
#' @param include_trunc If to include the truncation time in the data
#' @param n The number of samples near the truncation time
#' @param thresh Maximum distance from truncation time of the n samples
#' @export
truncate_data <- function(tree, truncation_time, include_trunc = F, n = 2, thresh = 0.01){
  # browser()
  tree_data <- phylodyn::summarize_phylo(tree)
  if (include_trunc) {
    # find the nearest sample greater than the truncation time
    nearest_times <- find_nearest_later_time(tree_data$samp_times, truncation_time, n = n, thresh = thresh)
    truncation_time <- min(nearest_times)
  }

  samps_keep <- tree_data$samp_times >= truncation_time
  tree_data$samp_times <- tree_data$samp_times[samps_keep]

  coals_keep <- tree_data$coal_times > min(tree_data$samp_times)
  tree_data$coal_times <- tree_data$coal_times[coals_keep]

  numcoals_removed <- sum(!coals_keep)
  numsamps_removed <- sum(tree_data$n_sampled[!samps_keep])

  tree_data$n_sampled <- tree_data$n_sampled[samps_keep]
  tree_data$n_sampled[1] <- tree_data$n_sampled[1] + numsamps_removed - numcoals_removed

  tree_data$samp_times <- tree_data$samp_times - truncation_time
  tree_data$coal_times <- tree_data$coal_times - truncation_time

  return(list(tree = tree_data, truncation_time = truncation_time))
}

#' Compute the date of the Monday of the week for truncation
#'
#' @param year The year of the truncation time
#' @param week The week of the truncation time
#' @export
compute_truncation_time <- function(year, week) {
  date <- lubridate::ymd(paste0(year, "-01-01"))
  lubridate::week(date) <- week
  return(lastmon(date))
}

#' Forecast starting at given time
#' @param tree A phylo object
#' @param last_time The last sampling time in the tree
#' @param week_start The week of first forecast timepoint
#' @param year_start The year of first forecast timepoint
#' @param formula The formula for BNPR forecast
#' @param label A character string as label for the time series
#' @param pred The number of weeks to forecast
#' @param pref Whether using preferential sampling
#' @param include_trunc Whether to include the truncation timepoint
#' @param n The number of samples near the truncation time
#' @param thresh Maximum distance from truncation time of the n samples
#' @return A INLA object for forecast
#' @export
forecast_starting <- function(tree, last_time, week_start, year_start,
                              formula, label, pred = 4, pref = TRUE, n = 1, include_trunc = T, thresh = 0.019) {
  #browser()
  truncation_time <- lubridate::decimal_date(compute_truncation_time(year_start, week_start))
  tree_trunc <- truncate_data(tree, last_time - truncation_time, n = n, include_trunc = include_trunc, thresh = thresh)
  forecast_res <- BNPR_forecast(tree_trunc$tree, last_time = last_time - tree_trunc$truncation_time,
                                formula = formula, pred = pred, pref = pref)
  forecast_res$truncation_time <- last_time - tree_trunc$truncation_time
  #print(forecast_res$truncation_time)
  forecast_df <- BNPR_to_df(forecast_res, label = label, last_time - tree_trunc$truncation_time)
  return(list(res = forecast_res, df = forecast_df, pred = pred))
}

#' Fit BNPR_PS with truncated data
#'
#' @param year_start The year of the truncation time
#' @param week_start The week of the truncation time
#' @param tree The phylo object
#' @param last_time The last sampling time
#' @param include_trunc Whether to include the truncation timepoint
#' @param n The number of samples near the truncation time
#' @param thresh Maximum distance from truncation time of the n samples
#' @param ... arguments to BNPR_PS
#' @export
truncate_BNPR_PS <- function(year_start, week_start, tree, last_time, n = 1, include_trunc = T, thresh = 0.019, ...) {
 # browser()
  truncation_time <- lubridate::decimal_date(compute_truncation_time(year_start, week_start))
  tree_trunc <- truncate_data(tree, last_time - truncation_time, n = n, include_trunc = include_trunc, thresh = thresh)
  BNPR_PS_trunc <- phylodyn::BNPR_PS(tree_trunc$tree, ...)
  return(list(bnpr = BNPR_PS_trunc, trunc_time = last_time - tree_trunc$truncation_time))
}
