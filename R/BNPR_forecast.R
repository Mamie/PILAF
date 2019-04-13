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
  result$effpop <- exp(-result$result$summary.random$time$`0.5quant`)
  result$effpopmean <- exp(-result$result$summary.random$time$mean)
  result$effpop975 <- exp(-result$result$summary.random$time$`0.025quant`)
  result$effpop025 <- exp(-result$result$summary.random$time$`0.975quant`)
  result$summary <- with(result$result$summary.random$time,
                         data.frame(time = ID, mean = exp(-mean), sd = sd * exp(-mean),
                                    quant0.025 = exp(-`0.975quant`), quant0.5 = exp(-`0.5quant`),
                                    quant0.975 = exp(-`0.025quant`)))

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
create_grid <- function(samp_times, coal_times, last_time, pred){
  grid_week <- list()
  samp_times_real <- last_time - samp_times
  coal_times_real <- last_time - coal_times
  min_time <- min(coal_times_real)
  max_time <- max(samp_times_real)
  min_date <- as.character(lubridate::date_decimal(min_time))
  max_date <- as.character(lubridate::date_decimal(max_time))
  min_mon_date <- lastmon(as.Date(min_date))
  max_mon_date <- lastmon(as.Date(max_date))

  train_time <- seq(min_mon_date, max_mon_date, by = 7)
  train_week <- rev(as.numeric(strftime(train_time, format = "%V")))
  train_grid_0 <- rev(last_time - lubridate::decimal_date(time2grid(train_time)))
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

  grid_week <- create_grid(samp_times, coal_times, last_time, pred = pred)
  grid <- grid_week$train$grid
  week <- grid_week$train$week
  time_pred <- grid_week$test$time
  week_pred <- grid_week$test$week
  # grid <- seq(min(samp_times), max(coal_times), length.out = lengthout + 1)

  if (is.null(n_sampled))
    n_sampled <- rep(1, length(samp_times))
  coal_data <- phylodyn:::coal_stats(grid = grid, samp_times = samp_times,
                          n_sampled = n_sampled, coal_times = coal_times)
  coal_data <- with(coal_data, phylodyn:::condense_stats(time = time,
                                                event = event, E = E))
  hyper <- list(prec = list(param = c(prec_alpha, prec_beta)))
  if (!use_samp) {
    data <- with(coal_data, data.frame(y = event, time = time,
                                       E_log = E_log))
    data <- cbind(data, week = week)
    if (pred != 0) {
      data <- rbind(data, data.frame(y = rep(NA, pred),
                                     time = time_pred,
                                     week = week_pred,
                                     E_log = rep(NA, pred)))
    }

    family <- "poisson"
  }
  else if (use_samp) {
    if (events_only)
      samp_data <- phylodyn:::samp_stats(grid = grid, samp_times = samp_times)
    else samp_data <- phylodyn:::samp_stats(grid = grid, samp_times = samp_times,
                                 n_sampled = n_sampled)
    data <- phylodyn:::joint_stats(coal_data = coal_data, samp_data = samp_data)
    if (is.null(fns)) {
      formula <- Y ~ -1 + beta0 + f(time, model = "rw1",
                                    hyper = hyper, constr = FALSE) + f(time2, w,
                                                                       copy = "time", fixed = FALSE, param = c(0, beta1_prec))
    }
    else {
      vals <- NULL
      bins <- sum(data$beta0 == 0)
      for (fni in fns) {
        if (log_fns)
          vals <- cbind(vals, c(rep(0, bins), log(fni(samp_data$time))))
        else vals <- cbind(vals, c(rep(0, bins), fni(samp_data$time)))
      }
      data$fn <- vals
      formula <- Y ~ -1 + beta0 + fn + f(time, model = "rw1",
                                         hyper = hyper, constr = FALSE) + f(time2, w,
                                                                            copy = "time", fixed = FALSE, param = c(0, beta1_prec))
    }
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
  return(list(result = mod, data = data, grid = grid, x = coal_data$time))
}
