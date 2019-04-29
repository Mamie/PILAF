context("Test grid computation")

library(PILAF)

test_that("Test for Mondays", {
  samp_times <- lubridate::decimal_date(as.Date(c("2019-04-29", "2019-04-22")))
  coal_times <- lubridate::decimal_date(as.Date(c("2019-04-15")))
  last_time <- max(samp_times)
  samp_times <- last_time - samp_times
  coal_times <- last_time - coal_times
  res <- create_grid(samp_times, coal_times, last_time, 1)
  expect_equal(res$train$week, c(18, 17, 16))
  expect_equal(res$test$week, 19)
  expect_equal(res$train$grid, c(samp_times, coal_times))
  expect_equal(res$test$time, lubridate::decimal_date(as.Date("2019-04-29")) -
                 lubridate::decimal_date(as.Date("2019-05-06")))
})

test_that("Test for other dates in a week", {
  samp_times <- lubridate::decimal_date(as.Date(c("2019-04-30", "2019-04-25")))
  coal_times <- lubridate::decimal_date(as.Date(c("2019-04-15")))
  last_time <- max(samp_times)
  samp_times <- last_time - samp_times
  coal_times <- last_time - coal_times
  mondays <- lubridate::decimal_date(as.Date(c("2019-04-29", "2019-04-22", "2019-04-15")))
  res <- create_grid(samp_times, coal_times, last_time, 1)
  expect_equal(res$train$week, c(18, 17, 16))
  expect_equal(res$test$week, 19)
  expect_equal(res$train$grid, last_time - mondays)
  expect_equal(res$test$time, last_time - lubridate::decimal_date(as.Date("2019-05-06")))
})

test_that("Test creating grid for given timepoint of interest", {
  test <- c(1, 2, 3, 4)
  grid <- time2grid(test)
  expect_equal(grid, c(0.5, 1.5, 2.5, 3.5, 4.5))
})

