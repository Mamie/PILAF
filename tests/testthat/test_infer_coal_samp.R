context("Test grid computation")

library(PILAF)

test_that("Simulation code and forecasting runs", {
  samp_times <- lubridate::decimal_date(as.Date(c("2019-01-01", "2018-12-18")))
  coal_times <- lubridate::decimal_date(as.Date(c("2017-12-12")))
  last_time <- max(samp_times)
  samp_times <- last_time - samp_times
  coal_times <- last_time - coal_times
  res <- create_grid(samp_times, coal_times, last_time, 1)

})

test_that("A single simulation iteration does not give NA in ILI", {
  n = 1
  lim = c(0, 100)
  flu.Ne = PILAF:::YearlyLogisticTraj
  flu.sampNum = 100
  ILI.sampNum = 1000
  sim.data = PILAF::SimulateILISampCoalCounts(lim, flu.Ne, flu.sampNum, ILI.sampNum)
})
