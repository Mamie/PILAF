context("Test simulation of coalescent, sampling, ILI counts")

library(PILAF)

test_that("Simulation code and forecasting runs", {
  n = 10
  lim = c(0, 100)
  flu.Ne = PILAF:::YearlyLogisticTraj
  flu.sampNum = 100
  ILI.sampNum = 1000
  sim.data = PILAF::SimulateILISampCoalCountsN(n, lim, flu.Ne, flu.sampNum, ILI.sampNum)
  plot(sim.data)
  train.time = seq(0, 90)
  test.time = seq(91, 100)
  sim.train = sim.data[sim.data$time %in% train.time,]
  sim.test = sim.data[!sim.data$time %in% train.time,]
  sim.fitted = forecast(sim.train, test.time)
  AE = evaluate(sim.fitted, sim.test)
  # plot(sim.fitted, sim.train)
  expect_equal_to_reference(sim.fitted, "test_simulation.rds")
})

test_that("A single simulation iteration does not give NA in ILI", {
  n = 1
  lim = c(0, 100)
  flu.Ne = PILAF:::YearlyLogisticTraj
  flu.sampNum = 100
  ILI.sampNum = 1000
  sim.data = PILAF::SimulateILISampCoalCounts(lim, flu.Ne, flu.sampNum, ILI.sampNum)
})
