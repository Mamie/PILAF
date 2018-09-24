context("Test simulation of coalescent, sampling, ILI counts")

library(PILAF)

test_that("Simulation code and forecasting runs", {
  n = 1
  lim = c(0, 81)
  flu.Ne = PILAF:::YearlyLogisticTraj
  flu.sampNum = 30
  ILI.sampNum = 1000
  sim.data = PILAF::SimulateILISampCoalCountsN(n, lim, flu.Ne, flu.sampNum, ILI.sampNum)
  plot(sim.data)
  train.idx = 1:(nrow(sim.data)-5)
  sim.train = sim.data[train.idx,]
  sim.test = sim.data[-train.idx,]
  sim.fitted = forecast(sim.train, sim.data$time[nrow(sim.data)-(4:0)])

  expect_equal_to_reference(res, "test_simulation.rds")
})
