context("Test simulation of coalescent, sampling, ILI counts")

library(PILAF)

test_that("Simulation code runs", {
  n = 2
  lim = c(0, 51)
  flu.Ne = PILAF:::YearlyLogisticTraj
  flu.sampNum = 300
  ILI.sampNum = 1000
  res = PILAF::SimulateILISampCoalCountsN(n, lim, flu.Ne, flu.sampNum, ILI.sampNum)
  expect_equal_to_reference(res, "test_simulation.rds")
})
