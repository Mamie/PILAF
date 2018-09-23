#' Poisson Process Simulation
#'
#' Simulation of Poisson process given intensity.
#'
#' @param traj A function of time.
#' @param c A numeric as proportionality constant on traj for intensity.
#' @param lim A tuple that specifies the start and end time.
#' @return A numeric vector of event times.
#' @export
#' @seealso PILAF::ComputePropCnst
PoissonTimeSampling = function(traj, c, lim) {
  upperbound = 1.1 * max(c * traj(seq(lim[1], lim[2], by=1)))
  time = 0
  times = 0
  while (time < lim[2]) {
    y = rexp(1, upperbound)
    if (runif(1) < c * traj(time + y)/upperbound) {
      times = c(times, time)
    }
    time = time + y
  }
  return(times)
}

#' A Logistic Trajectory with Period 52
#'
#' A logistic trajectory (based on `phylodyn::logistic_traj`) with period 52.
#'
#' @inheritParams phylodyn::logistic_traj
YearlyLogisticTraj = function(...) {
  phylodyn::logistic_traj(..., period=52)
}

#' Compute Proportionality Constant
#'
#' Compute proportionality constant to produce a correct expected number of
#' samples
#' @param n Number of samples.
#' @param traj A function of time (intensity).
#' @param lim A numeric tuple of start and end time.
#' @return A numeric scalar as proportionality constant.
#' @export
ComputePropCnst = function(n, traj, lim) {
  n/stats::integrate(traj, lim[1], lim[2])$value
}


#' Simulate Coalescent Counts
#'
#' Simulate coalescent counts given sampling events and intensity.
#'
#' @param grid A numeric vector of grid.
#' @param samp.times A numeric vector of sampling times.
#' @param samp.num A numeric vector of number of samples per time.
#' @param traj A function of time.
#' @return A data frame of coalescent counts.
#' @export
#' @seealso phylodyn::coalsim
SimulateCoalCounts = function(grid, samp.times, samp.num, traj) {
  coal.events = phylodyn::coalsim(samp_times=samp.times,
                                  n_sampled=samp.num,
                                  traj=traj)
  coal.counts = phylodyn:::coal_stats(grid,
                                      coal.events$samp_times,
                                      coal.events$coal_times,
                                      coal.events$n_sampled)
  coal.counts = phylodyn:::condense_stats(coal.counts$time,
                                          coal.counts$event,
                                          coal.counts$E)
  return(coal.counts)
}

#' Convert Ordered Time as Time to Last Time Point
#'
#' Convert ordered time as time to last timepoint.
#'
#' @param time A numeric vector of time.
#' @return A list consisting of a numeric vector of time with respect to last
#' time point and last time point.
AsBackwardTime = function(time) {
  last.time = max(time)
  backward.time = last.time - time
  return(list(time=backward.time, last.time=last.time))
}

#' Convert Time to Last Timepoint as Normal Time
#'
#' Convert time to last timepoint as normal time.
#'
#' @param time A numeric vector of time.
#' @param last.time Actual time of last timepoint
#' @return A vector of time in normal order
AsForwardTime = function(time, last.time) {
  forward.time = last.time - time
  return(forward.time)
}


#' Simulate ILI, Sampling, and Coalescent counts
#'
#' Simulate flu sampling, coalescence and ILI counts given effective population
#' size trajectory and number of expected samples.
#'
#' @param lim A numeric tuple of start and end times.
#' @param flu.Ne A function of time (for flu effective popluation size).
#' @param flu.sampNum A numeric scalar of expected number of flu samples in lim.
#' @param ILI.sampNum A numeric scalar of expected number of sampling events.
#' @param i An optional numeric scalar indicating iteration of simulation
#' in lim
#' @return A list containing coalescent, sampling and ILI counts
#' @export
SimulateILISampCoalCounts = function(lim, flu.Ne, flu.sampNum, ILI.sampNum, i=1) {
  flu.c = PILAF::ComputePropCnst(flu.sampNum, flu.Ne, lim)
  flu.sampTimes = PILAF::PoissonTimeSampling(flu.Ne, flu.c, lim)
  flu.nsamp = rep(1, length(flu.sampTimes))

  ILI.c = PILAF:::ComputePropCnst(ILI.sampNum, flu.Ne, lim)
  ILI.sampTimes = PILAF::PoissonTimeSampling(flu.Ne, ILI.c, lim)
  ILI.nsamp = rep(1, length(ILI.sampTimes))

  lim.backward = PILAF:::AsBackwardTime(lim)
  grid = seq(min(lim.backward$time) - 0.5, max(lim.backward$time) + 0.5, by=1)
  coal = PILAF::SimulateCoalCounts(grid, flu.sampTimes, flu.nsamp, flu.Ne)
  samp = phylodyn:::samp_stats(grid, flu.sampTimes, flu.nsamp)
  ILI = phylodyn:::samp_stats(grid, ILI.sampTimes, ILI.nsamp)
  time.forward = PILAF:::AsForwardTime(coal$time, lim.backward$last.time)
  pilaf = PILAF::PILAF(time=time.forward, coal=coal$event,
                     samp=samp$count, ILI=ILI$count,
                     coal.E=coal$E, samp.E=samp$E,
                     ILI.E=ILI$E, iter=i)
  return(pilaf)
}

#' Simulate ILI, Sampling, and Coalescent counts N Times
#'
#' Simulate ILI, sampling and coalescent event counts for n Times
#'
#' @param n Number of simulations.
#' @param lim A numeric tuple of start and end times.
#' @param flu.Ne A function of time (for flu effective popluation size).
#' @param flu.sampNum A numeric scalar of expected number of flu samples in lim.
#' @param ILI.sampNum A numeric scalar of expected number of sampling events.
#' @export
SimulateILISampCoalCountsN = function(n, lim, flu.Ne, flu.sampNum, ILI.sampNum, seed=1) {
  set.seed(seed)
  sim.all = c()
  p <- dplyr::progress_estimated(n)
  for (i in seq(n)) {
    sim = SimulateILISampCoalCounts(lim, flu.Ne, flu.sampNum, ILI.sampNum, i=i)
    sim.all = rbind(sim.all, sim)
    p$pause(0.1)$tick()$print()
  }
  class(sim.all) = "PILAF"
  return(sim.all)
}
