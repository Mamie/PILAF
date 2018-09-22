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


#' Simulate ILI, Sampling, and Coalescent counts
#'
#' Simulate flu sampling, coalescence and ILI counts given effective population
#' size trajectory and number of expected samples.
#'
#' @param lim A numeric tuple of start and end times.
#' @param flu.Ne A function of time (for flu effective popluation size).
#' @param flu.sampNum A numeric scalar of expected number of flu samples in lim.
#' @param ILI.sampNum A numeric scalar of expected number of sampling events.
#' in lim
#' @return A list containing coalescent, sampling and ILI counts
#' @export
SimulateILISampCoalCounts = function(lim, flu.Ne, flu.sampNum, ILI.sampNum) {
  flu.c = PILAF::ComputePropCnst(flu.sampNum, flu.Ne, lim)
  flu.sampTimes = PILAF::PoissonTimeSampling(flu.Ne, flu.c, lim)
  flu.nsamp = rep(1, length(flu.sampTimes))

  ILI.c = PILAF:::ComputePropCnst(ILI.sampNum, flu.Ne, lim)
  ILI.sampTimes = PILAF::PoissonTimeSampling(flu.Ne, ILI.c, lim)
  ILI.nsamp = rep(1, length(ILI.sampTimes))

  grid = seq(lim[1] - 0.5, lim[2] + 0.5, by=1)
  coal.counts = PILAF::SimulateCoalCounts(grid, flu.sampTimes, flu.nsamp, flu.Ne)
  samp.counts = phylodyn:::samp_stats(grid, flu.sampTimes, flu.nsamp)
  ILI.counts = phylodyn:::samp_stats(grid, ILI.sampTimes, ILI.nsamp)
  return(list(coal=coal.counts, samp=samp.counts, ILI=ILI.counts))
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
    sim = SimulateILISampCoalCounts(lim, flu.Ne, flu.sampNum, ILI.sampNum)
    sim = data.frame(time=sim$coal$time, coal=sim$coal$event,
                     samp=sim$samp$count, ILI=sim$ILI$count,
                     coal.E=sim$coal$E, samp.E=sim$samp$E,
                     ILI.E=sim$ILI$E, iter=i)
    sim.all = rbind(sim.all, sim)
    p$pause(0.1)$tick()$print()
  }
  return(sim.all)
}

#' Plot Simulations
#'
#' Plot simulated counts.
#'
#' @param sim A data frame from PILAF::SimulateILISampCoalCountsN
#' @export
#' @import ggplot2
#' @import magrittr
PlotSimulations = function(sim) {
  sim = with(sim,
       data.frame(time = time, coalescent=coal, sampling=samp,ILI=ILI,iter=iter))
  sim = tidyr::gather(sim, type, counts, -c(time, iter))
  ggplot(data=sim) +
    geom_line(aes(x=time, y=counts, group=iter, color=iter, alpha=0.4), size=0.1) +
    facet_wrap(~type, scales='free', ncol=1) +
    theme_classic() +
    theme(axis.ticks.x=element_blank(),
          legend.position='none',
          strip.background=element_blank())
}
