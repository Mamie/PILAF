flu.Ne = PILAF:::YearlyLogisticTraj
lim = c(0, 52) # start and end time
grid = seq(lim[1] - 0.5, lim[2] + 0.5, by=1)

flu.sampNum = 100
ILI.sampNum = 1000

res = PILAF:::SimulateILISampCoalCounts(lim, flu.Ne, flu.sampNum, ILI.sampNum)

library(ggplot2)

data.frame(time = res$coal$time,
           `coalescent counts`=res$coal$event,
           `sampling counts`=res$samp$count,
           `ILI counts`=res$ILI$count) %>%
  tidyr::gather(type, counts, -time) %>%
  ggplot(data=.) +
    geom_line(aes(x=time, y=counts, group=type, color=type), size=0.1) +
    facet_wrap(~type, scales='free', ncol=1) +
    theme_classic() +
    theme(axis.title=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position='none')
