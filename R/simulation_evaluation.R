RAE <- function(x, y) {
  abs(x - y)/y
}

coverage <- function(l, u, y) {
  mapply(function(lower, upper, truth) {
    if (truth <= upper & truth >= lower) return(1)
    else return(0)
  }, l, u, y)
}

RW <- function(l, u, y) {
  (u - l)/y
}

plot_MRE <- function(forecast_res, truth, color) {
  forecast_res <- select(forecast_res, iter, time, mean, type)
  truth <- select(truth, iter, time, ILI)
  left_join(forecast_res, truth, by = c('time', 'iter')) %>%
    mutate(RAE = RAE(mean, ILI)) %>%
    group_by(time, type) %>%
    summarize(mean = mean(RAE), quant0.25 = quantile(RAE, 0.25), quant0.75 = quantile(RAE, 0.75)) %>%
    ggplot(data = ., aes(x = time, y = mean, color = type)) +
      geom_linerange(aes(x = time, ymin = quant0.25, ymax = quant0.75), size = 0.2, position = position_dodge(width=0.5)) +
      geom_point(size = 3, alpha = 0.9, position = position_dodge(width=0.5)) +
      geom_line(size = 0.3) +
      scale_color_manual(values = color) +
      theme_classic() +
      ylab('RAE') +
      scale_y_continuous(limits = c(0, 1)) +
      scale_x_reverse(breaks = c(-2, -4, -6, -8, -10)) +
      theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
            legend.position = 'top')
}


plot_coverage <- function(forecast_res, truth, color) {
  forecast_res <- forecast_res %>%
    select(iter, time, quant0.025, quant0.975, type)
  truth <- truth %>%
    select(iter, time, ILI)
  left_join(forecast_res, truth, by = c('time', 'iter')) %>%
    mutate(coverage = coverage(quant0.025, quant0.975, ILI)) %>%
    group_by(time, type) %>%
    summarize(mean = mean(coverage)) %>%
    ggplot(data = ., aes(x = time, y = mean, color = type)) +
    geom_point(size = 3, alpha = 0.9, position = position_dodge(width=0.5)) +
    geom_line(size = 0.3) +
    theme_classic() +
    ylab('Coverage') +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_reverse(breaks = c(-2, -4, -6, -8, -10)) +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
    scale_color_manual(values = color) +
    theme(legend.position = 'none')
}

plot_MRW <- function(forecast_res, truth, color) {
  forecast_res <- forecast_res %>%
    select(iter, time, quant0.025, quant0.975, type)
  truth <- truth %>%
    select(iter, time, ILI)
  left_join(forecast_res, truth, by = c('time', 'iter')) %>%
    mutate(RW = RW(quant0.025, quant0.975, ILI)) %>%
    group_by(time, type) %>%
    summarize(mean = mean(RW), quant0.25 = quantile(RW, 0.25), quant0.75 = quantile(RW, 0.75)) %>%
    ggplot(data = ., aes(x = time, y = mean, color = type)) +
    geom_line(size = 0.3) +
    geom_linerange(aes(x = time, ymin = quant0.25, ymax = quant0.75), size = 0.2, position = position_dodge(width=0.5)) +
    geom_point(size = 3, alpha = 0.9, position = position_dodge(width=0.5)) +
    theme_classic() +
    ylab('RW') +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_reverse(breaks = c(-2, -4, -6, -8, -10)) +
    scale_color_manual(values = color) +
    theme(legend.position = 'none')
}
