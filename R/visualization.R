#' plot real CDC time series data
#'
#' @params data A data frame containing year, week and total ILI counts
#' @param value The name of the variable that contains ILI counts
#' @export
plot_traj_heatmap <- function(data, value) {
  data_wide <- data %>%
    mutate(season = ifelse(week >= 30, paste0(year, '-', year + 1),
                           paste0(year - 1, '-', year))) %>%
    select(-year) %>%
    tidyr::spread_('week', value)
  data_wide_mat <- data.matrix(data_wide %>% select(-season)) %>% .[ , c(seq(30, 53), seq(1, 29))]
  rownames(data_wide_mat) <- data_wide$season
  iheatmapr::main_heatmap(data_wide_mat, name = value, colors = 'Blues') %>%
    add_row_labels() %>%
    add_col_labels() %>%
    add_row_summary(summary_function = 'mean')
}


plot_LOSO_performance <- function(performances) {
  perf <- performances %>%
    tidyr::gather(metrics, performance, -c(season, task, model))

  plot_metric <- function(data) {
    ggplot(data) +
      geom_hline(aes(yintercept = task), size = 0.1, color = 'gray') +
      geom_point(aes(x = performance, y = task, color = model), size = 0.5) +
      scale_y_reverse() +
      facet_wrap( ~ season, ncol = 6) +
      theme_classic() +
      theme(strip.background = element_blank(),
            legend.position = 'bottom',
            axis.line = element_blank()) +
      scale_color_manual(values = c("#2F408E", "#E5801C")) +
      scale_x_log10()
  }

  p_MRE <- plot_metric(perf %>%
                         filter(metrics == 'MRE')) + xlab("Mean relative error")

  p_MRW <- plot_metric(perf %>%
                         filter(metrics == 'MRW')) + xlab("Mean relative width")

  p_MCV <- plot_metric(perf %>%
                         filter(metrics == 'MCV')) + xlab("Mean coverage")
  return(list(p_MRE = p_MRE, p_MRW = p_MRW, p_MCV = p_MCV))
}


plot_trajectory <- function(df, estimated_root_height = NULL) {
  xmin <- c()
  xmax <- c()
  year <- floor(min(df$Time))
  while (year < max(df$Time)) {
    xmin <- c(xmin, year + 0.9)
    xmax <- c(xmax, year + 1.2)
    year <- year + 1
  }
  winter_seasons <- data.frame(xmin = xmin, xmax = xmax) %>%
    filter(xmin < max(df$Time) & xmax > min(df$Time))
  p <- ggplot(data = df) +
    geom_rect(data = winter_seasons, aes(xmin = xmin, xmax = xmax, ymin = min(df$Lower), ymax = max(df$Upper)), fill = 'steelblue', alpha = 0.1) +
    geom_line(aes(x = Time, y = Median), size = 1) +
    geom_line(aes(x = Time, y = Upper), size = 0.5) +
    geom_line(aes(x = Time, y = Lower), size = 0.5) +
    scale_y_log10() +
    theme_classic() +
    ylab(expression(N["e"])) +
    theme(axis.text.x = element_text(size = 13), axis.title = element_blank())

  if(!is.null(estimated_root_height)) {
    p <- p + geom_vline(aes(xintercept = estimated_root_height), linetype = "dotted")
  }
  return(p)
}

plot_MCC_Ne <- function(MCC, root_time, last_time, Ne, breaks) {
  par(mar = c(2, 3.8, 0, 1))
  layout(matrix(c(1, 2), ncol = 1), c(4, 4), c(1.2, 1))
  plot(ape::ladderize(MCC), show.tip.label=F, show.node.label=F)
  ape::axisPhylo(root.time = root_time, backward = F)
  plot.new()
  vps <- baseViewports()
  pushViewport(vps$figure)
  vp1 <-plotViewport(c(1,0.1,0,0.5))
  p <- plot_trajectory(Ne) +
    scale_x_continuous(limits = c(root_time, last_time), breaks = breaks)
  print(p, vp = vp1)
}
