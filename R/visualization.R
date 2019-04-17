#' plot real CDC time series data
#'
#' @param data A data frame containing year, week and total ILI counts
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


plot_trajectory <- function(df, ylimits, breaks, root_time) {
  xmin <- c()
  xmax <- c()
  year <- floor(min(df$Time))
  while (year < max(df$Time)) {
    xmin <- c(xmin, year + 0.77)
    xmax <- c(xmax, year + 1.40)
    year <- year + 1
  }
  winter_seasons <- data.frame(xmin = xmin, xmax = xmax) %>%
    mutate(ymin = min(ylimits), ymax = max(ylimits)) %>%
    filter(xmin < max(df$Time) & xmax > min(df$Time))
  print(winter_seasons)
  p <- ggplot(data = df) +
    geom_rect(data = winter_seasons, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = 'steelblue', alpha = 0.2) +
    geom_ribbon(aes(x = Time, ymin = Lower, ymax = Upper), fill = "lightgray", alpha = 0.9) +
    geom_line(aes(x = Time, y = Median)) +
    scale_y_log10(limits = ylimits) +
    theme_classic() +
    ylab(expression(N["e"])) +
    facet_wrap(~Label, ncol = 1) +
    theme(axis.text.x = element_text(size = 13),
          axis.title.x = element_blank(),
          strip.background = element_blank()) +
    scale_x_continuous(limits = c(root_time, max(winter_seasons$xmax, df$Time)),
                       breaks = breaks)
  return(p)
}


#' Plot MCC tree and skygrid effective population size
#' @param MCC A nexus object
#' @param root_time The actual time for the root of the tree
#' @param last_time The last time to be plotted
#' @param Ne A list of data frame containing Ne summary statistics
#' @export
plot_MCC_Ne <- function(MCC, root_time, last_time, Ne, breaks, main = NULL,
                        ratio = c(1.2, 1), mask_range = NULL, ylimits = NULL) {
  par(mar = c(2, 3.8, 2, 1))
  layout(matrix(c(1, 2), ncol = 1), c(4, 4), ratio)

  plot(ape::ladderize(MCC), show.tip.label = F, show.node.label = F, main = main)
  ape::axisPhylo(root.time = root_time, backward = F)
  plot.new()
  vps <- gridBase::baseViewports()
  grid::pushViewport(vps$figure)
  vp1 <- grid::plotViewport(c(1,0.1,0,0.5))
  if (!is.null(mask_range)) Ne <- Ne[!(Ne$Time < max(mask_range) & Ne$Time > min(mask_range)),]
  if (is.null(ylimits)) ylimits <- c(min(Ne$Lower), max(Ne$Upper))
  p <- plot_trajectory(Ne, ylimits, breaks, root_time)
  print(p, vp = vp1)
  dev.new()
  dev.off()
}

#' Convert BNPR object summary statistics to data frame
#' @param BNPR_obj A BNPR object
#' @param label Label to be displayed on the plot
#' @param last_time Last sampling time
#' @export
BNPR_to_df <- function(BNPR_obj, label, last_time) {
  with(BNPR_obj,
       data.frame(Time = x,
                  Mean = effpopmean,
                  Median = effpop,
                  Upper = effpop975,
                  Lower = effpop025,
                  Label = label,
                  stringsAsFactors = F)) %>%
    mutate(Time = last_time - Time)
}
