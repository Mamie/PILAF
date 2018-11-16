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
  return(list(p_MRE = p_MRE, p_MRW = p_MRW))
}
