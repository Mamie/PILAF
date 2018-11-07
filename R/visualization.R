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
