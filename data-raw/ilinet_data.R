library(cdcfluview)
library(tidyverse)
walk(c("national", "hhs", "census", "state"), ~{

  ili_df <- ilinet(region = .x)

  print(glimpse(ili_df))

  ggplot(ili_df, aes(week_start, unweighted_ili, group=region, color=region)) +
    geom_line() +
    viridis::scale_color_viridis(discrete=TRUE) +
    labs(x=NULL, y="Unweighted ILI", title=ili_df$region_type[1]) +
    theme_ipsum_rc(grid="XY") +
    theme(legend.position = "none") -> gg

  print(gg)

})
