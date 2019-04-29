library(phylodyn)
library(dplyr)
library(ggplot2)
library(PILAF)

H3N2_tree <- ape::read.nexus(system.file("extdata/H3N2/MCC.tree",
                                         package = "PILAF"))
H3N2_last_time <- 2019.0739726027398
H3N2_root_time <- H3N2_last_time - max(phylodyn::summarize_phylo(H3N2_tree)$coal_times)
pred <- 20
BNPR_PS_formula <- Y ~ -1 + beta0 +
  f(time, model="rw1",hyper = hyper, constr = FALSE) +
  f(time2, w, copy="time", fixed=FALSE, param=c(0, beta1_prec))
H3N2_Ne_forecast <- BNPR_forecast(H3N2_tree, last_time = H3N2_last_time, formula = BNPR_PS_formula, pred = pred, pref = T)

tests <- purrr::map(seq(100, 500, by = 100), ~BNPR_PS(H3N2_tree, lengthout = .x))
outs <- purrr::imap(tests, ~data.frame(time = H3N2_last_time- .x$summary$time,
                                      effpop = .x$effpop,
                                      n = .y))
outs <- dplyr::bind_rows(outs)
outs$n <- as.factor(outs$n)
ggplot(data = outs) +
  geom_line(aes(x = time, y = effpop, color = n, group = n)) +
  theme_classic() +
  geom_line(data = data.frame(time = H3N2_last_time - H3N2_Ne_forecast$result$summary.random$time$ID,
                               effpop = H3N2_Ne_forecast$effpop),
             aes(x = time, y = effpop))

# forecast starting debug
year_start <- 2015
week_start <- 45
forecast_2015week46 <- forecast_starting(H3N2_tree, last_time = H3N2_last_time, week_start, year_start, label = "H3N2 BNPR PS", formula = BNPR_PS_formula, pred = 5)

truncation_time <- lubridate::decimal_date(compute_truncation_time(year_start, week_start))
tree_trunc <- truncate_data(H3N2_tree, H3N2_last_time - truncation_time)
BNPR_PS_trunc <- BNPR_PS(tree_trunc)

plot(forecast_2015week46$df$Time, forecast_2015week46$res$effpop,type="l",col="black")
points(truncation_time - BNPR_PS_trunc$summary$time, BNPR_PS_trunc$effpop, type="l", col="red")
