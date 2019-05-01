library(phylodyn)
library(dplyr)
library(ggplot2)
library(PILAF)

H3N2_tree <- ape::read.nexus(system.file("extdata/H3N2/MCC.tree",
                                         package = "PILAF"))
H3N2_last_time <- 2019.0739726027398
H3N2_root_time <- H3N2_last_time - max(phylodyn::summarize_phylo(H3N2_tree)$coal_times)
pred <- 20
rw1_formula <- Y ~ -1 + beta0 +
  f(time, model="rw1",hyper = hyper, constr = FALSE) +
  f(time2, w, copy="time", fixed=FALSE, param=c(0, beta1_prec))
ar2_formula <- Y ~ -1 + beta0 +
  f(time, model="ar", order = 2, hyper = hyper, constr = FALSE) +
  f(time2, w, copy="time", fixed=FALSE, param=c(0, beta1_prec))

H3N2_Ne_forecast <- BNPR_forecast(H3N2_tree, last_time = H3N2_last_time, formula = rw1_formula, pred = pred, pref = T)

# test BNPR PS with lengthout = 100, 200, ..., 500
lengthout <- seq(100, 500, by = 100)
tests <- purrr::map(lengthout, ~BNPR_PS(H3N2_tree, lengthout = .x))
outs <- purrr::map2(tests, lengthout, ~data.frame(time = H3N2_last_time- .x$summary$time,
                                      effpop = .x$effpop,
                                      n = .y))
outs <- dplyr::bind_rows(outs)
outs$n <- as.factor(outs$n)
# black line is BNPR PS forecast and colored lines are BNPR PS with different lengthout parameters
# BNPR PS forecast is consistent with BNPR PS in fitting
ggplot(data = outs) +
  geom_line(aes(x = time, y = effpop, color = n, group = n)) +
  theme_classic() +
  geom_line(data = data.frame(time = H3N2_last_time - H3N2_Ne_forecast$result$summary.random$time$ID,
                               effpop = H3N2_Ne_forecast$effpop),
             aes(x = time, y = effpop))

# forecast starting debug
year_start <- 2015
week_start <- 45
forecast_2015week46_rw1 <- forecast_starting(H3N2_tree, last_time = H3N2_last_time, week_start, year_start, label = "H3N2 BNPR PS", formula = rw1_formula, pred = 5)

truncation_time <- lubridate::decimal_date(compute_truncation_time(year_start, week_start))
tree_trunc <- truncate_data(H3N2_tree, H3N2_last_time - truncation_time, include_trunc = T)
BNPR_PS_trunc <- BNPR_PS(tree_trunc$tree, lengthout = 190)

# BNPR PS with truncated data matches with BNPR PS forecast in fitting
plot(forecast_2015week46_rw1$df$Time, forecast_2015week46_rw1$res$effpop,type="l",col="black")
points(truncation_time - BNPR_PS_trunc$summary$time, BNPR_PS_trunc$effpop, type="l", col="red")

forecast_2015week46_ar2 <- forecast_starting(H3N2_tree, last_time = H3N2_last_time, week_start, year_start, label = "H3N2 BNPR PS", formula = ar2_formula, pred = 5)

plot(forecast_2015week46_rw1$df$Time, forecast_2015week46_rw1$res$effpop, type="l",col="black", xlab = "Time", ylab = "Ne")
points(truncation_time - BNPR_PS_trunc$summary$time, BNPR_PS_trunc$effpop, type="l", col="red")
points(forecast_2015week46_ar2$df$Time, forecast_2015week46_ar2$res$effpop, type="l",col="blue")
abline(v = forecast_2015week46_ar2$df$Time[6], col = "gray")

# H1N1 truncate BNPR PS debug
test <- truncate_BNPR_PS(2014, 45, H1N1_tree, H1N1_last_time) # sampling times 2014 - 2014.268, coal times 2013.206 - 2014.265
BNPR_to_df(test$bnpr, "H3N2_Ne", test$trunc_time)$Time # 2013.211 - 2014.263
trunc_tree <- truncate_data(H3N2_tree, H1N1_last_time - test$trunc_time)
range(test$trunc_time - c(trunc_tree$samp_times, trunc_tree$coal_times))
range(test$trunc_time - test$bnpr$result$summary.random$time$ID)

