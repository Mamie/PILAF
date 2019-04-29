library(phylodyn)
library(dplyr)
library(ggplot2)
library(PILAF)

H3N2_tree <- ape::read.nexus(system.file("extdata/H3N2/MCC.tree",
                                         package = "PILAF"))
H3N2_last_time <- 2019.0739726027398
H3N2_root_time <- H3N2_last_time - max(phylodyn::summarize_phylo(H3N2_tree)$coal_times)
ar2crw1_formula <- y ~ -1 + f(time, model = "ar", order = 2) +
  f(week, model = "rw1", cyclic = T, constr = F)
pred <- 20
H3N2_Ne_forecast <- BNPR_forecast(H3N2_tree, last_time = H3N2_last_time, formula = ar2crw1_formula, pred = pred)

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

