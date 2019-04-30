context("Test data truncation")

library(PILAF)

test_that("Test for truncation", {
  H3N2_tree <- ape::read.nexus(system.file("extdata/H3N2/MCC.tree",
                                           package = "PILAF"))
  H3N2_phylo <- phylodyn::summarize_phylo(H3N2_tree)
  truncation_time <- .08
  H3N2_truncated <- truncate_data(H3N2_tree, truncation_time, include_trunc = T)
  truncation_time <- H3N2_truncated$truncation_time
  H3N2_truncated <- H3N2_truncated$tree
  keep_idx <- which(H3N2_phylo$samp_times >= truncation_time)
  expect_equal(H3N2_truncated$samp_times[1], H3N2_phylo$samp_times[min(keep_idx)] - truncation_time)
  expect_true(H3N2_truncated$coal_times[1] > H3N2_phylo$samp_times[min(keep_idx)])
  expect_equal(sum(H3N2_truncated$n_sampled), 1 + length(H3N2_truncated$coal_times))
})

test_that("Test for computation of truncation time", {
  date <- as.character(compute_truncation_time(2019, 16))
  expect_equal(date, "2019-04-15")
})

test_that("BNPR truncation", {
  H1N1_tree <- ape::read.nexus(system.file("extdata/H1N1/MCC.tree",
                                           package = "PILAF"))
  H1N1_last_time <- 2019.0931506849315
  test <- truncate_BNPR_PS(2014, 45, H1N1_tree, H1N1_last_time, include_trunc = F) # sampling times 2014 - 2014.268, coal times 2013.206 - 2014.265
  # truncation time at 2014.838
  fit_range <- range(BNPR_to_df(test$bnpr, "H3N2_Ne", test$trunc_time)$Time)
  expect_true(abs(fit_range[1] - 2013.206) < 0.1)
  expect_true(abs(fit_range[2] - 2014.268) < 0.1)
})
