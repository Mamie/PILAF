context("Test data truncation")

library(PILAF)

test_that("Test for truncation", {
  H3N2_tree <- ape::read.nexus(system.file("extdata/H3N2/MCC.tree",
                                           package = "PILAF"))
  H3N2_phylo <- phylodyn::summarize_phylo(H3N2_tree)
  truncation_time <- .08
  H3N2_truncated <- truncate_data(H3N2_tree, truncation_time)
  keep_idx <- which(H3N2_phylo$samp_times >= truncation_time)
  expect_equal(H3N2_truncated$samp_times[1], H3N2_phylo$samp_times[min(keep_idx)] - truncation_time)
  expect_true(H3N2_truncated$coal_times[1] > H3N2_phylo$samp_times[min(keep_idx)])
  expect_true(sum(H3N2_truncated$n_sampled), 1 + length(H3N2_truncated$coal_times))
})

test_that("Test for computation of truncation time", {
  date <- as.character(compute_truncation_time(2019, 16))
  expect_equal(date, "2019-04-15")
})
