library(sbmt)
data(la_byhour_edgelist)

K = 3

fit_dir_dc0 = sbmt(la_byhour_edgelist,  degreeCorrect = 0, directed = TRUE, klPerNetwork = 2, maxComms = K, seed = 1)
fit_undir_dc0 = sbmt(la_byhour_edgelist,  degreeCorrect = 0, directed = FALSE, klPerNetwork = 2, maxComms = K, seed = 1)

fit_dir_dc2 = sbmt(la_byhour_edgelist,  degreeCorrect = 2, directed = TRUE, klPerNetwork = 2, maxComms = K, seed = 1)
fit_undir_dc2 = sbmt(la_byhour_edgelist,  degreeCorrect = 2, directed = FALSE, klPerNetwork = 2, maxComms = K, seed = 1)

fit_dir_dc3 = sbmt(la_byhour_edgelist,  degreeCorrect = 3, directed = TRUE, klPerNetwork = 2, maxComms = K, seed = 1)
fit_undir_dc3 = sbmt(la_byhour_edgelist,  degreeCorrect = 3, directed = FALSE, klPerNetwork = 2, maxComms = K, seed = 1)

test_that("test time-dependent and not time-dependent code gives same results when T = 1", {
  expect_equal(fit_dir_dc0$theta, NA)
  expect_equal(fit_undir_dc0$theta, NA)
  
  expect_equal(apply(fit_dir_dc2$theta, 2, aggregate, by = list(fit_dir_dc2$FoundComms), sum)[[1]]$x, rep(1,K)) #directed degree corrrection case
  expect_equal(apply(fit_dir_dc2$theta, 2, aggregate, by = list(fit_dir_dc2$FoundComms), sum)[[2]]$x, rep(1,K))
  expect_equal(aggregate(fit_undir_dc2$theta, by = list(fit_undir_dc2$FoundComms), sum)$x, rep(1,K))
  
  expect_equal(aggregate(fit_dir_dc3$theta, by = list(fit_dir_dc3$FoundComms), sum)$x, rep(1,K))
  expect_equal(aggregate(fit_undir_dc3$theta, by = list(fit_undir_dc3$FoundComms), sum)$x, rep(1,K))
})


