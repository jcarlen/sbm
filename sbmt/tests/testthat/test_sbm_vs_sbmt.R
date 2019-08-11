library(sbmt)
data(zach)

test_that("test time-dependent and not time-dependent code gives same results when T = 1", {
  expect_equal(sbm(zach, seed = 1)$HighestScore, sbmt(list(zach), seed = 1)$HighestScore)
  expect_equal(sbm(zach, directed = T, seed = 1)$HighestScore, sbmt(list(zach), directed = T, seed = 1)$HighestScore)
  expect_equal(sbm(zach, degreeCorrect = 1, seed = 1)$HighestScore,  sbmt(list(zach), degreeCorrect = 1, seed = 1)$HighestScore)
  expect_equal(sbm(zach, directed = T, degreeCorrect = 1, seed = 1)$HighestScore, 
               sbmt(list(zach), directed = T, degreeCorrect = 1, seed = 1)$HighestScore)
  expect_equal(sbm(zach, directed = T, degreeCorrect = 1, seed = 1)$HighestScore, 
               sbmt(list(zach), directed = T, degreeCorrect = 2, seed = 1)$HighestScore)
})

