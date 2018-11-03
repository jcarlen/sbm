library(sbm)
data(zach)
data(zach_weighted)

test_that("sbm works on undirected zachary karate club with binary edges", {
  expect_equal( sum(sbm(zach)$FoundComms) %in% c(5, 29) , TRUE)
  expect_equal( sum(sbm(zach, degreeCorrect = T)$FoundComms), 17)
  expect_equal( round(sbm(zach, degreeCorrect = T, directed = F)$HighestScore), -370)
})

test_that("sbm works on directed, not degree corrected zachary karate club", {
  expect_equal( sum(sbm(zach_weighted, directed = T)$FoundComms) %in% c(14, 20) , TRUE)
})

test_that("sbm works on zachary karate club with weighted edges", {
  expect_equal( sum(sbm(zach_weighted)$FoundComms) %in% c(10, 24) , TRUE)
})

test_that("sbm doesn't crash with weird input", {
  expect_equal( attributes(sbm(zach_weighted, degreeCorrect = 2.3))$degreeCorrect , FALSE)
  expect_equal( attributes(sbm(zach_weighted, degreeCorrect = TRUE))$degreeCorrect , TRUE)
})
