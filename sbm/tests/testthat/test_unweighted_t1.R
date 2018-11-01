library(sbm)
data(zach)
data(zach_weighted)

test_that("sbm works on zachary karate club with binary edges", {
  expect_equal( sum(sbm(zach)$FoundComms) %in% c(5, 29) , TRUE)
  expect_equal( sum(sbm(zach, degreeCorrect = T)$FoundComms), 17)
  expect_equal( round(sbm(zach, degreeCorrect = T, directed = F)$HighestScore), -370)
})

test_that("sbm works on zachary karate club with weighted edges", {
  expect_equal( sum(sbm(zach_weighted)$FoundComms)==24 | sum(sbm(zach_weighted)$FoundComms)==10 , TRUE)
})

test_that("sbm doesn't crash with weird input", {
  expect_equal( attributes(sbm(zach_weighted, degreeCorrect = 2.3))$degreeCorrect , TRUE)
})
