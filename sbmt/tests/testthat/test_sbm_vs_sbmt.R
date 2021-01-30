library(sbmt)
data("zach")
data("la_byhour_edgelist")

test_that("test time-dependent and not time-dependent code gives same results when T = 1, directed", {
  expect_equal(sbm(zach, seed = 1)$HighestScore, sbmt(list(zach), seed = 1)$HighestScore)
  expect_equal(sbm(zach, directed = T, seed = 1)$HighestScore, sbmt(list(zach), directed = T, seed = 1)$HighestScore)
  expect_equal(sbm(zach, degreeCorrect = 1, seed = 1)$HighestScore,  sbmt(list(zach), degreeCorrect = 1, seed = 1)$HighestScore)
  expect_equal(sbm(zach, directed = T, degreeCorrect = 1, seed = 1)$HighestScore, 
               sbmt(list(zach), directed = T, degreeCorrect = 1, seed = 1)$HighestScore)
  expect_equal(sbm(zach, directed = T, degreeCorrect = 1, seed = 1)$HighestScore, 
               sbmt(list(zach), directed = T, degreeCorrect = 2, seed = 1)$HighestScore)
})

test_that("test time-dependent and not time-dependent code gives same results when T = 1, undirected", {
  
  # not degree corrected
  tmp.sbmt = sbmt(edgelist.time = list(la_byhour_edgelist[[18]]), maxComms = 3, degreeCorrect = 0, directed = F, 
       klPerNetwork = 10, tolerance = 1e-4, seed = 1, seedComms = NULL)
  
  tmp.sbm = sbm(edgelist = la_byhour_edgelist[[18]], maxComms = 3, degreeCorrect = 0, directed = F, 
      klPerNetwork = 10, tolerance = 1e-4, seed = 1, seedComms = NULL)
  
  expect_true(identical(tmp.sbm$EdgeMatrix, tmp.sbmt$EdgeMatrix[[1]]))
  expect_true(identical(tmp.sbm$FoundComms, tmp.sbmt$FoundComms))

  # degree corrected
  tmp.sbmt = sbmt(edgelist.time = list(la_byhour_edgelist[[18]]),
                  maxComms = 3, degreeCorrect = 1, directed = F,
                  klPerNetwork = 10, tolerance = 1e-4, seed = 1, seedComms = NULL)
  
  tmp.sbm = sbm(edgelist = la_byhour_edgelist[[18]],
                maxComms = 3, degreeCorrect = 1, directed = F,
                klPerNetwork = 10, tolerance = 1e-4, seed = 1, seedComms = NULL)
  
  expect_true(identical(tmp.sbm$EdgeMatrix, tmp.sbmt$EdgeMatrix[[1]]))
  expect_true(identical(tmp.sbm$FoundComms, tmp.sbmt$FoundComms))
  
  # but a bug when we use degreeCorrect = 2 in sbmt ?? 
  # This should be the same as degreeCorrect = 1 when there is only one time slice in the network
  # see https://github.com/jcarlen/sbm/issues/3
  # 
  # tmp.sbmt = sbmt(edgelist.time = list(la_byhour_edgelist[[18]]),
  #                 maxComms = 3, degreeCorrect = 2, directed = F,
  #                 klPerNetwork = 10, tolerance = 1e-4, seed = 1, seedComms = NULL)
  # 
  # expect_true(identical(tmp.sbm$EdgeMatrix, tmp.sbmt$EdgeMatrix[[1]]))
  # expect_true(identical(tmp.sbm$FoundComms, tmp.sbmt$FoundComms))
  
})
