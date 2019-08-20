library(sbmt)
data(la_byhour_edgelist)
data(zach)

# edgelist_to_adj ----

test_that("edgelist to adjacency, as.array=TRUE, undireced", {
  
    #selfedges FALSE ----
  
      #on diagonal 
  expect_equal(edgelist_to_adj(la_byhour_edgelist, selfedges = FALSE, as.array = T, directed = FALSE)[9,9,],
               setNames(rep(0,24), as.character(1:24)))
  
      #off diagonal
  expect_equal(edgelist_to_adj(la_byhour_edgelist, selfedges = FALSE, as.array = T, directed = FALSE)[1,2,],
               setNames(sapply(la_byhour_edgelist, function(x) {
                    sum(x[,3][(x[,1]=="3008" & x[,2]=="3031") | (x[,1]=="3031" & x[,2]=="3008")])
                 }), as.character(1:24)))
  
    #selfedges TRUE ----
      #on diagonal
  
  expect_equal(edgelist_to_adj(la_byhour_edgelist, selfedges = T, as.array = T, directed = FALSE)[9,9,],
               setNames(sapply(la_byhour_edgelist, function(x) {
                 sum(x[,3][(x[,1]=="3038" & x[,2]=="3038")])
                 }), as.character(1:24)))
       
     #off diagonal        
  expect_equal(edgelist_to_adj(la_byhour_edgelist, selfedges = T, as.array = T, directed = FALSE)[1,2,],
               setNames(sapply(la_byhour_edgelist, function(x) {
                 sum(x[,3][(x[,1]=="3008" & x[,2]=="3031") | (x[,1]=="3031" & x[,2]=="3008")])
               }), as.character(1:24)))
  
               
})

test_that("edgelist to adjacency, as.array=TRUE, directed", {

    #selfedges FALSE ----
  
    #on diagonal 
  expect_equal(edgelist_to_adj(la_byhour_edgelist, selfedges = FALSE, as.array = T, directed = TRUE)[9,9,],
               setNames(rep(0,24), as.character(1:24)))
  
    #off diagonal
  expect_equal(edgelist_to_adj(la_byhour_edgelist, selfedges = FALSE, as.array = TRUE, directed = TRUE)[1,2,],
               setNames(sapply(la_byhour_edgelist, function(x) {
                 sum(x[,3][(x[,1]=="3008" & x[,2]=="3031")])}), as.character(1:24)))
  
    #selfedges TRUE ----
    #on diagonal
  
  expect_equal(edgelist_to_adj(la_byhour_edgelist, selfedges = TRUE, as.array = TRUE, directed = TRUE)[9,9,],
               setNames(sapply(la_byhour_edgelist, function(x) {
                 sum(x[,3][(x[,1]=="3038" & x[,2]=="3038")])
               }), as.character(1:24)))
  
   #off diagonal        
  expect_equal(edgelist_to_adj(la_byhour_edgelist, selfedges = TRUE, as.array = TRUE, directed = TRUE)[1,2,],
               setNames(sapply(la_byhour_edgelist, function(x) {
                 sum(x[,3][(x[,1]=="3008" & x[,2]=="3031")])}), as.character(1:24)))
  
  
})

test_that("edgelist to adjacency, as.array=FALSE, undirected", {
  
    #selfedges FALSE ----
  
  #on diagonal 
  expect_equal(sum(sapply(edgelist_to_adj(la_byhour_edgelist, selfedges = FALSE, as.array = FALSE, directed = FALSE), diag)), 
                   0)
  
  #off diagonal
  expect_equal(edgelist_to_adj(la_byhour_edgelist, selfedges = FALSE, as.array = FALSE, directed = FALSE)[[11]][2,5],
               sum((la_byhour_edgelist[[11]][,3])[
                 (la_byhour_edgelist[[11]][,1]=="3031" & la_byhour_edgelist[[11]][,2]=="3005") | 
                   (la_byhour_edgelist[[11]][,1]=="3005" & la_byhour_edgelist[[11]][,2]=="3031")
                 ]))

  
    #selfedges TRUE ----
  
  #on diagonal
  
  expect_equal(lapply(edgelist_to_adj(la_byhour_edgelist, selfedges = TRUE, as.array = FALSE, directed = FALSE), diag)[[12]][2], 
               setNames((la_byhour_edgelist[[12]][,3])[
                 (la_byhour_edgelist[[12]][,1]=="3031" & la_byhour_edgelist[[12]][,2]=="3031")
                 ], "3031"))
  
  #off diagonal        
  expect_equal(edgelist_to_adj(la_byhour_edgelist, selfedges = TRUE, as.array = FALSE, directed = FALSE)[[11]][2,5],
               sum((la_byhour_edgelist[[11]][,3])[
                 (la_byhour_edgelist[[11]][,1]=="3005" & la_byhour_edgelist[[11]][,2]=="3031") | 
                   (la_byhour_edgelist[[11]][,1]=="3031" & la_byhour_edgelist[[11]][,2]=="3005")
                 ]))
  
})

test_that("edgelist to adjacency, as.array=FALSE, directed", {
   
    #selfedges FALSE ----
  
  #on diagonal 
  expect_equal(sum(sapply(edgelist_to_adj(la_byhour_edgelist, selfedges = FALSE, as.array = FALSE, directed = TRUE), diag)), 
               0)
  
  #off diagonal
  expect_equal(edgelist_to_adj(la_byhour_edgelist, selfedges = FALSE, as.array = FALSE, directed = TRUE)[[11]][2,5],
               (la_byhour_edgelist[[11]][,3])[
                 (la_byhour_edgelist[[11]][,1]=="3031" & la_byhour_edgelist[[11]][,2]=="3005")])
  
  
    #selfedges TRUE ----
  #on diagonal
  
  expect_equal(lapply(edgelist_to_adj(la_byhour_edgelist, selfedges = TRUE, as.array = FALSE, directed = TRUE), diag)[[12]][2], 
               setNames((la_byhour_edgelist[[12]][,3])[
                 (la_byhour_edgelist[[12]][,1]=="3031" & la_byhour_edgelist[[12]][,2]=="3031")
                 ], "3031"))
  
  #off diagonal        
  expect_equal(edgelist_to_adj(la_byhour_edgelist, selfedges = TRUE, as.array = FALSE, directed = TRUE)[[11]][2,5],
               (la_byhour_edgelist[[11]][,3])[
                 (la_byhour_edgelist[[11]][,1]=="3031" & la_byhour_edgelist[[11]][,2]=="3005")])
  
  
  
})

# tdd_sbm_llik - add more tests ----

test_that("tdd_sbm_llik - likelihood for discrete model", {

  #directed
  A = edgelist_to_adj(la_byhour_edgelist, as.array = TRUE, directed = TRUE, selfedges = TRUE)
  model1 = sbmt(la_byhour_edgelist, degreeCorrect = 3, directed = TRUE, klPerNetwork = 2, maxComms = 3)
  
  expect_gt(tdd_sbm_llik(A, model1$FoundComms, model1$EdgeMatrix, directed = TRUE, selfEdges = FALSE), 
            tdd_sbm_llik(A, model1$FoundComms, model1$EdgeMatrix, directed = TRUE, selfEdges = TRUE))
  
  
  #undirected
  A = edgelist_to_adj(la_byhour_edgelist, as.array = TRUE, directed = FALSE, selfedges = TRUE)
  model1 = sbmt(la_byhour_edgelist, degreeCorrect = 3, directed = FALSE, klPerNetwork = 2, maxComms = 3)
  
  expect_gt(tdd_sbm_llik(A, model1$FoundComms, model1$EdgeMatrix, directed = FALSE, selfEdges = FALSE), 
            tdd_sbm_llik(A, model1$FoundComms, model1$EdgeMatrix, directed = FALSE, selfEdges = TRUE))
  
})

# tdmm_sbm_llik - add more tests ----

# Test likelihood accounting is right: "HighestScore" should change by ~ same amount as tdd_sbm_llik 
# from initialization to final (note reported InitialScore gives varying levels of accuracy)

test_that("Likelihood accounting, T = 1", {
  A = edgelist_to_adj(list(zach), directed = T)
  model1 = sbmt(list(zach), directed= T, degreeCorrect = 3, klPerNetwork = 1, seed = 1)
  init =  rep_len(c(0,1), length(model1$FoundComms)); names(init) = names(model1$FoundComms)
  
  tmp = data.frame(zach, g1 = init[as.character(zach[,1])], g2 = init[as.character(zach[,2])])
  init.EdgeMatrix = matrix(aggregate(rep(1, nrow(tmp)), by = list(tmp$g1,tmp$g2), sum)$x, nrow = 2, ncol = 2)
  
  model.init = sbmt(list(zach), directed= T, degreeCorrect = 3, klPerNetwork = 1, seedComms = init, seed = 1)
  #These should be approximately equal:
  x1 = model.init$HighestScore - (-447.463) #InitialScore
  x2 = tdd_sbm_llik(A, model.init$FoundComms, model.init$EdgeMatrix, directed = T) - tdd_sbm_llik(A, init, init.EdgeMatrix, directed = T)
  expect_true(abs(x1 - x2) < .001)
})

test_that("Likelihood accounting, T = 24", {
  A = edgelist_to_adj(la_byhour_edgelist, directed = T)
  model2 = sbmt(la_byhour_edgelist, directed= T, degreeCorrect = 3, klPerNetwork = 1, seed = 1)
  init2 =  rep_len(c(0,1), length(model2$FoundComms)); names(init2) = names(model2$FoundComms)
  init.EdgeMatrix = lapply(la_byhour_edgelist, function(data.t) {
    tmp = data.frame(data.t, g1 = init2[as.character(data.t[,1])], g2 = init2[as.character(data.t[,2])])
    tmp2 = aggregate(data.t$x, by = list(tmp$g1,tmp$g2), sum)$x
    init.EdgeMatrix = matrix(tmp2, nrow = 2, ncol = 2)
  })
  model.init = sbmt(la_byhour_edgelist, directed= T, degreeCorrect = 3, klPerNetwork = 1, seedComms = init2, seed = 1)
  #These should be approximately equal:
  x1 = model.init$HighestScore - (-429267) #InitialScore - note it's rounded, so less precision below
  x2 = tdd_sbm_llik(A, model.init$FoundComms, model.init$EdgeMatrix, directed = T) -  tdd_sbm_llik(A, init2, init.EdgeMatrix, directed = T)
  expect_true(abs(x1 - x2) < 1)
})



