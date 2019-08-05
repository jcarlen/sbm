library(sbmt)
data(la_byhour_edgelist)

# edgelist_to_adj

test_that("edgelist to adjacency, as.array=TRUE, undireced", {
  
    #selfedges FALSE ----
  
      #on diagonal 
  expect_equal(edgelist_to_adj(la_byhour_edgelist, selfedges = F, as.array = T, directed = F)[9,9,],
               setNames(rep(0,24), as.character(1:24)))
  
      #off diagonal
  expect_equal(edgelist_to_adj(la_byhour_edgelist, selfedges = F, as.array = T, directed = F)[1,2,],
               setNames(sapply(la_byhour_edgelist, function(x) {
                    sum(x[,3][(x[,1]=="3008" & x[,2]=="3031") | (x[,1]=="3031" & x[,2]=="3008")])
                 }), as.character(1:24)))
  
    #selfedges TRUE ----
      #on diagonal
  
  expect_equal(edgelist_to_adj(la_byhour_edgelist, selfedges = T, as.array = T, directed = F)[9,9,],
               setNames(sapply(la_byhour_edgelist, function(x) {
                 sum(x[,3][(x[,1]=="3038" & x[,2]=="3038")])
                 }), as.character(1:24)))
       
     #off diagonal        
  expect_equal(edgelist_to_adj(la_byhour_edgelist, selfedges = T, as.array = T, directed = F)[1,2,],
               setNames(sapply(la_byhour_edgelist, function(x) {
                 sum(x[,3][(x[,1]=="3008" & x[,2]=="3031") | (x[,1]=="3031" & x[,2]=="3008")])
               }), as.character(1:24)))
  
               
})

test_that("edgelist to adjacency, as.array=TRUE, directed", {

    #selfedges FALSE ----
  
    #on diagonal 
  expect_equal(edgelist_to_adj(la_byhour_edgelist, selfedges = F, as.array = T, directed = T)[9,9,],
               setNames(rep(0,24), as.character(1:24)))
  
    #off diagonal
  expect_equal(edgelist_to_adj(la_byhour_edgelist, selfedges = F, as.array = T, directed = T)[1,2,],
               setNames(sapply(la_byhour_edgelist, function(x) {
                 sum(x[,3][(x[,1]=="3008" & x[,2]=="3031")])}), as.character(1:24)))
  
    #selfedges TRUE ----
    #on diagonal
  
  expect_equal(edgelist_to_adj(la_byhour_edgelist, selfedges = T, as.array = T, directed = T)[9,9,],
               setNames(sapply(la_byhour_edgelist, function(x) {
                 sum(x[,3][(x[,1]=="3038" & x[,2]=="3038")])
               }), as.character(1:24)))
  
   #off diagonal        
  expect_equal(edgelist_to_adj(la_byhour_edgelist, selfedges = T, as.array = T, directed = T)[1,2,],
               setNames(sapply(la_byhour_edgelist, function(x) {
                 sum(x[,3][(x[,1]=="3008" & x[,2]=="3031")])}), as.character(1:24)))
  
  
})

test_that("edgelist to adjacency, as.array=FALSE, undirected", {
  
    #selfedges FALSE ----
  
  #on diagonal 
  expect_equal(sum(sapply(edgelist_to_adj(la_byhour_edgelist, selfedges = F, as.array = F, directed = F), diag)), 
                   0)
  
  #off diagonal
  expect_equal(edgelist_to_adj(la_byhour_edgelist, selfedges = F, as.array = F, directed = F)[[11]][2,5],
               sum((la_byhour_edgelist[[11]][,3])[
                 (la_byhour_edgelist[[11]][,1]=="3031" & la_byhour_edgelist[[11]][,2]=="3005") | 
                   (la_byhour_edgelist[[11]][,1]=="3005" & la_byhour_edgelist[[11]][,2]=="3031")
                 ]))

  
    #selfedges TRUE ----
  
  #on diagonal
  
  expect_equal(lapply(edgelist_to_adj(la_byhour_edgelist, selfedges = T, as.array = F, directed = F), diag)[[12]][2], 
               setNames((la_byhour_edgelist[[12]][,3])[
                 (la_byhour_edgelist[[12]][,1]=="3031" & la_byhour_edgelist[[12]][,2]=="3031")
                 ], "3031"))
  
  #off diagonal        
  expect_equal(edgelist_to_adj(la_byhour_edgelist, selfedges = T, as.array = F, directed = F)[[11]][2,5],
               sum((la_byhour_edgelist[[11]][,3])[
                 (la_byhour_edgelist[[11]][,1]=="3005" & la_byhour_edgelist[[11]][,2]=="3031") | 
                   (la_byhour_edgelist[[11]][,1]=="3031" & la_byhour_edgelist[[11]][,2]=="3005")
                 ]))
  
})

test_that("edgelist to adjacency, as.array=FALSE, directed", {
   
    #selfedges FALSE ----
  
  #on diagonal 
  expect_equal(sum(sapply(edgelist_to_adj(la_byhour_edgelist, selfedges = F, as.array = F, directed = T), diag)), 
               0)
  
  #off diagonal
  expect_equal(edgelist_to_adj(la_byhour_edgelist, selfedges = F, as.array = F, directed = T)[[11]][2,5],
               (la_byhour_edgelist[[11]][,3])[
                 (la_byhour_edgelist[[11]][,1]=="3031" & la_byhour_edgelist[[11]][,2]=="3005")])
  
  
    #selfedges TRUE ----
  #on diagonal
  
  expect_equal(lapply(edgelist_to_adj(la_byhour_edgelist, selfedges = T, as.array = F, directed = T), diag)[[12]][2], 
               setNames((la_byhour_edgelist[[12]][,3])[
                 (la_byhour_edgelist[[12]][,1]=="3031" & la_byhour_edgelist[[12]][,2]=="3031")
                 ], "3031"))
  
  #off diagonal        
  expect_equal(edgelist_to_adj(la_byhour_edgelist, selfedges = T, as.array = F, directed = T)[[11]][2,5],
               (la_byhour_edgelist[[11]][,3])[
                 (la_byhour_edgelist[[11]][,1]=="3031" & la_byhour_edgelist[[11]][,2]=="3005")])
  
  
  
})

# tdd_sbm_llik
# add more tests here

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


# tdmm_sbm_llik