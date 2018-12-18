library(sbm)
data(la_byhour_edgelist)

test_that("sbm doesn't fail on uncorrected, weighted, time-dependent network", {
    
    set.seed(1)
    #undirected and directed - no degree correction
    expect_equal(
        5,
        length( sbmt(la_byhour_edgelist,  degreeCorrect = F, directed = F, klPerNetwork = 5, maxComms = 2))
    )
    
    expect_equal(
        5,
        length( sbmt(la_byhour_edgelist,  degreeCorrect = F, directed = T, klPerNetwork = 5, maxComms = 2))
    )
})

test_that("sbm doesn't fail on corrected, undirected, weighted, time-dependent network", {
    
    #undirected and directed - degree correction level 1
    expect_equal(
        5,
        length(sbmt(la_byhour_edgelist,  degreeCorrect = 1, directed = F, klPerNetwork = 5, maxComms = 2))
    )
    
    expect_equal(
        5,
        length(sbmt(la_byhour_edgelist,  degreeCorrect = T, directed = T, klPerNetwork = 5, maxComms = 2))
    )
})
