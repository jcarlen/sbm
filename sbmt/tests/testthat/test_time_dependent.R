library(sbm)
data(la_byhour_edgelist)

test_that("sbmt on uncorrected, weighted, time-dependent network", {
    
    set.seed(1)
    #undirected and directed - no degree correction
    expect_equal(
        7,
        length( sbmt(la_byhour_edgelist,  degreeCorrect = F, directed = F, klPerNetwork = 5, maxComms = 2))
    )
    
    expect_equal(
        7,
        length( sbmt(la_byhour_edgelist,  degreeCorrect = F, directed = T, klPerNetwork = 5, maxComms = 2))
    )
})

test_that("sbmt on time-slice corrected, weighted, time-dependent network", {
    
    #undirected and directed - degree correction level 1
    expect_equal(
        7,
        length(sbmt(la_byhour_edgelist,  degreeCorrect = 1, directed = F, klPerNetwork = 5, maxComms = 2))
    )
    
    expect_equal(
        7,
        length(sbmt(la_byhour_edgelist,  degreeCorrect = T, directed = T, klPerNetwork = 5, maxComms = 2))
    )
})

test_that("sbmt on time-independent corrected, directed, weighted, time-dependent network", {
    
    #undirected and directed - degree correction level 2
    expect_equal(
    7,
    length(sbmt(la_byhour_edgelist,  degreeCorrect = 2, directed = T, klPerNetwork = 5, maxComms = 2))
    )
    
    expect_equal(
    7,
    length(sbmt(la_byhour_edgelist,  degreeCorrect = 3, directed = T, klPerNetwork = 5, maxComms = 2))
    )

})
