#' Zachary's Karate Club network (used for testing).
#'
#' The two columns list the edges by their associated nodes
#'
#' @format A data frame with 78 rows and 2 variables:
#' \describe{
#'   \item{X1}{first node involved in an edge}
#'   \item{X2}{second node involved in an edge}
#' }
#' 
#' @source \url{
#' Krivitsky P (2019). _ergm.count: Fit, Simulate and Diagnose
#' Exponential-Family Models for Networks with Count Edges_. The Statnet
#' Project (<URL: https://statnet.org>). R package version 3.4.0, <URL:
#  https://CRAN.R-project.org/package=ergm.count>
#' }
#' 
#' @references \{
#' 
#' Zachary, WW (1977). An Information Flow Model for Conflict and Fission in Small Groups.
#' Journal of Anthropological Research, 33(4), 452-473.
#' Sociomatrix in machine-readable format was retrieved from http://vlado.fmf.uni-lj.si/pub/
#  networks/data/ucinet/ucidata.htm.
#' 
#' Zachary, WW (1977). An Information Flow Model for Conflict and Fission in Small Groups.
#' Journal of Anthropological Research, 33(4), 452-473.
#' 
#' }
#'@details 
#' library(ergm.count)
#' data(zach)
#' zach = data.frame(as.edgelist(zach))
#' save(zach, file="../sbm/data/zach.rda")
"zach"

#' An (arbitrarily) weighted version of the Zachary's Karate Club network used for testing.
#'
#' The first two columns list the edges by their associated nodes, and
#' the third lists the (NOT meaningfully assigned) weights of the edges.
#'
#' @format A data frame with 78 rows and 3 variables:
#' \describe{
#'   \item{X1}{first node involved in an edge}
#'   \item{X2}{second node involved in an edge}
#'   \item{count}{weight of the edge}
#'   ...
#' }
#' @source \url{
#' Krivitsky P (2019). _ergm.count: Fit, Simulate and Diagnose
#' Exponential-Family Models for Networks with Count Edges_. The Statnet
#' Project (<URL: https://statnet.org>). R package version 3.4.0, <URL:
#  https://CRAN.R-project.org/package=ergm.count>
#' }
#' 
#' @references \{
#' 
#' Zachary, WW (1977). An Information Flow Model for Conflict and Fission in Small Groups.
#' Journal of Anthropological Research, 33(4), 452-473.
#' Sociomatrix in machine-readable format was retrieved from http://vlado.fmf.uni-lj.si/pub/
#  networks/data/ucinet/ucidata.htm.
#' 
#' Zachary, WW (1977). An Information Flow Model for Conflict and Fission in Small Groups.
#' Journal of Anthropological Research, 33(4), 452-473.
#' 
#' }
#'@details 
#' library(ergm.count)
#' data(zach)
#' zach = data.frame(as.edgelist(zach))
#' save(zach, file="../sbm/data/zach.rda")
#' zach$count = (1:13) #recycle
#' zach_weighted = zach
#' save(zach_weighted, file="../data/zach_weighted.rda")
"zach_weighted"


