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
#' @source
#' Krivitsky P (2019). _ergm.count: Fit, Simulate and Diagnose
#' Exponential-Family Models for Networks with Count Edges_. The Statnet
#' Project (<URL: https://statnet.org>). R package version 3.4.0, <URL:
#   <https://CRAN.R-project.org/package=ergm.count>
#' 
#' 
#' @references
#' Zachary, WW (1977). An Information Flow Model for Conflict and Fission in Small Groups.
#' Journal of Anthropological Research, 33(4), 452-473.
#' Sociomatrix in machine-readable format was retrieved from http://vlado.fmf.uni-lj.si/pub/
#  networks/data/ucinet/ucidata.htm.
#' 
#' Zachary, WW (1977). An Information Flow Model for Conflict and Fission in Small Groups.
#' Journal of Anthropological Research, 33(4), 452-473.
#' 
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
#' 
#' @source
#' Krivitsky P (2019). _ergm.count: Fit, Simulate and Diagnose
#' Exponential-Family Models for Networks with Count Edges_. The Statnet
#' Project (<URL: https://statnet.org>). R package version 3.4.0, <URL:
#  <https://CRAN.R-project.org/package=ergm.count>
#' 
#' 
#' @references 
#' Zachary, WW (1977). An Information Flow Model for Conflict and Fission in Small Groups.
#' Journal of Anthropological Research, 33(4), 452-473.
#' Sociomatrix in machine-readable format was retrieved from http://vlado.fmf.uni-lj.si/pub/
#  networks/data/ucinet/ucidata.htm.
#' 
#' Zachary, WW (1977). An Information Flow Model for Conflict and Fission in Small Groups.
#' Journal of Anthropological Research, 33(4), 452-473.
#' 
#'@details 
#'   library(ergm.count)
#'   data(zach)
#'   zach = data.frame(as.edgelist(zach))
#'   zach$count = (1:13) #recycle
#'   zach_weighted = zach
#'   save(zach_weighted, file="../data/zach_weighted.rda")
"zach_weighted"

#' A 
#'
#' The first two columns list the edges by their associated nodes, and
#' the third lists the (NOT meaningfully assigned) weights of the edges.
#'
#' @format A length-24 list of data frames for each our of the data, 
#' #' e.g. the first item in the list describes trips originating between midnight and 1am.
#' #' Each data frame has three columns
#' \describe{
#'   \item{start_station_id}{starting station of trip by system ID}
#'   \item{end_station_id}{ending station of trip by system IDe}
#'   \item{count}{Number of trips from start_station_id to end_station_id starting in the given hour}
#' }
#' 
#' @source
#' Metro Bike Share Trip Data
#'  LA Metro, 2017 <https://bikeshare.metro.net/about/data>
#'  Download url: <https://bikeshare.metro.net/wp-content/uploads/2017/01/Metro_trips_Q4_2016.zip> 
#'  last checked 2019-05-05, provides version last modified 2018-09-17
#'
"la_byhour_edgelist"
