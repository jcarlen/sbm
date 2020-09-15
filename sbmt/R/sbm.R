## sbm.R: Running extended Karrar SBM code via R
#
#' Estimate parameters of a static SBM by Kernighan-Lin algorithm. Allows for directed networks and degree correction.
#
#' @param edgelist Network in edgelist format. The first two columns designate edges "from" and "to", and the third, if present, is the count for each edge.
#' @param maxComms maximum number of communities (blocks) represented in the network
#' @param degreeCorrect whether to apply degree correction.:
#' #' \itemize{
#'   \item 0 = no degree correction
#'   \item 1 = the degree correction of Karrer and Newman - inherits directedness of graph
#' }
#' @param directed whether the network is directed and off-diagonal block parameters may be assymetrical
#' @param klPerNetwork this is the number of KL runs on a network
#' @param tolerance stopping criteria for KL. Prevents loops due to numerical errors.
#' @param seed a random seet set for reproducibility.
#' @param seedComms user-supplied starting values for KL runs. They will be converted to integer levels numbered starting at 0
#
#'@return FoundComms A vector of node labels with estimated block assignments.
#'@return EdgeMatrix block to block edges corresponding to FoundComms
#'@return HighestScore Highest score found by algorithm runs
#'@return llik unnormalized log-likelihood of result as calculated by tdd_sbm_llik
#'@return directed Whether input network was considered directed
#'@return degreeCorrect type of degreeCorrection used in fitting
#'@return klPerNetwork number of runs of KL algorithm used in fitting
#'@return tolerance tolerance value used for fitting
#'@return init seed = set for fitting, if supplied;  seedComms = initial communities (blocks) as fed to sbmtFit (conversion to zero-min integer vector may have occurred)
#'
#' 
sbm <- function(edgelist, maxComms = 2, degreeCorrect = F, directed = FALSE,
                klPerNetwork = 50, tolerance = 1e-4, seed = NULL, seedComms = NULL)

{
    # Argument Checks
    if (is.null(dim(edgelist)) || ncol(edgelist) < 2) {
        stop("edgelist must have at least two columns, three if weights are supplied")
    }
    A = edgelist_to_adj(list(edgelist), directed = directed) # store array from original
    if (is.list(edgelist)) {edgelist = as.matrix(edgelist)} #convert from data frame to matrix if necessary
    
    # Remove potential attributes
    edgelist = matrix(as.numeric(as.vector(edgelist)), ncol = ncol(edgelist))
    
    # Relevel IDs
    edgelist1 = as.factor(c(edgelist[,1], edgelist[,2]))
    N = length(unique(edgelist1)); link.nodes = 1:N; names(link.nodes) = levels(edgelist1)
    levels(edgelist1) = 1:N
    edgelist1 = matrix(as.numeric(edgelist1), ncol = 2)
    
    # Format weights
    if (ncol(edgelist) == 2) {
        weights = rep(1, nrow(edgelist))
        weighted = FALSE
    }
    
    if (ncol(edgelist) == 3) {
        weights = edgelist[,3]
        weighted = TRUE
    }
    
    
    # Put directed arg in logical terms
    if (directed == 0 || directed == 1) {directed = as.logical(directed)}
    if (!is.logical(directed)) {
        stop("Directed should be FALSE or TRUE (0 or 1 OK)")
    }
    
    if(degreeCorrect %in% c(0,1)) {degreeCorrect = as.numeric(degreeCorrect)} else {
        warning("Degree correct should be 0 or 1 for non-time-dependent network. Defaulted to 0 (no correction)")
        degreeCorrect = 0
        # relic of Karrer and Newman code to have degreeCorrect be 0|1
    }
    
    # set seed
    if (!is.null(seed)) set.seed(seed)
    
    # check seedComms if given, reformat if necessary
    if (!is.null(seedComms)) {
      if (length(seedComms) != N | length(unique(seedComms)) > maxComms) {
        stop("seedComms should have length == #nodes and unique values shouldn't exceed maxComms")
      }
      
      # if not integer format, convert to cpp appropriate (0,1,2,...) style
      if (!is.integer(seedComms)) {
        seedComms = as.numeric(as.factor(seedComms))-1
      } else {
        # convert to cpp range (0,1,...)
        if (min(seedComms) != 0) {
          seedComms = seedComms - min(seedComms)
        }
      } 
      # check cpp range (0,1,2,...)
      if (max(seedComms) > (maxComms - 1))  stop("seedComm range exceeds maxComms")
    } else 
    {
      seedComms = 0
    }
    
    ## Make the call...
    Results <- sbmFit(edgelist1-1, maxComms, degreeCorrect, directed, klPerNetwork, weights, tolerance, seedComms)
    
    cat("\nResults!\n")
    
    # re-level blocks in order of appearance so that label switching doesn't affect output
    Results$FoundComms = as.factor(Results$FoundComms)
    levels(Results$FoundComms) = levels(Results$FoundComms)[rank(sapply(0:(maxComms-1), function(x) {which(Results$FoundComms==x)[1]}))]
    tmp.levels = as.numeric(levels(Results$FoundComms))+1
    Results$FoundComms = as.numeric(as.character(Results$FoundComms))
    names(Results$FoundComms) = names(link.nodes) # Return found community membership in order of ID
    
    Results$EdgeMatrix = matrix(Results$EdgeMatrix, maxComms, maxComms, byrow = T)
    Results$EdgeMatrix = Results$EdgeMatrix[order(tmp.levels), order(tmp.levels)] #align with re-leveled blocks
    
    Results$llik = "Not yet implemented for sbm" #tdd_sbm_llik(A, Results$FoundComms, Results$EdgeMatrix, directed = directed)
    Results$degreeCorrect = degreeCorrect
    Results$directed = directed
    Results$klPerNetwork = klPerNetwork
    Results$tolerance = tolerance
    Results$init = list(seed = seed, seedComms = seedComms)
    
    Results
}


