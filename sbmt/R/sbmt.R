## sbmt.R: Time-dependent extended SBM in R/cpp, based on Karrar & Newman cpp code
##
## License? (You should have received a copy of the GNU General Public License
## along with RcppExamples.  If not, see <http://www.gnu.org/licenses/>.)

## To do: Check that seed can be set for reproducibility
##        Add functionality to set initial/seed communities
##        Allow network input and function converts to edgelist
##        Allow additional entry types and time-value attribute in edgelist.time
##        More checks for proper data entry, otherwise can cause R error
##        Update la edgelist in sbm package to agree with paper (la_byhour)
##        Implement undirected degree correction for static sbm on directed network?

#' Estimate parameters of the TDD-SBM (time-dependent discrete stochastic block model) by Kernighan-Lin algorithm.
#' Note that although the motivation for this function and package is analyzing time-sliced networks, it can handle multilayer networks more generally as long as the primary input is a list of networks in edgelist format. (Slice parameters are independent given block membership.) 
#
#' @param edgelist.time A (time) series of networks represented as a list of edgelists. 
#' Assumes all edgelist slices have the same names and number of columns. The first two columns designate edges "from" and "to", and the third, if present, is the count for each edge.
#' @param klPerNetwork this is the number of KL runs on a network
#' @param maxComms maximum number of communities (blocks) represented in the network
#' @param degreeCorrect whether to apply degree correction:
#' \itemize{
#'   \item 0 = no degree correction
#'   \item 1 = the degree correction of Karrer and Newman at every time step - inherits directedness of graph
#'   \item 2 = time-independent degree correction over time-dependent data - inherits directedness of graph
#'   \item 3 = time-independent degree correction over time-dependent data - undirected regardless
#'}
#' @param directed whether the network is directed and off-diagonal block parameters may be assymetrical
#'@param tolerance stopping criteria for KL. Prevents loops due to numerical errors.
#'@param seed a random seet set for reproducibility.
#'@param seedComms user-supplied starting values for KL runs. They will be converted to integer levels numbered starting at 0.
#'To avoid confusion, they must have names equal to the unique nodes in the network. E.g., If nodes are numbered numerically from 1 then seedComm names are "1", "2",... 
#
#'@return FoundComms A vector of node labels with estimated block assignments.
#'@return EdgeMatrix time-specific block-to-block edges corresponding to FoundComms
#'@return theta nodal degree correction parameters. NA if degreeCorrect == 0; not yet implemented for degreeCorrect == 1; 
#'N x 2 data frame (one column for each direction) if network is directed and degreeCorrect == 2; and a length-N numeric vector otherwise.
#'@return llik unnormalized log-likelihood of result as calculated by tdd_sbm_llik
#'@return directed Whether input network was considered directed
#'@return degreeCorrect type of degreeCorrection used in fitting
#'@return klPerNetwork number of runs of KL algorithm used in fitting
#'@return tolerance tolerance value used for fitting
#'@return init seed = set for fitting, if supplied;  seedComms = initial communities (blocks) as fed to sbmtFit (conversion to zero-min integer vector may have occurred)
#' 
#' 
sbmt <- function(edgelist.time, maxComms = 2, degreeCorrect = 0, directed = F,
                 klPerNetwork = 50, tolerance = 1e-4, seed = NULL, seedComms = NULL) {
    
    # check it's the right type of data/list format ----
    if (!is.list(edgelist.time) | is.null(edgelist.time[[1]])) {
        stop("edgelist.time should be a list where each entry is an edgelist of at least two columns, three if weights are supplied")
    }
        
    # extract unique nodes and re-level them to 1:N ----
    nodes = sort(unique(as.character(unlist(sapply(1:length(edgelist.time), function(x) {unlist(edgelist.time[[x]][,1:2])})))))
    N = length(nodes)
    link.nodes = 0:(N-1); names(link.nodes) = nodes
    ncol.min = min(sapply(edgelist.time, ncol))
    # reformulate input ----
    A.time = edgelist_to_adj(edgelist.time, directed = directed) # store array from original
    edgelist.time = lapply(edgelist.time, function(edgelist) {
            # Argument Checks ----
            if (is.null(dim(edgelist)) || ncol(edgelist) < 2) {
                stop("edgelist must have at least two columns, three if weights are supplied")
            }
            
            if (ncol(edgelist) > 3 ) {
                edgelist = edgelist[,1:ncol.min]
                warning("currently only edglist format is supported with first two columns indicating from & to nodes and third column is weights, if present. Truncating to minimum column number.")
            }
            
            if (is.list(edgelist)) {edgelist = as.matrix(edgelist)} #convert from data frame to matrix if necessary
            
            # Remove potential attributes
            #edgelist = matrix(as.numeric(as.vector(edgelist)), ncol = ncol(edgelist))
            
            # Relevel IDs
            edgelist1 = cbind(link.nodes[as.character(edgelist[,1])], link.nodes[as.character(edgelist[,2])])
            
            
            if (ncol(edgelist) == 2) {
              edgelist1 = cbind(edgelist1, rep(1, nrow(edgelist)))
            }
            
            if (ncol(edgelist) == 3) {
              edgelist1 = cbind(edgelist1, edgelist[,3])
            }
            
            edgelist1 = matrix(as.numeric(edgelist1), ncol = 3)
            
            edgelist1
        })

    # put directed arg in logical terms
    if (directed == 0 || directed == 1) {directed = as.logical(directed)}
    if (!is.logical(directed)) {
        stop("directed should be FALSE or TRUE (0 or 1 OK)")
    }
    
    degreeCorrect = as.integer(degreeCorrect) #to handle T/F input
    if (is.na(degreeCorrect)) stop("degreeCorrect should be a logical or an integer in the specific range, with 0 for no degree correction")
    if (degreeCorrect == 3 & directed == FALSE) {degreeCorrect = 2} #Degree correct three only applied to directed graphs, otherwise graph and degree correction are undirected --> degree correct == 2 case.

    # check seedComms if given, reformat if necessary ----
    if (!is.null(seedComms)) {
      # check number of seed blocks
      if (length(unique(seedComms)) > maxComms) {
        stop("unique values of seedComms shouldn't exceed maxComms") 
      }
      # to avoid confusions, seedComms, if given, should have names matching the node names (as characters). they will be ordered accordingly.
      if (!identical(sort(names(seedComms)), sort(nodes))) {
        stop("to avoid confusion, seedComms should have names equal to node set (as characters)")
      } else { #sort accordingly
        seedComms = seedComms[order(names(seedComms))]
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
    } else {
      seedComms = 0
    }
    
    # set seed ----
    if (!is.null(seed)) set.seed(seed)
    
    # MAKE THE CALL ----
    Results <- sbmtFit(edgelist.time, maxComms, directed, klPerNetwork, degreeCorrect, N, tolerance, seedComms)

    # reformulate output ----
    
    # reformat best matrices
    Time = length(edgelist.time)
    
    #   Results$FoundComms ----
    # re-level blocks in order of appearance so that label switching doesn't affect output
    Results$FoundComms = as.factor(Results$FoundComms)
    levels(Results$FoundComms) = levels(Results$FoundComms)[rank(sapply(0:(maxComms-1), function(x) {which(Results$FoundComms==x)[1]}))]
    tmp.levels = as.numeric(levels(Results$FoundComms))+1
    Results$FoundComms = as.numeric(as.character(Results$FoundComms))
    
    # return found community membership with ID as name
    names(Results$FoundComms) = names(link.nodes)
    
    # reorder FoundComms by ids if they're numeric (otherwise stays character/alphabetical sorted)
    if (sum(is.na((as.numeric(names(Results$FoundComms)))))==0) {Results$FoundComms = Results$FoundComms[order(as.numeric(names(Results$FoundComms)))]}
    
    #   Results$EdgeMatrix ----
    # reformed edge matrices
    Results$EdgeMatrix= sapply(1:Time, function(x) {
      matrix(Results$EdgeMatrix[ (maxComms^2 * (x-1) + 1) : (maxComms^2 * x) ], nrow = maxComms, ncol = maxComms, byrow = T)
    }, simplify = F,  USE.NAMES = F)
    Results$EdgeMatrix = lapply(Results$EdgeMatrix, function(x) {x[order(tmp.levels), order(tmp.levels)]}) #align with re-leveled blocks
    #   Results$theta ----
    # extract theta (degree correction parameters)
    if (degreeCorrect == 0 ) {
      Results$theta = NA
    }
    if (degreeCorrect == 1 ) {
      Results$theta = "theta not yet implemented"
    }
    # directed degree correction parameters
    if (degreeCorrect == 2 & directed) {
      node_in_degree = apply(A.time, 2, sum)
      node_out_degree = apply(A.time, 1, sum)
      role_in_sum = aggregate(node_in_degree, by = list(Results$FoundComms), sum)$x[Results$FoundComms + 1]
      role_out_sum = aggregate(node_out_degree, by = list(Results$FoundComms), sum)$x[Results$FoundComms + 1]
      Results$theta = cbind(theta_in = node_in_degree/role_in_sum, theta_out = node_out_degree/role_out_sum)
    }
    # undirected degree correction parameters
    if (degreeCorrect == 3  | degreeCorrect == 2 & !directed) {
        if (!directed) {
          node_degree = (apply(A.time, 1, sum) + apply(A.time, 2, sum))/2
        }
        if (directed) {
          node_degree = apply(A.time, 1, sum) + apply(A.time, 2, sum)
        }
        node_degree = node_degree[names(Results$FoundComms)] #reorder to agree with Results$FoundComms
        role_sum = aggregate(node_degree, by = list(Results$FoundComms), sum)$x[Results$FoundComms + 1]
        Results$theta = node_degree/role_sum
        # check: aggregate(theta, by = list(Results$FoundComms), sum)
    }
    
    #   Results$ all other stuff ----
    Results$llik = ifelse(degreeCorrect == 3, tdd_sbm_llik(A.time, roles = Results$FoundComms, omega = Results$EdgeMatrix, directed = directed), "llik currently only implemented for degreeCorrect = 3" )
    Results$degreeCorrect = degreeCorrect
    Results$directed = directed
    Results$klPerNetwork = klPerNetwork
    Results$tolerance = tolerance
    Results$init = list(seed = seed, seedComms = seedComms)
    
    # set method for plot
    class(Results) <- "sbmt"
    
    Results
}


# The overall process via sbmt/:

    # We have three sets of parameters
    # g - cluster membership
        # repeatedly initialize randomly and use KL to optimize
        # proposed change in g -> update thetas
    # theta - degree correction by node
    # omega - cluster to cluster parameter (mean for poisson)


    # Initialize cluster assignments  - check index? ####

    # Initialize thetas (closed-form)####
  
    # Calculate omegas (closed form| g and thetas) ####
  
      # Update assignments (KL) ####

      # Update thetas ####
  
      # Update omegas -> ###
