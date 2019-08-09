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
#
#' @param edgelist.time Time-sliced network to model, formatted as a list of networks in edgelist format
#' @param klPerNetwork this is the number of KL runs on a network
#' @param maxComms maximum number of communities represented in the network
#' @param seedComms user-supplied starting values for KL runs. NOT YET IMPLEMENTED.
#' @param directed whether the network is directed and off-diagonal block parameters may be assymetrical
#' @param degreeCorrect whether to apply degree correction.:
#' \itemize{
#'   \item 0 = no degree correction
#'   \item 1 = the degree correction of Karrer and Newman at every time step
#'   \item 2 = time-independent degree correction over time-dependent data - inherets directedness of graph
#'   \item 3 = time-independent degree correction over time-dependent data - undirected regardless
#'}
#'@param tolerance stopping criteria for KL. Prevents loops due to numerical errors.
#'@param seed a random seet set for reproducibility. (Implemented?)
#
#'@return FoundComms is a vector of node labels with estimated block assignments.
#'@return EdgeMatrix is 
sbmt <- function(edgelist.time, maxComms = 2, degreeCorrect = 0, directed = F, klPerNetwork = 50, tolerance = 1e-4, seedComms = NULL, seed = NULL) {
    
        #check it's the right type of data/list format
        if (!is.list(edgelist.time) | is.null(edgelist.time[[1]])) {
            stop("edgelist.time should be a list where each entry is an edgelist of at least two columns, three if weights are supplied")
        }
        
        #extract unique nodes and re-level them to 1:N
        nodes = sort(unique(as.character(unlist(sapply(1:length(edgelist.time), function(x) {unlist(edgelist.time[[x]][,1:2])})))))
        N = length(nodes)
        link.nodes = 0:(N-1); names(link.nodes) = nodes

        ncol.min = min(sapply(edgelist.time, ncol))
        
        edgelist.time = lapply(edgelist.time, function(edgelist) {
            # Argument Checks
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

    # Put directed arg in logical terms
    if (directed == 0 || directed == 1) {directed = as.logical(directed)}
    if (!is.logical(directed)) {
        stop("directed should be FALSE or TRUE (0 or 1 OK)")
    }

    degreeCorrect = as.integer(degreeCorrect) #to handle T/F input
    if (is.na(degreeCorrect)) stop("degreeCorrect should be a logical or an integer in the specific range, with 0 for no degree correction")
    if (degreeCorrect == 3 & directed == FALSE) {degreeCorrect = 2} #Degree Correct three used when graph is directed to induce non-directed degree correction.
    #seed?

    ## Make the call...

    Results <- sbmtFit(edgelist.time, maxComms, directed, klPerNetwork, degreeCorrect, N, tolerance)

    #Reformat best matrices
    Time = length(edgelist.time)
    
    # have levels numbered in order of appearance so that lable switching doesn't affect output
    Results$FoundComms = as.factor(Results$FoundComms)
    levels(Results$FoundComms) = levels(Results$FoundComms)[rank(sapply(0:(maxComms-1), function(x) {which(Results$FoundComms==x)[1]}))]
    tmp.levels = as.numeric(levels(Results$FoundComms))+1
    Results$FoundComms = as.numeric(as.character(Results$FoundComms))
    # Return found community membership in order of ID
    names(Results$FoundComms) = names(link.nodes)
    
    Results$EdgeMatrix= sapply(1:Time, function(x) {
      matrix(Results$EdgeMatrix[ (maxComms^2 * (x-1) + 1) : (maxComms^2 * x) ], nrow = maxComms, ncol = maxComms, byrow = T)
    }, simplify = F,  USE.NAMES = F)
    Results$EdgeMatrix = lapply(Results$EdgeMatrix, function(x) {x[order(tmp.levels), order(tmp.levels)]})
    
    Results$directed = directed
    Results$klPerNetwork = klPerNetwork
    Results$degreeCorrect = degreeCorrect
    Results$tolerance = tolerance
    
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
