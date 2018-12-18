## sbmt.R: Time-dependent extended SBM in R/cpp, based on Karrar & Newman cpp code
##
## License? (You should have received a copy of the GNU General Public License
## along with RcppExamples.  If not, see <http://www.gnu.org/licenses/>.)

## To do: Check that seed can be set for reproducibility
##        Add functionality to set initial/seed communities
##        Allow network input and function converts to edgelist
##        Allow additional entry types and time-value attribute in edgelist.time
##        More checks for proper data entry, otherwise can cause R error

#' @KLPerNetwork this is the number of KL runs on a network
#' @maxComms maximum number of communities represented in the network
#' @seedComms user-supplied starting values for KL runs. NOT YET IMPLEMENTED.
# '@directed whether the network is directed and off-diagonal block parameters may be assymetrical
# '@degreeCorrect whether to apply degree correction. 0 = no degree correction; 1 = the degree correction of Karrer and Newman, 2 = time-independent degree correction over time-dependent data
# '@seed a random seet set for reproducibility. (Implemented?)



sbmt <- function(edgelist.time, klPerNetwork = 50, maxComms = 2, seedComms = NULL, directed = F, degreeCorrect = 0, seed = NULL) {
    
        #check it's the right type of data/list format
        if (!is.list(edgelist.time) | is.null(edgelist.time[[1]])) {
            stop("edgelist.time should be a list where each entry is an edgelist of at least two columns, three if weights are supplied")
        }
        
        #extract unique nodes and re-level them to 1:N
        nodes = sort(unique(as.character(unlist(sapply(1:length(edgelist.time), function(x) {unlist(edgelist.time[[x]][,1:2])})))))
        N = length(nodes)
        link.nodes = 1:N; names(link.nodes) = nodes

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
            
            # Format weights
            if (ncol(edgelist) == 2) {edgelist1 = cbind(edgelist1, rep(1, nrow(edgelist)))}
                #weighted = FALSE
            
            if (ncol(edgelist) == 3) {edgelist1 = cbind(edgelist1, edgelist[,3])}
                #weighted = TRUE
         
            edgelist1 = matrix(as.numeric(edgelist1)-1, ncol = 3)
            
            edgelist1
        })
        
        # Create total matrix; remove diag?
        # edgelist.total = lapply(edgelist.time,
        # function (x){
        #    tmp = igraph::graph.data.frame(x, vertices = 0:(N-1))
        #    igraph::as_adjacency_matrix(tmp, sparse = F, attr = "V3")
        #    })
        
        # edgelist.total = Reduce('+', edgelist.total)
        # edgelist.total = reshape2::melt(edgelist.total)

    # Put directed arg in logical terms
    if (directed == 0 || directed == 1) {directed = as.logical(directed)}
    if (!is.logical(directed)) {
        stop("directed should be FALSE or TRUE (0 or 1 OK)")
    }

    degreeCorrect = as.integer(degreeCorrect) #to handle T/F input
    if (is.na(degreeCorrect)) stop("degreeCorrect should be a logical or an integer in the specific range, with 0 for no degree correction")
    #seed?

    ## Make the call...

    Results <- sbmtFit(edgelist.time,
    #as.matrix(edgelist.total[,1:2], ncol = 2), as.vector(edgelist.total[,3]),
    maxComms, directed, klPerNetwork, degreeCorrect, N)

    #Reformat best matrices
    T = length(edgelist.time)
    Results$TimeMatrices = sapply(1:T, function(x) {
        matrix(Results$EdgeMatrix[ (maxComms^2 * (x-1) + 1) : (maxComms^2 * x) ], nrow = maxComms, ncol = maxComms)
    }, simplify = F,  USE.NAMES = F)
    Results$link.nodes = link.nodes
    
    Results
}



# The overall process vis sbmt/:

    # We have three sets of parameters
    # g - cluster membership
        # repeatedly initialize randomly and use KL to optimize
        # proposed change in g -> update thetas
    # theta - degree correction by node
    # omega - cluster to cluster parameter (mean for poisson)


    # Initialize cluster assignments  - check index? ####

    # Initialize thetas (closed-form)####
  
    # Calculate omegas (closed from| g and thetas) ####
  
      # Update assignments (KL) ####

      # Update thetas ####
  
      # Update omegas -> ###
