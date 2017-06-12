## KL.R: call CPP code for Kernighan Lin Algorithm, developed by Karrer and Newman
##
## License? (You should have received a copy of the GNU General Public License
## along with RcppExamples.  If not, see <http://www.gnu.org/licenses/>.)

## To do: Check that seed can be set for reproducibility
##        Update ComputeProposal function to allow for time-dependent data
##        Allow network input and function converts to edgelist
##        Allow additional entry types and time-value attribute in edgelist.time

#' @maxComms maximum number of communities represented in the network
# '@degreeCorrect whether to apply the degree correction of Karrer and Newman
# '@directed whether the network is directed and off-diagonal block parameters may be assymetrical
#' @KLPerNetwork this is the number of KL runs on a network



#R Function
# need to pass AdjList objects to cpp

KLt <- function(edgelist.time, klPerNetwork = 50, maxComms = 2, seedComms = NULL, directed = F, degreeCorrect = F, seed = NULL) {
    
    if (!is.list(edgelist.time)) {
        stop("edgelist.time should be a list where each entry is an edgelist of at least two columns, three if weights are supplied")
    } else {
        
        #extract unique nodes and re-level them to 1:N
        nodes = as.factor(unique(unlist(sapply(1:length(edgelist.time), function(x) {unlist(edgelist.time[[x]][,1:2])}))))
        N = length(nodes)
        i = 1
        
        for (edgelist in edgelist.time) {
            # Argument Checks
            if (is.null(dim(edgelist)) || ncol(edgelist) < 2) {
                stop("edgelist must have at least two columns, three if weights are supplied")
            }
            if (is.list(edgelist)) {edgelist = as.matrix(edgelist)} #convert from data frame to matrix if necessary
            
            # Remove potential attributes
            edgelist = matrix(as.numeric(as.vector(edgelist)), ncol = ncol(edgelist))
            
            # Relevel IDs
            edgelist1 = as.factor(c(edgelist[,1], edgelist[,2]))
            levels(edgelist1) = 1:N
            edgelist1 = matrix(as.numeric(edgelist1), ncol = 2) - 1 #for 0 indexing in cpp
            
            # Format weights
            if (ncol(edgelist) == 2) {
                edgelist1 = cbind(edgelist1, rep(1, nrow(edgelist)))
                #weighted = FALSE
            }
            
            if (ncol(edgelist) == 3) {
                edgelist1 = cbind(edgelist1, edgelist[,3])
                #weighted = TRUE
            }
            
            edgelist.time[[i]] = edgelist1
            i = i + 1
        }
    
    # Create total matrix; remove diag?
    edgelist.total = lapply(edgelist.time,
    function (x){
        tmp = igraph::graph.data.frame(x, vertices = 0:(N-1))
        igraph::as_adjacency_matrix(tmp, sparse = F, attr = "V3")
        })
        
    edgelist.total = Reduce('+', edgelist.total)
    edgelist.total = reshape2::melt(edgelist.total)
    
    }
    
    # Put directed arg in logical terms
    if (directed == 0 || directed == 1) {directed = as.logical(directed)}
    if (!is.logical(directed)) {
        stop("Directed should be FALSE or TRUE (0 or 1 OK)")
    }
    
    degreeCorrect = as.numeric(degreeCorrect) # relic of Karrer and Newman code to have degreeCorrect be 0|1
    
    #seed?
    
    ## Make the call...
    
    Results <- RunKLt(edgelist.time, as.matrix(edgelist.total[,1:2], ncol = 2), as.vector(edgelist.total[,3]), maxComms, directed, klPerNetwork, degreeCorrect, N)

}
