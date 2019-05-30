## sbm.R: Running extended Karrar SBM code via R
##
## License? (You should have received a copy of the GNU General Public License
## along with RcppExamples.  If not, see <http://www.gnu.org/licenses/>.)

## To do: Check that seed can be set for reproducibility
## Make seedComms usable for non-random initialization

#' Estimate parameters of a static SBM by Kernighan-Lin algorithm. Allows for directed networks and degree correction.
#
# '@maxComms maximum number of communities represented in the network
# '@degreeCorrect whether to apply the degree correction of Karrer and Newman
# '@directed whether the network is directed and off-diagonal block parameters may be assymetrical
# '@KLPerNetwork this is the number of KL runs on a network
# '@tolerance stopping criteria for KL. Prevents loops due to numerical errors.

sbm <- function(edgelist, maxComms = 2, degreeCorrect = F, directed = FALSE, klPerNetwork = 50, tolerance = 1e-4, seedComms = NULL, seed = NULL)

{
    # Argument Checks
    if (is.null(dim(edgelist)) || ncol(edgelist) < 2) {
        stop("edgelist must have at least two columns, three if weights are supplied")
    }
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
    
    ## Make the call...
    Results <- sbmFit(edgelist1-1, maxComms, degreeCorrect, directed, klPerNetwork, weights, tolerance)
    
    cat("\nResults!\n")
    
    names(Results$FoundComms) = names(link.nodes) # Return found community membership in order of ID
    Results$EdgeMatrix = matrix(Results$EdgeMatrix, maxComms, maxComms, byrow = T)
    Results$degreeCorrect = degreeCorrect
    Results$directed = directed
    Results$klPerNetwork = klPerNetwork
    Results$weighted = weighted
    Results$tolerance = tolerance
    
    Results
}
