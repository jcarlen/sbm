## sbm.R: Running extended Karrar SBM code via R
##
## License? (You should have received a copy of the GNU General Public License
## along with RcppExamples.  If not, see <http://www.gnu.org/licenses/>.)

#' @KLPerNetwork this is the number of KL runs on a network

sbm <- function(edgelist, maxComms = 2, degreeCorrect = 0, directed = FALSE, klPerNetwork = 50, seedComms = NULL, seed = NULL)

{
    
    if (is.null(dim(edgelist)) || ncol(edgelist) < 2) {
        stop("edgelist must have at least two columns, three if weights are supplied")
    }
    
    if (ncol(edgelist) == 3) {
        weights = edgelist[,3]
        edgelist = edgelist[,1:2]
    }
    
    if (ncol(edgelist ==2)) {
        weights = rep(1, nrow(edgelist))
    }
 
    #seed?
    
    ## Make the call...
    Results <- sbmFit(edgelist-1, maxComms, degreeCorrect, directed, klPerNetwork, weights)
    
    cat("\nResults!\n")
    Results$EdgeMatrix = matrix(Results$EdgeMatrix, maxComms, maxComms)
    attr(Results, "degreeCorrect") <- degreeCorrect == 1
    attr(Results, "directed") <- directed
    attr(Results, "klPerNetwork") <- klPerNetwork
    
    Results
}

