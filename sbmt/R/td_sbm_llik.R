# Functions to evaluate log likelihood for TDD-SBM and TDMM-SBM 
# Used to compare output of discrete and mixed models, and different estimation methods
# Although sbmt returns a "highest scores" it's an un-normalized llik that's not comparable to mixed model results
# Jane Carlen
# TO DO: add unirected options for these functions
#        add tests for adj_to_edgelist in test_likelihood (at least for now)

# ------------------------------- Helpers ----------------------------------------------

#' Convert a multilayer edgelist to a multilayer adjacency matrix. 
#' 
#' Convert a multilayer edgelist to a multilayer adjacency matrix (as a list of adjacency matrices or N x N x T array depending on the `as.array` parameter). 
#' Used in pre-processing for likelihood functions.
#' 
#' @param edgelist.time A (time) series of networks represented as a list of edgelists. 
#' Assumes all edgelist slices have the same names and number of columns. The first two columns designate edges "from" and "to", and the third, if present, is the count (or more generally the weight) for that edge.
#' @param selfEdges If true, include self-edges in converted adjacency matrix. If false, diagonal of adjaceny matrix is zero.
#' @param as.array If true, return an N x N x Time array instead of a list of adjacency matrices.
#' @param directed Are the edges in the edgelist.time directed?

edgelist_to_adj <- function(edgelist.time, selfEdges = TRUE, as.array = TRUE, directed = TRUE) {
  Time = length(edgelist.time)
  E = do.call("rbind", edgelist.time)
  S = unique(c(as.character(E[,1]), as.character(E[,2])))
  N = length(S)
  Nodes = 1:N; names(Nodes) = S
  adjlist = lapply(1:Time, function(t) {
    #add count if not present
    if (ncol(edgelist.time[[t]]) == 2) {edgelist.time[[t]][,3] = 1}
    names(edgelist.time[[t]]) = c("from", "to", "count")
    # aggregate by "from" "to" pairs
    edgelist.time[[t]] = stats::aggregate(count ~ from + to, data = edgelist.time[[t]], sum)
    # convert to adjacency matrix
    adjlist.t = matrix(0,N,N, dimnames = list(S, S));
    adjlist.t[ cbind(   Nodes[as.character(edgelist.time[[t]][,1])], 
                        Nodes[as.character(edgelist.time[[t]][,2])] ) ] = edgelist.time[[t]][,3]
    # remove self-edges?
    if (selfEdges == FALSE) {diag(adjlist.t) = 0}
    if (directed == FALSE) {
      adjlist.t = adjlist.t+t(adjlist.t)
      diag(adjlist.t) = diag(adjlist.t)/2
    }
    
    return(adjlist.t)
  })
  
  if (as.array == TRUE) {
    adjlist = array(unlist(adjlist), dim = c(N,N,Time), dimnames = list(S, S, 1:Time))
  }
  
 return(adjlist)
  
}

#' Convert a multilayer adjacenty matrix to a multilayer edgelist.
#' 
#' Convert representation of time-sliced network as N x N x T array to a length-T list of edglists for each time period. (Handle NA?)
#' Resulting edgelists have three columns, "from", "to", and "count" (which is more generally the edge weight).
#' 
#' @param A is a (time) series of network data represented as a N x N x Time array (each slice represented as an adjacency matrix).
#' @param directed Are the edges in the edgelist directed?
#' @param selfEdges If true, include self-edges in output edgelists. If false, remove. Note tdsbm methods allow/include selfedges.
#' @param removeZeros Probably want to remove zeros for efficiency, but maybe not if edgelist-length consistency over time is desired

adj_to_edgelist <- function(A, directed = FALSE, selfEdges = TRUE, removeZeros = TRUE) {
  discrete_edge_list = apply(A, 3, function(x) {
    N = dim(A)[1]
    indices = expand.grid(1:N,1:N)
    edge_list = data.frame(cbind(indices, as.vector(x)))
    names(edge_list) = c("from", "to", "count")
    if (!directed) { edge_list = edge_list[edge_list$from <= edge_list$to,] }
    if (!selfEdges) { edge_list = edge_list[edge_list$from != edge_list$to,] }
    if (removeZeros) {edge_list = edge_list[edge_list$count>0, ]}
    return(edge_list)
  })
  return(discrete_edge_list)
}

#' Calculate the un-normalized td-sbm log-likelihood for a single network slice.
#' 
#' Calculate the un-normalized td-sbm log-likelihood for a single network slice given mu as input (works for discrete or continuous). Helper for tdd_sbm_llik and tdmm_sbm_llik
#' 
#' @param A_t is a single-time network, represented as a N x N adjacency matrix
#' @param mu mu is a N x N matrix of edge expected values at this time slice
#' @param directed Are the edges in the edgelist directed?
#' @param selfEdges is whether to sum over self-edge indices in the likelihood calculation, or exclude them

td_sbm_llik_t <- function(A_t, mu, directed = TRUE, selfEdges = TRUE) {
  if (!selfEdges) {
    diag(mu) = 1
  }
  term0 = A_t * log(mu)
  #print(A_t[!is.finite(term0)]) #<- should be empty (only return numeric(0))
  term0[!is.finite(term0)] = 0 # note we let 0*log(0) = 0
  term1 = ifelse(directed, sum(term0), sum(term0[upper.tri(term0, diag = TRUE)]))
  if (!selfEdges) {
    diag(mu) = 0
  }
  term2 = ifelse(directed, sum(mu), sum(mu[upper.tri(mu, diag = TRUE)]))
  #term3 = sum(log(factorial(A_t)))
  return(term1 - term2)
}

#' Calculate number of parameters for tdd-sbm
#
#' @param N number of nodes
#' @param K number of blocks
#' @param Time number of time slices
#' @param directed Is the network directed?

tdd_n_param <- function(N, K, Time, directed = TRUE) {
  ifelse(directed, 2*N - K + Time*K^2, 2*N - K + Time*(K*(K-1)/2))
}

#' Calculate number of parameters for tdmm-sbm
#
#' @param N number of nodes
#' @param K number of blocks
#' @param Time number of time slices
#' @param directed Is the network directed?

tdmm_n_param <- function(N, K, Time, directed = TRUE) {
  ifelse(directed, K*N - K + Time*K^2, K*N - K + Time*(K*(K-1)/2))
}

# ---------------------------------- Functions (discrete membership TDD-SBM) -----------------------------------------

#' Calculate un-normalized log-likelihood for a TDD-SBM
#' 
#' Calculate un-normalized log-likelihood for a TDD-SBM (time-dependent discrete-membership stochastic block model).
#' This is not the same as the "highest score" returned by sbmt, which is unnormalized in a different way
#' 
#' @param A (time) series of network data represented  as N x N x Time array.
#' @param roles is a length-N list of estimated block assignment for each node. If not 0-indexed will substract min value so roles are 0,1,2,...
#' @param omega is a Time x K x K array describing block-to-block traffic at each time period OR Time-length list of K x K matrices.
#' @param degreeCorrect degree correction type:
#' \itemize{
#'   \item 0 = no degree correction
#'   \item 1 = the degree correction of Karrer and Newman at every time step - inherits directedness of graph
#'   \item 2 = time-independent degree correction over time-dependent data - inherits directedness of graph
#'   \item 3 = time-independent degree correction over time-dependent data - undirected regardless
#'}
#' @param directed Is the network directed?
#' @param selfEdges is whether to sum over self-edge indices in the likelihood calculation, or exclude them
#' @examples
#' data("la_byhour_edgelist")
#' A = edgelist_to_adj(la_byhour_edgelist, as.array = TRUE)
#' model1 = sbmt(la_byhour_edgelist,  degreeCorrect = 3, 
#'   directed = TRUE, klPerNetwork = 2, maxComms = 3)
#' tdd_sbm_llik(A, model1$FoundComms, model1$EdgeMatrix)

tdd_sbm_llik <- function(A, roles, omega, degreeCorrect = 3, directed = TRUE, selfEdges = TRUE) {
  
  # arg checks
  if (min(roles)!=0) {
    roles = roles - min(roles)
    warning(
      paste("roles should be 0-indexed. Re-indexed to min 0 (by subtracting min value).\n new roles: ", paste(roles, collapse = " "), "\n")
    )
  
  }
  if (length(intersect(class(omega), c("list", "array")))==0) {stop("omega should be a Time x K x K array or Time-length list of K x K matrices")}
  if (degreeCorrect %in% c(1,2)) {stop("tdd_sbm_llik not yet implement for degree correction types 1 or 2")}
  
  
  N = dim(A)[1]
  Time = dim(A)[3]
  K = ifelse("list" %in% class(omega), nrow(omega[[1]]), dim(omega)[1])

  # add node names 1:N if none
  if(is.null(names(roles))) {roles = setNames(roles, 1:N)} 
     
  omega = array(unlist(omega), dim = c(K,K,Time)) #make array if not already
  roles = roles[rownames(A)] # put into same order as A if not already
  
  A_total = apply(A, c(1,2), sum)
  degree_total = colSums(A_total) + rowSums(A_total)
  
  #may need to catch more cases her, e.g. factor roles 
  role_degree = data.frame(id = names(roles), role = roles + 1, degree_total = degree_total[names(roles)])
 
  if (degreeCorrect == 3) { #undirected time-independent degree correction
    role_degree$role_sum = aggregate(role_degree$degree_total, by = list(role_degree$role), sum)$x[role_degree$role]
    role_degree$frac = 1/table(roles)[role_degree$role]
    theta = role_degree$degree_total/role_degree$role_sum
  } 
  if (degreeCorrect == 0) { 
    theta = 1/table(roles)[role_degree$role]  #undirected time-independent degree correction
  }

  lik = sum(
          sapply(1:Time, function(t) {
            mu_t = diag(theta) %*% matrix( (omega[,,t])[as.matrix(expand.grid(role_degree$role, role_degree$role))], N, N) %*% diag(theta)
            td_sbm_llik_t(A[,,t], mu_t, directed, selfEdges)
          })
        ) 
  
  return(lik)
}

# ------------------------------- Functions (mixed-membership TDMM-SBM) ----------------------------------------------

#' Calculate un-normalized log-likelihood for a TDMM-SBM (time-dependent mixed-membership stochastic block model)
#
#' @param A is a A (time) series of network data represented  as N x N x Time array.
#' @param C is a N x K matrix of mixed group membership whose columns sum to 1. If no rownames (indicating corresponding nodes) they will be assigned 1:N.
#' @param omega is a K x K x Time array describing block-to-block traffic at each time period, or Time-length list of K x K matrices, or a K^2 x Time array.
#' All types will be internally convereted to a Time x K x K array 
#' @param directed Is the network directed?
#' @param selfEdges is whether to sum over self-edge indices in the likelihood calculation, or exclude them
#' @examples
#' # Uses output from python implementation (slighly reformatted)
#' data("la_byhour_edgelist", "la_mixed_roles_2", "la_mixed_omega_2")
#' A = edgelist_to_adj(la_byhour_edgelist, as.array = TRUE)
#' tdmm_sbm_llik(A, la_mixed_roles_2, la_mixed_omega_2)

tdmm_sbm_llik <- function(A, C, omega, directed = TRUE, selfEdges = TRUE) {
  
  
  N = dim(A)[1]
  Time = dim(A)[3]
  K = ncol(C)
  # add node names 1:N if none
  if(is.null(rownames(C))) {rownames(C) = 1:N} 
  
  if ( !identical(sort(rownames(C)), sort(rownames(A))) ) {stop("A and C should have matching rownames")}
  C = as.matrix(C[rownames(A),])
  
  #check omega
  if (! (  identical(as.numeric(dim(omega)), as.numeric(c(K, K, Time))) | 
           identical(as.numeric(dim(omega)), as.numeric(c(K^2, Time))) | 
           (length(omega) ==Time && identical( as.numeric(dim(omega[[1]])), as.numeric(c(K,K)))) ) ) {
    warning("omega should be a K x K x Time array, K^2 x Time array, or Time-length list of K x K matrices")
  }
  omega = array(unlist(omega), dim = c(K,K,Time))
  
  lik = sum(
    sapply(1:Time, function(t) {
      mu_t = C %*% omega[,,t] %*% t(C)
      td_sbm_llik_t(A[,,t], mu_t, directed, selfEdges)
    })
  ) 
  return(lik)
}



