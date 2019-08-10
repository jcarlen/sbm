# Functions to evaluate log likelihood for TDD-SBM and TDMM-SBM 
# Used to compare output ofdiscrete and mixed models, and different estimation methods
# Although sbmt returns a "highest scores" it's an un-normalized llik that's not comparable to mixed model results
# Jane Carlen
# TO DO: add unirected options for these functions

# ------------------------------- Helpers ----------------------------------------------

#' Convert time-sliced edgelist to time-sliced adjacency matrix (helper for likelihood functions)
#
#' @param edgelist A (time) series of network data represented as a list of edgelists. Assumes all edgelist slices have same names and number of columns, first two columns designate "from" and "to" for edges and third, if present, is count.
#' @param selfedges If true, include self-edges in edgelist in converted adjacency matrix. If false, diagonal of adjaceny matrix is zero.
#' @param as.array If true, return an N x N x Time array instead of a list of adjacency matrices.
#' @param directed Are the edges in the edgelist directed?
#'
edgelist_to_adj <- function(edgelist, selfedges = TRUE, as.array = TRUE, directed = TRUE) {
  Time = length(edgelist)
  E = do.call("rbind", edgelist)
  S = unique(c(as.character(E[,1]), as.character(E[,2])))
  N = length(S)
  Nodes = 1:N; names(Nodes) = S
  adjlist = lapply(1:Time, function(t) {
    #add count if not present
    if (ncol(edgelist[[t]]) == 2) {edgelist[[t]][,3] = 1}
    names(edgelist[[t]]) = c("from", "to", "count")
    # aggregate by "from" "to" pairs
    edgelist[[t]] = stats::aggregate(count ~ from + to, data = edgelist[[t]], sum)
    # convert to adjacency matrix
    adjlist.t = matrix(0,N,N, dimnames = list(S, S));
    adjlist.t[ cbind(   Nodes[as.character(edgelist[[t]][,1])], 
                        Nodes[as.character(edgelist[[t]][,2])] ) ] = edgelist[[t]][,3]
    # remove self-edges?
    if (selfedges == FALSE) {diag(adjlist.t) = 0}
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

#' Calculate the un-normalized log-likelihood for a single timeslice given mu as input (works for discrete or continuous). Helper for tdd_sbm_llik and tdmm_sbm_llik
#' #
#' @param A_t is a single-time network, represented as a N x N adjacency matrix
#' @param mu mu is a N x N matrix of edge expected values at this time slice
#' @param directed Are the edges in the edgelist directed?
#' @param selfEdges is whether to sum over self-edge indices in the likelihood calculation, or exclude them
#' 
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

#' Calculate number of parameters for discrete model
#
#' @param N number of nodes
#' @param K number of blocks
#' @param Time number of time slices
#' @param directed Is the network directed?
tdd_n_param <- function(N, K, Time, directed = TRUE) {
  ifelse(directed, 2*N - K + Time*K^2, 2*N - K + Time*(K*(K-1)/2))
}

#' Calculate number of parameters for mixed model
#
#' @param N number of nodes
#' @param K number of blocks
#' @param Time number of time slices
#' @param directed Is the network directed?
tdmm_n_param <- function(N, K, Time, directed = TRUE) {
  ifelse(directed, K*N - K + Time*K^2, K*N - K + Time*(K*(K-1)/2))
}

# ---------------------------------- Functions (discrete membership TDD-SBM) -----------------------------------------

#' Calculate un-normalized log-likelihood for a TDD-SBM (time-dependent discrete-membership stochastic block model).
#' This is not the same as the "highest score" returned by sbmt, which is unnormalized in a different way
#' 
#
#' @param A (time) series of network data represented  as N x N x Time array.
#' @param roles is a length-N list of estimated block assignment for each node
#' @param omega is a Time x K x K array describing block-to-block traffic at each time period or Time-length list of K x K matrices.
#' @param directed Is the network directed?
#' @param selfEdges is whether to sum over self-edge indices in the likelihood calculation, or exclude them
#' @examples
#' data("la_byhour_edgelist")
#' A = edgelist_to_adj(la_byhour_edgelist, as.array = TRUE)
#' model1 = sbmt(la_byhour_edgelist,  degreeCorrect = 3, directed = TRUE, klPerNetwork = 2, maxComms = 3)
#' tdd_sbm_llik(A, model1$FoundComms, model1$EdgeMatrix)
#' 
tdd_sbm_llik <- function(A, roles, omega, directed = TRUE, selfEdges = TRUE) {
  
  N = dim(A)[1]
  Time = dim(A)[3]
  K = length(unique(roles))
  
  omega = array(unlist(omega), dim = c(K,K,Time)) #make array if not already
  roles = roles[rownames(A)] # put into same order as A if not already
  
  A_total = apply(A, c(1,2), sum)
  degree_total = colSums(A_total) + rowSums(A_total)
  
  #may need to catch more cases her, e.g. factor roles 
  role_degree = data.frame(id = names(roles), role = roles + 1, degree_total = degree_total[names(roles)])
  role_degree$role_sum = aggregate(role_degree$degree_total, by = list(role_degree$role), sum)$x[role_degree$role]
  theta = role_degree$degree_total/role_degree$role_sum
  
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
#' @param C is a N x K matrix of mixed group membership whose columns sum to 1
#' @param omega is a Time x K x K array describing block-to-block traffic at each time period or Time-length list of K x K matrices.
#' @param directed Is the network directed?
#' @param selfEdges is whether to sum over self-edge indices in the likelihood calculation, or exclude them
#' @examples
#' # Uses output from python implementation (slighly reformatted)
#' data("la_byhour_edgelist", "la_mixed_roles_2", "la_mixed_omega_2")
#' A = edgelist_to_adj(la_byhour_edgelist, as.array = TRUE)
#' tdmm_sbm_llik(A, la_mixed_roles_2, la_mixed_omega_2)
#' 

tdmm_sbm_llik <- function(A, C, omega, directed = TRUE, selfEdges = TRUE) {
  
  N = dim(A)[1]
  Time = dim(A)[3]
  K = ncol(C)
  if ( !identical(sort(rownames(C)), sort(rownames(A))) ) {stop("A and C should have matching rownames")}
  C = as.matrix(C[rownames(A),])
  omega = array(unlist(omega), dim = c(K,K,Time)) #time period starts as a row, ends up as [,,t]
  
  lik = sum(
    sapply(1:Time, function(t) {
      mu_t = C %*% omega[,,t] %*% t(C)
      td_sbm_llik_t(A[,,t], mu_t, directed, selfEdges)
    })
  ) 
  return(lik)
}



