# Functions to evaluate log likelihood for TDD-SBM and TDMM-SBM 
# Used to compare output ofdiscrete and mixed models, and different estimation methods
# Although sbmt returns a "highest scores" it's an un-normalized llik that's not comparable to mixed model results
# Jane Carlen
# TO DO: add unirected options for these functions

# ------------------------------- Helper (convert time-sliced edgelist to array) ----------------------------------------------

#' @edgelist A (time) series of network data represented as a list of edgelists. Assumes all edgelist slices have same names and number of columns, first two columns designate "from" and "to" for edges and third, if present, is count.
#' @selfects If true, allow non-zero diagonal of converted adjacency matrices. If false zeros out the diagonal.
#' @as.array If true, return an N x N x Time array instead of a list of adjacency matrices.

edeglist_to_adj <- function(edgelist, selfedges = FALSE, as.array = TRUE, directed = TRUE) {
  Time = length(edgelist)
  E = do.call("rbind", edgelist)
  S = unique(c(as.character(E[,1]), as.character(E[,2])))
  N = length(S)
  Nodes = setNames(1:N, S)
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
    return(adjlist.t)
    
  })
  
  if (as.array == TRUE) {
    A = array(unlist(adjlist), dim = c(N,N,Time), dimnames = list(S, S, 1:Time))
    return(A)
  } else {return(adjlist)}  
}


# Single timeslice likelihood given mu as input (works for discrete or continuous)
#' @A_t is a single time network, represented as a N x N adjacency mtrix
#' 
td_sbm_llik_t <- function(A_t, mu, K = 2) {
  diag(mu) = 1
  term0 = A_t * log(mu)
  #print(A_t[!is.finite(term0)]) #<- should be empty (only return numeric(0))
  term0[!is.finite(term0)] = 0 # note we're letter 0*log(0) = 0
  term1 = sum(term0)
  diag(mu) = 0
  term2 = sum(mu)
  #term3 = sum(log(factorial(A_t)))
  return(term1 - term2)
}

# ------------------------------- Functions (mixed-membership TDMM-SBM) ----------------------------------------------

#' @C is a N x K matrix of mixed group membership whose columsn sum to 1
#' @A is a A (time) series of network data represented  as N x N x Time array.
#' @omega is a Time x K x K array describing block-to-block traffic at each time period
#' @K is the number of blocks
#' 

tdmm_sbm_llik <- function(A, C, omega, K = ncol(C), N = dim(A)[1], Time = dim(A)[3]) {
  lik = sum(
    sapply(1:Time, function(t) {
      w_t = matrix(omega[,,t], K, K)
      mu_t = C %*% w_t %*% t(C)
      td_sbm_llik_t(A[,,t], mu_t, K = K)
    })
  ) 
  return(lik)
}

# ---------------------------------- Functions (discrete membership TDD-SBM) -----------------------------------------

#' @A is a A (time) series of network data represented  as N x N x Time array.

tdd_sbm_llik <- function(A, roles, omega,  K = length(unique(roles)), N = dim(A)[1], Time = dim(A)[3]) {
  
  A_total = apply(A, c(1,2), sum)
  degree_total = colSums(A_total) + rowSums(A_total)
  
  role_degree = left_join(data.frame(id = names(roles), role = roles + 1), #may need to catch more cases her, e.g. factor roles
                          data.frame(id = names(degree_total), degree_total), by = "id")
  
  theta = (role_degree %>% group_by(role) %>% mutate(theta = degree_total/sum(degree_total)))$theta
  
  lik = sum(
          sapply(1:Time, function(t) { #small constant added to omega prevents degenerate calculation
            mu_t = diag(theta) %*% matrix( (omega[,,t])[as.matrix(expand.grid(role_degree$role, role_degree$role))], N, N) %*% diag(theta)
            td_sbm_llik_t(A[,,t], mu_t, K)
          })
        ) 
  return(lik)
}

# --------------------------------------- Examples--------------------------------------------------------------
# 

# data("la_byhour_edgelist")
# data1 = la_byhour_edgelist
# A = edeglist_to_adj(data1, as.array = TRUE)
# Time = dim(A)[3]
#  
# # discrete ----
# 
# model1 =  sbmt(la_byhour_edgelist,  degreeCorrect = 3, directed = T, klPerNetwork = 5, maxComms = 3)
# roles = model1$FoundComms[rownames(A)]
# K = length(unique(roles))
# omega = array(unlist(model1$TimeMatrices), dim = c(K,K,Time))
# tdd_sbm_llik(A, roles[rownames(A)], omega, K = length(unique(roles)), N = dim(A)[1], Time = dim(A)[3])
# 
# # mixed ----
# # from python output
# omega1 = ny_hm_continuous.omega.3
# roles1 = ny_hm_continuous.roles.3
# omega1 = ny_hm_continuous.omega.2
# roles1 = ny_hm_continuous.roles.2
#  
# K = sqrt(ncol(omega1))
# omega = array(unlist(t(omega1)), dim = c(K,K,Time))
# C = roles1[  match(rownames(A[,,1]), as.character(roles1[,1])),  ]
# rownames(C) = C[,1]; C = as.matrix(C[,2:ncol(C)])
#  
# # calculate likelihood
# tdmm_sbm_llik(A, C, omega)

 
# # check omega  = 
# array(unlist(
#    lapply(1:24, function(t) {
#      tmp = left_join(expand.grid(0:(K-1),0:(K-1)),
#        aggregate(data1[[t]]$x, by = list(roles[as.character(data1[[t]][,1])], 
#                                        roles[as.character(data1[[t]][,2])]), sum),
#        by = c("Var1" = "Group.1", "Var2" = "Group.2"))
#       matrix(tmp[,3], K, K)
#    })),
#   dim = c(K,K,Time))
# 
# omega[is.na(omega)] = 0
