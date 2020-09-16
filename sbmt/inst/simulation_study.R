# TO DO
# add generate_multilayer_array to package to facilitate simulation?
# Note degree correcton can also lead to a more parsimonious and interpretable model representation where there is degree heterogeneity
# because a unique class is not needed for each degree-activity leve
# ---------------------------------------------------------------------------------------------------------------
# new functions ----

# generate N x N x T array (time-sliced adjacency matrices) based on TDD-SBM (type == "discrete) or TDMM-SBM (type = "mixed") model
# defaults to no degree correction (all dc_factors == 1)
# for type "discrete" roles are a length-N vector of block assignments. For type "mixed" roles are a G x N matrix of assignment weights (and dc_factor is ignored)
# (add this as a simulation function to abmt package?)
generate_multilayer_array <- function(N, Time, roles, omega, dc_factors = rep(1, N), type = "discrete") {
  # checks
  if (type == "discrete" & !identical(dc_factors, rep(1, N))) {
   if (! identical(aggregate(dc_factors ~ roles, FUN = "sum")[,2], rep(1, length(unique(roles))))) {
     warning("degree correction factors not normalized to sum 1 by group")
   }
  }
  
  edge_array = array(0, dim = c(N, N, Time))
  for (i in 1:N) {
    role_i = roles[i]
    for (j in 1:N) {
      role_j = roles[j]
      for (time in 1:Time) {
          if (type == "discrete") { ijt = rpois(dc_factors[i]*dc_factors[j]*omega[role_i, role_j, time], n= 1) }
          if (type == "mixed") { ijt = rpois(roles[, i]%*% omega[, , time] %*% t(roles[, j]), n= 1) }
          edge_array[i, j, time] = ijt
      } 
    }
  }
  return(edge_array)
}

# libraries ----
# devtools::install_github("jcarlen/sbm", subdir = "sbmt") 
library(sbmt)
library(fossil) #for adj rand index

# ---------------------------------------------------------------------------------------------------------------
# parameters  ----

n_roles = 2
Time = 16 #use a power of two for compatibility with ppsbm hist method
N = 30

a = 10
b = 5
ymax =  max(2*a, 2*b)

#roles
roles_discrete = rep(1:n_roles, length.out = N)

# curves ----
par(mfrow = c(n_roles, n_roles))
curve(b*sin(x*pi/Time) + b, 0, Time, ylim = c(0, ymax))
curve(a*sin(x*2*pi/Time)+a, 0, Time, ylim = c(0, ymax))
curve(-a*sin(x*2*pi/Time)+a, 0, Time, ylim = c(0, ymax))
curve(b*sin(x*pi/Time)+ b, 0, Time, ylim = c(0, ymax))

# show overlapping
par(mfrow = c(1,1))
curve(a*sin(x*2*pi/Time)+a, 0, Time, ylim = c(0, ymax))
curve(-a*sin(x*2*pi/Time)+a, 0, Time, col = "red", add = TRUE)
curve(b*sin(x*pi/Time) + b, 0, Time, col = "blue", add = TRUE)

x = seq(.5, Time)

omega_11 = b*sin(x*pi/Time) + b
omega_12 =  a*sin(x*2*pi/Time)+a
omega_21 = -a*sin(x*2*pi/Time)+a
omega_22 = omega_11

omega = array(rbind(omega_11, omega_21, omega_12, omega_22), dim = c(n_roles, n_roles, Time)) #left most index moves fastest

# for degree corrected ----

dc_factor = seq(0,1,length.out = n_roles+1)[-1]
#dc_factor = seq(0,1,length.out = n_roles+5)[-1]

dc_factors = rep(dc_factor, each = round(N/n_roles))[1:N]
#apply sum to 1 constraint
f.tmp = function(v) {v/sum(v)}
dc_factors = as.vector(aggregate(dc_factors ~ roles_discrete, FUN = "f.tmp")[,-1])
# check 
identical(aggregate(dc_factors ~ roles_discrete, FUN = "sum")[,2], rep(1, n_roles))
  
# adjust omega after apply sum 1 constraint to degree-correction factors
dc_omega = omega*array(table(roles_discrete) %*% t(table(roles_discrete)), dim = c(2,2,16))

# ---------------------------------------------------------------------------------------------------------------
# no degree correction case ----
N_sim = 10
degree_correct = 0
role_results = 1:N_sim

  # # checks 
  # discrete_edge_array = generate_multilayer_array(N, Time, roles_discrete, omega)
  # dim(discrete_edge_array)
  # par(mfrow = c(n_roles, n_roles))
  # plot(apply(discrete_edge_array[roles_discrete==1,roles_discrete==1,],3,mean), type = "l")
  # plot(apply(discrete_edge_array[roles_discrete==1,roles_discrete==2,],3,mean), type = "l")
  # plot(apply(discrete_edge_array[roles_discrete==2,roles_discrete==1,],3,mean), type = "l")
  # plot(apply(discrete_edge_array[roles_discrete==2,roles_discrete==2,],3,mean), type = "l")
  # image(discrete_edge_array[(1:N)[order((1:N)%%n_roles)], (1:N)[order((1:N)%%n_roles)], 1])
  # image(discrete_edge_array[(1:N)[order((1:N)%%n_roles)], (1:N)[order((1:N)%%n_roles)], 7])
  # image(discrete_edge_array[(1:N)[order((1:N)%%n_roles)], (1:N)[order((1:N)%%n_roles)], 13])
  # image(discrete_edge_array[(1:N)[order((1:N)%%n_roles)], (1:N)[order((1:N)%%n_roles)],19])

# - fit sbmt ----

set.seed(1)
sbmt_ari = 1:N_sim

for (s in 1:N_sim) {
  
  discrete_edge_array = generate_multilayer_array(N, Time, roles_discrete, omega)
  discrete_edge_list = adj_to_edgelist(discrete_edge_array, directed = TRUE, selfEdges = TRUE)
  discrete_sbmt = sbmt(discrete_edge_list, maxComms = 2, degreeCorrect = 0, directed = TRUE, klPerNetwork = 10)
  plot(discrete_sbmt)
  sbmt_ari[s] = adj.rand.index(discrete_sbmt$FoundComms[order(as.numeric(names(discrete_sbmt$FoundComms)))], roles_discrete)
}

# - results ----

# role detection
sbmt_ari

# omega detection
plot(discrete_sbmt)
# overall curve distances?

# ---------------------------------------------------------------------------------------------------------------
# degree correction case  ----

set.seed(1)
dc_ari = 1:N_sim #evaluate with adjusted rand index

# - fit sbmt ----

for (s in 1:N_sim) {
  
  dc_discrete_edge_array = generate_multilayer_array(N, Time, roles_discrete, dc_omega, dc_factors)
  dc_discrete_edge_list = adj_to_edgelist(dc_discrete_edge_array, directed = TRUE, selfEdges = TRUE)
  dc_discrete_sbmt = sbmt(dc_discrete_edge_list, maxComms = 2, degreeCorrect = 3, directed = TRUE, klPerNetwork = 10)
  # to compare likelihood
  # dc_discrete_sbmt2 = sbmt(dc_discrete_edge_list, maxComms = n_roles*2, degreeCorrect = 3, directed = TRUE, klPerNetwork = 10)
  # diff in llik vs. diff in param
  # (dc_discrete_sbmt2$llik - dc_discrete_sbmt$llik)/(tdd_n_param(N, n_roles*2, Time) - tdd_n_param(N, n_roles, Time)) 
  plot(dc_discrete_sbmt)
  dc_ari[s] = adj.rand.index(dc_discrete_sbmt$FoundComms[order(as.numeric(names(dc_discrete_sbmt$FoundComms)))], roles_discrete)
}


# - results ----

# role detection
dc_ari

# omega detection
plot(dc_discrete_sbmt)
# overall curve distances?

# ---------------------------------------------------------------------------------------------------------------
# ppsbm ----
library(ppsbm)

# no degree correction case. their model works as expected ----
# - fit ppsbm ----

# Use the "hist" method because agrees more closely with out discrete time slices and requires little data manipulation

Nijk = sapply(adj_to_edgelist(discrete_edge_array, directed = TRUE, selfEdges = FALSE, removeZeros = FALSE), "[[", 3); dim(Nijk)
discrete_ppsbm = mainVEM(list(Nijk=Nijk,Time=Time), N, Qmin = 1, Qmax = 4, directed=TRUE, 
                         method='hist', d_part=5, n_perturb=10, n_random=0)

# - results ----

# number of blocks selected
selected_Q = modelSelection_Q(list(Nijk=Nijk,Time=Time), N, Qmin = 1, Qmax = 4, directed = TRUE, sparse = FALSE, discrete_ppsbm)$Qbest
selected_Q
selected_Q == n_roles #should equal n_roles

# role detection
apply(discrete_ppsbm[[selected_Q]]$tau, 2, which.max)
adj.rand.index(apply(discrete_ppsbm[[selected_Q]]$tau, 2, which.max), roles_discrete)

# omegas
par(mfrow = c(n_roles, n_roles))
apply(exp(discrete_ppsbm[[selected_Q]]$logintensities.ql), 1, plot, type = "l")

# degree correction case ----
# - fit ppsbm ----

dc_Nijk = sapply(adj_to_edgelist(dc_discrete_edge_array, directed = TRUE, selfEdges = FALSE, removeZeros = FALSE), "[[", 3); dim(dc_Nijk)
dc_discrete_ppsbm = mainVEM(list(Nijk=dc_Nijk,Time=Time), N, Qmin = 1, Qmax = 6, directed=TRUE, 
                            method='hist', d_part=5, n_perturb=10, n_random=0)
# - results ----

# number of blocks selected
selected_Q = modelSelection_Q(list(Nijk=dc_Nijk,Time=Time), N, Qmin = 1, Qmax = 6, directed = TRUE, sparse = FALSE, dc_discrete_ppsbm)$Qbest
selected_Q
selected_Q == n_roles

# role detection

# Wants a seperate class for each degree-correcton level
apply(dc_discrete_ppsbm[[selected_Q]]$tau, 2, which.max)
adj.rand.index(apply(dc_discrete_ppsbm[[selected_Q]]$tau, 2, which.max), roles_discrete)

# omegas
par(mfrow = c(selected_Q, selected_Q)); par(mai = rep(.5, 4))
apply(exp(dc_discrete_ppsbm[[selected_Q]]$logintensities.ql), 1, plot, type = "l", col = "blue")

# with true number of groups? gets it right
apply(dc_discrete_ppsbm[[n_roles]]$tau, 2, which.max)
adj.rand.index(apply(dc_discrete_ppsbm[[n_roles]]$tau, 2, which.max), roles_discrete)
s
par(mfrow = c(n_roles, n_roles)); par(mai = rep(.5, 4))
apply(exp(dc_discrete_ppsbm[[n_roles]]$logintensities.ql), 1, plot, type = "l", col = "blue")

# Bike example shows how degree correction ib model -> group statins with similar behavior across activity levels
# ---------------------------------------------------------------------------------------------------------------
# mixed-membership example ----

roles_mixed = matrix(c(rep(.5, 2*N/3), rep(c(0,1), N/3), rep(c(1,0), N/3)), nrow = n_roles, ncol = N)

