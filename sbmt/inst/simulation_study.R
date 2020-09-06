# TO DO
# add generate_multilayer_array and dgelist as function to package 
# to enable simulation and complement edgelist_to_adj
# ---------------------------------------------------------------------------------------------------------------
# new functions ----

# add this as a simulation function to package?
# defaults to no degree correction (all have multilplier 1)
generate_multilayer_array <- function(N, Time, roles, omega, dc_factor = rep(1, N)) {
  discrete_edge_array = array(0, dim = c(N, N, Time))
  for (i in 1:N) {
    role_i = roles[i]
    for (j in 1:N) {
      role_j = roles[j]
      for (time in 1:Time) {
        ijt = rpois(dc_factor[i]*dc_factor[j]*omega[role_i, role_j, time], n= 1)
        discrete_edge_array[i, j, time] = ijt
      }
    }
  }
  return(discrete_edge_array)
}

# convert N x N x T edgelist array to length T list of edges at each time period (includes edges valued 0, ignores NA)  
adj_to_edgelist <- function(edge_array, directed = FALSE, selfEdge = FALSE) {
  discrete_edge_list = apply(edge_array, 3, function(x) {
    indices = data.frame(which(is.finite(x), arr.ind = TRUE))
    names(indices) = c("from", "to")
    if (!directed) { indices = indices[indices$from <= indices$to,] }
    if (!selfEdge) { edge_list = indices[indices$from != indices$to,] }
    edge_list = data.frame(indices, count = as.vector(x))
  })
}


# ---------------------------------------------------------------------------------------------------------------
# libraries ----
library(sbmt)

# parameters  ----

n_roles = 2
Time = 24 
N = 30

a = 10
b = 5
ymax =  max(2*a, 2*b)

# curves ----
par(mfrow = c(n_roles, n_roles))
curve(b*sin(x*pi/Time) + b, 0, Time, ylim = c(0, ymax))
curve(a*sin(x*2*pi/Time)+a, 0, Time, ylim = c(0, ymax))
curve(-a*sin(x*2*pi/Time)+a, 0, Time, ylim = c(0, ymax))
curve(b*sin(x*pi/Time)+ b, 0, Time, ylim = c(0, ymax))

# show overlapping
curve(a*sin(x*2*pi/Time)+a, 0, Time, ylim = c(0, ymax))
curve(-a*sin(x*2*pi/Time)+a, 0, Time, col = "red", add = TRUE)
curve(b*sin(x*pi/Time) + b, 0, Time, col = "blue", add = TRUE)

x = seq(.5, Time)

omega_11 = b*sin(x*pi/Time) + b
omega_12 =  a*sin(x*2*pi/Time)+a
omega_21 = -a*sin(x*2*pi/Time)+a
omega_22 = omega_11

omega = array(rbind(omega_11, omega_21, omega_12, omega_22), dim = c(n_roles, n_roles, Time)) #left most index moves fastest
    
# simulate roles and edges ----      

# no degree correction case ----
N_sim = 10
roles_discrete = rep(1:n_roles, length.out = N)
degree_correct = 0
role_results = 1:N_sim

  # # checks ----
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

# sbmt ----

set.seed(1)
for (s in 1:N_sim) {
  
  discrete_edge_array = generate_multilayer_array(N, Time, roles_discrete, omega)
  discrete_edge_list = adj_to_edgelist(discrete_edge_array, directed = TRUE, selfEdge = FALSE)
  discrete_sbmt = sbmt(discrete_edge_list, maxComms = 2, degreeCorrect = 0, directed = TRUE, klPerNetwork = 3)
  role_results[s] = N - sum(diag(table(discrete_sbmt$FoundComms[order(as.numeric(names(discrete_sbmt$FoundComms)))], roles_discrete)))
}

role_results

#omega results
par(mfrow = c(n_roles, n_redgeoles))
sapply(tmp$EdgeMatrix, as.vector) %>% apply(1, plot, type = "l")

# Degree correct case  ----

roles_discrete
dc_factors = rep(c(.5,.1), each = round(N/n_roles))[1:N]

# ppsbm ----
  
install.packages("ppsbm")
libray

l
