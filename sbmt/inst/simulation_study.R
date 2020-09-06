library(sbmt)

# simulate data ----

n_roles = 2
Time = 24 
par(mfrow = c(2,2))
N = 30

a = 10
b = 5
ymax =  max(2*a, 2*b)

#curves
par(mfrow = c(2,2))
curve(b*sin(x*pi/Time) + b, 0, Time, ylim = c(0, ymax))
curve(a*sin(x*2*pi/Time)+a, 0, Time, ylim = c(0, ymax))
curve(-a*sin(x*2*pi/Time)+a, 0, Time, ylim = c(0, ymax))
curve(b*sin(x*pi/Time)+ b, 0, Time, ylim = c(0, ymax))

x = seq(.5, Time)

omega_11 = b*sin(x*pi/Time) + b
omega_12 =  a*sin(x*2*pi/Time)+a
omega_21 = -a*sin(x*2*pi/Time)+a
omega_22 = omega_11

# show overlapping
curve(omega_12, 0, Time, ylim = c(0, ymax))
curve(omega_21, 0, Time, col = "red", add = TRUE)
curve(omega_11, 0, Time, col = "blue", add = TRUE)

omega = array(rbind(omega_11, omega_21, omega_12, omega_22), dim = c(n_roles, n_roles, Time)) #left most index moves fastest
              
roles_discrete = rep(1:n_roles, length.out = N)
set.seed(1)

discrete_edge_array = array(0, dim = c(N, N, Time))

for (i in 1:N) {
  role_i = roles_discrete[i]
  for (j in 1:N) {
    role_j = roles_discrete[j]
    for (time in 1:Time) {
      ijt = rpois(omega[role_i, role_j, time], n= 1)
      discrete_edge_array[i, j, time] = ijt
    }
  }
}
  
  # checks
  dim(discrete_edge_array)
  plot(apply(discrete_edge_array[roles_discrete==1,roles_discrete==1,],3,mean), type = "l") 
  plot(apply(discrete_edge_array[roles_discrete==1,roles_discrete==2,],3,mean), type = "l") 
  plot(apply(discrete_edge_array[roles_discrete==2,roles_discrete==1,],3,mean), type = "l") 
  plot(apply(discrete_edge_array[roles_discrete==2,roles_discrete==2,],3,mean), type = "l") 
  image(discrete_edge_array[(1:N)[order((1:N)%%n_roles)], (1:N)[order((1:N)%%n_roles)], 1])
  image(discrete_edge_array[(1:N)[order((1:N)%%n_roles)], (1:N)[order((1:N)%%n_roles)], 7])
  image(discrete_edge_array[(1:N)[order((1:N)%%n_roles)], (1:N)[order((1:N)%%n_roles)], 13])
  image(discrete_edge_array[(1:N)[order((1:N)%%n_roles)], (1:N)[order((1:N)%%n_roles)],19])
  
# sbmt ----

tmp = apply(discrete_edge_array, 3, function(x) {data.frame(cbind(which(is.finite(x), arr.ind = TRUE)), count = as.vector(x))})
  
# ppsbm ----
  
install.packages("ppsbm")
libray

l