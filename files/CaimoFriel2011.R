
#### Simulate and estimate ERGM using a Bayesian approach (Caimo and Friel, 2011)
# toy example

library(ergm)
library(coda)
library(igraph)

N <- 150
theta1 <- -2.5 

x0 <- matrix(0, nrow = N, ncol = N) # initialize with an empty graph
iter <- 100000


t=1
for (t in 1:iter){
  # choose an edge
  i_loop <- sample(c(1:N), 1)
  j_loop <- sample(c(1:N), 1)
  while (i_loop == j_loop) { # to ensure that I don't pick an element of W where i=j since such an element is set to 0
    i_loop <- sample(c(1:N), 1)
    j_loop <- sample(c(1:N), 1)
  }
  
  # Metropolis-Hastings Algorithm
  if (t == 1){
    x_loop_old <- x0
    x_loop_candidate <- x_loop_old
    x_loop_candidate[i_loop,j_loop] <- ifelse(x_loop_old[i_loop,j_loop] == 0, 1, 0)
    x_loop_candidate[j_loop,i_loop] = x_loop_candidate[i_loop,j_loop] # since undirected network
    log_theta1_changestats <- 0.5*theta1*(sum(x_loop_candidate) - sum(x_loop_old)) # 0.5 since don't want to count edge twice
    log_thetavec_changestats <- log_theta1_changestats # only the number of edges in this case
    log.u_lambda <- log(runif(1))
    if (as.numeric(log.u_lambda) < as.numeric(log_thetavec_changestats)){
      x_loop_updated <- x_loop_candidate
    }
    else{
      x_loop_updated <- x_loop_old
    }
  }
  
  if (t > 1){
    x_loop_old <- x_loop_updated
    x_loop_candidate <- x_loop_old
    x_loop_candidate[i_loop,j_loop] <- ifelse(x_loop_old[i_loop,j_loop] == 0, 1, 0)
    x_loop_candidate[j_loop,i_loop] = x_loop_candidate[i_loop,j_loop] # since undirected network
    log_theta1_changestats <- 0.5*theta1*(sum(x_loop_candidate) - sum(x_loop_old)) # 0.5 since don't want to count edge twice
    log_thetavec_changestats <- log_theta1_changestats # only the number of edges in this case
    log.u_lambda <- log(runif(1))
    if (as.numeric(log.u_lambda) < as.numeric(log_thetavec_changestats)){
      x_loop_updated <- x_loop_candidate
    }
    else{
      x_loop_updated <- x_loop_old
    }
  }
}

x_observed_network <- x_loop_updated
sum(x_observed_network)
class(x_observed_network)
isSymmetric(x_observed_network)
networkgraph <- graph_from_adjacency_matrix(x_observed_network, mode = "upper")
plot(networkgraph)

rm(x_loop_updated)
rm(log_theta1_changestats)
rm(log_thetavec_changestats)



########################################################################################################
########################################################################################################
########################################################################################################
# recover theta1
# reference: pg. 149 of Exponential Random Graph Models for Social Networks: Theories, Methods and Applications


burn_in <- 0.2
ndraw <- 5000
aux_iter <- 100
theta_store <- matrix(NA, nrow = ndraw, ncol = 1)
accept_reject_store = matrix(NA, nrow = ndraw, ncol = 1)
theta_store[,1] = 0 # initialize with 0
accept_reject_store[,1] = 0

iteration = 2
for (iteration in 2:ndraw){
  theta_current <- theta_store[iteration-1,]
  theta_proposal <- rnorm(1, mean = theta_current, sd = 1)
  
  for (t in 1:aux_iter){
    i_loop <- sample(c(1:N), 1)
    j_loop <- sample(c(1:N), 1)
    while (i_loop == j_loop) { # to ensure that I don't pick an element of W where i=j since such an element is set to 0
      i_loop <- sample(c(1:N), 1)
      j_loop <- sample(c(1:N), 1)
    }
    
    # Metropolis-Hastings Algorithm
    if (t == 1){
      x_loop_old <- x_observed_network
      x_loop_candidate <- x_loop_old
      x_loop_candidate[i_loop,j_loop] <- ifelse(x_loop_old[i_loop,j_loop] == 0, 1, 0)
      x_loop_candidate[j_loop,i_loop] = x_loop_candidate[i_loop,j_loop] # since undirected network
      log_theta1_changestats <- 0.5*theta_proposal*(sum(x_loop_candidate) - sum(x_loop_old)) # 0.5 since don't want to count edge twice
      log_thetavec_changestats <- log_theta1_changestats # only the number of edges in this case
      log.u_lambda <- log(runif(1))
      if (as.numeric(log.u_lambda) < as.numeric(log_thetavec_changestats)){
        x_loop_proposed <- x_loop_candidate
      }
      else{
        x_loop_proposed <- x_loop_old
      }
    }
    
    if (t > 1){
      x_loop_old <- x_loop_proposed
      x_loop_candidate <- x_loop_old
      x_loop_candidate[i_loop,j_loop] <- ifelse(x_loop_old[i_loop,j_loop] == 0, 1, 0)
      x_loop_candidate[j_loop,i_loop] = x_loop_candidate[i_loop,j_loop] # since undirected network
      log_theta1_changestats <- 0.5*theta_proposal*(sum(x_loop_candidate) - sum(x_loop_old)) # 0.5 since don't want to count edge twice
      log_thetavec_changestats <- log_theta1_changestats # only the number of edges in this case
      log.u_lambda <- log(runif(1))
      if (as.numeric(log.u_lambda) < as.numeric(log_thetavec_changestats)){
        x_loop_proposed <- x_loop_candidate
      }
      else{
        x_loop_proposed <- x_loop_old
      }
    }
  }
  
  log_theta1_stats <- 0.5*(theta_current-theta_proposal)*(sum(x_loop_proposed) - sum(x_observed_network)) 
  log_thetavec_stats <- log_theta1_stats # only the number of edges in this case
  if (as.numeric(log.u_lambda) < as.numeric(log_thetavec_stats)){
    theta_updated <- theta_proposal
    accept_reject_store[iteration,] = 1
  }
  else{
    theta_updated <- theta_current
    accept_reject_store[iteration,] = 0
  }
  
  theta_store[iteration,] = theta_updated
  
}

mean(theta_store[c((burn_in*ndraw+1):ndraw)])
theta_store_mcmc <- as.mcmc(theta_store)
traceplot(theta_store_mcmc)
plot(density(theta_store[c((burn_in*ndraw+1):ndraw)]))
mean(accept_reject_store)
