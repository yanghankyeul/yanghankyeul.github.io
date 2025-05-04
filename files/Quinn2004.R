# Bayesian Mixed Factor Analysis
# Implementation of Quinn (2004)

rm(list=ls())


library(dplyr)
library(MASS)
library(tmvtnorm)
library(MCMCpack)
library(truncnorm)
library(rstudioapi)

setwd(dirname(getActiveDocumentContext()$path))

data(PErisk)

# X* is of dimension N by J
# Phi is of dimension N by K
# Lambda is of dimension J by K

# assuming that K = 2,
# if both columns of Lambda are free for certain row j, Phi_jemptycircle is the entire Phi

constrained_variables_list <- list()
constrained_variables_list[[1]] <- c("courts", "-")



N <- nrow(PErisk) # no. of countries
J <- 5 # 2 observed continuous indicators and 3 ordinal variables
K <- 2 # 1 dimension


burn_in <- 0.2
iter <- 1000

X <- PErisk[,c(2:6)]
class(X[,2])



### Rescale and level ordered factors as to start from 0

X_scaled <- X
X_scaled[,2] <- scale(X_scaled[,2])
X_scaled[,5] <- scale(X_scaled[,5])
X <- X_scaled
X[,2] <- as.numeric(X[,2])
X[,5] <- as.numeric(X[,5])

for (j in 1:J){
  if("factor" %in% class(X[,j]) & min(X[,j]) != 1){
    X[,j] = factor(X[,j], levels = c(sort(as.numeric(as.character(unique(X[,j]))))), labels = c(1:length(unique(X[,j]))), ordered = TRUE)
  }
}




I_Kminus1 <- diag((K-1))
a_0j <- 0.001
b_0j <- 0.001
Lambda_prior_mean <- 0
Lambda_prior_precision <- 0.25





#### Create storage matrices and arrays
X_star_store <- array(data = NA, dim = c(N, J, iter))
Phi_store <- array(data = NA, dim = c(N, K, iter))
Lambda_store <- array(data = NA, dim = c(J, K, iter))
Psi_store <- matrix(NA, nrow = iter, ncol = J) # need to store only Psi_jj

for (j in 1:J){
  if ("factor" %in% class(X[,j])){
    assign(paste0("gamma_j",j,"_store"), matrix(NA, nrow = iter, ncol = (length(unique(X[,j]))+1)))
    assign(paste0("gamma_j",j,"_candidate_acceptreject_store"), matrix(NA, nrow = iter, ncol = 1))
  }
}




###########################################################################################################
###########################################################################################################



### Assign constraints, free rows and columns
lambda2_1_constrained01 <- 1
lambda2_2_constrained01 <- 0
lambda2_3_constrained01 <- 0
lambda2_4_constrained01 <- 0
lambda2_5_constrained01 <- 0
lambda2_1_positiveconstrained01 <- 0



# determine free and fixed indices for each row of J by K Lambda
for (j in 1:J){
  if ("numeric" %in% class(X[,j])){
    assign(paste0("free_indices_row",j), c(2))
  }
  
  if ("factor" %in% class(X[,j])) {
    assign(paste0("free_indices_row",j), c(1,2))
  }
  
  if ("numeric" %in% class(X[,j])){
    assign(paste0("fixed_indices_row",j), c(1))
  }
  
  if ("factor" %in% class(X[,j])) {
    assign(paste0("fixed_indices_row",j), NA)
  }
}

free_indices_row1
free_indices_row2
free_indices_row3
free_indices_row4
free_indices_row5

fixed_indices_row1
fixed_indices_row2
fixed_indices_row3
fixed_indices_row4
fixed_indices_row5


###########################################################################################################
###########################################################################################################


### Initialize
# initialize x_star
X_star_store[,,1] = 0


# initialize Lambda
Lambda_store[,,1] = 0

Lambda_store[1,2,1] = -1 # the second column initialized to -1 for negatively constrained (in this case, courts)
for (j in 1:J){
  if ("numeric" %in% class(X[,j])){
    Lambda_store[j,1,] = 0 # the first column constrained to 0 for continuous variables
  }
}

# initialize and constrain Phi
Phi_store[,,1] = 0
Phi_store[,1,] = 1 # the first column constraine to 1


for (j in 1:J){
  if ("factor" %in% class(X[,j])) {
    Psi_store[,j] = 1 # constrain Psi_jj to 1 for ordinal variables
  }
  Psi_store[1,j] = ifelse(is.na(Psi_store[1,j]), rchisq(1, df = 1), Psi_store[1,j]) # initialize Psi_jj to a draw from chi-squared for continuous variables
}


for (j in 1:J){
  if ("factor" %in% class(X[,j])){
    matrix_loop <- get(paste0("gamma_j",j,"_store"))
    matrix_loop[1,] = sort(rchisq(ncol(matrix_loop), df = 1))
    matrix_loop[,1] = -Inf
    matrix_loop[,2] = 0
    matrix_loop[,ncol(matrix_loop)] = Inf
    assign(paste0("gamma_j",j,"_store"), matrix_loop)
  }
}

for (j in 1:J){
  if ("factor" %in% class(X[,j])){
    gamma_j_candidate_acceptreject_store_loop <- get(paste0("gamma_j",j,"_candidate_acceptreject_store"))
    gamma_j_candidate_acceptreject_store_loop[1,] = 0
    assign(paste0("gamma_j",j,"_candidate_acceptreject_store"), gamma_j_candidate_acceptreject_store_loop)
  }
}


for (j in 1:J){
  if ("factor" %in% class(X[,j])){
    assign(paste0("t2j_",j), 0.25/length(unique(X[,j])))
  }
}

t2j_3 <- 0.25/length(unique(X[,3]))
t2j_4 <- 0.25/length(unique(X[,4]))


i=1
j=1
t=2



for (t in 2:iter){
  # update x_star_ij
  for (i in 1:N){
    for (j in 1:J){
      if ("numeric" %in% class(X[,j])){
        X_star_store[i,j,t] =  X[i,j]
      }
      else{
        mean_loop <- t(matrix(Lambda_store[j,,t-1])) %*% matrix(Phi_store[i,,t-1])
        X_star_store[i,j,t] = rtruncnorm(1, a = get(paste0("gamma_j",j,"_store"))[t-1,(as.numeric(as.character(X[i,j])))], b = get(paste0("gamma_j",j,"_store"))[t-1,(as.numeric(as.character(X[i,j])))+1], mean = mean_loop, sd = 1)
        # draw from Normal truncated to the interval(gamma_j(x_{ij}-1), gamma_j(x_{ij})]
        # for the first observation (i = 1) for j = 3, the category is 2 (=x_{ij}).
        # if we look at gamma_j3_store which consists of seven columns because of six categories,
        # then the first column is for c_0 = -Inf,
        # the second column is for c_1 = 0 and so forth.
        # meaning the interval(gamma_j(x_{ij}-1), gamma_j(x_{ij})] should be coded as "as.numeric(as.character(X[i,j]))"th column and "as.numeric(as.character(X[i,j]+1))"
        # this is to take into account of the 0th column at the very left for the gamma_j_store with -Inf
      }
    }
  }
  
  # update Phi_i - just need to do the second column
  for (i in 1:N){
    x_star_iloop <- matrix(X_star_store[i,,t]) # x*_i
    Psi_jj_matrix <- diag(Psi_store[t-1,])
    Lambda_2toK_matrix <- matrix(Lambda_store[,c(2:K),t-1])
    Lambda_column1_vector <- matrix(Lambda_store[,1,t-1])
    Phi_i_variance_loop <- solve(I_Kminus1 + t(Lambda_2toK_matrix) %*% solve(Psi_jj_matrix) %*% Lambda_2toK_matrix)
    Phi_i_mean_loop <- Phi_i_variance_loop %*% (t(Lambda_2toK_matrix) %*% solve(Psi_jj_matrix) %*% (x_star_iloop - Lambda_column1_vector))
    Phi_store[i,2,t] = rnorm(1, mean = Phi_i_mean_loop, sd = sqrt(Phi_i_variance_loop))
  }
  
  # update Lambda
  for (j in 1:J){
    if (is.na(get(paste0("fixed_indices_row",j)))){
      Phi_fixed_jloop <- matrix(0, nrow = N)
      lambda_fixed_jloop <- matrix(0, nrow = length(get(paste0("fixed_indices_row",j))), ncol = 1)
    }
    if (!is.na(get(paste0("fixed_indices_row",j)))) {
      Phi_fixed_jloop <- matrix(Phi_store[c(1:N),get(paste0("fixed_indices_row",j)),t])
      lambda_fixed_jloop <- matrix(Lambda_store[j,get(paste0("fixed_indices_row",j)),t-1], ncol = 1) # dimension Q by 1
    }
    Psi_jj_matrix <- diag(Psi_store[t-1,])
    Phi_free_jloop <- matrix(Phi_store[c(1:N),get(paste0("free_indices_row",j)),t], ncol = length(get(paste0("free_indices_row",j)))) # N by F: at least one column of Phi is free
    x_star_jloop <- matrix(X_star_store[c(1:N),j,t])
    Lambda_j_variance_loop <- solve(diag(length(get(paste0("free_indices_row",j))))*Lambda_prior_precision + solve(Psi_jj_matrix)[j,j]*(t(Phi_free_jloop)%*%Phi_free_jloop))
    Lambda_j_mean_loop <- Lambda_j_variance_loop %*% ((diag(length(get(paste0("free_indices_row",j))))*Lambda_prior_precision)%*%matrix(0, nrow = length(get(paste0("free_indices_row",j)))) + solve(Psi_jj_matrix)[j,j]*t(Phi_free_jloop)%*%(x_star_jloop - Phi_fixed_jloop%*%lambda_fixed_jloop))
    if ((get(paste0("lambda2_",j,"_constrained01")) == 1)) {
      if (get(paste0("lambda2_",j,"_positiveconstrained01")) == 1){
        Lambda_store[j,get(paste0("free_indices_row",j)),t] = rtmvnorm(1, mean = c(Lambda_j_mean_loop), sigma = Lambda_j_variance_loop, lower = c(-Inf, 0), upper = c(Inf, Inf), algorithm = "gibbs")
      }
      if (get(paste0("lambda2_",j,"_positiveconstrained01")) == 0){
        Lambda_store[j,get(paste0("free_indices_row",j)),t] = rtmvnorm(1, mean = c(Lambda_j_mean_loop), sigma = Lambda_j_variance_loop, lower = c(-Inf, -Inf), upper = c(Inf,0), algorithm = "gibbs")
      }
    }
    if ((get(paste0("lambda2_",j,"_constrained01")) == 0)) {
      Lambda_store[j,get(paste0("free_indices_row",j)),t] = mvrnorm(1, mu = Lambda_j_mean_loop, Sigma = Lambda_j_variance_loop)
    }
  }

  
  # update Psi_jj
  for (j in 1:J){
    if ("numeric" %in% class(X[,j])){
      x_star_jloop <- matrix(X_star_store[c(1:N),j,t])
      Phi_loop <- matrix(Phi_store[c(1:N),,t], ncol = 2)
      lambda_jloop <- matrix(Lambda_store[j,,t], ncol = 1)
      prod <- t(x_star_jloop - Phi_loop %*% lambda_jloop) %*% (x_star_jloop - Phi_loop %*% lambda_jloop)
      Psi_store[t,j] = rinvgamma(1, shape = (a_0j + N)/2, scale = (b_0j + prod)/2)
    }
  }
  
  
  # update gamma_j
  for (j in 1:J){
    if ("factor" %in% class(X[,j]) & length(unique(X[,j])) > 2){
      cutpoint_numerator_list <- list()
      cutpoint_denominator_list <- list()
      obs_numerator_list <- list()
      obs_denominator_list <- list()
      
      gamma_j_cutpoint_candidate_length <- length(unique(X[,j])) - 2 # no. of cutpoints to be drawn
      # There are six categories for j = 4. In this case, need the cutpoints for 2,3,4,5 and thus, gamma_j_cutpoint_candidate_length is 4
      C_j <- as.numeric(as.character(max(unique(X[,j])))) # no. of categories
      #t2j <- 1.5
      
      gamma_j_cutpoint_candidate_store_temp <- matrix(NA, nrow = 1, ncol = (gamma_j_cutpoint_candidate_length+3))
      gamma_j_cutpoint_candidate_store_temp[1,1] = -Inf # gamma_j0 = -Inf
      gamma_j_cutpoint_candidate_store_temp[1,2] = 0 # gamma_j1 = 0
      gamma_j_cutpoint_candidate_store_temp[1,ncol(gamma_j_cutpoint_candidate_store_temp)] = Inf # gamma_jC = Inf
      
      
      # Summary: wrt indexing columns, index by adding 1 when using gamma_j_store and the candidate_store
      for (c in 2:(C_j-1)){ # coded as c starting from c = 2
        if (c == 2){ # need to check the upper bound b
          gamma_j_cutpoint_candidate_store_temp[1,c+1] = rtruncnorm(1, a = 0, b = get(paste0("gamma_j",j,"_store"))[t-1,c+2], mean = get(paste0("gamma_j",j,"_store"))[t-1,c+1], sd = sqrt(get(paste0("t2j_",j))))
          # draw from Normal truncated to (0, \gamma_{j(c+1)}) with mean \gamma_{jc}
          # because for c = 2, \gamma^(can)_{j(c-1)} = \gamma_{j1} = 0 by definition
          # coded as a = 0 and b with (c+2)th column from the gamma_j_store
        }
        else {
          gamma_j_cutpoint_candidate_store_temp[1,c+1] = rtruncnorm(1, a = gamma_j_cutpoint_candidate_store_temp[1,c], b = get(paste0("gamma_j",j,"_store"))[t-1,c+2], mean = get(paste0("gamma_j",j,"_store"))[t-1,c+1], sd = sqrt(get(paste0("t2j_",j))))
          # draw from Normal truncated to (\gamma^(can)_{j(c-1)}, \gamma_{j(c+1)}) with mean \gamma_{jc}
          # eg/ if c = 3, then lower bound is the \gamma^(can)_{j(2)} which is the previous candidate draw for c = 2 which is the second column of gamma_j_cutpoint_candidate_store_temp
          # while the upper bound and the mean follow the same pattern as above
        }
      }
      
      # calculating cutpoint numerator and denominator
      for (c in 2:(C_j-1)){
        if (c == 2){
          cutpoint_numerator_loop <- pnorm((get(paste0("gamma_j",j,"_store"))[t-1,c+2] - get(paste0("gamma_j",j,"_store"))[t-1,c+1])/sqrt(get(paste0("t2j_",j)))) - pnorm((0 - get(paste0("gamma_j",j,"_store"))[t-1,c+1])/sqrt(get(paste0("t2j_",j))))
          # Phi((\gamma_{c+1} - \gamma_{c})/t_j) - Phi((\gamma^{can}_{c-1} - \gamma_{c})/t_j)
          # = Phi((\gamma_{c+1} - \gamma_{c})/t_j) - Phi((0 - \gamma_{c})/t_j) for c = 2
          cutpoint_numerator_list[[(c)]] <- cutpoint_numerator_loop
          
          cutpoint_denominator_loop <- pnorm((gamma_j_cutpoint_candidate_store_temp[1,c+2] - gamma_j_cutpoint_candidate_store_temp[1,c+1])/sqrt(get(paste0("t2j_",j)))) - pnorm((0 - gamma_j_cutpoint_candidate_store_temp[1,c+1])/sqrt(get(paste0("t2j_",j))))
          # Phi((\gamma^{can}_{c+1} - \gamma^{can}_{c})/t_j) - Phi((\gamma_{c-1} - \gamma^{can}_{c})/t_j)
          cutpoint_denominator_list[[c]] <- cutpoint_denominator_loop
        }
        
        else{
          cutpoint_numerator_loop <- pnorm((get(paste0("gamma_j",j,"_store"))[t-1,c+2] - get(paste0("gamma_j",j,"_store"))[t-1,c+1])/sqrt(get(paste0("t2j_",j)))) - pnorm((gamma_j_cutpoint_candidate_store_temp[1,c] - get(paste0("gamma_j",j,"_store"))[t-1,c+1])/sqrt(get(paste0("t2j_",j))))
          # Phi((\gamma_{c+1} - \gamma_{c})/t_j) - Phi((\gamma^{can}_{c-1} - \gamma_{c})/t_j)
          cutpoint_numerator_list[[c]] <- cutpoint_numerator_loop
          
          cutpoint_denominator_loop <- pnorm((gamma_j_cutpoint_candidate_store_temp[1,c+2] - gamma_j_cutpoint_candidate_store_temp[1,c+1])/sqrt(get(paste0("t2j_",j)))) - pnorm((get(paste0("gamma_j",j,"_store"))[t-1,c] - gamma_j_cutpoint_candidate_store_temp[1,c+1])/sqrt(get(paste0("t2j_",j))))
          # Phi((\gamma^{can}_{c+1} - \gamma^{can}_{c})/t_j) - Phi((\gamma_{c-1} - \gamma^{can}_{c})/t_j)
          cutpoint_denominator_list[[c]] <- cutpoint_denominator_loop
        }
      }
      
      # calculating observation numerator and denominator
      for (i in 1:N){
        mean_loop <- t(matrix(Lambda_store[j,,t])) %*% matrix(Phi_store[i,,t])
        obs_numerator_loop <- pnorm(gamma_j_cutpoint_candidate_store_temp[1,(as.numeric(as.character(X[i,j])))+1] - mean_loop) - pnorm(gamma_j_cutpoint_candidate_store_temp[1,(as.numeric(as.character(X[i,j])))] - mean_loop)
        obs_numerator_list[[i]] <- obs_numerator_loop
        
        obs_denominator_loop <- pnorm(get(paste0("gamma_j",j,"_store"))[t-1,(as.numeric(as.character(X[i,j])))+1] - mean_loop) - pnorm(get(paste0("gamma_j",j,"_store"))[t-1,(as.numeric(as.character(X[i,j])))] - mean_loop)
        obs_denominator_list[[i]] <- obs_denominator_loop
      }
      
      alpha <- (prod(unlist(obs_numerator_list))/prod(unlist(obs_denominator_list))) * (prod(unlist(cutpoint_numerator_list))/prod(unlist(cutpoint_denominator_list)))
      unif_draw <- runif(1)
      

      if(unif_draw < alpha){
        gamma_jstore_loop <- get(paste0("gamma_j",j,"_store"))
        gamma_jstore_loop[t,] = gamma_j_cutpoint_candidate_store_temp
        assign(paste0("gamma_j",j,"_store"), gamma_jstore_loop)
        
        gamma_j_candidate_acceptreject_store_loop <- get(paste0("gamma_j",j,"_candidate_acceptreject_store"))
        gamma_j_candidate_acceptreject_store_loop[t,1] = 1
        assign(paste0("gamma_j",j,"_candidate_acceptreject_store"), gamma_j_candidate_acceptreject_store_loop)
      }
      else{
        gamma_jstore_loop <- get(paste0("gamma_j",j,"_store"))
        gamma_jstore_loop[t,] = gamma_jstore_loop[t-1,]
        assign(paste0("gamma_j",j,"_store"), gamma_jstore_loop)
        
        gamma_j_candidate_acceptreject_store_loop <- get(paste0("gamma_j",j,"_candidate_acceptreject_store"))
        gamma_j_candidate_acceptreject_store_loop[t,1] = 0
        assign(paste0("gamma_j",j,"_candidate_acceptreject_store"), gamma_j_candidate_acceptreject_store_loop)
      }
    }
  }
  
  
  cat(t, paste0("of ", iter,"\r")) 
  flush.console()
  
}


mean(Lambda_store[1,1,c((burn_in*iter):iter)]) # -0.041
mean(Lambda_store[2,1,c((burn_in*iter):iter)]) # 0
mean(Lambda_store[3,1,c((burn_in*iter):iter)]) # 3.517
mean(Lambda_store[4,1,c((burn_in*iter):iter)]) # 3.146
mean(Lambda_store[5,1,c((burn_in*iter):iter)]) # 0

mean(Lambda_store[1,2,c((burn_in*iter):iter)]) # -2.930
mean(Lambda_store[2,2,c((burn_in*iter):iter)]) # 0.750
mean(Lambda_store[3,2,c((burn_in*iter):iter)]) # -1.963
mean(Lambda_store[4,2,c((burn_in*iter):iter)]) # -2.278
mean(Lambda_store[5,2,c((burn_in*iter):iter)]) # -0.721


mean(Psi_store[,1])
mean(Psi_store[,2])
mean(Psi_store[,3])
mean(Psi_store[,4])
mean(Psi_store[,5])

plot(density(Phi_store[6,2,])) # Bolivia
Phi_store[N,2,] # Zimbabwe


#mean(gamma_j1_candidate_acceptreject_store[,1])
mean(gamma_j3_candidate_acceptreject_store[,1])
mean(gamma_j4_candidate_acceptreject_store[,1])

head(gamma_j3_store)
head(gamma_j4_store)
tail(gamma_j3_store)
tail(gamma_j4_store)



# library(MCMCpack)
# data(PErisk)
# post.samp <- MCMCmixfactanal(~ courts+barb2+prsexp2+prscorr2+gdpw2,
#                              factors=1, data=PErisk,
#                              lambda.constraints = list(courts=list(2, "-")),
#                              burin=100, mcmc=2000, thin=1,
#                              verbose=TRUE, L0=0.25,
#                              store.lambda=TRUE, store.scores=TRUE, tune=0.25)
# summary(post.samp)
# View(post.samp)
