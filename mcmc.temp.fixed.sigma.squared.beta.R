##
## MCMC pallette algorithm
##
## John Tipton
##
## Created 10.20.2013

##
## Model: y_t = K_t %*% X %*% B_t + epsilon_t
##
##        B_t ~ N(mu_B, Sigma_B) 
##
##        mu_B ~ N(mu_0, Sigma_0)  
##
## X_t = f(Y_tau) where the Y's are the PRISM predictions for years tau where the PRISM exits and f() is a PCA transformation
##
## y_t are the fort observations for year t
##
## K_t is the selection matrix that ties observation locations for the fort data to the PRISM data
##
## Sigma_B is the matrix of eigenvalues or eigenvectors? for shrinkage. Something like: sigma^2_B * Lambda^(-1) for the eigenvector(value) matrix Lambda
##
## mu_0 and Sigma_0 are hyperparameters
##


##
## y is a list of t years of observational data
## 
##

pallette.mcmc <- function(y, X, mu.0, Sigma.0, sigma.squared.beta, n.mcmc, alpha.epsilon, beta.epsilon, alpha.beta, beta.beta, mu.beta, KT, loc.id, na.rows){

  ##
  ## Libraries and Subroutines
  ##

  make.epsilon.vector <- function(s, beta){
    y[, s][loc.id[[s]]] - KT[[s]] %*% X %*% beta[, s]
  }

#  make.sum.sigma.beta <- function(beta, mu.beta){
#    temp <- vector(length = t)
#    for(s in 1:t){
#    	  temp[s] <- t(beta[, s] - mu.beta) %*% Lambda.inv %*% (beta[, s] - mu.beta)
#    }
#    return(sum(temp))
#  }
  
  ##
  ## Initialize parameters
  ## 
  
  X.pca <- prcomp(X, center = FALSE, retx = TRUE)
  lambda <- X.pca$sdev^2
  X <- X.pca$x[, 1:num.pca]
  t <- dim(y)[2]
  if(is.null(dim(X)) == TRUE){
    ncells <- length(X)} else {ncells <- dim(X)[1]
  }

  ncells
  
  tau <- ifelse(is.null(dim(X)) == TRUE, 1, dim(X)[2])
  
  nt <- c()
  
  for(s in 1:t){
  	nt[s] <- length(which(is.na(y[, s]) == FALSE))
  }
  
  nt.sum <- sum(nt)
  sigma.squared.epsilon <- rigamma(1, alpha.epsilon, beta.epsilon)
  Sigma.epsilon <- list()
  Sigma.epsilon.inv <- list()
  
  for(s in 1:t){
    Sigma.epsilon[[s]] <- sigma.squared.epsilon * diag(nt[s])
    Sigma.epsilon.inv[[s]] <- solve(Sigma.epsilon[[s]])
  } 
  
  Lambda <- diag(lambda[1:num.pca], nrow  = length(lambda[1:num.pca])) ## eigenvalues??
  Lambda.determinant <- prod(X.pca$sdev[1:num.pca]^2)
  Lambda.inv <- solve(Lambda)
  
#  sigma.squared.beta <- rigamma(1, alpha.beta, beta.beta)
  Sigma.beta <- sigma.squared.beta * Lambda
  Sigma.beta.inv <- solve(Sigma.beta)
  Sigma.0.inv <- solve(Sigma.0)
  beta <- matrix(0, nrow = tau, ncol = t)
  n.burn <- floor(n.mcmc / 10)  
  
  ##
  ## Initialize Storage
  ##
  
  beta.save <- array(dim = c(tau, t, n.mcmc))
  sigma.squared.beta.save <- vector(length = n.mcmc)
  sigma.squared.epsilon.save <- vector(length = n.mcmc)
  mu.beta.save <- matrix(NA, nrow = tau, ncol = n.mcmc)
  fort.raster <- matrix(0, nrow = ncells, ncol = t)
  
  ##
  ## Begin MCMC loop
  ##
  for(k in 1:n.mcmc){
  	if(k %% 10 == 0) cat(" ", k)

  	##
  	## Sample Beta
  	##
  	
  	for(s in 1:t){
      devs <- rnorm(tau)
      beta.A.chol <- chol(t(X) %*% t(KT[[s]]) %*% Sigma.epsilon.inv[[s]] %*% KT[[s]] %*% X + Sigma.beta.inv)
      beta.b <- t(X) %*% t(KT[[s]]) %*% Sigma.epsilon.inv[[s]] %*% y[, s][loc.id[[s]]] + Sigma.beta.inv %*% mu.beta
      beta[, s] <- backsolve(beta.A.chol, backsolve(beta.A.chol, beta.b, transpose = TRUE) + devs)
  	}
  	
  	##
  	## Sample mu.beta
  	##
  	
  	devs <- rnorm(tau)
  	mu.beta.A.chol <- chol(t * Sigma.beta.inv + Sigma.0.inv)
  	mu.beta.b <- apply(Sigma.beta.inv %*% beta, 1, sum) + Sigma.0.inv %*% mu.0
  	mu.beta <- backsolve(mu.beta.A.chol, backsolve(mu.beta.A.chol, mu.beta.b, transpose = TRUE) + devs)
  
  ##
  ## Sample sigma.squared.beta <- this is treated as fixed
  ##
  
#  sigma.squared.beta <- rigamma(1, alpha.beta + t * Lambda.determinant / 2, 1 / 2 * make.sum.sigma.beta(beta, mu.beta) + beta.beta)  
#  Sigma.beta <- sigma.squared.beta * Lambda
#  Sigma.beta <- solve(Sigma.beta)
  
    ##
    ## Sample sigma.squared.epsilon
    ##
  
    epsilon.vector <- unlist(lapply(1:t, make.epsilon.vector, beta = beta))
    sigma.squared.epsilon <-rigamma(1, 1 / 2 * nt.sum + alpha.epsilon, 1 / 2 * t(epsilon.vector) %*% epsilon.vector + beta.epsilon)
    
    for(s in 1:t){
      Sigma.epsilon[[s]] <- sigma.squared.epsilon * diag(nt[s])
#      Sigma.epsilon.inv[[s]] <- solve(Sigma.epsilon[[s]])
      Sigma.epsilon.inv[[s]] <- 1 / sigma.squared.epsilon * diag(nt[s])
  }
  
  ##
  ## Simulate random field
  ##
  
  if(k > n.burn){
  	for(s in 1:t){
  	  if(dim(beta)[1] == 1){
  	    fort.raster[, s] <- fort.raster[, s] + X * beta[s] / (n.mcmc - n.burn)
  	  } else {
  	    fort.raster[, s] <- fort.raster[, s] + X %*% beta[, s] / (n.mcmc - n.burn)
  	  }
  	}
  }
  
  ##
  ## Save variables
  ## 
  
  beta.save[, , k] <- beta
  sigma.squared.beta.save[k] <- sigma.squared.beta
  sigma.squared.epsilon.save[k] <- sigma.squared.epsilon
  mu.beta.save[, k] <- mu.beta  	
  }
  
##
## Write output
##
  
list(beta.save = beta.save, sigma.squared.beta.save = sigma.squared.beta.save, sigma.squared.epsilon.save = sigma.squared.epsilon.save, mu.beta.save = mu.beta.save, n.mcmc = n.mcmc, fort.raster = fort.raster)
}


