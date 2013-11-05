library(fields)
library(raster)
library(rgdal)
library(pscl, verbose = FALSE)

setwd("/Volumes/Data Drive 1/Dropbox/PalEON Working Folder/Temp")
setwd("~/Dropbox/PalEON Working Folder/Temp")

## Load pre-processed fort data
load(file = "Temperature.RData")

source("mcmc.temp.fixed.sigma.squared.beta.R")

##
## Initialize MCMC parameters
##

num.pca <- 3 # number of principal components to use
mu.0 <- rep(0, num.pca)
Sigma.0 <- diag(num.pca)
mu.beta <- rep(0, num.pca)
sigma.squared.beta <- 0.188 / 100  # posterior mcmc mean is 0.188 from other MCMC model fits
n.mcmc <- 1000
alpha.epsilon <- 10
beta.epsilon <- 10
alpha.beta <- 10
beta.beta <- 10

##
## Fit model using MCMC
##

start <- Sys.time()
mcmc.out <- pallette.mcmc(y.celsius, X.celsius, mu.0, Sigma.0, sigma.squared.beta, n.mcmc, alpha.epsilon, beta.epsilon, alpha.beta, beta.beta, mu.beta, KT, loc.id.unique, na.rows)
finish <- Sys.time() - start ## 100 iterations about 1 hour
finish

## Put MCMC mean surface in raster format
for(s in 1:73){
  values(fort.raster[[s]])[ -na.rows] <- mcmc.out$fort.raster[,s]
}

## Plot MCMC mean surface
#library(animation)
#saveGIF(
for(s in 1:73){
  image(fort.raster[[s]], main = paste0("Reconstructed temperature surface for year ", s + 1819), col = topo.colors(100))
  image.plot(matrix(values(fort.raster[[s]]), nrow = 132, ncol = 240), add = TRUE, col = topo.colors(100))
#  points(year[[s]][, 1:2])
}#, "Climate.Reconstruction.full.mcmc.shrinkage.gif") 


## Posterior mean of beta regression parameters
apply(mcmc.out$beta.save, 1, mean)

##
## Plot of Posterior quantities
##

## Plot each Beta vector through time (temporal evolution of beta)
dim(mcmc.out$beta.save)
for(i in 1:num.pca){
  matplot(t(mcmc.out$beta.save[i, , ]) , type = 'l', main = paste0("trace plot of beta", i, " for each of 73 years"))
}

## Plot of each beta 
for(i in 1:73){
  matplot(t(mcmc.out$beta.save[, i, ]) , type = 'l', main = paste0("trace plot of all betas for year ", i + 1819), xlab = "mcmc iterations")
}

plot(mcmc.out$sigma.squared.epsilon.save, type = 'l', main = "trace plot of sigma.squared.epsilon")

matplot(t(mcmc.out$mu.beta.save), type = 'l', main = "trace plots for mu.beta")

