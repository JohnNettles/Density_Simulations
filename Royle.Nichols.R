detection_binary <- (detection_matrix > 0) *1


# Bundle data to give to JAGS ---------------------------------------------

n_sites = dim(detection_binary)[1] #number of sites
n = dim(detection_binary)[1] * dim(detection_binary)[2] #number of sites times number of samples
n_bins = dim(detection_binary)[2]


jags.dat <- list(n_sites = n_sites, n_bins = n_bins, D = detection_binary)



# Poisson process ---------------------------------------------------------

cat(file = "model.txt", "
model {

    #Priors
    lambda ~ dunif(0,10)
    p ~ dunif(0,1)
  
    #Likelihood
    for (i in 1:n_sites) {
    #Latent abundance (Poisson)
    N[i] ~ dpois(lambda)
    Dens[i] <- N[i]/(2.62*2.62*3.14159)
  
      #Detection process (binomial)
      for (j in 1:n_bins) {
      D[i,j] ~ dbin(p, N[i])
      }
    }

    #Derived quantities
    totalN <- sum(N[])
    AvgDens <- mean(Dens[])


    # Assess model fit using Chisquare discrepancy
      for (i in 1:n_sites) {
        for (j in 1:n_bins) {
        
        #Compute fit statistics for observed data
        eval[i,j] <- p*N[i]
        E[i,j] <- pow((D[i,j]-eval[i,j]),2)/(eval[i,j] + 0.5) #Error
        
        #Generate replicate data and compute fit statistics
        D.new[i,j] ~ dbin(p,N[i])
        E.new[i,j] <- pow((D.new[i,j]-eval[i,j]),2)/(eval[i,j] + 0.5)
        }
      }
      
      fit <- sum(E[,])
      fit.new <- sum(E.new[,])
      
}    
    ")


##### INITIAL VALUES #####
# Inits function
# use the maximum count at each site as a first guess of what N might be 
# add 1 to avoid zeros
Nst <- apply(detection_binary, 1, max) + 1 #max count of observed for each site plus 1 (so no zeros to avoid JAGS issues)
inits <- function(){list(N = Nst, lambda=runif(1,1,5), 
                         p=runif(1, 0.1, 0.9))}

##### PARAMETERS TO MONITOR #####
params <- c("N", "totalN", "Dens", "AvgDens", "lambda", "p",
            "fit","fit.new")

##### MCMC DIMENSIONS #####
nc = 3
nb =3000
ni = 10000
nt =10
n.iter = ni + nb

##### RUN THE MODEL IN JAGS #####
starttime = Sys.time()
jags.fit =  jags(data = jags.dat, inits = inits, parameters.to.save = params, model.file = "model.txt", 
                 n.thin = nt, n.chains = nc, n.burnin = nb, n.iter = n.iter, parallel = F)
Sys.time() - starttime



# Negative Binomial -------------------------------------------------------

cat(file = "model.txt", "
model {

    #Priors
    lambda ~ dunif(0,10)
    p ~ dunif(0,1)
    r ~ dunif(0,50) #overdispersion parameter
  
    theta.RN <- r/(r+lambda)
  
    #Likelihood
    for (i in 1:n_sites) {
      #Latent abundance (Negative Binomial)
      N[i] ~ dnegbin(theta.RN,r)
      Dens[i] <- N[i]/(2.62*2.62*3.14159)
  
      #Detection process (binomial)
      for (j in 1:n_bins) {
      D[i,j] ~ dbin(p, N[i])
      }
    }

    #Derived quantities
    totalN <- sum(N[])
    AvgDens <- mean(Dens[])


    # Assess model fit using Chisquare discrepancy
      for (i in 1:n_sites) {
        for (j in 1:n_bins) {
        
        #Compute fit statistics for observed data
        eval[i,j] <- p*N[i]
        E[i,j] <- pow((D[i,j]-eval[i,j]),2)/(eval[i,j] + 0.5) #Error
        
        #Generate replicate data and compute fit statistics
        D.new[i,j] ~ dbin(p,N[i])
        E.new[i,j] <- pow((D.new[i,j]-eval[i,j]),2)/(eval[i,j] + 0.5)
        }
      }
      
      fit <- sum(E[,])
      fit.new <- sum(E.new[,])
      
}    
    ")


##### INITIAL VALUES #####
# Inits function
# use the maximum count at each site as a first guess of what N might be 
# add 1 to avoid zeros
Nst <- apply(detection_binary, 1, max) + 1 #max count of observed for each site plus 1 (so no zeros to avoid JAGS issues)
inits <- function(){list(N = Nst, lambda=runif(1,1,5), 
                         p=runif(1, 0.1, 0.9))}

##### PARAMETERS TO MONITOR #####
params <- c("N", "totalN", "Dens", "AvgDens", "lambda", "p",
            "fit","fit.new")

##### MCMC DIMENSIONS #####
nc = 3
nb =3000
ni = 10000
nt =10
n.iter = ni + nb

##### RUN THE MODEL IN JAGS #####
starttime = Sys.time()
jags.fit =  jags(data = jags.dat, inits = inits, parameters.to.save = params, model.file = "model.txt", 
                 n.thin = nt, n.chains = nc, n.burnin = nb, n.iter = n.iter, parallel = F)
Sys.time() - starttime


# Convergence Diagnostics -------------------------------------------------

# extract summary statistics
out = jags.fit$summary

hist(out[,"Rhat"])

# visualize trace and posterior plots
traceplot(jags.fit) 

# view table of JAGS model output - can be copied and pasted to excel
jags.View(jags.fit, digits=3)

# posterior predictive check
plot(jags.fit$sims.list$fit, jags.fit$sims.list$fit.new, main = "", xlab = "Discrepancy
measure for actual data set", ylab = "Discrepancy measure for perfect data sets")
abline(0,1, lwd = 2, col = "red")

# calculate Bayesian p-value
mean(jags.fit$sims.list$fit.new > jags.fit$sims.list$fit) #for each iteration, is fit.new greater than fit

# check betas
whiskerplot(jags.fit, parameters="N", quantiles=c(0.025,0.975), zeroline=TRUE)

