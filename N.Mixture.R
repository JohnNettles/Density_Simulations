# load packages
library(jagsUI)
library(dplyr)
library(lubridate)
library(tidyr)



# Create timestamps and format data ---------------------------------------

#Create made up start time and add 1 minute for each simulated step
time.start <- as.POSIXct("2025-01-01 12:00:00")
time.end <- time.start + minutes(n_steps-1)

time.stamps <- seq(time.start, time.end, by= "1 min")


#Create a list of vectors, each with the detection history for a camera
vector_list <- list()
for (i in 1:n_cams){
  name <- paste("Camera",i,sep="_")
  vector_list[[name]] <- captures[i,]
}

#Format data into a data frame and then add a timestamp column
data <- as.data.frame(vector_list)
data$Timestamps <- c(time.stamps)

#Change to long format (so one row for every detection)
data_long <- data %>%
  pivot_longer(cols = -Timestamps, names_to = "Camera", values_to = "Detection") %>%
  filter(Detection == 1) %>%  # Keep only detections
  arrange(Camera, Timestamps)  # Ensure correct ordering



# Group and filter detections ---------------------------------------------

#Function to remove detections that are within 30 minutes of the previous detection
filter_detections <- function(df) {
  df <- arrange(df, Timestamps) #sort timestamps
  
  if (nrow(df) == 0) return(df[0,]) #what to do for detectors with no detections (empty case; return empty df)
  if (nrow(df) == 1) return(df) #If only one detection, keep it
  
  keep <- rep(FALSE, nrow(df)) #Initialize a vector of F so all are the same length
    keep[1] <- TRUE #always keep the first detection
  
  for (i in 2:nrow(df)) { #(Length of times would be the number of minutes until the next detection)
    last_kept_index <- max(which(keep))
    if (difftime(df$Timestamps[i],df$Timestamps[last_kept_index], units='mins') >= 30) {
      keep[i] <- TRUE
    }
    #max(which(keep)) finds the most recent (maximum) index of keep that was true
    #want to subtract the timestamp of the previous detection from that of the current detection
    #If that distance is greater than 30 minutes, then we include it
  }
  return(df[keep,])
}

#Apply function to each camera and count the number of detection events at each camera during a 5-day window
filtered_counts <- data_long %>%
  group_by(Camera) %>%
  group_modify(~ filter_detections(.x)) %>%
  ungroup()  %>% # Convert list output from `group_map()` back into a data frame
  mutate(Time_Bin = floor_date(Timestamps, unit = '1 day')) %>%
  group_by(Camera, Time_Bin) %>%
  summarise(Detection_Count = n(), .groups='drop') #Count detections per bin


#Convert the table to wide format and back to a data frame
counts <- data.frame(
  filtered_counts %>%
  pivot_wider(
    names_from = Time_Bin,
    values_from = Detection_Count,
    values_fill = list(Detection_Count = 0) #Fill missing counts with 0
  )
)


#Create vector of names for time bins
names <- rep(NA,ncol(counts))
names[1] <- "Camera"

for (i in 2:ncol(counts)){
  name <- paste("TimeBin",i-1,sep="_")
  names[i] <- name
}

#Edit column/row names
colnames(counts) <- c(names)
rownames(counts) <- c(counts[,1])
counts <- counts[,-1]

#Convert to matrix
detection_matrix <- as.matrix(counts)
dm <- c(detection_matrix)






# Bundle data to send to JAGS ---------------------------------------------

R = dim(detection_matrix)[1] #number of sites
n = dim(detection_matrix)[1] * dim(detection_matrix)[2] #number of sites times number of samples
n_bins = dim(detection_matrix)[2]

# also need a site index and covariates that have the same length 
# as the variable C
site = 1:R                # ‘Short’ version of site covariate; list of sites
site.p <- rep(site, n_bins)    # ‘Long’ version of site covariate; repeat for each rep
cbind(dm, site.p)  # Check that all went right

jags.dat <- list(R = R, n = n, C = dm, site.p = site.p)



# Poisson process ---------------------------------------------------------

##### SPECIFY MODEL CODE #####
cat(file="model.txt", "
model {

    # Priors
    alpha.lam ~ dunif(-10,10)
    alpha.p ~ dunif(-10,10)
    
    # Likelihood
    # Biological model for true abundance
    # Loop over R sites
    for (i in 1:R) {                        
      N[i] ~ dpois(lambda[i])
      Dens[i] <- N[i]/(2.62*2.62*3.14159)
      log(lambda[i]) <- alpha.lam
    }
    
    # Observation model for replicated counts
    # Loop over all n observations; i is 1:112 here
    for (i in 1:n){
      C[i] ~ dbin(p[i],N[site.p[i]]) #site.p is just 1 to R for each rep (length of 112)
      lp[i] <- alpha.p  #logit(p) follows linear relationship
      p[i] <- exp(lp[i])/(1+exp(lp[i])) #logit link
    }
    
    # Derived quantities
    # Estimate total population size across all sites
    totalN <- sum(N[]) #sum of all site-specific abundances
    AvgDens <- mean(Dens[])
    
    # Assess model fit using Chisquare discrepancy
    for (i in 1:n) {
    
      #Compute fit statistics for observed data
      eval[i] <- p[i]*N[site.p[i]]
      E[i] <- pow((C[i]-eval[i]),2)/(eval[i] + 0.5) #Error
      
      #Generate replicate data and compute fit statistics
      C.new[i] ~ dbin(p[i],N[site.p[i]])
      E.new[i] <- pow((C.new[i]-eval[i]),2)/(eval[i] + 0.5)
    }
    
    fit <- sum(E[])
    fit.new <- sum(E.new[])
    
    
}
 ")

##### INITIAL VALUES #####
# Inits function
# use the maximum count at each site as a first guess of what N might be 
# add 1 to avoid zeros
Nst <- apply(detection_matrix, 1, max) + 1 #max count of observed for each site plus 1 (so no zeros to avoid JAGS issues)
inits <- function(){list(N = Nst, alpha.lam=rnorm(1, 0, 1), 
                         alpha.p=rnorm(1, 0, 1))}

##### PARAMETERS TO MONITOR #####
params <- c("N", "totalN", "Dens", "AvgDens", "alpha.lam", "alpha.p",
            "fit","fit.new")

##### MCMC DIMENSIONS #####
nc = 3
nb =20000
ni = 50000
nt =50
n.iter = ni + nb

##### RUN THE MODEL IN JAGS #####
starttime = Sys.time()
jags.fit =  jags(data = jags.dat, inits = inits, parameters.to.save = params, model.file = "model.txt", 
                 n.thin = nt, n.chains = nc, n.burnin = nb, n.iter = n.iter, parallel = F)
Sys.time() - starttime



# Negative Binomial -------------------------------------------------------

##### SPECIFY MODEL CODE #####
cat(file="model.text","
model {

    #Priors
    alpha.lam ~ dunif(-10,10)

    alpha.p ~ dunif(-10,10)

    #Likelihood
    # Biological model for true abundance
    # Loop over R sites
        for (i in 1:R) {
          N[i] ~ dnegbin(p.nb[i],r)
          Dens[i] <- N[i]/(2.62*2.62*3.14159)
          p.nb[i] <- r/(r+lambda[i])
          log(lambda[i]) <- alpha.lam
        }

    # Observation model for replicated counts
    # Loop over all n observations; i is 1:600 here
        for (j in 1:n){
          C[j] ~ dbin(p[j],N[site.p[j]])
          lp[j] <- alpha.p  #logit(p) follows linear relationship
          p[j] <- exp(lp[j])/(1+exp(lp[j])) #logit link
     }
        
    # Derived quantities
    # Estimate total population size across all sites
        totalN <- sum(N[]) #sum of all site-specific abundances
        AvgDens <- mean(Dens[])
            
    # Assess model fit using Chisquare discrepancy
    for (i in 1:n) {
    
      #Compute fit statistics for observed data
      eval[i] <- p[i]*N[site.p[i]]
      E[i] <- pow((C[i]-eval[i]),2)/(eval[i] + 0.5) #Error
      
      #Generate replicate data and compute fit statistics
      C.new[i] ~ dbin(p[i],N[site.p[i]])
      E.new[i] <- pow((C.new[i]-eval[i]),2)/(eval[i] + 0.5)
    }
    
    fit <- sum(E[])
    fit.new <- sum(E.new[])
  }      
    
    ")


##### INITIAL VALUES #####
# Inits function
# use the maximum count at each site as a first guess of what N might be 
# add 1 to avoid zeros
Nst <- apply(detection_matrix, 1, max) + 1 #max count of observed for each site plus 1 (so no zeros to avoid JAGS issues)
inits <- function(){list(N = Nst, alpha.lam=rnorm(1, 0, 1), 
                         alpha.p=rnorm(1, 0, 1))}

##### PARAMETERS TO MONITOR #####
params <- c("N","Dens", "totalN","AvgDens", "alpha.lam", "alpha.p","fit","fit.new")

##### MCMC DIMENSIONS #####
nc = 3
nb =3000
ni = 6000
nt =50
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
