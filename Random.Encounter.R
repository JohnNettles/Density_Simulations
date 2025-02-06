
Dens <- (y/t) * (pi/(v*r*(2+theta)))

#r = radius of detection zone (detection distance)
#theta = angle of detection zone
#v = animal speed of movement
#y = # of detections 
#t = total sampled time


# Define variables --------------------------------------------------------

r.EDD <- r.EDD
theta.EDA <- theta.EDA
v <- 0.5 #converted to meters per minute
t <- n_steps


# Format data -------------------------------------------------------------

#Create an empty results matrix
counts <- rep(NA,n_cams)

#need to sum the detections for each row (camera)
for (i in 1:n_cams) {
  counts[i] <- sum(captures[i,])
  
}

#create y vector
y <- counts


# Calculate density for each camera ---------------------------------------
#answer is animals per square meter
#multiply by area of study space (in square meters) to get abundance

REM <- function(y,t,v,r,theta,n_cams) {
  
  dens <- rep(NA,n_cams)
  
  for(i in 1:n_cams) {
    dens[i] <- (y[i]/t) * (pi/(v*r.EDD[i]*(2+theta.EDA[i])))
  }
  
  return(dens)
}

REM.res <- REM(y=y,t=t,v=v,r=r.EDD,theta=theta.EDA,n_cams=n_cams)

abundance <- mean(REM.res) * A 



# Bootstrap for variance --------------------------------------------------

boot.fn <- function(data,nboot) {
  counts <- rep(NA,n_cams)
  boot.rem <- rep(NA,nboot)
  
  for (i in 1:nboot) {
    #Resample detections
    boot.caps <- matrix(sample(data,n_cams*n_steps,replace=T),nrow=n_cams) 
    
    #need to sum the detections for each row (camera)
    counts <- rowSums(boot.caps)
    
    #Calculate the average density for each bootstrap iteration
    boot.rem[i] <- mean(REM(y = counts, t=t, v=v, r=r.EDD,theta=theta.EDA, n_cams=n_cams)) 
  }
  
  return(boot.rem)
}
  
#Run the function
estD.boot <- boot.fn(captures, nboot = 100)

#Convert to abundance
meanN_boot <- mean(estD.boot)*A
meanN_boot
se_estN <- sd(estD.boot)*A
se_estN
CI_estN <- quantile(estD.boot*A, c(0.025, 0.975))
CI_estN

