
Dens <- (y/t) * (pi/(v*r*(2+theta)))

#r = radius of detection zone (detection distance)
#theta = angle of detection zone
#v = animal speed of movement
#y = # of detections 
#t = total sampled time


# Define variables --------------------------------------------------------

r <- r
theta <- theta
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

density <- rep(NA, n_cams)

for(i in 1:n_cams) {
  
  density[i] <- (y[i]/t) * (pi/(v*r*(2+theta)))
}

abundance <- mean(density) * 1000000 * (sa.width/1000)*(sa.length/1000)

# Bootstrap for variance

for (i in 1:10000) {
  
  bootstrap <- sample(density,length(density),replace=T)
  bs.sd[i] <- sd(bootstrap)
}




P <- A/a # Number of possible cameras (grid cells) in study area

boot.fn <- function(data, nboot, P) {
  estN <- rep(NA, nboot)
  for (i in 1:nboot){
    boot.dat <- list(count = matrix(sample(data$count, data$ncam*data$nocc, replace = T), 
                                    nrow = data$ncam),
                     ncam = data$ncam,
                     nocc = data$nocc) 
    estN[i] <- IS.fn(boot.dat, P)
  }
  return(estN)
}

estN.boot <- boot.fn(data.list, nboot = 100, P = P)

se_estN <- sd(estN.boot)
se_estN
CI_estN <- quantile(estN.boot, c(0.025, 0.975))
CI_estN

