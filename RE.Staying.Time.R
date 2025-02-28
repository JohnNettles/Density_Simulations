dens <- (Y * st)/(a*t)

#y = # of detections
#st = staying time
#a = area of detection zone
#t = length of research period


# Define variables --------------------------------------------------------

r.EDD <- r.EDD
theta.EDA <- theta.EDA
t <- n_steps
a <- pi * r.EDD * theta.EDA/360 # Area of a each camera in square meters


# Compute staying time ----------------------------------------------------

count_ones <- function(row) {
  rle_result <- rle(row)
  ones_positions <- which(row == 1)
  
  counts <- rep(NA, length(row))
  
  start_index <- 1
  for (i in seq_along(rle_result$values)) {
    if (rle_result$values[i] == 1) {
      count <- rle_result$lengths[i] 
      counts[start_index] <- count
    }
    start_index <- start_index + rle_result$lengths[i]
  }
  return(counts)
}

counts_matrix <- t(apply(captures,1,count_ones))
st <- mean(counts_matrix, na.rm=T)

# Get counts -------------------------------------------------------------

#Create an empty results matrix
counts <- rep(NA,n_cams)

#need to sum the detections for each row (camera)
for (i in 1:n_cams) {
  counts[i] <- sum(captures[i,])
}

#create y vector
y <- counts


# Calculate density for each camera --------------------------------------
#answer is animals per square meter
#multiply by area of study space (in square meters) to get abundance

REST <- function(y,st,a,t,n_cams) {
  
  dens <- rep(NA,n_cams)
  
  for(i in 1:n_cams) {
    dens[i] <- (y[i] * st)/(a[i]*t)
  }
  
  return(dens)
}

REST.res <- REST(y=y,st=st,t=t,a=a,n_cams=n_cams)

abundance <- mean(REST.res) * A 


# Bootstrap for variance --------------------------------------------------

boot.fn <- function(data,nboot) {
  counts <- rep(NA,n_cams)
  boot.rest <- rep(NA,nboot)
  
  for (i in 1:nboot) {
    #Resample detections
    boot.caps <- matrix(sample(data,n_cams*n_steps,replace=T),nrow=n_cams) 
    
    #need to sum the detections for each row (camera)
    counts <- rowSums(boot.caps)
    
    #Calculate the average density for each bootstrap iteration
    boot.rest[i] <- mean(REST(y = counts, st=st, t=t, a=a, n_cams=n_cams)) 
  }
  
  return(boot.rest)
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
