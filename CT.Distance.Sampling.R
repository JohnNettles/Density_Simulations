library(tidyverse)
library(Distance)

#Total counts
counts <- sum(captures)


# Sampling effort ---------------------------------------------------------

#Temporal effort
t.e <- (n_steps*60)/2 #length of sampling period (in seconds)/time between snapshot moments

#Spatial effort (detection area)
a <- a

#total sampling effort
samp.e <- t.e*a



# Using distance to model detection ----------------------------------------

sampled_secs <- seq(2,n_steps*60,by=2) #This is the list of seconds to use
total_samples <- n_individ * length(sampled_secs)
sampled_x <- numeric(total_samples)
sampled_y <- numeric(total_samples)

index <- 1

for (i in 1:n_individ){
  print(i)
  
  x_seq <- matrix(NA, nrow=n_steps, ncol=120)
  y_seq <- matrix(NA, nrow=n_steps, ncol=120)
  
  for (j in 1:n_steps){
    #Create a vector of x locations between the jth location and the next
    x_seq[j,] <- seq(locations[i,1,j],locations[i,1,j+1],length.out=120)
    y_seq[j,] <- seq(locations[i,2,j],locations[i,2,j+1],length.out=120)
  }
  
  #Sample 60 if those locations, one for each second between steps
  #With replacement because they could stay in the same spot for multiple seconds
  sampled_indices <- sort(sample(1:120,60,replace=T))
  sampled_x_locs <- x_seq[,sampled_indices]
  sampled_y_locs <- y_seq[,sampled_indices]
    
  #Store in matrix
  x_locs[i,] <- as.vector(t(sampled_x_locs))
  y_locs[i,] <- as.vector(t(sampled_y_locs))
    

  #Only include the specified sample moments (i.e., every other second). Combine all points in one long vector
  sampled_x[index:(index+length(sampled_secs)-1)] <- x_locs[i,sampled_secs]
  sampled_y[index:(index+length(sampled_secs)-1)] <- y_locs[i,sampled_secs]
  
  index <- index + length(sampled_secs)
}  


#For each camera, create a matrix with the coordinates of the polygon. See if each point falls within polygon

#Pre-allocate memory
in.view <- matrix(FALSE, nrow=n_cams, ncol=length(sampled_x))
in.view_x <- vector(mode="list",length=n_cams)
in.view_y <- vector(mode="list",length=n_cams)

#Convert detection points to an sf object
points <- st_as_sf(data.frame(x=sampled_x,y=sampled_y),coords=c('x','y'))

for(k in 1:n_cams) {
  print(k)
  #create polygons
  camera_poly <- st_as_sf(st_sfc(sectors_polygons[[k]]))
    
  #Each row of in.view will be a list of 0's and 1's for whether or not that point is in view of camera k.
  in.view[k,] <- st_within(points,camera_poly, sparse=F)
  
  #If the camera is in view, extract the x value and put it
  in.view_x[[k]] <- sampled_x[in.view[k,]]
  in.view_y[[k]] <- sampled_y[in.view[k,]]
}


#Distance to camera
dist <- vector()

for (l in 1:n_cams) {
  x <- in.view_x[[l]] 
  y <- in.view_y[[l]]
  
  x.dist <- x - cameras$X[l]
  y.dist <- y - cameras$Y[l]
  
  dist.new <- sqrt(x.dist^2 + y.dist^2)
  
  dist <- c(distance,dist.new)
}

#Want to remove some of the detections at greater distances to represent missed detections
dist.close <- subset(dist, dist <= 6)
dist.lowmid <- subset(dist, dist > 6 & dist <= 9)
dist.highmid <- subset(dist, dist > 9 & dist <= 12)
dist.far <- subset(dist, dist > 12)

distance <- c(
  dist.close,
  sample(dist.lowmid, size=round(length(dist.lowmid)*0.95)),
  sample(dist.highmid, size=round(length(dist.highmid)*0.85)),
  sample(dist.far, size=round(length(dist.far)*0.7))
)


#Detection function
data <- data.frame(distance)
detect <- ds(data,key='hr',adjustment=NULL,transect='point')
plot(detect)

d.prob.fun <- function(d,b_hat,sigma_hat) {
  d.prob <- as.numeric(1 - exp(-(d/sigma_hat)^(-b_hat)))
  return(d.prob)
}

data <- data.frame(distance)
detect <- ds(data,key='hr',adjustment=NULL,transect='point')
plot(detect)


#This stuff plots the predicted values but isn't needed
b_hat <- detect$ddf$par[1] #shape parameter (b)
sigma_hat <- detect$ddf$par[2] #scale parameter (sigma)


distances <- seq(0,max(data$distance),length.out=100)
detection_probs <- d.prob.fun(d=distances, b_hat=b_hat,sigma_hat=sigma_hat)
data.new <- data.frame(distance=distances, detection_prob=detection_probs)

ggplot(data.new, aes(x = distance, y = detection_prob)) +
  geom_line(color = "blue", lwd = 1)


# Calculate density -------------------------------------------------------
s <- summary(detect)
p_hat <- s$ds$average.p


CT.DS <- function(n,theta,w,t,p_hat){
  d_hat.m2 <- (2*n)/(theta*(w^2)*t*p_hat)
  d_hat.km2 <- d_hat.m2 * 1000000
  return(d_hat.km2)
}
  
b <- CT.DS(n=length(distance), theta=mean(theta.EDA), w=mean(r.EDD), t=t.e, p_hat=p_hat) 
  

# Bootstrap for variance --------------------------------------------------

boot.fn <- function(data,nboot) {
  boot.ctds <- rep(NA,nboot)
  
  for (i in 1:nboot) {
    #Resample detection distances
    boot.dist <- sample(distance,length(distance),replace=T)
    boot.data <- data.frame(boot.dist)
    colnames(boot.data) <- c("distance")
    
    #need to calculate detection probability
    invisible(d <- ds(boot.data,key='hr',adjustment=NULL,transect='point'))
    d <- summary(d)
    boot.p <- d$ds$average.p
    
    #Calculate the average density for each bootstrap iteration
    boot.ctds[i] <- mean(CT.DS(n = length(distance), theta=mean(theta.EDA), w=mean(r.EDD), t=t.e, p_hat=boot.p)) 
  }
  
  return(boot.ctds)
}

#Run the function
estD.boot <- boot.fn(distance, nboot = 100)

#Convert to abundance
meanN_boot <- mean(estD.boot)
meanN_boot
se_estN <- sd(estD.boot)
se_estN
CI_estN <- quantile(estD.boot, c(0.025, 0.975))
CI_estN  
  
  
