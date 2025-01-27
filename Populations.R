set.seed(14)


# Sample from home range --------------------------------------------------

n_individ <- 10         #Population size
n_steps <- 100          #Number of time steps
space_size <- 100       #Size of 2D region
max_speed <- 10          #Max speed in m/minute
HR_sd_max <- 20         #Maximum home range sd
HR_sd_min <- 5          #Minimum home range sd


#Uniform distribution to sample standard deviation (home range size)
HR_sd <- runif(n_individ, HR_sd_min, HR_sd_max) 

#Sample home range centers for each individual
HR_center <- data.frame(
  ID = 1:n_individ,
  x = runif(n_individ, 0, space_size),
  y = runif(n_individ, 0, space_size)
)

#Matrix to store positions over time
locations <- array(NA, dim=c(n_individ, 2, n_steps+1))
locations[,,1] <- as.matrix(HR_center[,c("x","y")])


#Simulate movement
for (t in 2:(n_steps+1)){
  for (i in 1:n_individ){
    
  #For every individual, sample a random point from home range
  proposed_x <- rnorm(1, mean=HR_center$x[i], sd=HR_sd[i])  
  proposed_y <- rnorm(1, mean=HR_center$y[i], sd=HR_sd[i])
  
  #Compute distance
  dx <- proposed_x - locations[i,1,t-1]
  dy <- proposed_y - locations[i,2,t-1]
  step_distance <- sqrt(dx^2 + dy^2)
  
  #If too far, move to max distance
  over_limit <- step_distance > max_speed
  dx[over_limit] <- dx[over_limit]/step_distance[over_limit]*max_speed
  dy[over_limit] <- dy[over_limit]/step_distance[over_limit]*max_speed
  
  #Define new positions
  locations[i,1,t] <- locations[i,1,t-1]+dx
  locations[i,2,t] <- locations[i,2,t-1]+dy

  #Keep within bounds
  locations[i,1,t] <- pmax(0,pmin(space_size, locations[i,1,t]))
  locations[i,2,t] <- pmax(0,pmin(space_size, locations[i,2,t])) 

  }
  }


# Recursive ---------------------------------------------------------------

library(circular)

n_individ <- 10               #Population size
n_steps <- 1000                #Number of time steps
space_size <- 100             #Size of 2D region
max_speed <- 4               #Max speed
cor.angle <- 0.6              #Correlation between successive angles
cor.speed <- 0.4              #Correlation between successive speeds
acceleration_strength <- 0.2  #Strength of acceleration towards center


#Sample home range centers for each individual
HR_center <- data.frame(
  ID = 1:n_individ,
  x = runif(n_individ, 0, space_size),
  y = runif(n_individ, 0, space_size)
)

#Sample starting speed and angle
speeds<- runif(n_individ, 0, max_speed)
angles <- runif(n_individ, 0, 2*pi)

#Matrix to store positions over time
locations <- array(NA, dim=c(n_individ, 2, n_steps+1))
locations[,,1] <- as.matrix(HR_center[,c("x","y")])
steps <- array(NA, dim=c(n_individ,2,n_steps+1))
steps[,1,1] <- speeds
steps[,2,1] <- angles

#Simulate movement
for (t in 2:(n_steps+1)){
  for (i in 1:n_individ){

#Create random variation in speed. If its greater than max, set to max. If its less than 0, set to 0.
noise_speed <- rnorm(1, mean=0, sd=sqrt(1-cor.speed^2)*max_speed)
steps[i,1,t] <- cor.speed * steps[i,1,t-1] + noise_speed
steps[i,1,t] <- ifelse(steps[i,1,t] >= max_speed, max_speed, ifelse(steps[i,1,t] <= 0, 0, steps[i,1,t]))

#Calculate angle to the center of home range
dx_to_center <- HR_center$x[i] - locations[i,1,t-1]
dy_to_center <- HR_center$y[i] - locations[i,2,t-1]
angle_to_center <- atan2(dy_to_center, dx_to_center)

#Sample random angle using von Mises distribution
new_angle <- circular::rvonmises(1, mu=circular(steps[i,2,t-1]), kappa=10*cor.angle)
new_angle <- as.numeric(new_angle)

#Angle plus bias toward the home range center
mean_angle <- atan2(sin(new_angle) + acceleration_strength * sin(angle_to_center),
                    cos(new_angle) + acceleration_strength * cos(angle_to_center))

#Updated direction
steps[i,2,t] <- mean_angle

#Update velocities
movement_x <- steps[i,1,t] * cos(steps[i,2,t])
movement_y <- steps[i,1,t] * sin(steps[i,2,t])

#Update positions
locations[i,1,t] <- locations[i,1,t-1] + movement_x
locations[i,2,t] <- locations[i,2,t-1] + movement_y

#Keep within bounds
locations[i,1,t] <- pmax(0,pmin(space_size, locations[i,1,t]))
locations[i,2,t] <- pmax(0,pmin(space_size, locations[i,2,t])) 

angles[i] <- mean_angle
}
}



#################################################################################################
#Plot results

library(ggplot2)

movement_data <- data.frame(
  ID = rep(1:n_individ, n_steps + 1),
  Time = rep(0:n_steps, each = n_individ),
  x = as.vector(locations[, 1, ]),
  y = as.vector(locations[, 2, ])
)

ggplot(movement_data[movement_data$ID <= 3, ], aes(x = x, y = y, color = factor(ID))) +
  geom_path() +
  geom_point() +
  labs(title = "Movement of Individuals in 2D Space", x = "X Position", y = "Y Position") +
  theme_minimal()
