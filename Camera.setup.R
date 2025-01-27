library(spsurvey)
library(sf)


# GRTS sampling -----------------------------------------------------------

#Create a vector of all possible x and y values
sa.width <- space_size*0.75 #don't sample the entire region, could also specify a non-square grid
sa.length <- space_size*0.75 #don't sample the entire region, could also specify a non-square grid

x <- rep(seq(1,sa.length, by=1),sa.width)
y <- rep(seq(1,sa.width, by=1),each=sa.length)

#Put these into a data frame
coords <- data.frame(x,y)

#Convert them to a spatial data frame
points_sf <- st_as_sf(coords,coords=c('x','y'))

#Set the number of cameras (sampling points)
n_cams <- 15

#Conduct GRTS sample
grts_sample <- grts(points_sf, n_base = n_cams, stratum_var = NULL, seltype = 'equal',projcrs_check = F)

#Create data frame with x/y values and ID for the sample results
cameras <- data.frame(
  ID = seq(1,n_cams,by=1),
  X = grts_sample$sites_base$X,
  Y = grts_sample$sites_base$Y
)

#Plot results
plot(coords)
points(x=cameras$X,y=cameras$Y,col='red',pch=16, cex=0.5)



# Camera Viewsheds --------------------------------------------------------

# Specifications
r <- 15 # camera maximum distance in meters
theta <- 42 # camera viewshed in degrees
a <- pi * r^2 * theta/360 # Area of a single camera in square meters
A <- sa.width*sa.length # area of study area in square meters

# Sector arc details
theta_start <- 0*(pi/180) #Start all angles at 0, convert to radians
theta_end <- theta*(pi/180) #End the angles at theta, convert to radians
n_points <- 100 #Number of points used to create arc


#Create a list with parameter values for each sector
sectors_params <- list()

for (i in 1:n_cams){
  sectors_params[[i]] <- list(x=cameras$X[i],y=cameras$Y[i], r=r, theta_start=theta_start, theta_end=theta_end)
  
}

#Define a function to create each sector based on parameter values  
create_sector <- function(x,y,r,theta_start,theta_end,n_points=100) {
  
  #Generate arc points
  angles <- seq(theta_start, theta_end, length.out=n_points)
  arc_x <- x + r*cos(angles)
  arc_y <- y + r*sin(angles)
  
  #Combine points into polygon
  polygon_coords <- rbind(
    c(x, y),
    cbind(arc_x, arc_y),
    c(x, y)
  )
  
  #Create a spatial polygon
  sector_polygon <- st_polygon(list(polygon_coords))
  #sector_sf <- st_sfc(sector_polygon)
  
}  


#Create a list of polygons by applying the create sector function to each value of the sectors_params list  
sectors_polygons <- lapply(sectors_params, function(params) {
  create_sector(params$x,params$y,params$r,params$theta_start,params$theta_end)
}) 


#Create matrix to store capture history
captures <- matrix(NA,nrow=n_cams,ncol=n_steps)



# Create capture histories ------------------------------------------------

#For each camera, create a matrix with the coordinates of the polygon
for(i in 1:n_cams) {
  points <- st_coordinates(sectors_polygons[[i]])
  
  #For each time step, determine if there were any points within that polygon 
  for(z in 1:n_steps) {
    
    #This produces a vector of 0's and 1's, one value for each individual
    c.z <- point.in.polygon(locations[,1,z],locations[,2,z],points[,1],points[,2])
    
    #To create our unmarked capture history: If there is at least one individual detected
    #during that occasion, then that cell is marked as a 1. 
    #So a 1 in row 4, column 12 means that the animal was detected at camera 4 during time step 12
    captures[i,z] <- ifelse(sum(c.z) >= 1,1,0)
    
  }
}

