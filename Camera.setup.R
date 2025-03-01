library(spsurvey)
library(sf)
library(sp)


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
n_cams <- 25

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
r <- 15 # camera maximum distance in meters, if constant
r.EDD <- runif(n_cams,min=6,max=15) #Simulate effective detection distance
theta <- 42 # camera viewshed in degrees, if constant
theta.EDA <- runif(n_cams,min=20, max=42) #Simulate effective detection angle
a <- pi * r.EDD^2 * theta.EDA/360 # Area of a each camera in square meters
A <- sa.width*sa.length # area of study area in square meters

# Sector arc details
theta_start <- (90-0.5*theta.EDA)*(pi/180) #Center all sectors at 90 degrees, convert to radians
theta_end <- (90+0.5*theta.EDA)*(pi/180) #Center all sectors at 90 degrees, convert to radians
n_points <- 100 #Number of points used to create arc


#Create a list with parameter values for each sector
sectors_params <- lapply(1:n_cams, function(i) {
  list(x = cameras$X[i],
       y = cameras$Y[i],
       r = r.EDD[i],
       theta_start = theta_start[i],
       theta_end = theta_end[i])
})


#Define a function to create each sector based on parameter values  
create_sector <- function(x,y,r,theta_start,theta_end,n_points=100) {
  
    
  #Generate arc points
  angles <- seq(theta_start, theta_end, length.out=n_points)
  arc_x <- x + r*cos(angles)
  arc_y <- y + r*sin(angles)
  
  #Combine points into polygon
  polygon_coords <- rbind(
    c(x, y), #camera position
    cbind(arc_x, arc_y),
    c(x, y)
  )

  #Return a spatial polygon
  st_polygon(list(polygon_coords))
}  


#Create a list of polygons by applying the create sector function to each value of the sectors_params list  
sectors_polygons <- lapply(sectors_params, function(params) {
  create_sector(params$x,params$y,params$r,params$theta_start,params$theta_end)
}) 

sectors_sf <- st_sfc(sectors_polygons)
plot(sectors_sf)

# Create capture histories ------------------------------------------------

#Create matrix to store capture history
captures <- matrix(NA,nrow=n_cams,ncol=n_steps)

#For each camera, create a matrix with the coordinates of the polygon
for(i in 1:n_cams) {
  points <- st_coordinates(sectors_polygons[[i]])
  
  #For each time step, determine if there were any points within that polygon 
  for(z in 1:n_steps) {
    
    #This produces a vector of 0's and 1's, one value for each individual
    c.z <- point.in.polygon(locations[,1,z],locations[,2,z],points[,1],points[,2])
    
    #To create our unmarked capture history: Count the number of detections at each camera during that occasion
    #So a 3 in row 4, column 12 means that 3 animals were detected at camera 4 during time step 12
    captures[i,z] <- ifelse(sum(c.z) >= 1,sum(c.z),0)
    
  }
}
