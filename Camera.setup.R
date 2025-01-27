library(spsurvey)
library(sf)

#Create a vector of all possible x and y values
x <- rep(seq(1,space_size*.75, by=1),space_size*.75)
y <- rep(seq(1,space_size*.75, by=1),each=space_size*.75)

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
