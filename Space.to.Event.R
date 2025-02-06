library(dplyr)

# Define variables --------------------------------------------------------

# Space to event on each occasion (columns)
n_cams <- n_cams
sampled_minute <- trunc(runif(1,0,60))
sample_order <- sample(c(1:n_cams),n_cams,replace=F) #randomly order the cameras

# Sample for one minute each hour, starting with a random value, length = n_occ
# Assumes animal is constantly available for capture
sampled_times <- seq(sampled_minute,n_steps,by=60)
n_occ_ste <- length(sampled_times)

# Format data -------------------------------------------------------------

#Use just the sampled times
sample <- data.frame(captures[,c(sampled_times)])
sample$camera <- seq(1,n_cams)
sample_order <- paste0(sample_order)
sample <- sample %>% arrange(factor(camera,levels=sample_order))


#Create an empty results matrix
data <- rep(NA,n_occ_ste)


#Loop over cameras and sampling occasions
  for(i in 1:n_occ_ste){
      
      sub_sample <- sample[,i]  
      STE <- match(TRUE, sub_sample >= 1) #Count the number of rows until at least one capture
    
    data[i] <- STE
  }


data <- list(toevent=matrix(data*mean(a), ncol = n_occ_ste),
            censor=n_cams*mean(a), 
            A=A)



# STE functions -----------------------------------------------------------

# Estimate abundance from STE
STE.estN.fn <- function(data){
  opt <- optim(log(1/mean(data$toevent, na.rm = T)), 
               exp.logl.fn, 
               data = data, 
               control = list(fnscale = -1),
               hessian = T)
  # Estimate of lambda
  estlam <- exp(opt$par)
  
  # estlam is average density per m2
  estN <- estlam * data$A
  
  # Delta method for variance
  varB <- -ginv(opt$hessian)
  form <- sprintf("~ %f * exp(x1)", data$A)
  SE_N <- deltamethod(g = as.formula(form), mean = opt$par, cov = varB, ses = T)
  
  return(list(estN = estN,
              SE_N = SE_N,
              LCI = estN - SE_N * 1.96,
              UCI  = estN + SE_N * 1.96) )
}

# Estimate abundance with Space-to-Event
STE.estN.fn(data)
