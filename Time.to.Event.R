library(msm)
library(MASS)


# Sampling details --------------------------------------------------------

n_cams <- n_cams
censor <- 4 # K (number of sampling periods per occasion)
mvmt_speed <- 50 #Time it takes individual to cross the entire field of view (in minutes)
occ_length <- censor * mvmt_speed # (length of each sampling occasion in minutes)
n_occ_tte <- 0.5 * n_steps/occ_length # (number of sampling occasions (using only every other potential sampling period))

# sample every other potential occasion. Each occasion includes a number of sampling periods equal to 'censor'
# sampling period length is determined by the amount of time it takes an animal to cross the widest part of the FOV
# Assumes animal is constantly available for capture
sampled_times <- vector()

for (i in 1:n_occ_tte) {
  
  start <- 1 + occ_length*2*(i-1)
  end <- start + (mvmt_speed * censor)
  x <- seq(from=start, to=end)
  sampled_times <- c(st,x)
  
}

#Use just the sampled times
sample <- captures[,c(sampled_times)]

#Create an empty results matrix
data <- matrix(NA,nrow=n_cams,ncol=n_occ_tte)



# Format data -------------------------------------------------------------

#Loop over cameras and sampling occasions
for(b in 1:n_cams) {
  for(i in 1:n_occ_tte){
    first_period <- (censor*i)-(censor-1) #First column to use
    last_period <- (censor*i) #Last column to use
    
    sub_sample <- sample[b,first_period:last_period,drop=F] #Keep only the sampling periods that fall within this occasion
    TTE <- match(TRUE, sub_sample >= 1) #Count the number of columns until a 1 (capture)
    
    data[b,i] <- TTE
  }
}


#Store data as a list to input into TTE function
data.list <- list(toevent = matrix(data, ncol=n_occ_tte),
             censor = censor,
             A=A,
             mean_a=a)



# TTE functions -----------------------------------------------------------

# Exponential log likelihood
exp.logl.fn <- function(data, param){
  # param: beta parameter for lambda
  lambda <- exp(param) 
  logL <- 0
  for(i in 1:nrow(data$toevent)) {
    for (j in 1:ncol(data$toevent)) {
      if(!is.na(data$toevent[i, j])) {  #If the value is not NA (meaning there was a detection)
        tmp <- dexp(data$toevent[i, j], lambda) #Probability density function for the time to event
      } else {
        tmp <- pexp(data$censor, lambda, lower.tail = F) #Cumulative density for greater than our censor
      }
      logL <- logL + log(tmp) #Calculates log likelihood (including right-censored data) by summing all values of log(p(x))
    }
  }
  return(logL)
}


# Estimate abundance with Time-to-Event
TTE.estN.fn <- function(data){
  opt <- optim(par=log(1/mean(data$toevent, na.rm = T)), #initial value is given by 1/sample mean
               fn=exp.logl.fn, #Maximizing for our log likelihood function
               data = data, 
               control = list(fnscale = -1), #Turns it into a maximization problem
               hessian = T)
  # Estimate of lambda
  estlam <- exp(opt$par)
  
  # If estlam is average abundance in a camera area
  P <- data$A/data$mean_a # Number of possible cameras (grid cells) in study area
  estN <- estlam * P
  
  # Variance with msm::deltamethod
  varB <- -ginv(opt$hessian)
  form <- sprintf("~ %f * exp(x1)", P)
  SE_N <- msm::deltamethod(g = as.formula(form), mean = opt$par, cov = varB, ses = T)
  
  return(list(estN = estN,
              SE_N = SE_N,
              LCI = estN - SE_N * 1.96,
              UCI  = estN + SE_N * 1.96) )
}

#Calculate estimated N
TTE.estN.fn(data.list)
