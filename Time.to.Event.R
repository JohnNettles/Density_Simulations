

# Sampling details
n_cams <- 15
n_occ <- 20
censor <- 4 # K (number of sampling periods per occasion)
mvmt_speed <- 60 #Time it takes individual to cross the entire field of view (in minutes)

sampled_times <- sort(c(seq(1,n_steps,8),seq(2,n_steps,8),seq(3,n_steps,8),seq(4,n_steps,8)))
#Will need to trim the capture history to just include the sampled times rather than all time.





# Example dataset

#would need to define number of occasions and number of sampling periods
#If individual is inside a, then it would be a 1, otherwise would be a 0


dat.tte <- list(toevent = matrix(c(NA, NA, 1, NA, NA, NA, NA, 4, NA, 2, 1, NA, 
                                   NA, NA, NA, NA, NA, 3, NA, 2), ncol = nocc),
                censor = censor,
                A = A,
                mean_a = a)
dat.tte




# Exponential log likelihood
exp.logl.fn <- function(data, param){
  # param: beta parameter for lambda
  lambda <- exp(param) 
  logL <- 0
  for(i in 1:nrow(data$toevent)) {
    for (j in 1:ncol(data$toevent)) {
      if(!is.na(data$toevent[i, j])) {
        tmp <- dexp(data$toevent[i, j], lambda)
      } else {
        tmp <- pexp(data$censor, lambda, lower.tail = F)
      }
      logL <- logL + log(tmp)
    }
  }
  return(logL)
}

# Estimate abundance with Time-to-Event
TTE.estN.fn <- function(data){
  opt <- optim(log(1/mean(data$toevent, na.rm = T)), 
               exp.logl.fn, 
               data = data, 
               control = list(fnscale = -1),
               hessian = T)
  # Estimate of lambda
  estlam <- exp(opt$par)
  
  # If estlam is average abundance in a camera area
  P <- data$A/data$mean_a # Number of possible cameras (grid cells) in study area
  estN <- estlam * P
  
  # Variance with msm::deltamethod
  varB <- -ginv(opt$hessian)
  form <- sprintf("~ %f * exp(x1)", P)
  SE_N <- deltamethod(g = as.formula(form), mean = opt$par, cov = varB, ses = T)
  
  return(list(estN = estN,
              SE_N = SE_N,
              LCI = estN - SE_N * 1.96,
              UCI  = estN + SE_N * 1.96) )
}

TTE.estN.fn(dat.tte)

## Space-to-Event
# Space to event on each occasion (columns)
ncam <- 8
nocc <- 12
dat.ste <- list(toevent = matrix(c(NA, NA, 1, NA, NA, NA, NA, 4, NA, 8, 1, NA) * a, 
                                 ncol = nocc),
                censor = ncam * a,
                A = A)
dat.ste