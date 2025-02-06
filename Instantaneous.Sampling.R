
# Define variables --------------------------------------------------------

P <- A/a # Number of possible cameras (grid cells) in study area

# Count of animals at each location (rows) and occasion (columns)
# on the hour, every hour
sampled_times <- seq(60,n_steps,by=60)
n_occ_is <- length(sampled_times)



# Format data -------------------------------------------------------------

data <- captures[,sampled_times]

data.list <- list(count = matrix(data, ncol = n_occ_is),
                ncam = n_cams,
                nocc = n_occ_is )




# IS functions ------------------------------------------------------------

# Estimate abundance
IS.fn <- function(data, P) {
  Ni <- apply(data$count, 1, sum) #sum across all rows of data$count (counts at each location across all occasions)
  tot <- sum(Ni/a) #summation of n/a's for each camera
  estD <- 1/(data$ncam * data$nocc) * tot
  estN <- estD * A
  return(estN)
}

IS.fn(data.list, P)


# Bootstrap for variance
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

estN.boot <- boot.fn(data.list, nboot = 1000, P = P)

mean(estN.boot)
se_estN <- sd(estN.boot)
se_estN
CI_estN <- quantile(estN.boot, c(0.025, 0.975))
CI_estN


