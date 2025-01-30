P <- A/a # Number of possible cameras (grid cells) in study area

# Count of animals at each location (rows) and occasion (columns)
# on the hour, every hour
sampled_times <- seq(60,n_steps,by=60)
n_occ_is <- length(sampled_times)


data <- captures[,sampled_times]

data.list <- list(count = matrix(sample, ncol = n_occ_is),
                ncam = n_cams,
                nocc = n_occ_is )


# Estimate abundance
IS.fn <- function(data, P) {
  Nj <- apply(data$count, 2, sum)
  tot.ct <- sum(Nj)
  estN <- P / (data$ncam * data$nocc) * tot.ct
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

estN.boot <- boot.fn(data.list, nboot = 100, P = P)

se_estN <- sd(estN.boot)
se_estN
CI_estN <- quantile(estN.boot, c(0.025, 0.975))
CI_estN


