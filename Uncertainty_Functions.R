#-------------------------------------------------------------------------------------------------------
#
#     This script analyzes four approaches of uncertainty assessment in allometric equations  
#     Approach one: Slope and intercept sampling based on variance-covariance matrix 
#     Approach two: bootstrapped (sampling with replacement)
#     Approach three: analytical approach (or CI and PI approach) (based on Yanai et al. (2010))
#     Approach four: Bayesian approach adapted from the paper Rejou-Mechain et al. (2017) and the R package "Biomass"
#     
#-------------------------------------------------------------------------------------------------------

# Input data documentary for the function "Uncertain_Slope_inter"
# Mex_mod 1st column: DBH in cm (for developing allometric equation)
# Mex_mod 2nd column: biomass in kg (for developing allometric equation)
# X_Var_caseStudy = log(DBH_case_study), DBH_case_study is the DBH values in your chosen case study sites
# n_iter: the number of iterations
# Output data documentary for the function "Uncertain_Slope_inter": mean AGB value per tree, each value is for one iteration


#### Approach one: Slope and intercept sampling based on variance-covariance matrix---------------------
Uncertain_Slope_inter <-  function(Mex_mod, X_Var_caseStudy, n_iter){
  
  
  ## Obtain the original slope, intercept, and sigma values from model fitting
  BIO_log <- log(Mex_mod[, 2])   # log(Biomass)
  DBH_log <- log(Mex_mod[, 1])   # log(DBH)
  
  mod <- lm(BIO_log ~ DBH_log)
  # summary.lm(mod)
  intercept <- as.numeric(coef(mod)[1])
  slope <- as.numeric(coef(mod)[2])
  MSE <- sigma(mod)^2
  
  n_row <- length(Mex_mod[, 2])
  n_casestudy <- length(X_Var_caseStudy)
  
  ## Prepare slope and intercept, and create a matrix (1st column is intercept, 2nd column is slope)
  X <- data.frame(Intercept = 1, Slope = DBH_log)
  X <- data.matrix(X)
  
  ## Calculate variance-covariance matrix, MSE = 0.2137
  V_C <- MSE * solve(t(X) %*% X)
  
  ## Randomly simulate intercept and slope values, intercept = -0.46213; slope = 1.65874
  cof_mod4 <- c(intercept, slope)
  cof <- rmvnorm(n_iter, mean=cof_mod4, sigma= V_C)
  
  
  ## Calculate AGB
  cof <- data.frame(cof)
  for (i in seq(1:n_iter)){
    # cof$out[[i]] <- exp(cof[i, 1] + cof[i, 2]*X_Var_caseStudy + 0.2137/2)     # with the Baskerville correction
    cof$out1[[i]] <- exp(cof[i, 1] + cof[i, 2]*X_Var_caseStudy)                 # without the Baskerville correction
  }
  
  ## Pre-calculated plot density (cm/ha) based on tree size and nested plot design
  if (n_row ==48){
    PD =6557.6953
  } else if (n_row ==245){
    PD = 322.1505
  } else if (n_row ==93){
    PD = 1527.778
  }
  
  
  ## Convert the AGB from kg to kg per tree
  # AGB_Tot_V_C <- map_dbl(cof$out, sum)/n_casestudy
  AGB_Tot_V_C1 <- (map_dbl(cof$out1, sum)*PD)/n_casestudy
  
  
  return(list(AGB_Tot_V_C1))
  
}
Uncertain_Slope_inter_Para <-  function(Para_Mod, X_Var_caseStudy, n_iter){
  
  
  n_casestudy <- length(X_Var_caseStudy)
  
  ## Prepare slope and intercept, and create a matrix (1st column is intercept, 2nd column is slope)
  X_Var <- Para_Mod[, 1]^2 * Para_Mod[, 2]
  
  X <- data.frame(Intercept = 1, Slope = X_Var)
  X <- data.matrix(X)
  
  ## Calculate variance-covariance matrix, MSE = 0.2137
  V_C <- 0.3070696^2 * solve(t(X) %*% X)
  
  ## Randomly simulate intercept and slope values, intercept = -3.06835; slope = 2.68143
  cof_mod4 <- c(0.02076, 0.00003161)
  cof <- rmvnorm(n_iter, mean=cof_mod4, sigma= V_C)
  
  
  ## Calculate AGB
  cof <- data.frame(cof)
  for (i in seq(1:n_iter)){
    # cof$out[[i]] <- exp(cof[i, 1] + cof[i, 2]*X_Var_caseStudy + 0.2137/2)     # with the Baskerville correction
    cof$out1[[i]] <- cof[i, 1] + cof[i, 2]*X_Var_caseStudy                 # without the Baskerville correction
  }
  
  
  ## Convert the AGB from kg to kg per tree
  # AGB_Tot_V_C <- map_dbl(cof$out, sum)/n_casestudy
  AGB_Tot_V_C1 <- (map_dbl(cof$out1, sum)*953.75)/n_casestudy
  
  
  return(list(AGB_Tot_V_C1))
  
}


#### Approach two: bootstrapped (sampling with replacement)------------------------------------------------------
Uncertain_Boot <- function(Mex_mod,  X_Var_caseStudy, n_iter){
  
  n_casestudy <- length(X_Var_caseStudy)
  
  Y_X <- data.frame(Y = log(Mex_mod[, 2]), X = log(Mex_mod[, 1]))
  
  n_row <- length(Mex_mod[, 2])
  row_index <- replicate(n_iter, sample(n_row, n_row, replace = T))
  
  Y_X_Out <- vector("list", length = n_iter)
  for (i in seq(1:n_iter)){
    Y_X_Out[[i]] <- Y_X[row_index[, i], ]
  }
  
  
  ## Model fitting
  AGB_model <- function(df) {
    lm(Y ~ X, data = df)
  }
  models <- map(Y_X_Out, AGB_model)
  
  
  ## Extract intercept, slope, and sigma from models
  intercept <- map(models, coef) %>%
    map_dbl(1)
  
  slope <- map(models, coef) %>%
    map_dbl(2)
  
  sigma <- map_dbl(models, sigma) 
  
  data4 <- data.frame(intercept, slope, sigma)
  
  
  
  ## Calculate AGB 
  for (i in seq(1:n_iter)){
    #data4$out[[i]] <- exp(data4$intercept[i] + data4$slope[i]*X_Var_caseStudy + (data4$sigma[i])^2/2)
    data4$out1[[i]] <- exp(data4$intercept[i] + data4$slope[i]*X_Var_caseStudy)
  }
  
  
  ## Pre-calculated plot density (cm/ha) based on tree size and nested plot design
  if (n_row ==48){
    PD =6557.6953
  } else if (n_row ==245){
    PD = 322.1505
  } else if (n_row ==93){
    PD = 1527.778
  }
  
  #AGB_Tot_boot <- map_dbl(data4$out, sum)/n_casestudy
  AGB_Tot_boot1 <- (map_dbl(data4$out1, sum)*PD)/n_casestudy
  
  
  return(list(AGB_Tot_boot1))
}
Uncertain_Boot_Para <- function(Para_Mod,  X_Var_caseStudy, n_iter){
  
  n_casestudy <- length(X_Var_caseStudy)
  
  X_Var <- Para_Mod[, 1]^2 * Para_Mod[, 2]
  Y_X <- data.frame(Y = Para_Mod[, 3], X_Var)
  
  row_index <- replicate(n_iter, sample(655, 655, replace = T))
  
  Y_X_Out <- vector("list", length = n_iter)
  for (i in seq(1:n_iter)){
    Y_X_Out[[i]] <- Y_X[row_index[, i], ]
  }
  
  
  ## Model fitting
  AGB_model <- function(df) {
    lm(Y ~ X_Var, data = df)
  }
  models <- map(Y_X_Out, AGB_model)
  
  
  ## Extract intercept, slope, and sigma from models
  intercept <- map(models, coef) %>%
    map_dbl(1)
  
  slope <- map(models, coef) %>%
    map_dbl(2)
  
  sigma <- map_dbl(models, sigma) 
  
  data4 <- data.frame(intercept, slope, sigma)
  
  
  ## Calculate AGB 
  for (i in seq(1:n_iter)){
    #data4$out[[i]] <- exp(data4$intercept[i] + data4$slope[i]*X_Var_caseStudy + (data4$sigma[i])^2/2)
    data4$out1[[i]] <- data4$intercept[i] + data4$slope[i]*X_Var_caseStudy
  }
  
  #AGB_Tot_boot <- map_dbl(data4$out, sum)/n_casestudy
  AGB_Tot_boot1 <- (map_dbl(data4$out1, sum)*953.75)/n_casestudy
  
  
  return(list(AGB_Tot_boot1))
}



#### Approach three: analytical approach (or CI and PI approach)---------------------------------------
Uncertain_Analytical <- function(Mex_mod, X_Var_caseStudy, n_iter){
  
  n_casestudy <- length(X_Var_caseStudy)
  n_row <- length(Mex_mod[, 2])
  
  ## Obtain the original slope, intercept, and sigma values from model fitting
  BIO_log <- log(Mex_mod[, 2])   # log(Biomass)
  DBH_log <- log(Mex_mod[, 1])   # log(DBH)
  
  mod <- lm(BIO_log ~ DBH_log)
  # summary.lm(mod)
  intercept <- as.numeric(coef(mod)[1])
  slope <- as.numeric(coef(mod)[2])
  MSE <- sigma(mod)^2
  
  
  AGB_CI <- matrix(0, nrow=n_casestudy, ncol = n_iter)
  AGB_PI <- matrix(0, nrow=n_casestudy, ncol = n_iter)
  
  denominator <- sum((DBH_log - mean(DBH_log))^2)
  

  # The original approach to sample sigma, sample all sigma values and store them in matrix
  sigma_s_CI <- rnorm(n = n_iter, mean = 0, sd = sigma(mod))
  sigma_s_PI <- data.frame(replicate(n_casestudy, rnorm(n = n_iter, mean = 0, sd = sigma(mod))))
  
  
  for (j in seq(1: n_iter)){     # iteration
    for (i in seq(1:n_casestudy)){         # number of trees
      numerator <- (X_Var_caseStudy[i] - mean(DBH_log))^2
      
      CI <-  sigma_s_CI[j] * sqrt((1/n_row) + (numerator/denominator))
      AGB_CI[i, j] <- exp(intercept + slope * X_Var_caseStudy[i] + CI)
      
      PI <- sigma_s_PI[j, i] * sqrt(1 + (1/n_row) + (numerator/denominator))
      AGB_PI[i, j] <- exp(intercept + slope * X_Var_caseStudy[i] + PI)
    }
  }
  
  ## Pre-calculated plot density (cm/ha) based on tree size and nested plot design
  if (n_row ==48){
    PD =6557.6953
  } else if (n_row ==245){
    PD = 322.1505
  } else if (n_row ==93){
    PD = 1527.778
  }
  
  # Convert kg to kg per tree
  AGB_CI_out <- (colSums(AGB_CI)*PD)/n_casestudy      
  AGB_PI_out <- (colSums(AGB_PI)*PD)/n_casestudy   
  
  return(list(AGB_CI_out, AGB_PI_out))
  
}
Uncertain_Analytical_Para <- function(Para_Mod, X_Var_caseStudy, n_iter){
  
  n_casestudy <- length(X_Var_caseStudy)
  
  AGB_CI <- matrix(0, nrow=n_casestudy, ncol = n_iter)
  AGB_PI <- matrix(0, nrow=n_casestudy, ncol = n_iter)
  
  X_Var <- Para_Mod[, 1]^2 * Para_Mod[, 2]
  denominator <- sum((X_Var - mean(X_Var))^2)
  
  # The original approach to sample sigma, sample all sigma values and store them in matrix
  sigma_s_CI <- rnorm(n = n_iter, mean = 0, sd = 0.3070696)
  sigma_s_PI <- data.frame(replicate(n_casestudy, rnorm(n = n_iter, mean = 0, sd = 0.3070696)))
  
  
  for (j in seq(1: n_iter)){     # iteration
    for (i in seq(1:n_casestudy)){         # number of trees
      numerator <- (X_Var_caseStudy[i] - mean(X_Var))^2
      
      CI <-  sigma_s_CI[j] * sqrt((1/655) + (numerator/denominator))
      AGB_CI[i, j] <- 0.02076 + 0.00003161 * X_Var_caseStudy[i] + CI
      
      PI <- sigma_s_PI[j, i] * sqrt(1 + (1/655) + (numerator/denominator))
      AGB_PI[i, j] <- 0.02076 + 0.00003161 * X_Var_caseStudy[i] + PI
    }
  }
  
  AGB_CI_out <- (colSums(AGB_CI)*953.75)/n_casestudy      
  AGB_PI_out <- (colSums(AGB_PI)*953.75)/n_casestudy   
  
  return(list(AGB_CI_out, AGB_PI_out))
  
}


#### Approach four: Bayesian approach adapted from the paper Rejou-Mechain et al (2017) and the R package "Biomass"------

## Functions to run the Bayesian approach
likelihood <- function(y, X, param){
  theta = param[-length (param)]
  sd = param[length(param)]
  pred = X %*% theta
  singlelikelihoods = dnorm (y, mean = pred, sd = sd, log = T)
  sumll = sum(singlelikelihoods)
  return(sumll)
}
prior <- function(param){
  theta = param[-length (param)]
  sd = param[length(param)]
  thetaprior = dnorm(theta, mean=0, sd= 10, log = T)
  sdprior = dunif (sd, min= 0, max= 30, log = T)
  return(sum(thetaprior)+sdprior)
}
posterior <- function(y, X, param){
  return ( likelihood(y, X, param) + prior (param))
}
proposalfunction <- function(param, sdprop){
  return(rnorm ( length(param) , mean = param, sd= sdprop))
}
run_metropolis_MCMC <- function(y, X, startvalue, iterations, sdprop){
  # chain : row i = set of parameters at iteration i
  chain = array (dim = c (iterations+1, length(startvalue)))
  colnames (chain) <- names (startvalue)
  lp = array (dim= iterations+1)
  chain[ 1,] = startvalue
  for (i in 1:iterations){
    proposal = proposalfunction (chain[i,], sdprop)
    probab = exp(posterior (y, X, proposal) - posterior(y, X, chain[i,]))
    # probab ranges between 0 and 1
    
    # probab > 1 <=> new proposal is more likely considering the data
    if ( runif (1) < probab){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,]
    }
    lp[i +1 ] = likelihood (y, X, chain[i+1,])
  }
  return(cbind (chain, lp))
}


## Function to perform uncertainty analysis based on the calculated intercept and slope 
Uncertain_Baye <- function(Mex_mod, X_Var_caseStudy, param4){
  
  n_iter <- nrow(param4)
  n_casestudy <- length(X_Var_caseStudy)
  out <- vector("numeric", length = n_iter)
  
  ## Pre-calculated plot density (cm/ha) based on tree size and nested plot design
  n_row <- length(Mex_mod[, 2])
  if (n_row ==48){
    PD =6557.6953
  } else if (n_row ==245){
    PD = 322.1505
  } else if (n_row ==93){
    PD = 1527.778
  }
  
  for (i in seq(1:n_iter)){
    out[[i]] <- (sum(exp(param4$intercept[i] + param4$logagbt[i]*X_Var_caseStudy))*PD)/n_casestudy
  }
  
  return(out)
}
Uncertain_Baye_Para <- function(X_Var_caseStudy, param4){
  
  n_iter <- nrow(param4)
  n_casestudy <- length(X_Var_caseStudy)
  out <- vector("numeric", length = n_iter)
  
  for (i in seq(1:n_iter)){
    out[[i]] <- (sum(param4$intercept[i] + param4$logagbt[i]*X_Var_caseStudy)*953.75)/n_casestudy
  }
  
  return(out)
}






