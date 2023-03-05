#### Load package
library(mvtnorm)
library(tidyverse)


#### Applying four uncetainty methods for the Yucatán---------------------------------------------------
# output mean value for Yucatan is kg/ha (divided by 1000 to converted to Mg/ha)

#### (1) the slope-intercept method
Yuca_SI_10 <- Uncertain_Slope_inter(Yuca_Mod, Yuca_Tree10_log, 1000)
Yuca_SI_100 <- Uncertain_Slope_inter(Yuca_Mod, Yuca_Tree100_log, 1000)
Yuca_SI_1000 <- Uncertain_Slope_inter(Yuca_Mod, Yuca_Tree1000_log, 1000)
Yuca_SI_10000 <- Uncertain_Slope_inter(Yuca_Mod, Yuca_Tree10000_log, 1000)

# Report CV values as indicator of uncertainty (Unit: %)
sd(Yuca_SI_10[[1]])/mean(Yuca_SI_10[[1]]) *100     
sd(Yuca_SI_100[[1]])/mean(Yuca_SI_100[[1]]) *100     
sd(Yuca_SI_1000[[1]])/mean(Yuca_SI_1000[[1]]) *100     
sd(Yuca_SI_10000[[1]])/mean(Yuca_SI_10000[[1]]) *100     



#### (2) Bootstrapped sampling
Yuca_Boot_10 <- Uncertain_Boot(Yuca_Mod, Yuca_Tree10_log, 1000)
Yuca_Boot_100 <- Uncertain_Boot(Yuca_Mod, Yuca_Tree100_log, 1000)
Yuca_Boot_1000 <- Uncertain_Boot(Yuca_Mod, Yuca_Tree1000_log, 1000)
Yuca_Boot_10000 <- Uncertain_Boot(Yuca_Mod, Yuca_Tree10000_log, 1000)

# Report CV values as indicator of uncertainty (Unit: %)
sd(Yuca_Boot_10[[1]])/mean(Yuca_Boot_10[[1]]) *100     
sd(Yuca_Boot_100[[1]])/mean(Yuca_Boot_100[[1]]) *100     
sd(Yuca_Boot_1000[[1]])/mean(Yuca_Boot_1000[[1]]) *100     
sd(Yuca_Boot_10000[[1]])/mean(Yuca_Boot_10000[[1]]) *100     



#### (3) Analytical approach
Yuca_Analy_10 <- Uncertain_Analytical(Yuca_Mod, Yuca_Tree10_log, 100)
Yuca_Analy_100 <- Uncertain_Analytical(Yuca_Mod, Yuca_Tree100_log, 100)
Yuca_Analy_1000 <- Uncertain_Analytical(Yuca_Mod, Yuca_Tree1000_log, 100)
Yuca_Analy_10000 <- Uncertain_Analytical(Yuca_Mod, Yuca_Tree10000_log, 100)

# Report CV values as indicator of uncertainty (Unit: %)
# Yuca_Analy_Number[[1]] is for CI; Yuca_Analy_Number[[2]] is for PI;
sd(Yuca_Analy_10[[1]])/mean(Yuca_Analy_10[[1]]) *100     
sd(Yuca_Analy_10[[2]])/mean(Yuca_Analy_10[[2]]) *100     

sd(Yuca_Analy_100[[1]])/mean(Yuca_Analy_100[[1]]) *100     
sd(Yuca_Analy_100[[2]])/mean(Yuca_Analy_100[[2]]) *100     

sd(Yuca_Analy_1000[[1]])/mean(Yuca_Analy_1000[[1]]) *100     
sd(Yuca_Analy_1000[[2]])/mean(Yuca_Analy_1000[[2]]) *100     

sd(Yuca_Analy_10000[[1]])/mean(Yuca_Analy_10000[[1]]) *100     
sd(Yuca_Analy_10000[[2]])/mean(Yuca_Analy_10000[[2]]) *100     


#### (4) Bayesian approach
# Calculate intercept, slope and sigma from the Bayesian approach
# The example below produce 1100000 set of [intercept, slope and sigma], discard the first 100000, and then sample 1 out of 100
Y4 <- log(Yuca_Mod[, 2])
X4 = as.matrix(data.frame(intercept=1, logagbt=log(Yuca_Mod[, 1])))
start4 = c(0,0,1)
names(start4) <- c(colnames (X4), "sd")

chain4 <- run_metropolis_MCMC(Y=Y4, X=X4, startvalue = start4, iterations = 1100000, sdprop = c(0.01, 0.002, 0.01))
burnIn <- 100000
param4 <- data.frame(chain4[seq(burnIn +1, dim(chain4)[1 ], by= 100), ])


## Applying the calculated intercept, slope and sigma to tree10, tree100, tree1000, tree10000
case_study <- list(Yuca_Tree10_log, Yuca_Tree100_log, Yuca_Tree1000_log, Yuca_Tree10000_log)

out_case_study <- vector("list", length = 4)
out_CV <- vector("numeric", length = 4)

for (i in seq(1:4)){
  out_case_study[[i]] <- Uncertain_Baye(Yuca_Mod, case_study[[i]], param4)
  out_CV[[i]] <- sd(out_case_study[[i]]) / mean(out_case_study[[i]])
}

# Report CV values as indicator of uncertainty (Unit: %)
out_CV*100   






#### Applying four uncetainty methods for the Hawaii---------------------------------------------------
# output mean value for Hawaii is kg/ha (divided by 1000 to converted to Mg/ha)

#### (1) the slope-intercept method
Hawa_SI_10 <- Uncertain_Slope_inter(Hawa_Mod, Hawa_Tree10_log, 1000)
Hawa_SI_100 <- Uncertain_Slope_inter(Hawa_Mod, Hawa_Tree100_log, 1000)
Hawa_SI_1000 <- Uncertain_Slope_inter(Hawa_Mod, Hawa_Tree1000_log, 1000)
Hawa_SI_10000 <- Uncertain_Slope_inter(Hawa_Mod, Hawa_Tree10000_log, 1000)

# Report CV values as indicator of uncertainty (Unit: %)
sd(Hawa_SI_10[[1]])/mean(Hawa_SI_10[[1]]) *100     
sd(Hawa_SI_100[[1]])/mean(Hawa_SI_100[[1]]) *100     
sd(Hawa_SI_1000[[1]])/mean(Hawa_SI_1000[[1]]) *100     
sd(Hawa_SI_10000[[1]])/mean(Hawa_SI_10000[[1]]) *100     



#### (2) Bootstrapped sampling
Hawa_boot_10 <- Uncertain_Boot(Hawa_Mod,  Hawa_Tree10_log, 1000)
Hawa_boot_100 <- Uncertain_Boot(Hawa_Mod, Hawa_Tree100_log, 1000)
Hawa_boot_1000 <- Uncertain_Boot(Hawa_Mod, Hawa_Tree1000_log, 1000)
Hawa_boot_10000 <- Uncertain_Boot(Hawa_Mod, Hawa_Tree10000_log, 1000)

# Report CV values as indicator of uncertainty (Unit: %)
sd(Hawa_boot_10[[1]])/mean(Hawa_boot_10[[1]]) *100     
sd(Hawa_boot_100[[1]])/mean(Hawa_boot_100[[1]]) *100     
sd(Hawa_boot_1000[[1]])/mean(Hawa_boot_1000[[1]]) *100    
sd(Hawa_boot_10000[[1]])/mean(Hawa_boot_10000[[1]]) *100     


#### (3) Analytical approach
Hawa_Analy_10 <- Uncertain_Analytical(Hawa_Mod, Hawa_Tree10_log, 100)
Hawa_Analy_100 <- Uncertain_Analytical(Hawa_Mod, Hawa_Tree100_log, 100)
Hawa_Analy_1000 <- Uncertain_Analytical(Hawa_Mod, Hawa_Tree1000_log, 100)
Hawa_Analy_10000 <- Uncertain_Analytical(Hawa_Mod, Hawa_Tree10000_log, 100)


# Report CV values as indicator of uncertainty (Unit: %)
# Hawa_Analy_Number[[1]] is for CI; Hawa_Analy_Number[[2]] is for PI;
sd(Hawa_Analy_10[[1]])/mean(Hawa_Analy_10[[1]]) *100     
sd(Hawa_Analy_10[[2]])/mean(Hawa_Analy_10[[2]]) *100     

sd(Hawa_Analy_100[[1]])/mean(Hawa_Analy_100[[1]]) *100     
sd(Hawa_Analy_100[[2]])/mean(Hawa_Analy_100[[2]]) *100     

sd(Hawa_Analy_1000[[1]])/mean(Hawa_Analy_1000[[1]]) *100     
sd(Hawa_Analy_1000[[2]])/mean(Hawa_Analy_1000[[2]]) *100     

sd(Hawa_Analy_10000[[1]])/mean(Hawa_Analy_10000[[1]]) *100     
sd(Hawa_Analy_10000[[2]])/mean(Hawa_Analy_10000[[2]]) *100     



#### (4) Bayesian approach
# Calculate intercept, slope and sigma from the Bayesian approach
# The example below produce 1100000 set of [intercept, slope and sigma], discard the first 100000, and then sample 1 out of 100
Y4 <- log(Hawa_Mod[, 2])
X4 = as.matrix(data.frame(intercept=1, logagbt=log(Hawa_Mod[, 1])))
start4 = c(0,0,1)
names(start4) <- c(colnames (X4), "sd")

chain4 <- run_metropolis_MCMC(y=Y4, X=X4, startvalue = start4, iterations = 1100000, sdprop = c(0.01, 0.002, 0.01))
burnIn <- 100000
param4 <- data.frame(chain4[seq(burnIn +1, dim(chain4)[1 ], by= 100), ])


## Applying the calculated intercept, slope and sigma to tree10, tree100, tree1000, tree10000
case_study <- list(Hawa_Tree10_log, Hawa_Tree100_log, Hawa_Tree1000_log, Hawa_Tree10000_log)

out_case_study <- vector("list", length = 4)
out_CV <- vector("numeric", length = 4)

for (i in seq(1:4)){
  out_case_study[[i]] <- Uncertain_Baye(Hawa_Mod, case_study[[i]], param4)
  out_CV[[i]] <- sd(out_case_study[[i]]) / mean(out_case_study[[i]])
}

# Report CV values as indicator of uncertainty (Unit: %)
out_CV*100   







#### Applying four uncetainty methods for the Chaco---------------------------------------------------
# output mean value for Chaco is Mg/ha 

#### (1) the slope-intercept method
Chaco_SI_10 <- Uncertain_Slope_inter(Chaco_Mod, Chaco_Tree10_log, 1000)
Chaco_SI_100 <- Uncertain_Slope_inter(Chaco_Mod, Chaco_Tree100_log, 1000)
Chaco_SI_1000 <- Uncertain_Slope_inter(Chaco_Mod, Chaco_Tree1000_log, 1000)
Chaco_SI_10000 <- Uncertain_Slope_inter(Chaco_Mod, Chaco_Tree10000_log, 1000)

# Report CV values as indicator of uncertainty (Unit: %)
sd(Chaco_SI_10[[1]])/mean(Chaco_SI_10[[1]]) *100     
sd(Chaco_SI_100[[1]])/mean(Chaco_SI_100[[1]]) *100     
sd(Chaco_SI_1000[[1]])/mean(Chaco_SI_1000[[1]]) *100     
sd(Chaco_SI_10000[[1]])/mean(Chaco_SI_10000[[1]]) *100     




#### (2) Bootstrapped sampling
Chaco_Boot_10 <- Uncertain_Boot(Chaco_Mod, Chaco_Tree10_log, 1000)
Chaco_Boot_100 <- Uncertain_Boot(Chaco_Mod, Chaco_Tree100_log, 1000)
Chaco_Boot_1000 <- Uncertain_Boot(Chaco_Mod, Chaco_Tree1000_log, 1000)
Chaco_Boot_10000 <- Uncertain_Boot(Chaco_Mod, Chaco_Tree10000_log, 1000)

# Report CV values as indicator of uncertainty (Unit: %)
sd(Chaco_Boot_10[[1]])/mean(Chaco_Boot_10[[1]]) *100     
sd(Chaco_Boot_100[[1]])/mean(Chaco_Boot_100[[1]]) *100     
sd(Chaco_Boot_1000[[1]])/mean(Chaco_Boot_1000[[1]]) *100     
sd(Chaco_Boot_10000[[1]])/mean(Chaco_Boot_10000[[1]]) *100     



#### (3) Analytical approach
Chaco_Analy_10 <- Uncertain_Analytical(Chaco_Mod, Chaco_Tree10_log, 100)
Chaco_Analy_100 <- Uncertain_Analytical(Chaco_Mod, Chaco_Tree100_log, 100)
Chaco_Analy_1000 <- Uncertain_Analytical(Chaco_Mod, Chaco_Tree1000_log, 100)
Chaco_Analy_10000 <- Uncertain_Analytical(Chaco_Mod, Chaco_Tree10000_log, 100)

# Report CV values as indicator of uncertainty (Unit: %)
# Chaco_Analy_Number[[1]] is for CI; Chaco_Analy_Number[[2]] is for PI;
sd(Chaco_Analy_10[[1]])/mean(Chaco_Analy_10[[1]]) *100     
sd(Chaco_Analy_10[[2]])/mean(Chaco_Analy_10[[2]]) *100     

sd(Chaco_Analy_100[[1]])/mean(Chaco_Analy_100[[1]]) *100     
sd(Chaco_Analy_100[[2]])/mean(Chaco_Analy_100[[2]]) *100     

sd(Chaco_Analy_1000[[1]])/mean(Chaco_Analy_1000[[1]]) *100     
sd(Chaco_Analy_1000[[2]])/mean(Chaco_Analy_1000[[2]]) *100     

sd(Chaco_Analy_10000[[1]])/mean(Chaco_Analy_10000[[1]]) *100     
sd(Chaco_Analy_10000[[2]])/mean(Chaco_Analy_10000[[2]]) *100  


#### (4) Bayesian approach
# Calculate intercept, slope and sigma from the Bayesian approach
# The example below produce 1100000 set of [intercept, slope and sigma], discard the first 100000, and then sample 1 out of 100
Y4 <- log(Chaco_Mod[, 2])
X4 = as.matrix(data.frame(intercept=1, logagbt=log(Chaco_Mod[, 1])))
start4 = c(0,0,1)
names(start4) <- c(colnames (X4), "sd")

chain4 <- run_metropolis_MCMC(y=Y4, X=X4, startvalue = start4, iterations = 1100000, sdprop = c(0.01, 0.002, 0.01))
burnIn <- 100000
param4 <- data.frame(chain4[seq(burnIn +1, dim(chain4)[1 ], by= 100), ])


## Applying the calculated intercept, slope and sigma to tree10, tree100, tree1000, tree10000
case_study <- list(Chaco_Tree10_log, Chaco_Tree100_log, Chaco_Tree1000_log, Chaco_Tree10000_log)

out_case_study <- vector("list", length = 4)
out_CV <- vector("numeric", length = 4)

for (i in seq(1:4)){
  out_case_study[[i]] <- Uncertain_Baye(Chaco_Mod, case_study[[i]], param4)
  out_CV[[i]] <- sd(out_case_study[[i]]) / mean(out_case_study[[i]])
}

# Report CV values as indicator of uncertainty (Unit: %)
out_CV*100   



#### Applying four uncetainty methods for the Paranaense---------------------------------------------------
# output mean value for Chaco is m3/ha 

#### (1) the slope-intercept method
Para_SI_10 <- Uncertain_Slope_inter_Para(Para_Mod, Para_Tree10_log, 1000)
Para_SI_100 <- Uncertain_Slope_inter_Para(Para_Mod, Para_Tree100_log, 1000)
Para_SI_1000 <- Uncertain_Slope_inter_Para(Para_Mod, Para_Tree1000_log, 1000)
Para_SI_10000 <- Uncertain_Slope_inter_Para(Para_Mod, Para_Tree10000_log, 1000)

# Report CV values as indicator of uncertainty (Unit: %)
sd(Para_SI_10[[1]])/mean(Para_SI_10[[1]]) *100         
sd(Para_SI_100[[1]])/mean(Para_SI_100[[1]]) *100       
sd(Para_SI_1000[[1]])/mean(Para_SI_1000[[1]]) *100     
sd(Para_SI_10000[[1]])/mean(Para_SI_10000[[1]]) *100     


#### (2) Bootstrapped sampling
Para_Boot_10 <- Uncertain_Boot_Para(Para_Mod, Para_Tree10_log, 1000)
Para_Boot_100 <- Uncertain_Boot_Para(Para_Mod, Para_Tree100_log, 1000)
Para_Boot_1000 <- Uncertain_Boot_Para(Para_Mod, Para_Tree1000_log, 1000)
Para_Boot_10000 <- Uncertain_Boot_Para(Para_Mod, Para_Tree10000_log, 1000)

# Report CV values as indicator of uncertainty (Unit: %)
sd(Para_Boot_10[[1]])/mean(Para_Boot_10[[1]]) *100     
sd(Para_Boot_100[[1]])/mean(Para_Boot_100[[1]]) *100     
sd(Para_Boot_1000[[1]])/mean(Para_Boot_1000[[1]]) *100     
sd(Para_Boot_10000[[1]])/mean(Para_Boot_10000[[1]]) *100 



#### (3) Analytical approach
Para_Analy_10 <- Uncertain_Analytical_Para(Para_Mod, Para_Tree10_log, 100)
Para_Analy_100 <- Uncertain_Analytical_Para(Para_Mod, Para_Tree100_log, 100)
Para_Analy_1000 <- Uncertain_Analytical_Para(Para_Mod, Para_Tree1000_log, 100)
Para_Analy_10000 <- Uncertain_Analytical_Para(Para_Mod, Para_Tree10000_log, 100)

# Report CV values as indicator of uncertainty (Unit: %)
# Para_Analy_Number[[1]] is for CI; Para_Analy_Number[[2]] is for PI;
sd(Para_Analy_10[[1]])/mean(Para_Analy_10[[1]]) *100     
sd(Para_Analy_10[[2]])/mean(Para_Analy_10[[2]]) *100     

sd(Para_Analy_100[[1]])/mean(Para_Analy_100[[1]]) *100     
sd(Para_Analy_100[[2]])/mean(Para_Analy_100[[2]]) *100     

sd(Para_Analy_1000[[1]])/mean(Para_Analy_1000[[1]]) *100     
sd(Para_Analy_1000[[2]])/mean(Para_Analy_1000[[2]]) *100    


sd(Para_Analy_10000[[1]])/mean(Para_Analy_10000[[1]]) *100     
sd(Para_Analy_10000[[2]])/mean(Para_Analy_10000[[2]]) *100  



#### (4) Bayesian approach
# Calculate intercept, slope and sigma from the Bayesian approach
# The example below produce 1100000 set of [intercept, slope and sigma], discard the first 100000, and then sample 1 out of 100
Y4 <- Para_Mod[, 3]
X_Var <- Para_Mod[, 1]^2 * Para_Mod[, 2]
var4 = data.frame(intercept=1, logagbt=X_Var)
X4 = as.matrix(var4)
start4 = c(0,0,1)
names(start4) <- c(colnames (X4), "sd")

chain4 <- run_metropolis_MCMC(y=Y4, X=X4, startvalue = start4, iterations = 1100000, sdprop = c(0.01, 0.002, 0.01))
burnIn <- 100000
param4 <- data.frame(chain4[seq(burnIn +1, dim(chain4)[1 ], by= 100), ])


## Applying the calculated intercept, slope and sigma to tree10, tree100, tree1000, tree10000
case_study <- list(Para_Tree10_log, Para_Tree100_log, Para_Tree1000_log, Para_Tree10000_log)

out_case_study <- vector("list", length = 4)
out_CV <- vector("numeric", length = 4)

for (i in seq(1:4)){
  out_case_study[[i]] <- Uncertain_Baye_Para(case_study[[i]], param4)
  out_CV[[i]] <- sd(out_case_study[[i]]) / mean(out_case_study[[i]])
}

# Report CV values as indicator of uncertainty (Unit: %)
out_CV*100   


