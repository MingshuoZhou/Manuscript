#decane first try, alpha = v*Tr + 0.62
library(HetCalibrate)
#save(exp_data, file = "exp_data.RData")

# model: x is input, cpara is calibration parameters
f.sim <- function(x, cpara) {
  return(c(x * cpara+1.0))
}
# derivative of the model function (with repective to the parameter)
df.sim <- function(x, cpara) {
  return(c(x))
}

# setting for lower and upper bounds of parameters
cpara_min <- -20
cpara_max <- 10
# setting for the inital guess
cpara_init.vt <- c(0.3,0.5,0.6,0.7,0.8)

# load experimental data
load("exp_data1.RData")
X <- exp_data1[,1]
Z <- exp_data1[,2]

# model fit
model <- vector("list", length(cpara_init.vt))
jj <- 0
for(cpara.init in cpara_init.vt){
  jj <- jj + 1
  model[[jj]] <- mleHomCalibrate(X = X, Z = Z, cpara_max = cpara_max, cpara_min = cpara_min,
                                 init = list("cpara" = cpara.init),
                                 covtype = "Gaussian", orthogonal = TRUE, f.sim = f.sim, df.sim = df.sim)
  
}
llmax.index <- which.max(sapply(model, function(x) x$ll))
model <- model[[llmax.index]]

# calibration estimate
cpara.est <- model$cpara
cat("parameter estimate:", cpara.est, "\n")

# predictions at test data
xgrid <- seq(min(X), max(X), length.out = 1000)
xgrid <- matrix(xgrid, ncol = 1)
predictions <- predict(x = xgrid, object =  model)

## Display mean predictive surface
plot(X, Z, xlab = "Tr", ylab = "alpha")
lines(xgrid, predictions$mean, col = 'red', lwd = 2)
lines(xgrid, f.sim(xgrid, cpara.est), col = 4, lty = 2, lwd = 2)
## Display 95% prediction intervals
lines(xgrid, qnorm(0.025, predictions$mean, sqrt(predictions$sd2 + predictions$nugs)),
      col = 3, lty = 3, lwd = 2)
lines(xgrid, qnorm(0.975, predictions$mean, sqrt(predictions$sd2 + predictions$nugs)),
      col = 3, lty = 3, lwd = 2)

## Display discrepancy
plot(xgrid, predictions$mean - f.sim(xgrid, cpara.est), xlab = "Tr", ylab = "discrepancy", type = "l", col = 2, lwd = 2)
abline(h=0, lty = 2)

