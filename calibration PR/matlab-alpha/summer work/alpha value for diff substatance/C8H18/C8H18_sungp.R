#save(exp_data, file = "exp_data.RData")

library(HetCalibrate)
f.sim <- function(x, cpara) {
  out <- cpara[1] + x*cpara[2]
  return(c(out))
}
df.sim <- function(x, cpara) {
  if(is.null(dim(x))){
    return(c(1,x))
  }else{
    return(cbind(1,x))
  }
}
# setting for lower and upper bounds of parameters
cpara_min <- c(-40,-5)
cpara_max <- c(40,15)
# setting for the inital guess
cpara_init.vt <- c(-2,-5)

# load experimental data
load("exp_data.RData")
X <- exp_data[,1]
Z <- exp_data[,2]

# model fit

model <- mleHomCalibrate(X = X, Z = Z, cpara_max = cpara_max, cpara_min = cpara_min,
                         init = list("cpara" = cpara_init.vt),
                         covtype = "Gaussian", orthogonal = TRUE, f.sim = f.sim, df.sim = df.sim)


# calibration estimate
cpara.est <- model$cpara
cat("parameter estimate:", cpara.est, "\n")

# predictions at test data
xgrid <- seq(min(X), max(X), length.out = 1000)
xgrid <- matrix(xgrid, ncol = 1)
predictions <- predict(x = xgrid, object =  model)

## Display mean predictive surface
plot(X, Z)
lines(xgrid, predictions$mean, col = 'red', lwd = 2)
lines(xgrid, f.sim(xgrid, cpara.est), col = 4, lty = 2, lwd = 2)
## Display 95% prediction intervals
lines(xgrid, qnorm(0.025, predictions$mean, sqrt(predictions$sd2 + predictions$nugs)),
      col = 3, lty = 3, lwd = 2)
lines(xgrid, qnorm(0.975, predictions$mean, sqrt(predictions$sd2 + predictions$nugs)),
      col = 3, lty = 3, lwd = 2)

## Display discrepancy
plot(xgrid, predictions$mean - f.sim(xgrid, cpara.est), xlab = "day", ylab = "discrepancy", type = "l", col = 2, lwd = 2)
abline(h=0, lty = 2)
