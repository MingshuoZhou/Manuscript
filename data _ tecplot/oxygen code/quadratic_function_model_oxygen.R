#save(exp_data1, file = "exp_data2.RData")
library(HetCalibrate)
f.sim <- function(x, cpara) {
  out <- cpara[1]*(x-1)^2 + cpara[2]
  return(c(out))
}
df.sim <- function(x, cpara) {
  if(is.null(dim(x))){
    return(c((x-1)^2,1))
  }else{
    return(cbind((x-1)^2,1))
  }
}
# setting for lower and upper bounds of parameters
cpara_min <- c(-10,-10)
cpara_max <- c(10,10)
# setting for the inital guess
cpara_init.mx <- as.matrix(expand.grid(c(-1,-0.8,-0.5,-0.3,-0.1),c(0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)))
colnames(cpara_init.mx) <- NULL

# load experimental data
load("exp_data2.RData")
X <- exp_data1[,1]
Z <- exp_data1[,2]

model <- vector("list",nrow(cpara_init.mx))
# model fit
for (jj in 1:nrow(cpara_init.mx)){
  cpara_init.vt <- cpara_init.mx[jj,]
  model[[jj]] <- mleHomCalibrate(X = X, Z = Z, cpara_max = cpara_max, cpara_min = cpara_min,
                                 init = list("cpara" = cpara_init.vt),
                                 covtype = "Gaussian", orthogonal = TRUE, f.sim = f.sim, df.sim = df.sim)
}
llmax.index<-which.max(sapply(model,function(x) x$ll))
model <-model[[llmax.index]]
print(model$cpara)


# calibration estimate
cpara.est <- model$cpara
cat("parameter estimate:", cpara.est, "\n")

# predictions at test data
xgrid <- seq(min(X), max(X), length.out = 1000)
xgrid <- matrix(xgrid, ncol = 1)
predictions <- predict(x = xgrid, object =  model)

## Display mean predictive surface
plot(X, Z,xlab = "Tr",ylab = "alpha")
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
