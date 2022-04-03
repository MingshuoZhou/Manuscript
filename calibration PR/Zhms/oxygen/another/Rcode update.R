#save(exp_data1, file = "exp_data2.RData")
library(HetCalibrate)
f.sim <- function(x, cpara) {
  
  if(is.null(dim(x))){
    x1 <- x[1]
    x2 <- x[2]
  }else{
    x1 <- x[,1]
    x2 <- x[,2]
  }  
  
  out <- cpara[1] + x1*cpara[2] + x2*cpara[3]
  return(c(out))
}
df.sim <- function(x, cpara) {
  if(is.null(dim(x))){
    return(c(1,x[,1],x[,2]))
  }else{
    return(cbind(1,x))
  }
}

# setting for lower and upper bounds of parameters
cpara_min <- rep(-10, 3)  
cpara_max <- rep(10, 3)
cpara_init.mx <- as.matrix(expand.grid(c(-2.5,0,2.5),c(-2.5,0,2.5),c(-2.5,0,2.5)))
colnames(cpara_init.mx) <- NULLl

# load experimental data
load("exp_data2.RData")
X <- as.matrix(exp_data1[,1:2])
Z <- exp_data1[,3]

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
xgrid <- cbind(seq(min(X[,1]), max(X[,1]), length.out = 82), P = rep(0.4717, 82))
predictions <- predict(x = xgrid, object =  model)

## output mean predictive surface
out <- predictions$mean 
write.csv(out, file = "predictmean.csv", row.names = FALSE)
## output upper bound
#out <- qnorm(0.975, predictions$mean, sqrt(predictions$sd2 + predictions$nugs))
#write.csv(out, file = "upperbound.csv", row.names = FALSE)
## output upper bound
#out <- qnorm(0.025, predictions$mean, sqrt(predictions$sd2 + predictions$nugs))
#write.csv(out, file = "lowerbound.csv", row.names = FALSE)
## output model 
#out <- f.sim(xgrid, cpara.est)
#write.csv(out, file = "modeloutput.csv", row.names = FALSE)