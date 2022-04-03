library(plgp)
library(randtoolbox)
# read functions
cov_ortho_gen <- function(X1, X2 = NULL, theta, f.sim, df.sim = NULL, cpara, #theta->gamma
                          inputBounds = NULL, MC.num = 500, nugget = 1e-8){
  
  if(is.null(inputBounds)) inputBounds <- apply(X1, 2, range)#get the maximum or minimum value of x1 in rows
  
  cov.ori <- covar.sep(X1, X2, d=-1/(4*log(theta)), g=0)
  
  if(ncol(X1) == 1){
    xi <- matrix(seq(0, 1, length.out = MC.num), ncol = ncol(X1))#seq ->linespace
  }else{
    xi <- sobol(MC.num, ncol(X1))
  }
  
  xi <- xi * matrix(rep(inputBounds[2,] - inputBounds[1,], each = nrow(xi)), ncol = ncol(xi)) + 
    matrix(rep(inputBounds[1,], each = nrow(xi)), ncol = ncol(xi)) 
  
  Wm <- covar.sep(xi, d=-1/(4*log(theta)), g=0)
  wx1 <- apply(X1, 1, function(x) covar.sep(matrix(x, ncol = ncol(X1)), xi, d=-1/(4*log(theta)), g=0))
  if(is.null(X2)) {
    wx2 <- wx1
  }else{
    wx2 <- apply(X2, 1, function(x) covar.sep(matrix(x, ncol = ncol(X2)), xi, d=-1/(4*log(theta)), g=0))
  }
  
  Fm <- df.sim(xi, cpara)
  if(length(cpara) == 1) Fm <- matrix(Fm, ncol = 1)
  
  rm.index <- apply(Fm, 2, function(x) all(x==0))
  if(all(rm.index)){
    cov.ortho <- cov.ori
  }else{
    if(length(cpara) == 1){
      FWF_inverse <- 1/(t(Fm) %*% Wm %*% Fm)
    }else{
      if(any(rm.index)) {
        Fm <- Fm[,!rm.index]
      }
      FWF <- t(Fm) %*% Wm %*% Fm
      FWF_inverse <- solve(FWF + diag(nugget * diag(FWF)))
    }
    cov.ortho <- cov.ori - t(wx1) %*% Fm %*% FWF_inverse %*% t(Fm) %*% wx2
  }
  
  return(cov.ortho)
}

nl <- function(para, X, Y) 
{
  d <- ncol(X)
  theta <- para[1:d]                             
  cpara <- para[(length(para)-d):length(para)]
  n <- length(Y)
  K <- cov_ortho_gen(X, NULL, theta, f.sim, df.sim, cpara) 
  Ki <- solve(K+diag(1e-4,nrow(K)))
  ldetK <- determinant(K+diag(1e-4,nrow(K)), logarithm=TRUE)$modulus
  ll <- - (n/2)*log(t(Y-f.sim(X,cpara)) %*% Ki %*% (Y-f.sim(X,cpara))) - (1/2)*ldetK
  return(-ll)
}

predict.func <- function(out, newdata){

  d <- ncol(X)
  print(d)
  para <- out$par
  theta <- para[1:d]                             
  cpara <- para[(length(para)-d):length(para)]
  print(cpara)
  
  KX <- cov_ortho_gen(newdata, X, theta, f.sim, df.sim, cpara) 
  KXX <- cov_ortho_gen(newdata, NULL, theta, f.sim, df.sim, cpara) 
  K <- cov_ortho_gen(X, NULL, theta, f.sim, df.sim, cpara) 
  Ki <- solve(K+diag(1e-4,nrow(K)))
  mup <- f.sim(newdata,cpara) + KX %*% Ki %*% (Z - f.sim(X,cpara))
  tau2hat <- drop(t(Z-f.sim(X,cpara)) %*% Ki %*% (Z-f.sim(X,cpara)) / nrow(X))
  Sigmap <- diag(tau2hat*(KXX + diag(1e-4, nrow(newdata)) - KX %*% Ki %*% t(KX)))
  
  return(list(mean=mup, sd=sqrt(pmax(Sigmap,0)), model=f.sim(newdata,cpara)))
}

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

# read data
data.df <- read.csv("raw alpha for exp_data.csv")
X <- as.matrix(data.df[,1:2])
Z <- data.df[,3]

# estimate parameters
out <- optim(c(0.5, 0.5, 0, 0, 0), nl, method="L-BFGS-B", lower=c(1e-6,1e-6,-5,-5,-5), upper=c(1-1e-6, 1-1e-6,5,5,5), X=X, Y=Z) 

# prediction
newdata.df <- read.csv("alpha.csv")
newdata <- as.matrix(newdata.df[,1:2])

predictions <- predict.func(out, newdata)

# output predictions
output.df <- data.frame(newdata, mean=predictions$mean, true = newdata.df[,3], model=predictions$model, upper=predictions$mean+qnorm(0.975)*predictions$sd,
                        lower=predictions$mean+qnorm(0.025)*predictions$sd)

write.csv(output.df, file = "predictions.csv", row.names = FALSE)


