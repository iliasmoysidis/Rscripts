# In this script we implement the coordinate descent
# method to find the solution of the lasso for linear regression


soft_thres = function(k, a) {
  max(0, 1-k/abs(a))*a
}


# Number of variables
p      = 10

# Number of samples
n      = 5

# Penalizing parameter
lambda = 1.3


beta   = rnorm(p)
beta_0 = rep(10, p)
X      = matrix(rnorm(n*p), nrow = n, ncol = p)
y      = rnorm(n)


tolerance = 1e-05
maxiter   = 100
error     = 1
iter      = 0


while (error > tolerance & iter < maxiter) {
  
  for (i in 1:p) {
    beta[i] = soft_thres(lambda/norm(X[,i], "2"), t(X[,i])%*%(y-X[,-i]%*%beta[-i])/t(X[,i])%*%X[,i])
  }
  
  error  = norm(beta-beta_0, "2")
  beta_0 = beta
  iter   = iter + 1
}