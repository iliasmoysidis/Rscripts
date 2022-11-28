model_function = function(x, beta){
  n = length(x)
  p = length(beta)
  
  h = rep(0, n)
  for (i in 1:n){
    for (j in 1:p){
      h[i] = h[i] + cos(j*beta[j]*x[i])
    }
  }
  
  return(h)
}


jacobian = function(x, beta){
  n = length(x)
  p = length(beta)
  
  J = matrix(nrow = n, ncol = p)
  for (i in 1:n){
    for (j in 1:p){
      J[i, j] = -i*x[i]*sin(i*beta[j]*x[i])
    }
  }
  
  return(J)
}



n = 20 # number of samples
p = 30 # number of variables


e = rnorm(n, mean = 0, sd = 0.01) # noise
x = rnorm(n, mean = 0, sd = 1) # samples
beta_true = rnorm(p, mean = 0, sd = 0.5) # covariates
z = model_function(x, beta_true) + e # observations


beta = beta_true+rnorm(p, mean = 0, sd = 0.1)
lambda = 0.000000002
max_iter = 300
Delta = c()



for (i in 1:max_iter){
  beta_old = beta
  
  J = jacobian(x, beta)
  y = z - model_function(x, beta) + J %*% beta
  beta = solve(t(J) %*% J + lambda * diag(1, nrow = p)) %*% t(J) %*% y
  
  Delta[i] = norm(x = beta - beta_old, type = "2")
}

plot(Delta)