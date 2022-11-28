jacob = function(beta, x){
  n = length(x)
  p = length(beta)
  J = matrix(nrow = n, ncol = p)
  
  for (i in 1:n){
    a = -x[i]/(beta[2]+x[i])
    b = beta[1]*x[i]/(beta[2]+x[i])^2
    J[i, ] = c(a, b)
  }
  
  return(J)
}

model_function = function(beta, x){
  n = length(x)
  A = c()
  
  for (i in 1:n){
    A[i] = beta[1]*x[i]/(beta[2]+x[i])
  }
  
  return(A)
}





n = 50
p = 2
beta_true = c(0.6, 1.5)




e = rnorm(n, mean = 0, sd = 0.1)
x = rnorm(n, mean = 0, sd = 1)
y = model_function(beta_true, x)+e





beta = c(0.5, 1.4) # The algorithm depends heavily on the initial values... :_(
Delta = c()
Theta = c()


for (i in 1:50){
  beta_old = beta
  
  J = jacob(beta, x)
  r = y - model_function(beta, x)
  beta = beta - solve(t(J) %*% J) %*% t(J) %*% r
  Delta[i] = norm(x = beta - beta_old, type = "2")
  Theta[i] = norm(x = y - model_function(beta, x), type = "2")
}

beta
plot(Delta)
plot(Theta)
