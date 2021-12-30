GMD = function(X, Q, R, tolerance = 1e-05, maxiter = 100){
  
  n = dim(X)[1]
  p = dim(X)[2]
  
  U = matrix(0, nrow = n, ncol = n)
  V = matrix(0, nrow = p, ncol = p)
  D = rep(0, min(n, p))
  
  X_0 = X
  
  for (j in 1:min(n, p)) {
    
    
    error = 1
    iter  = 0
    u_0   = rnorm(n)
    v_0   = rnorm(p)
    
    while (error > tolerance & iter < maxiter) {
      
      x      = X_0%*%R%*%v_0
      x_norm = as.numeric(sqrt(t(x)%*%Q%*%x))
      u      = x/x_norm
      
      y      = t(X_0)%*%Q%*%u
      y_norm = as.numeric(sqrt(t(y)%*%R%*%y))
      v      = y/y_norm
      
      error = max(norm(u-u_0, "2"), norm(v-v_0, "2"))
      iter  = iter+1
      u_0   = u
      v_0   = v
      
    }
    
    d   = as.numeric(t(u)%*%Q%*%X_0%*%R%*%v)
    X_0 = X_0-d*u%*%t(v)
    
    U[,j] = as.vector(u)
    V[,j] = as.vector(v)
    D[j]  = d
    
  }
  
  result = list("d" = D, "u" = U, "v" = V)
  return(result)
}