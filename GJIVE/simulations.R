library(pracma); library(Matrix); library(fields)
library(expm); library(matrixcalc)

euc_dist = function(x){
  
  # Each component of this matrix expresses its euclidian distance from the diagonal
  abs(matrix((1:x)-1, nrow = x, ncol = x, byrow = TRUE)-((1:x)-1))
  
}
nor_dist = function(x){
  
  euc_dist(x)/(1+euc_dist(x))
  
}
tax_dist = function(x){
  # The taxicab distance is the maximum distance from the diagonal
  # d(i,j) = max(i,j) if i != j and 0 if i = j
  
  M = matrix(0, nrow = x, ncol = x)
  
  # Set the components of the upper diagonal equal to their
  # column coordinate
  coordinates    = which(row(M)<col(M))
  M[coordinates] = col(M)[coordinates]
  
  # Because the coordinates are symmetrics, so is their
  # distance from the diagonal for the components of 
  # the lower diagonal
  M              = M+t(M)
  
  return(M)
}
log_dist = function(x){
  
  log(1+euc_dist(x))
  
}
row_cov  = function(n, rho){
  
  rho^euc_dist(n)
  
}
col_cov  = function(p, rho){
  
  K = length(p)
  M = list()
  
  for(k in 1:K){
    
    if (k%%3 == 0) {
      
      M[[k]] = rho^nor_dist(p[k])
      
    } else if (k%%3 == 1) {
      
      M[[k]] = rho^log_dist(p[k])
      
    } else {
      
      M[[k]] = rho^(tax_dist(p[k])/p[k])
      
    }
    
  }
  
  return(M)
}

conc_list    = function(A, B, C){
  
  K = length(A)
  D = c()
  
  for(k in 1:K){
    D = cbind(D, A[[k]]%*%B[[k]]%*%t(C[[k]]))
  }
  
  return(D)
}
obs_data_gen = function(Sigma, Delta, rank_joint, rank_indiv){
  
  n = dim(Sigma)[1]
  K = length(Delta)
  p = sapply(Delta, function(x) dim(x)[1])
  
  sqrt_Sigma = sqrtm(Sigma)
  sqrt_Delta = lapply(Delta, sqrtm)
  
  # Create a basis from which I will randomly draw
  # vectors in order to make joint and individual
  # structure orthogonal, as well as individual
  # structure orthogonal among the categories
  
  Base    = randortho(n, type = "orthonormal")
  
  U_joint = sqrt_Sigma%*%Base[ ,1:rank_joint]
  V_joint = bdiag(sqrt_Delta)%*%randortho(sum(p), type = "orthonormal")[ ,1:rank_joint]
  
  eigenvals = sort(sample(200:250, size = rank_joint), decreasing = T)
  D_joint   = diag(eigenvals, nrow = rank_joint)
  
  
  # Create individual structure
  
  U_indiv = list()
  D_indiv = list()
  V_indiv = list()
  
  # The coordinates that will give the number
  # of vectors from the common base to create
  # the individual structure
  
  a = cumsum(c(0,rank_indiv[-K]))+rank_joint+1
  b = cumsum(rank_indiv)+rank_joint
  
  for(k in 1:K){
    
    U_indiv[[k]] = sqrt_Sigma%*%Base[ ,a[k]:b[k]]
    V_indiv[[k]] = sqrt_Delta[[k]]%*%randortho(p[k], type = "orthonormal")[ ,1:rank_indiv[k]]
    
    eigenvals    = sort(sample(200:250, size = rank_indiv[k]), decreasing = T)
    D_indiv[[k]] = diag(eigenvals, nrow = rank_indiv[k])
    
  }
  
  
  
  # Generating the observed data
  
  # Matrix of joint structure
  J = U_joint%*%D_joint%*%t(V_joint)
  
  # Matrix of individual structure
  A = conc_list(U_indiv, D_indiv, V_indiv)
  
  # Noise
  E = sqrt_Sigma%*%matrix(rnorm(n*sum(p), mean = 0, sd = 1), nrow = n, ncol = sum(p))%*%bdiag(sqrt_Delta)
  
  # Observed data
  X = J+A+E
  
  
  
  
  result = list("U_joint"  = U_joint,
                "D_joint"  = D_joint,
                "V_joint"  = V_joint,
                "U_indiv"  = U_indiv,
                "D_indiv"  = D_indiv,
                "V_indiv"  = V_indiv,
                "Obs_data" = X,
                "Noise"    = E)
  return(result)
}


QR_norm = function(X, Q, R){
  sqrt(sum(diag(Q%*%X%*%R%*%t(X))))
}
GMD     = function(X, Q, R, r, tolerance = 1e-05, maxiter = 1000){
  
  n = dim(X)[1]
  p = dim(X)[2]
  
  U = matrix(0, nrow = n, ncol = r)
  V = matrix(0, nrow = p, ncol = r)
  D = matrix(0, nrow = r, ncol = r)
  
  X_0 = X
  
  for (j in 1:r) {
    
    
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
      u_0 = u
      v_0 = v
      
    }
    
    d   = as.numeric(t(u)%*%Q%*%X_0%*%R%*%v)
    X_0 = X_0-d*u%*%t(v)
    
    U[,j]  = as.vector(u)
    V[,j]  = as.vector(v)
    D[j,j] = d
    
  }
  
  result = list("u" = U, "v" = V, "d" = D)
  return(result)
}
deconc  = function(X, p){
  
  K           = length(p)
  par_sums    = c(0, cumsum(p))
  coordinates = cbind(par_sums[-(K+1)]+1, par_sums[-1])
  Y           = lapply(1:K, function(x) X[, coordinates[x,1]:coordinates[x,2]])
  return(Y)
  
}
GJIVE   = function(X, Q, R, r, s,tolerance = 1e-05, maxiter = 100){
  
  K = length(R)
  n = dim(Q)[1]
  p = sapply(R, function(x) dim(x)[1])
  
  J_0   = matrix(0, nrow = n, ncol = sum(p))
  A_0   = matrix(0, nrow = n, ncol = sum(p))
  error = 1
  iter  = 0
  
  H = list()
  G = list()
  W = list()
  
  a = c()
  b = c()
  f = c()
  while (error > tolerance & iter < maxiter) {
    
    fit = GMD(X-A_0, Q, bdiag(R), r)
    U   = fit$u
    D   = fit$d
    V   = fit$v
    
    J = U%*%D%*%t(V)
    Y = deconc(X-J, p)
    P = diag(nrow = n)-U%*%t(U)%*%Q
    
    for (i in 1:K) {
      
      P    = diag(nrow = n)-U%*%t(U)%*%Q
      temp = P%*%Y[[i]]
      
      
      fit    = GMD(temp, Q, R[[i]], s[i])
      H[[i]] = fit$u
      G[[i]] = fit$d
      W[[i]] = fit$v
      
    }
    
    
    A = conc_list(H, G, W)
    
    a = c(a, QR_norm(J-J_0, Q, bdiag(R)))
    b = c(b, sqrt(sum(mapply(function(a,b) QR_norm(a, Q, b)^2, deconc(A-A_0, p), R))))
    f = c(f, sqrt(sum(mapply(function(a,b) QR_norm(a, Q, b)^2, deconc(X-J-A, p), R))))
    
    error = max(QR_norm(J-J_0, Q, bdiag(R)), sqrt(sum(mapply(function(a,b) QR_norm(a, Q, b)^2, deconc(A-A_0, p), R))))
    iter  = iter+1
    J_0   = J
    A_0   = A
    
  }
  
  result = list("U_joint" = U, "D_joint" = D, "V_joint" = V,
                "U_indiv" = H, "D_indiv" = G, "V_indiv" = W,
                "QR error" = f)
  return(result)
}