rMVN_BM_Chol = function(n, mu, Sigma){
  
  # Determine if the dimension is even
  # because Box-Muller produces only 
  # pairs of normal variables
  
  p = length(mu)
  a = p %% 2
  s = p^(1-a)*(p+1)^a/2
  
  
  # Create an s-dimensional standard normal
  
  C = list()
  for (i in 1:s) {
    U      = matrix(runif(2*n), nrow = n, ncol = 2)
    Z1     = sqrt(-2*log(U[,1]))*cos(2*pi*U[,2])
    Z2     = sqrt(-2*log(U[,1]))*sin(2*pi*U[,2])
    C[[i]] = matrix(c(Z1,Z2), nrow = n, ncol = 2)
  }
  
  Z = do.call(cbind, C)
  if(p %% 2 == 1) {
    Z = Z[,-1]
  }
  
  
  # Z~N(0,I) => X=mu+LZ~N(mu,Sigma)
  
  L = t(chol(Sigma))
  X = mu+L%*%t(Z)
  
  return(X)
}