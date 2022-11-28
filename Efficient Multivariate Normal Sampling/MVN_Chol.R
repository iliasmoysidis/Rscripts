rMVN_Chol = function(n, mu, Sigma) {
  
  p = dim(Sigma)[1]
  
  # Sigma = L%*%t(L)
  
  L = t(chol(Sigma))
  Z = matrix(rnorm(n*p), nrow = n, ncol = p)
  
  # Z~N(0,I) => X=mu+LZ~N(mu,Sigma) 
  
  X = mu+L%*%t(Z)
  
  return(X)
}