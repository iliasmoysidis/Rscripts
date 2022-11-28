box_muller = function(n){
  
  # Generate n samples of two independent
  # uniform (0,1) variables
  
  U = matrix(runif(2*n), nrow = n, ncol = 2)
  
  # Use the Box-Muller transformation to
  # get n samples of two independent standard
  # normal random variables
  
  Z1 = sqrt(-2*log(U[,1]))*cos(2*pi*U[,2])
  Z2 = sqrt(-2*log(U[,1]))*sin(2*pi*U[,2])
  
  Z = matrix(c(Z1,Z2), nrow = n, ncol = 2)
  
  return(Z)
}

standard_normal_pd = function(n, p){
  
  # Determine if the dimension is even
  # because Box-Muller produces only 
  # pairs of normal variables
  
  a = p %% 2
  s = p^(1-a)*(p+1)^a/2
  
  
  # Using the Box-Muller function to sample
  # s-pairs of independent standard normals
  
  C = lapply(as.list(rep(n, s)), box_muller)
  Z = do.call(cbind, C)
  
  # If p is odd we have to get rid onle of 
  # the variables, because Box-Muller generates
  # only pairs
  
  if(p %% 2 == 1){
    Z = Z[,-1]
  }
  
  return(Z)
}