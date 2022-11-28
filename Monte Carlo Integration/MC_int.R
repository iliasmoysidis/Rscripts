# We will write a function for integrating e^-x
# from a to b using only sampling from a uniform
# (0,1) distribution and the law of large numbers

MC_int = function(n, a, b) {
  
  u = runif(n)
  y = (b-a)*u+a
  
  int = (b-a)*mean(exp(-y))
  
  return(int)
}


int = MC_int(n = 10000, a = 2, b = 3)

# Compare with the actual value of the integral

abs(int - exp(-a) + exp(-b))