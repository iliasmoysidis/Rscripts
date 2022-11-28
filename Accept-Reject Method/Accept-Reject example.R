# Use an exponential distribution with rate 
# 1 to sample from a standard normal

# First, we will use the accept-reject
# method to sample from a folded standard
# normal with the proposal exponential.
# The folded standard normal has density
# f(x) = 2/sqrt(2*pi)*exp(-x^2/2)

# sup_(x>=0)f(x)/g(x) = sqrt(2*exp(1)/pi)

f = function(x) {
  2/sqrt(2*pi)*exp(-x^2/2)
}
g = function(x) {
  exp(-x)
}

d = 1.32



n_samples = 1000
n_run     = 0
iter      = 0
x         = c()

while (n_run < n_samples) {
  
  y = rexp(1)
  u = runif(1)
  
  if (u <= f(y)/(d*g(y))) {
    x     = c(x, y)
    n_run = n_run + 1
  }
  
  iter = iter + 1
}

# If Z follows a folded normal, then
# W = Z with probability 0.5 and 
# W = -Z with probability 0.5. follows
# a standard normal. We use a simple,
# and intuitive method to sample from
# a bernouli with prob 0.5.

n_samples = 1000
iter      = 0
s         = c()

while (iter < n_samples) {
  
  u = runif(1)
  
  if (u <= 0.5) {
    s     = c(s, 1)
    
  } else {
    s = c(s,-1)
  }
  iter = iter + 1
}


# Finally, the standard normal samples are 
# given by

Z = s*x
hist(Z)

