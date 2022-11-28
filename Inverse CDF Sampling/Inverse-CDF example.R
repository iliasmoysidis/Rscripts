sigma = sqrt(2)
n = 10000

true_dens = function(x, sigma) {
  (x/sigma^2)*exp(-(x^2)/(2*sigma^2))
}

sim_samples = function(n, sigma) {
  u = runif(n)
  x = sqrt(-2*(sigma^2)*log(1-u))
  return(x)
}

y = sim_samples(n, sigma)
hist(y, prob = T)
curve(true_dens(x, sigma), from = min(y), to = max(y), col = "red", add = TRUE)