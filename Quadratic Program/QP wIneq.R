############################################################################
##########Quadratic Program with only inequality constraints################
############################################################################

obj_const = function(A, b, x) {
  n = dim(A)[1]
  
  c = rep(0, n)
  for (i in 1:n){
    c[i] = b[i]-A[i, ] %*% x
  }
  
  return(c)
}

#####################################################################
######################Problem Construction###########################
#####################################################################


# minimize 1/2 x'Hx + f'x
# subject to Ax<=b

p = 50 # number of variables
n = 100 # number of inequalities

# p needs to be less than n for this method to work

# H
X = matrix(rnorm(n = p^2, mean = 0, sd = 1), nrow = p, ncol = p)
H = X%*%t(X)

# f
f = rnorm(n = p, mean = 0, sd = 1)

# A,b
A = matrix(rnorm(n = n*p, mean = 0, sd = 1), nrow = n, ncol = p)
b = rnorm(n = n, mean = 0, sd = 1)

###################################################################
########################Barrier Method#############################
###################################################################




max_iter = 30
mu = 1
x = solve(t(A) %*% A) %*% t(A) %*% b # initial value of x must be inside the feasible set
lambda = rexp(n, rate = 1)
Delta = c()
Theta = c()

for (i in 1:max_iter) {
  Theta[i] = 0.5 * t(x) %*% H %*% x + t(f) %*% x
  x_old = x
  
  # construct the jacobian
  c = obj_const(A, b, x)
  J = rbind(cbind(H, t(A)), cbind(-diag(lambda) %*% A, diag(c)))
  
  # construct gradient
  G = rbind(H %*% x + f + t(A) %*% lambda, -diag(lambda) %*% c - mu * rep(1, n))
  
  # construct the update
  J_inv = solve(J)
  a = 1 # step-size
  temp = c(x, lambda)-a * J_inv %*% G
  
  # make sure the dual variable has only positive values
  while (sum(temp[(p+1):(p+n)] > 0) < n) {
    a = 0.5 * a
    temp = c(x, lambda)-a * J_inv %*% G
  }
  
  x = temp[1:p]
  lambda = temp[(p + 1):(p + n)]
  
  mu = 0.9 * mu # reduce the penalty at each step
  
  Delta[i] = norm(x = x_old - x, type = "2")
}

plot(Delta)
plot(Theta)