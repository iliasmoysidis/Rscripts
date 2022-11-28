# Number of samples
n = 10

# Number of variables
p = 5

# Standard deviation of noise
sigma = 2

# Generate the data
X    = matrix(rnorm(n*p), nrow = n, ncol = p)
beta = runif(p)

# Generate the noise
epsilon = rnorm(n, mean = 0, sd = 2)

# Generate the observed data
y = X%*%beta + epsilon

# Use the lm function
fit = lm(y ~ X + 0)

# Estimate the coefficients
beta_hat = solve(t(X)%*%X)%*%t(X)%*%y

# Residuals
res    = y - X%*%beta_hat

# Residual standard error
res_se = norm(res, "2")/sqrt(n-p)

# Estimate standard deviation
sigma_hat = (sqrt(n-p)/sqrt(n))*res_se

# Estimate coefficients' standard error
coef_se = res_se*sqrt(diag(solve(t(X)%*%X)))

# t-values
t_vals = beta_hat/coef_se

# p-values
p_vals = 2*(pt(abs(t_vals), lower.tail = F, df = n-p))




# F-statistic and its p-value
F_val = t(beta_hat)%*%t(X)%*%X%*%beta_hat/(p*res_se^2)

p_val = pf(F_val, df1 = p, df2 = n-p, lower.tail = F)