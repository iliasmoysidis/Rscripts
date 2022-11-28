################################
library(magic)

QR_hdec = function(A) {
  
  n = dim(A)[1]
  p = dim(A)[2]
  
  e = c(1, rep(0, n-1))
  Q = diag(nrow = n)
  t = min(n-1, p)
  X = A
  
  for (i in 1:t) {
    
    x = X[,1]
    a = -sign(x[1])*norm(x, "2")
    e = c(1, rep(0 , length(x)-1))
    u = x-a*e
    v = u/norm(u, "2")
    Q = Q%*%adiag(diag(nrow = i-1), diag(nrow = length(v))-2*v%*%t(v))
    
    X = as.matrix((Q%*%A)[-c(1:i),-c(1:i)])
  }
  
  R = t(Q)%*%A
  
  return(list("Q" = Q, "R" = R))
  
}



# Wikipedia example

A = matrix(c(12,6,-4,-51,167,24,4,-68,-41), nrow =3, ncol=3)

temp = QR_hdec(A)
Q = temp$Q
R = temp$R
A-Q%*%R
t(Q)%*%Q


# Random example

n = 3
p = 5

A = matrix(rnorm(n*p), nrow = n, ncol = p)

temp = QR_hdec(t(A))

Q = temp$Q
R = temp$R
t(A)-Q%*%R
t(Q)%*%Q