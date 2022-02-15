# QR decomposition algorithm

proj = function(a,u) {
  
  e = as.numeric(t(u)%*%a)
  f = norm(u, "2")^2
  
  projection = (e/f)*u
  
  return(projection)
}

QR_gsdec = function(A) {
  
  n = dim(A)[1]
  p = dim(A)[2]
  
  U = matrix(nrow = n, ncol = p)
  Q = matrix(nrow = n, ncol = p)
  
  U[,1] = A[,1]
  Q[,1] = U[,1]/norm(U[,1], "2")
  
  for (j in 2:p) {
    U[,j] = A[,j]
    for (i in 1:(j-1)){
      U[,j] = U[,j] - proj(A[,j], U[,i])
    }
    Q[,j] = U[,j]/norm(U[,j], "2")
  }
  
  R = t(Q)%*%A
  
  
  return(list("Q" = Q, "R" = R))
}

# Wikipedia example

A = matrix(c(12,6,-4,-51,167,24,4,-68,-41), nrow =3, ncol=3)

temp = QR_gsdec(A)

Q = temp$Q
R = temp$R
A-Q%*%R
t(Q)%*%Q


# Random example

n = 5
p = 7

A = matrix(rnorm(n*p), nrow = n, ncol = p)

temp = QR_gsdec(t(A))

Q = temp$Q
R = temp$R
t(A)-Q%*%R
t(Q)%*%Q