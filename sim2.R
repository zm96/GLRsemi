library(MASS)

Ker <- function(u){
  v <- numeric(length(u))
  if(abs(u)<=1){
    .75 * (1 - u^2)
  }else{
    0
  }
}

p <- 2
q <- 3
n <- 100
m <- 1000# 重复次数
h <- .166
#h <- .25
#h <- .375
theta <- 0.5

Plmest <- function(p, q, theta, n, m, h){
  sigma <- matrix(c(1, 2/3, 2/3, 2/3, 2/3, 1, 2/3, 2/3, 2/3, 2/3, 1, 2/3, 2/3, 2/3, 2/3, 1), 4, 4)
  # theta = 0
  #A <- matrix(c(0, 0, 1, 0, 0, 1), nrow=2, ncol=3)
  # theta <- .5
  A <- matrix(c(0, 0, 1, 0, 0, -2), nrow=2, ncol=3)
  
  beta <- matrix(nrow=m, ncol=q)
  PLR <- matrix(nrow=m, ncol=1)
  Wald <- matrix(nrow=m, ncol=1)
  
  set.seed(123)
  for (j in 1:m){
    u <- runif(n)     
    #Simulate from a Multivariate Normal Distribution
    data <- mvrnorm(n, rep(0, 4), sigma)
    x1 <- data[, 1]
    x2 <- data[, 2]
    z1 <- data[, 3]
    z2 <- data[, 4]
    z3 <- rbinom(n, 1, .4)
    X <- cbind(x1, x2) 
    Z <- cbind(z1, z2, z3)
    error <- rnorm(n)
    # error <- 2/3 * rnorm(n, 0, 1/2) + 1/3 * rnorm(n, 0, 2)
    Y <- sin(6 * pi * u) * x1 + sin(2 * pi * u) * x2 + 2 * z1 + 2 * theta * z2 + theta * z3 + error 
    dim(Y) <- c(n, 1)
    S <- matrix(nrow=n,ncol=n)
    for (i in 1:n){
      xi <- X[i,]
      dim(xi) <- c(1, p)
      ze <- rep(0, p)
      dim(ze) <- c(1, p)
      x_tem = cbind(xi, ze)
      w <- (u - rep(u[i], n)) / h
      dim(w) <- c(n, 1)
      W <- apply(w, 1, Ker)
      W <- diag(W) / h
      D <- cbind(X, X*(w%*%rep(1,p)))
      S[i,] <- x_tem %*% ginv(t(D) %*% W %*% D) %*% t(D) %*% W
    }
    I <- diag(rep(n))
    Y_hat <- (I - S) %*% Y
    Z_hat <- (I - S) %*% Z
    beta_hat <- ginv(t(Z_hat) %*% Z_hat) %*% t(Z_hat) %*% Y_hat
    M_hat <- S %*% (Y - Z %*% beta_hat)
    RSS1 <- sum((Y - M_hat - Z %*% beta_hat)^2)
    # ginv广义逆
    beta_hat0 <-beta_hat - ginv(t(Z_hat) %*% Z_hat) %*% t(A) %*% ginv(A %*% ginv(t(Z_hat) %*% Z_hat) %*% t(A)) %*% A %*% beta_hat 
    M_hat0 <- S %*% (Y - Z %*% beta_hat0)
    RSS0 <- sum((Y - M_hat0 - Z %*% beta_hat0)^2)
    T_n <- n / 2 * (RSS0 - RSS1) / RSS1 
    sigmah = RSS1 /n * ginv(t(Z_hat) %*% Z_hat) %*% t(Z_hat) %*% (I - S) %*% t(I - S) %*% Z_hat %*% ginv(t(Z_hat) %*% Z_hat)
    wh <- t(beta_hat) %*% t(A) %*% ginv(A %*% sigmah %*% t(A)) %*% A %*% beta_hat
    beta[j,] <- t(beta_hat)
    PLR[j,] <- 2 * T_n
    Wald[j,] <- wh
    res <- cbind(beta, PLR, Wald)
  }
  return(res)
}

res <- Plmest(p, q, theta, n, m, h)
beta <- res[,1:3]
PLR <- res[,4]
Wald <- res[,5]
apply(beta, 2, mean)
apply(beta, 2, sd)

## Q-Q plot for Chi^2 data against true theoretical distribution:
qqplot(qchisq(ppoints(1000), df = 2), PLR,
       main = expression("Q-Q plot for" ~~ {chi^2}[nu == 2]))
qqplot(qchisq(ppoints(1000), df = 2), Wald,
       main = expression("Q-Q plot for" ~~ {chi^2}[nu == 2]))

# power
t <- (seq(from=0, to=30, by=1)) / 100
pw <- matrix(nrow=length(t), ncol=2)
al <- .05
for (k in 1:length(t)){
  res <- Plmest(p, q, t[k], n, m, h)
  count1 <- 0
  count2 <- 0
  for (l in 1:m){
    ans1 <- pchisq(res[l, 4], 2, 0, FALSE)
    if (ans1 < al){count1 = count1 +  1}
    ans2 <- pchisq(res[l, 5], 2, 0, FALSE)
    if (ans2 < al){count2 = count2 +  1}
  }
  pw[k, 1] <- count1 / m
  pw[k, 2] <- count2 / m
  print(k)
}
PLRtest <- pw[,1]
Waldtest <- pw[,2]
plot(t, PLRtest, type="l")
plot(t, Waldtest, type="l")