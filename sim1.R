library(MASS)

Ker <- function(u){
  v <- numeric(length(u))
  if(abs(u)<=1){
    .75 * (1 - u^2)
  }else{
    0
  }
}

p <- 1
q <- 3
n <- 100
m <- 1000 # 重复次数
h <- 0.4167
#h <- 0.625
#h <- 0.9375
theta <- 0.5

Plmest <- function(p, q, theta, n, m, h){
  sigma <- matrix(c(1, .5, .5, .5, 1, .5, .5, .5, 1), 3, 3)
  # theta = 0
  # A <- matrix(c(0, 0, 1, 0, 0, 1), nrow=2, ncol=3)
  # theta <- .5
  A <- matrix(c(0, 0, 1, 0, 0, -2), nrow=2, ncol=3)

  beta <- matrix(nrow=m, ncol=q)
  PLR <- matrix(nrow=m, ncol=1)
  Wald <- matrix(nrow=m, ncol=1)
  
  set.seed(123)
  for (j in 1:m){
    #Simulate from a Multivariate Normal Distribution
    data <- mvrnorm(n, rep(0, 3), sigma)
    index <- abs(data[,1])<1.645
    x1 <- data[index, 2]
    x2 <- data[index, 3]
    x3 <- rbinom(n, 1, .4)[index]
    u <- data[index,1]
    num <- length(u)
    x <- rep(1, num)
    dim(x) <- c(num, p)
    X <- cbind(x1, x2, x3) # z
    error <- rnorm(n)[index]
    # error <- 2/3 * rnorm(n, 0, 1/2)[index] + 1/3 * rnorm(n, 0, 2)[index]
    Y <- sin(2*u) + 2 * x1 + 2 * theta * x2 + theta * x3 + error 
    dim(Y) <- c(num, 1)
    S <- matrix(nrow=num,ncol=num)
    for (i in 1:num){
      x_tem = cbind(x[i,], rep(0, p))
      w <- (u - rep(u[i], num)) / h
      dim(w) <- c(num, 1)
      W <- apply(w, 1, Ker)
      W <- diag(W) / h
      D <- cbind(x, x*w)
      S[i,] <- x_tem %*% ginv(t(D) %*% W %*% D) %*% t(D) %*% W
    }
    I <- diag(rep(num))
    Y_hat <- (I - S) %*% Y
    Z_hat <- (I - S) %*% X
    beta_hat <- ginv(t(Z_hat) %*% Z_hat) %*% t(Z_hat) %*% Y_hat
    M_hat <- S %*% (Y - X %*% beta_hat)
    RSS1 <- sum((Y - M_hat - X %*% beta_hat)^2)
    # ginv广义逆
    beta_hat0 <-beta_hat - ginv(t(Z_hat) %*% Z_hat) %*% t(A) %*% ginv(A %*% ginv(t(Z_hat) %*% Z_hat) %*% t(A)) %*% A %*% beta_hat 
    M_hat0 <- S %*% (Y - X %*% beta_hat0)
    RSS0 <- sum((Y - M_hat0 - X %*% beta_hat0)^2)
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


DBEest <- function(p, q, theta, n, m){
  sigma <- matrix(c(1, .5, .5, .5, 1, .5, .5, .5, 1), 3, 3)
  beta <- matrix(nrow=m, ncol=q)
  set.seed(123)
  for (j in 1:m){
    #Simulate from a Multivariate Normal Distribution
    data <- mvrnorm(n, rep(0, 3), sigma)
    x3 <- rbinom(n, 1, .4)
    sort.index <- sort.int(data[,1], index.return = TRUE)$ix
    u <- data[sort.index,1]
    x <- rep(1, n)
    X <- cbind(data[sort.index, 2:3], x3[sort.index]) # z
    error <- rnorm(n)[sort.index]
    # error1 <- 2/3 * rnorm(n, 0, 1/2)[sort.index] + 1/3 * rnorm(n, 0, 2)[sort.index]
    Y <- sin(2*u) + 2 * X[, 1] + 2 * theta * X[, 2] + theta * X[, 3] + error 
    diffY <- diff(Y)
    diffu <- diff(u)
    diffX <- diff(X)
    # diffY <- diff(Y)[seq(1,n-1,2)]
    # diffu <- diff(u)[seq(1,n-1,2)]
    # diffX <- diff(X)[seq(1,n-1,2),]
    Z <- cbind(rep(1, n/2), diffu, diffX)
    # normal least square
    W_hat <- ginv(t(Z) %*% Z) %*% t(Z) %*% diffY
    beta[j,] <- t(W_hat[3:5])
  }
  return(beta)
}
res <- DBEest(p, q, theta, n, m)
apply(res, 2, mean)
apply(res, 2, sd)

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