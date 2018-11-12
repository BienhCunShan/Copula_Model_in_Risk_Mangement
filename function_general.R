# to get the variance matrix we need to compute conditional t-normal distribution
trans <- function(j,q,sig,G){
  sig1 <- sig
  if(j == q){
    sig1[q,1:q] <- sig[q,1:q]#/G
    sig1[1:q,q] <- sig[1:q,q]#/G
    sig1[q,q] <- sig[q,q]#*G
    return(sig1)
  }
  for (a in 1:(q-j)) {
    for (b in 1:(q-j)) {
      sig1[a,b] <- sig[j+a,j+b]#*G
    }
  }
  for (a in 1:j) {
    for (b in 1:j) {
      sig1[q-j+a,q-j+b] <- sig[a,b]
    }
  }
  sig1[1:(q-j),(q-j+1):q] <- sig[(j+1):q,1:j]
  sig1[(q-j+1):q,1:(q-j)] <- sig[1:j,(j+1):q]
  return(sig1)#/G)
}

# to get the variance matrix we need to compute conditional t-normal distribution
a.fun <- function(z.sta,z.sim,j){
  a <- matrix(0,q)
  for (i in 1:(q-j)) {
    a[i] <- z.sim[i+j]
  }
  for (i in (q-j+1):(q)) {
    a[i] <- z.sta[i+j-q]
  }
  return(a[1:(q-1)])
}

# to compute parameters
condi.par <- function(sig,a){
  sig11 <- sig[1:(q-1),1:(q-1)]
  sig12 <- sig[1:(q-1),q]
  sig21 <- sig[q,1:(q-1)]
  sig22 <- sig[q,q]
  mu <- 0 + sig21%*%solve(sig11)%*%(a)
  sigma <- sig22 - sig21%*%solve(sig11)%*%sig12
  return(c(mu,sigma))
}

# to get each j component in Kernal
# low and up are numbers
prob <- function(j,q,sig,G,z.sta,z.sim,low,up){
  sig <- trans(j,q,sig,G)
  a <- a.fun(z.sta,z.sim,j)
  temp <- condi.par(sig,a)
  mean <- temp[1]
  sigma <- temp[2]
  p <- dtruncnorm(z.sta[j],mean = mean, sd = sqrt(sigma), a = low, b = up)
  # print(j)
  # print(z.sta[j])
  # print(mean)
  # print(sigma)
  # print(c(low,up))
  # print(p)
  return(p)
}

# to compute K
# Blow and Bup are vectors (for the same i-th)
K <- function(q,sig,G,z.sta,z.sim,Blow,Bup){
  temp.pro <- 1
  for (j in 1:q) {
    #print(j)
    #print(Blow[j])
    #print(prob(j,q,sig,G,z.sta,z.sim,Blow[j],Bup[j]))
    temp.pro <- temp.pro*prob(j,q,sig,G,z.sta,z.sim,Blow[j],Bup[j])
    #print(temp.pro)
  }
  return(temp.pro)
}

# part2.function <- function(G,q,sig,z.sim.meam,z.sim.all,BL,BU,t){
#   temp.sum <- 0
#   for (i in 1:G) {
#     temp.sum <- temp.sum + K(q,sig,G,z.sta = z.sim.meam,z.sim = z.sim[i,],Blow=BL[t,],Bup=BU[t,])
#   }
#   temp.sum <- temp.sum / G
#   return(temp.sum)
# }

all.function <- function(G,q,sig,BL,BU){
  temp.p1 <- matrix(0 ,ncol = 1, nrow = q)
  temp.p2 <- matrix(0 ,ncol = 1, nrow = q)
  for (t in 1:q) {
    B.u <- BU[t,]
    B.l <- BL[t,]
    set.seed(123)
    z.sim <-rtmvnorm(n=G, sigma=sig, upper=B.u, lower = B.l)
    z.sim.meam <- c(mean(z.sim[,1]),mean(z.sim[,2]),mean(z.sim[,3]))
    temp.p1[t] <- dmvnorm(z.sim.meam,sigma = sig)
    part.2 <- 0
    for (i in 1:G) {
      part.2 <- part.2 + K(q,sig,G,z.sta = z.sim.meam,z.sim = z.sim[i,],Blow=B.l,Bup=B.u)
    }
    temp.p2[t] <- part.2/G
  }
  est <- log(temp.p1) + log(temp.p2)
  pro <- 1
  for (t in 1:q) {
    pro <- pro * est[t]
  }
  return(pro)
}

each.function <- function(G,q,sig,BL,BU){
  temp.p1 <- matrix(0 ,ncol = 1, nrow = q)
  temp.p2 <- matrix(0 ,ncol = 1, nrow = q)
  for (t in 1:q) {
    B.u <- BU[t,]
    B.l <- BL[t,]
    set.seed(123)
    z.sim <-rtmvnorm(n=G, sigma=sig, upper=B.u, lower = B.l)
    z.sim.meam <- c(mean(z.sim[,1]),mean(z.sim[,2]),mean(z.sim[,3]))
    temp.p1[t] <- dmvnorm(z.sim.meam,sigma = sig)
    part.2 <- 0
    for (i in 1:G) {
      part.2 <- part.2 + K(q,sig,G,z.sta = z.sim.meam,z.sim = z.sim[i,],Blow=B.l,Bup=B.u)
    }
    temp.p2[t] <- part.2/G
  }
  est <- log(temp.p1) + log(temp.p2)
  return(est)
}

# function sig.fun
# to prepare for the dev and gra
sig.fun <- function(sig,i,j){
  x <- 2e-10
  sig[i,j] <- sig[i,j] + x/2
  sig[j,i] <- sig[j,i] + x/2
  return(sig)
}

# function dev
# to find derivative
# to prepare for the gra
dev.all <- function(sig,i,j,G,q,BL,BU){
  sig2 <- sig.fun(sig,i,j)
  y1 <- all.function(G,q,sig,BL,BU)
  y2 <- all.function(G,q,sig2,BL,BU)
  return((y2-y1)/x)
}

dev.each <- function(sig,i,j,G,q,BL,BU){
  sig2 <- sig.fun(sig,i,j)
  y1 <- each.function(G,q,sig,BL,BU)
  y2 <- each.function(G,q,sig2,BL,BU)
  return((y2-y1)/x)
}

# function gra
# to get gradient
gra.all <- function(sig,G,q,BL,BU){
  g <- c(0)
  for (i in 1:q) {
    for (j in (i+1):(q-1)) {
      g[length(g) + 1] <- dev.all(sig,i,j,G,q,BL,BU)
    }
  }
  return(g[,-1])
}

gra.each <- function(sig,G,n,q,BL,BU){
  temp <- list(rep(0,n))
  for (i in 1:q) {
    for (j in (i+1):(q-1)) {
      temp[length(temp) + 1] <- dev.each(sig,i,j,G,q,BL,BU)
    }
  }
  g <- temp[[-1]]
  vec <- matrix(0, nrow = (q^2-q)/2 , ncol = n)
  for (i in 1:length(g)) {
    for (j in 1:length(g[[i]])) {
      vec[i,j] <- g[[i]][j]
    } 
  }
  return(vec)
}

# function L
# when you input y, q, theta, you can get the lower boundary
L <- function(y,q,theta){
  a <- pbinom(y,size = q - 1,prob = theta)
  #print(a)
  b <- dbinom(y,size = q - 1,prob = theta)
  #print(b)
  #print(a-b)
  temp <- qnorm(abs(a - b))
  return(temp)
}

# function U
# when you input y, q, theta, you can get the upper boundary
U <- function(y,q,theta){
  a <- pbinom(y,size = q - 1,prob = theta)
  temp <- qnorm(a)
  return(temp)
}

dev.theta.all <- function(sig,G,Y,theta,i){
  n <- nrow(Y)
  q <- ncol(Y)
  BL <- Y
  BU <- Y
  BL2 <- Y
  BU2 <- Y
  theta2 <- theta
  x <- 2e-10
  theta2[i] <- theta[i] + x
  for (i in 1:n) {
    for (j in 1:q) {
      BL[i,j] <- L(Y[i,j],q,theta[j])
      BU[i,j] <- U(Y[i,j],q,theta[j])
      BL2[i,j] <- L(Y[i,j],q,theta2[j])
      BU2[i,j] <- U(Y[i,j],q,theta2[j])
    } 
  y1 <- all.function(G,q,sig,BL,BU)
  y2 <- all.function(G,q,sig,BL2,BU2)
  return((y2-y1)/x)
}

gra.theta.all <- function(sig,G,Y,theta){
  q <- ncol(Y)
  temp <- c(0)
  for (i in 1:q) {
    temp[length(temp)+1] <- dev.theta.all(sig,G,Y,theta,i)
  }
  return(temp[-1])
}

dev.theta.each <- function(sig,G,Y,theta,i){
  n <- nrow(Y)
  q <- ncol(Y)
  BL <- Y
  BU <- Y
  BL2 <- Y
  BU2 <- Y
  theta2 <- theta
  x <- 2e-10
  theta2[i] <- theta[i] + x
  for (i in 1:n) {
    for (j in 1:q) {
      BL[i,j] <- L(Y[i,j],q,theta[j])
      BU[i,j] <- U(Y[i,j],q,theta[j])
      BL2[i,j] <- L(Y[i,j],q,theta2[j])
      BU2[i,j] <- U(Y[i,j],q,theta2[j])
    } 
  y1 <- each.function(G,q,sig,BL,BU)
  y2 <- each.function(G,q,sig,BL2,BU2)
  return((y2-y1)/x)
}

gra.theta.each <- function(sig,G,Y,theta){
  q <- ncol(Y)
  n <- nrow(Y)
  temp <- list(rep(0,n))
  for (i in 1:q) {
    temp[length(temp)+1] <- dev.theta.each(sig,G,Y,theta,i)
  }
  g <- temp[[-1]]
  vec <- matrix(0, nrow = q , ncol = n)
  for (i in 1:length(g)) {
    for (j in 1:length(g[[i]])) {
      vec[i,j] <- g[[i]][j]
    } 
  }
  return(vec)
}

comb.all <- function(sig,G,Y,theta){
# gra.theta.all <- function(sig,G,Y,theta)
# gra.all <- function(sig,G,q,BL,BU)
  n <- nrow(Y)
  q <- ncol(Y)
  BL <- Y
  BU <- Y
  for (i in 1:n) {
    for (j in 1:q) {
      BL[i,j] <- L(Y[i,j],q,theta[j])
      BU[i,j] <- U(Y[i,j],q,theta[j])
    } 
  g1 <- gra.all(sig,G,q,BL,BU)
  g2 <- gra.theta.all(sig,G,Y,theta)
  g3 <- matrix(0,ncol = 1, nrow = nrow(g1) + nrow(g2))
  g3[1:nrow(g1)] <- g1
  g3[(nrow(g1)+1):nrow(g1) + nrow(g2)] <- g2
  return(g3)
}

comb.each <- function(sig,G,Y,theta){
# gra.theta.each <- function(sig,G,Y,theta)
# gra.each <- function(sig,G,n,q,BL,BU)
  n <- nrow(Y)
  q <- ncol(Y)
  BL <- Y
  BU <- Y
  for (i in 1:n) {
    for (j in 1:q) {
      BL[i,j] <- L(Y[i,j],q,theta[j])
      BU[i,j] <- U(Y[i,j],q,theta[j])
    } 
  g1 <- gra.each(sig,G,q,BL,BU)
  g2 <- gra.theta.each(sig,G,Y,theta)
  g3 <- matrix(0,ncol = n, nrow = nrow(g1) + nrow(g2))
  g3[,1:nrow(g1)] <- g1
  g3[,(nrow(g1)+1):nrow(g1) + nrow(g2)] <- g2
  return(g3)
}

compute.each <- function(sig,G,Y,theta){
  temp <- comb.each(sig,G,Y,theta)
  n <- nrow(temp)
  q <- ncol(temp)
  temp.sum <- matrix(0,n,n)
  for (i in 1:q) {
    temp.sum <- temp.sum + temp[,i]%*%t(temp[,i])
  }
  return(solve(temp.sum))
}


new.sig <- function(phi){
  # to do
  sig <- 
  return(sig)
}

new.theta <- function(phi){
  # to do
  theta <-
  return(theta)
}

new.phi <- function(theta,sig){
  # to do
  phi <-
  return(phi)
}

BHHH <- function(sig,G,Y, step){
  # n is the number of branches
  # q is the number of risk indexes
  n <- nrow(Y)
  q <- ncol(Y)
  theta <- matrix(0,nrow = 1, ncol = q)
  for (i in 1:q) {
    theta[i] <- mean(Y[,i])/q
  }
  BL <- Y
  BU <- Y
  for (i in 1:n) {
    for (j in 1:q) {
      BL[i,j] <- L(Y[i,j],q,theta[j])
      BU[i,j] <- U(Y[i,j],q,theta[j])
    }  
  start <- new.phi(theta,sig)
  new <- start
  old <- new + step * compute.each(sig,G,Y,theta) %*% comb.all(sig,G,Y,theta)
  while(abs(new - old) > e-6){
    new <- old
    sig <- new.sig(new)
    theta <- new.theta(new)
    old <- new + step * compute.each(sig,G,Y,theta) %*% comb.all(sig,G,Y,theta)
  }
  return(new)
}




