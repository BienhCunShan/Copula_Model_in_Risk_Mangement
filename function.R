# function L
# when you input y, q, theta, you can get the lower boundary
L <- function(y,theta){
  a <- pbinom(y,size = 2,prob = theta)
  #print(a)
  b <- dbinom(y,size = 2,prob = theta)
  #print(b)
  #print(a-b)
  temp <- qnorm(abs(a - b))
  return(temp)
}

# function U
# when you input y, q, theta, you can get the upper boundary
U <- function(y,theta){
  a <- pbinom(y,size = 2,prob = theta)
  temp <- qnorm(a)
  return(temp)
}


K <- function(G,z.sta,z.g,low,up,rou){
  f1 <- dtruncnorm(z.sta[1], mean = rou*z.g[2], sd = sqrt(1-rou^2),
   a = low[1], b = up[1])
  #print(f1)
  f2 <- dtruncnorm(z.sta[2], mean = rou*z.sta[1], sd = sqrt(1-rou^2),
   a = low[2], b = up[2])
  #print(f2)
  return(f1*f2)
}

# t = 1,2,...,n
est.each <- function(G,rou,theta,t,Y){
  n <- nrow(Y)
  q <- ncol(Y)
  BL <- Y
  BU <- Y
  for (i in 1:n) {
    for (j in 1:q) {
      BL[i,j] <- L(Y[i,j],theta[j])
      BU[i,j] <- U(Y[i,j],theta[j])
    }
  }
  B.u <- BU[t,]
  B.l <- BL[t,]
  sig <- diag(1,q,q)
  sig[1,2] <- rou
  sig[2,1] <- rou
  set.seed(123)
  z.sim <-rtmvnorm.gibbs(n=G, sigma=sig, upper=B.u, lower = B.l)
  z.sim.meam <- c(mean(z.sim[,1]),mean(z.sim[,2]))
  part.1 <- dmvnorm(z.sim.meam,sigma = sig)
  part.2 <- 0
  for (i in 1:G) {
    part.2 <- part.2 + K(G,z.sta = z.sim.meam,z.g = z.sim[i,],low = B.l,up = B.u,rou = rou)
  }
  part.2 <- part.2/G
  return(log(part.1) - log(part.2))
}

est.all <- function(G,rou,theta,Y){
  n <- nrow(Y)
  est <- 0
  for (i in 1:n) {
    est <- est + est.each(G,rou,theta,i,Y)
  }
  return(est)
}

gra.phi.all <- function(G,phi,Y){
  theta <- phi[1:2]
  rou <- phi[3]
  x <- 2e-10
  theta1 <- theta
  theta1[1] <- theta[1] + x
  theta2 <- theta
  theta2[2] <- theta[2] + x
  rou1 <- rou + x
  y <- est.all(G,rou,theta,Y)
  y1 <- est.all(G,rou,theta1,Y)
  y2 <- est.all(G,rou,theta2,Y)
  yr <- est.all(G,rou1,theta,Y)
  gra.phi <- phi
  gra.phi[1] <- (y1 - y)/x
  gra.phi[2] <- (y2 - y)/x
  gra.phi[3] <- (yr - y)/x
  return(gra.phi)
}

dev.theta1.each <- function(G,rou,theta,t,Y){
  x <- 2e-10
  theta1 <- theta
  theta1[1] <- theta[1] + x
  y <- est.each(G,rou,theta,t,Y)
  y1 <- est.each(G,rou,theta1,t,Y)
  return((y1 - y)/x)
}


dev.theta2.each <- function(G,rou,theta,t,Y){
  x <- 2e-10
  theta1 <- theta
  theta1[2] <- theta[2] + x
  y <- est.each(G,rou,theta,t,Y)
  y1 <- est.each(G,rou,theta1,t,Y)
  return((y1 - y)/x)
}

dev.rou.each <- function(G,rou,theta,t,Y){
  x <- 2e-10
  rou1 <- rou + x
  y <- est.each(G,rou,theta,t,Y)
  y1 <- est.each(G,rou1,theta,t,Y)
  return((y1 - y)/x)
}

gra.phi.each <- function(G,t,phi,Y){
  theta <- phi[1:2]
  rou <- phi[3]
  dtheta1 <- dev.theta1.each(G,rou,theta,t,Y)
  dtheta2 <- dev.theta2.each(G,rou,theta,t,Y)
  drou <- dev.rou.each(G,rou,theta,t,Y)
  gra.phi <- phi
  gra.phi[1] <- dtheta1
  gra.phi[2] <- dtheta2
  gra.phi[3] <- drou
  return(gra.phi)
}

B <- function(G,phi,Y){
  n <- nrow(Y)
  B <- matrix(0,3,3)
  for (i in 1:n) {
    B <- B + gra.phi.each(G,i,phi,Y)%*%t(gra.phi.each(G,i,phi,Y))
  }
  return(B)
}
dis <- function(a,b){
  c = a - b
  return(sum(c^2))
}

BHHH <- function(G,Y,step){
  n <- nrow(Y)
  q <- ncol(Y)
  theta <- matrix(0, ncol = 1, nrow = q)
  sig <- diag(1,q,q)
  for (i in 1:q) {
    theta[i] <- mean(Y[,i])/n
  }
  phi <- theta
  rou <- cor(Y[,1],Y[,2])
  phi[length(theta) + 1] <- rou
  sig[1,2] <- rou
  sig[2,1] <- rou

  new <- as.matrix(phi)
  old <- new - step * solve(B(G,phi,Y)) %*% gra.phi.all(G,phi,Y)
  #print(phi)
  #print(old)
  c = 1
  while(dis(old,new)<0.01){
    c = c + 1
    #print(c)
    new <- old
    phi <- new
    #print(phi)
    old <- new - step * solve(B(G,phi,Y)) %*% gra.phi.all(G,phi,Y)
    #print(old)
  }

  return(new)
}


