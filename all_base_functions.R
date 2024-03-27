### All required functions for using ESS
### Functions related to Wood and Chan algorithm of drawing samples
### Fast.LS function for generating large Gaussian vector prior 
### LS.KLE function for generating large Gaussian vector using Karhunen-Lo\`eve expansion (KLE) 
### MH for sampling from \nu and \ell
### Covariance matrix and design matrix (using basis function) are also defined
### And all related and dependant functions are here

### Required libraries:
library(FastGP) # for tinv and rcpp_rmvnorm function
library(rSPDE) # Matern cov fct
library(mvtnorm);library(MASS)
library(Rfast) # matrnorm function
library(quadprog) # for the MAP (quad optim pb)


####################################################
######### Mat\'ern family covariance kernels #######
####################################################
#Given a \nu (smoothness parameter of matern kernel) finding a value of 
# l (length-scale parameter) such that the correlation between the 
# maximum seperation is some small value, say 0.05

## Matern family cov fct with \nu smooth para & theta length corr
k <- function(h, nu, l){
  matern.covariance(h, sqrt(2 * nu) / l, nu = nu, sigma = 1)
}

# function for uniroot:
fl <- function(l, para){ 
  #para[1]=x, para[2]=y and para[3]=nu of MK : Matern kernel function;
  #para[4]=pre-specified value of the correlation
  a <- k(abs(para[1]-para[2]), para[3], l)
  return(a - para[4])
}

# function for estimating l:
l_est <- function(nu, range, val){
  # nu : smoothness; range : c(min, max) of the range of variable
  # val : pre-specified value of the correlation between the maximum seperation
  para <- c(range[1], range[2], nu, val)
  rl <- uniroot(f=fl, interval = c(0.000001, 100000), para)
  return(rl$root)
}

# Covariance matrix
covmat <- function(knot, nu, l){
  k(outer(knot, knot, '-'), nu = nu, l = l) 
}



### Define design matrix ###
### The basis functions ###

#####################################################
######## Maatouk & Bay2017 Basis functions ##########
#####################################################
## hat basis functions
h <- function(x){
  ifelse(x >= -1 & x <= 1, 1-abs(x), 0)
}
hi <- function(x, u, i){
  delta <- (max(u)-min(u))/(length(u)-1)
  h((x - u[i]) / delta)
}
## integral of hat basis functions phi
## phi functions
phi1 <- function(x, u){
  delta <- (max(u) - min(u)) / (length(u)-1)
  ifelse(x >= u[1] & x <= u[2], delta/2-((u[2]-x)*hi(x,u,1))/2, delta/2)
}
phi <- function(x, u, i){
  delta <- (max(u)-min(u))/(length(u)-1)
  ifelse (x <= u[i] - delta, 0,
         ifelse (x >= u[i] - delta & x <= u[i], (x-u[i]+delta)^2/(2*delta),
                ifelse (x >= u[i] & x<= u[i] + delta, delta -
                         (-x+u[i]+delta)^2/(2*delta), delta)))
}
phii_v <- function(x, u, i){
  if (i == 1){
    phi1(x,u)
  }
  else {
    phi(x, u, i)
  }
}

##############################################

##############################################
####### function of design matrix ############
##############################################
## fct design matrix (hat basis function)
fcth <- function(x, u, N){
  n <- length(x)
  h <- matrix(NA, nrow = n, ncol = N)
  for (j in 1 : N){
    h[, j] <- hi(x, u, i = j)
  }
  return(h)
}

## fct design matrix (phi basis function)
fctphi <- function(x,u,N){
  n <- length(x)
  phi <- matrix(NA, nrow = n, ncol = N)
  for(j in 1 : N){
    phi[, j] <- phii_v(x,u,i=j)
  }
  return(phi)
}
####################################################


# Order of the circulant matrix:
# minimum value of g and m so that G can be embedded into C
min_g <- function(knot){
  N <- length(knot)
  g <- ceiling(log(2*N,2))   #m=2^g and m>=2(n-1) : Wood & Chan notation; 
  #since we are going upto n and not stopping at (n-1), the condition is modified!
  return("g" = g)
}

# forming the circulant matrix:
circulant <- function(x){
  n <- length(x)
  mat <- matrix(0, n, n)
  for (j in 1:n) {
    mat[j, ] <- c(x[-(1:(n+1-j))], x[1:(n+1-j)])
  }
  return(mat)
}

# Function for forming the vector of circulant matrix:
circ_vec <- function(knot,g,nu,l,tausq){
  delta_N <- 1/(length(knot)-1)
  m <- 2**g
  cj <- integer()
  for(j in 1:m){
    if(j<=(m/2))
      cj[j]=(j-1)*delta_N
    else
      cj[j]=(m-(j-1))*delta_N
  }
  x <- (tausq*k(cj,nu,l))
  return(x)
}


# Function for finding a g such that C is nnd:
# without forming the circulant matrix and without computing eigen values:
C.eval <- function(knot,g,nu,l,tausq){
  vec <- circ_vec(knot,g,nu,l,tausq)
  val <- fft(vec) # eigenvalues will be real as the circulant matrix formed by the 
  # vector is by construction is symmetric!
  ev <- min(Re(val))
  return(list("vec" = vec, "min.eig.val" = ev))
}


nnd_C <- function(knot,g,nu,l,tausq){
  C.vec <- C.eval(knot,g,nu,l,tausq)$vec
  eval <- C.eval(knot,g,nu,l,tausq)$min.eig.val
  if(eval>0)
    return(list("cj" = C.vec,"g" = g))
  else{
    g <- g+1
    nnd_C(knot,g,nu,l,tausq)
  }
}

# computing the eigen values of C using FFT:
eigval <- function(knot,nu,l,tausq){
  g <- min_g(knot)
  c.j <- nnd_C(knot,g,nu,l,tausq)$cj
  lambda <- Re(fft(c.j))
  if(min(lambda)>0)
    return(lambda)
  else
    stop("nnd condition is NOT satisfied!!!")
}


#################################################################
########## Samples drawn using Wood and Chan Algorithm ##########
#################################################################
samp.WC <- function(knot,nu,l,tausq,sseedWC=1){
  N <- length(knot)
  lambda <- eigval(knot,nu,l,tausq)
  m <- length(lambda)
  samp.vec <- rep(0,N)
  set.seed(sseedWC)
  a <- rep(0,m)
  a[1] <- sqrt(lambda[1])*rnorm(1)/sqrt(m)
  a[(m/2)+1] <- sqrt(lambda[(m/2)+1])*rnorm(1)/sqrt(m)
  i <- sqrt(as.complex(-1))
  for(j in 2:(m/2)){
    uj <- rnorm(1); vj <- rnorm(1)
    a[j] <- (sqrt(lambda[j])*(uj + i*vj))/(sqrt(2*m))
    a[m+2-j] <- (sqrt(lambda[j])*(uj - i*vj))/(sqrt(2*m))
  }
  samp <- fft(a)
  samp.vec <- Re(samp[1:N])
  return(samp.vec)
}
#############################################


#############################################
########## Functions for using ESS ##########
#############################################
ESS <- function(beta,nu_ess,y,X,sigsq,eta,seeds=1){
  thetamin <- 0; 
  thetamax <- 2*pi;
  
  set.seed(seeds)
  u <- runif(1)
  logy <- loglik(y,X,sigsq,eta,beta) + log(u); 
  
  theta <- runif(1,thetamin,thetamax); 
  thetamin <- theta - 2*pi; 
  thetamax <- theta;
  betaprime <- beta*cos(theta) + nu_ess*sin(theta);
  
  while(loglik(y,X,sigsq,eta,betaprime) <= logy){
    if(theta < 0)
      thetamin <- theta
    else
      thetamax <- theta
    theta <- runif(1,thetamin,thetamax)
    betaprime <- beta*cos(theta) + nu_ess*sin(theta)
  }
  return(betaprime)       
}

ESS.dec <- function(beta,nu_ess,y,X,sigsq,eta,seeds=1){
  thetamin <- 0;
  thetamax <- 2*pi;
  set.seed(seeds)
  u <- runif(1)
  logy <- loglik2(y,X,sigsq,eta,beta) + log(u);
  
  theta <- runif(1,thetamin,thetamax);
  thetamin <- theta - 2*pi;
  thetamax <- theta;
  betaprime <- beta*cos(theta) + nu_ess*sin(theta);
  
  while(loglik2(y,X,sigsq,eta,betaprime) <= logy){
    if(theta < 0)
      thetamin <- theta
    else
      thetamax <- theta
    theta <- runif(1,thetamin,thetamax)
    betaprime <- beta*cos(theta) + nu_ess*sin(theta)
  }
  return(betaprime)
}




## Defining the loglik function to be used in ESS:
## loglik calculates the log of the likelihood:
## nondecreasing constraints
loglik <- function(y,X,sigsq,eta,beta){
  mu <- y-(X%*%beta)
  val <- -sum(log(1+exp(-eta*beta)))-
    # eta*sum(beta)-sum(log(1+exp(eta*beta)))-
    sum(mu^2)/(2*sigsq)
  return(val)
}

## nonincreasing constraints:
loglik2 <- function(y,X,sigsq,eta,beta){
  mu <- y-(X%*%beta)
  val <- -sum(log(1+exp(eta*beta)))-sum(mu^2)/(2*sigsq)
  return(val)
}


## MH algo for \nu and \ell of Matern kernel:
nu.MH2 = function(nu.in,l.in,tau.in,xi.in,knot,range.nu=c(0.5,2.5),range.l=c(0.1,1),sd.nu=0.05,sd.l=0.1,seed=1){
  Kmat = covmat(knot,nu.in,l.in)
  Linv = solve(chol(Kmat))#tinv(Kmat)#inv_chol(Kmat)#
  set.seed(seed)
  nu.cand = exp(log(nu.in)+rnorm(1,0,sd.nu))
  l.cand = exp(log(l.in)+rnorm(1,0,sd.l))
  dnu = dunif(nu.cand,range.nu[1],range.nu[2])
  dl = dunif(l.cand,range.l[1],range.l[2])
  if(dnu > 0 && dl > 0){
    Kcand = covmat(knot,nu.cand,l.cand)
    Linv.cand = solve(chol(Kcand))#tinv(Kcand)#inv_chol(Kcand)
    t1 = sum((t(Linv.cand)%*%xi.in)^2)#sum((t(xi.in)%*%Linv.cand%*%xi.in)^2)
    t2 = sum((t(Linv)%*%xi.in)^2)
    r = exp(sum(log(diag(Linv.cand)))-sum(log(diag(Linv)))-((t1 - t2)/(2*tau.in)))*(nu.cand/nu.in)*(l.cand/l.in)
    alpha = min(r,1)
  }
  else{
    alpha = 0
    Linv.cand = Linv
  }
  u = runif(1)
  nu.out = (u < alpha)*nu.cand + (u >= alpha)*nu.in
  l.out = (u < alpha)*l.cand + (u >= alpha)*l.in
  cnt = (u < alpha)*1 + (u >= alpha)*0
  L_inv = (u < alpha)*Linv.cand + (u >= alpha)*Linv
  return(list("nu" = nu.out,"l" = l.out,"count" = cnt,"L_inv"=L_inv))
}



#############################################
########### Fast Large-scale ################
#############################################

### This function generates a multivariate normal (MVN) distribution with a covariance matrix extracted from a Mat\'ern kernel (MK)
### M represents the number of subdomains (M>=1)
### N1 represents the number of equally-spaced points of the 1st subdomain
### u: grid vector, should be of length M times N1
### nu represents the smoothness parameter of the MK (nu>0)
### l represents the correlation length-scale parameter of MK (l>0)
### tol: tolerance (relative to largest variance) for numerical lack of a positive-definiteness problem

library(mvtnorm)
# library(FastGP) # when using `rcppeigen_get_chol' for Cholesky factorization

Fast.LS <- function(u,M,N1,nu,l,tausq,tol,sseedLS=1){
  if(missing(tol)){
    tol <- 1e-5
  }
  if(missing(M)){
    M=1
  }
  if(M==0)
    stop("M cannot be zero")
  if(N1==0)
    stop("N1 cannot be zero")
  if(length(u)!= M*N1)
    stop("The length of the vector u must be M times N1")
  u1 <- u[1:N1]
  u2 <- u[(N1+1):(2*N1)]
  Gamma11 <- tausq*covmat(u1,nu=nu,l=l)+tol*diag(N1)
  set.seed(sseedLS)
  if(M==1){
    return(as.vector(rcpp_rmvnorm(n=1,S=Gamma11,mu=rep(0,N1))))
    # return(as.vector(mvtnorm::rmvnorm(n=1,mean=rep(0,N1),Gamma11,method='chol')))
  }
  else 
    Gamma12 <- tausq*k(outer(u2,u1,'-'),nu=nu,l=l)
  Ktilde <- Gamma12%*%tinv(Gamma11) # coupling matrix
  ## slightly faster for computing L
  # L <- rcppeigen_get_chol(Gamma11-Gamma12%*%t(Ktilde))%*%
  #   solve(((rcppeigen_get_chol(Gamma11))))
  L <- t(chol(Gamma11-Gamma12%*%t(Ktilde)))%*%
    solve(t(chol(Gamma11)))
  eta <- rcpp_rmvnorm(n=M,S=Gamma11,mu=rep(0,N1))
  # mvtnorm::rmvnorm(n=M,mean=rep(0,N1),Gamma11,method='chol')
  etaT <- matrix(NA,M,N1)
  etaT[1,] <- eta[1,]
  for(i in 2 : M){
    etaT[i,] <- Ktilde%*%(etaT[i-1,])+L%*%(eta[i,])
  }
  return(as.vector(t(etaT)))
}


## Combination Fast Large-scale and FFT from WC
Fast.LS.WC <- function(u,M,N1,nu,l,tausq,tol,sseedWC=1){
  if(missing(tol)){
    tol <- 1e-6
  }
  if(missing(M)){
    M=1
  }
  if(M==0)
    stop("M cannot be zero")
  if(N1==0)
    stop("N1 cannot be zero")
  if(length(u)!= M*N1)
    stop("The length of the vector u must be M times N1")
  u1 <- u[1:N1]
  u2 <- u[(N1+1):(2*N1)]
  Gamma12 <- tausq*k(outer(u2,u1,'-'),nu=nu,l=l)
  Gamma11 <- tausq*covmat(u1,nu=nu,l=l)+tol*diag(N1)
  if(M==1){
    return(samp.WC(u1,nu=nu,l=l,tausq=tausq,sseedWC))
  }
  else 
    # K <- Gamma12
    Ktilde <- Gamma12%*%tinv(Gamma11)
  L <- t(chol(Gamma11-Gamma12%*%t(Ktilde)))%*%
    solve(t(chol(Gamma11)))
  eta <- matrix(NA,M,N1)
  eta[1,] <- samp.WC(u1,nu=nu,l=l,tausq=tausq,sseedWC)
  etaT <- matrix(NA,M,N1)
  etaT[1,] <- eta[1,]
  for(i in 2 : M){
    eta[i,] <- samp.WC(u1,nu=nu,l=l,tausq=tausq,sseedWC)
    etaT[i,] <- Ktilde%*%(etaT[i-1,])+L%*%(eta[i,])
  }
  return(as.vector(t(etaT)))
}
################################

### Fast Large-scale for more than one sample
Fast.LS_v <- function(nbsim,u,M,N1,nu,l,tausq,tol,sseedLS=1){
  if(missing(tol)){
    tol <- 1e-6
  }
  if(missing(M)){
    M=1
  }
  if(M==0)
    stop("M cannot be zero")
  if(N1==0)
    stop("N1 cannot be zero")
  if(length(u)!= M*N1)
    stop("The length of the vector u must be M times N1")
  u1 <- u[1:N1]
  u2 <- u[(N1+1):(2*N1)]
  Gamma11 <- tausq*covmat(u1,nu=nu,l=l)+tol*diag(N1)
  set.seed(sseedLS)
  if(M==1){
    return(as.vector(mvtnorm::rmvnorm(n=nbsim,mean=rep(0,N1),Gamma11,method='chol')))
  }
  else 
    Gamma12 <- tausq*k(outer(u2,u1,'-'),nu=nu,l=l)
  Ktilde <- Gamma12%*%tinv(Gamma11) # coupling matrix
  ## slightly faster for computing L
  # L <- rcppeigen_get_chol(Gamma11-Gamma12%*%t(Ktilde))%*%
  #   solve(((rcppeigen_get_chol(Gamma11))))
  L <- t(chol(Gamma11-Gamma12%*%t(Ktilde)))%*%
    solve(t(chol(Gamma11)))
  eta <- t(mvtnorm::rmvnorm(n=nbsim,mean=rep(0,N1),Gamma11,method='chol'))
  etaT <- list()
  etaT[[1]] <- eta
  for(i in 2 : M){
    eta <- t(mvtnorm::rmvnorm(n=nbsim,mean=rep(0,N1),Gamma11,method='chol'))
    etaT[[i]] <- Ktilde%*%(etaT[[i-1]])+L%*%(eta)
  }
  return(return(do.call(rbind,etaT)))
}
##################################################




################ LS KLE one sample ################
###################################################
LS.KLE <- function(u,N1,p,M,nu,l,tausq,tol,sseedLS=1){
  if(missing(tol)){
    tol <- 1e-8
  }
  if(missing(p)){
    p = N1
  }
  if(missing(M)){
    M=1
  }
  if(M==0)
    stop("M cannot be zero")
  if(N1==0)
    stop("N1 cannot be zero")
  if(length(u)!= M*N1)
    stop("The length of the vector u must be M times N1")
  u1 <- u[1:N1]
  u2 <- u[(N1+1):(2*N1)]
  Gamma11 <- tausq*covmat(u1,nu,l)+tol*diag(N1)
  set.seed(sseedLS)
  if(M==1){
    return(as.vector(rcpp_rmvnorm(n=1,S=Gamma11,mu=rep(0,N1))))
  }
  else 
    Gamma12 <- tausq*k(outer(u1,u2,'-'),nu,l)
  eig11 <- eigen(Gamma11)
  value11 <- eig11$values[1:p]
  vector11 <- eig11$vectors[,1:p]
  K12 <- (((t(vector11)%*%(Gamma12))%*%
             (vector11))/
            sqrt(tcrossprod(value11)))
  L12 <- t(chol(diag(p)-crossprod(K12)))
  eta <- matrnorm(M,p)
  etaT <- matrix(NA,M,p)
  f <- matrix(NA,M,N1)
  f[1,] <- vector11%*%(sqrt(value11)*eta[1,])  
  etaT[1,] <- eta[1,]
  for(i in 2 : M){
    etaT[i,] <- t(K12)%*%etaT[(i-1),]+L12%*%eta[i,]
    f[i,] <- vector11%*%(sqrt(value11)*etaT[i,])  
  }
  return(as.vector(t(f)))
}
####################################################


## end
