#---------------------------------------------------
####### MCMC samples using ESS and LS algo #########
#---------------------------------------------------

setwd("~/Documents/Recherche/Article Maatouk and Pallavi Ray/code Maatouk & Pallavi Ray/Official codes")
source('all_base_functions.R')
library(Matrix) # crossprod
library(tmg) # HMC sampler
library(restrictedMVN) ## Gibbs sampler

#---------------------------------------------------
########## nonnegativity constraints ###############
#---------------------------------------------------

#---------------------------------------------------
### Function for drawing posterior samples using ESS and LS with fixed hyperparameters:
### For nonnegative functions estimation 
#---------------------------------------------------
pos.LS.ESS <- function(y, x, N1, M, nu, l, eta, mcmc, brn, thin, tau.in, sig.in, xi.in,
                       tau.fix, sig.fix, xi.fix, sseed, verbose, return.plot, tol){
  # y: Response variable; x: vector to form design matrix X (n x N)
  # N1: number of knots first subdomain; M: nb of subdomain
  # nu: smoothness parameter of Matern; l: length-scale parameter of Matern
  # eta: parameter of the approximation function of the indicator functions
  # mcmc, brn, thin : mcmc samples, burning and thinning for MCMC
  # tau.in, sig.in, xi.in : initial values (supplied by user or use the default values)
  # tau.fix,sig.fix,xi.fix : if fixed values of the parameters are to use
  # verbose : logical; if TRUE, prints current status; default is TRUE
  # return.plot : logical; if true a plot of estimate with 95% CI is returned; default is TRUE
  # tol: tolerance for numerical stability
  
  # OUTPUT: Posterior samples on xi,tau,sig and fhat with posterior mean, 95% CI of fhat, mode of posterior distribution z_star
  
  if (length(y) != length(x))
    stop("y and x should be of same length!")
  y <- y[order(x)]
  x <- sort(x)
  n <- length(y)
  N <- N1 * M
  delta <- 1 / (N - 1)
  my_knots <- seq(from = 0, to = 1, by = delta)
  X <- fcth(x, u = my_knots, N)
  
  if (missing(nu))
    stop ("nu needs to be supplied")
  if (nu == 0)
    stop("nu cannot be zero")
  if (!missing(l) && l == 0)
    stop ("l cannot be zero")
  if (missing(return.plot))
    return.plot <- TRUE
  if (!missing(sseed))
    set.seed(sseed)
  if (missing(sseed))
    set.seed(Sys.Date())
  if (missing(verbose))
    verbose <- TRUE
  if (missing(l))
    l <- l_est(nu = nu, range = range(my_knots), val = 0.05)
  
  # prior covariance K:
  K <- covmat(my_knots, nu, l)
  # prior precision:
  K_inv <- tinv(K)
  
  if (missing(tol))
    tol <- 1e-6
  if (missing(eta))
    eta <- 50
  if (missing(mcmc))
    mcmc <- 5000
  if (missing(brn))
    brn <- 1000
  if (missing(thin))
    thin <- 1
  em <- mcmc + brn
  ef <- mcmc / thin
  
  if (!missing(tau.fix))
    tau.in <- tau.fix
  if (!missing(sig.fix))
    sig.in <- sig.fix
  if (!missing(xi.fix))
    xi.in <- xi.fix
  
  if (missing(tau.fix) && missing(tau.in))
    tau.in <- 1
  if (missing(sig.fix) && missing(sig.in))
    sig.in <- 1
  if (missing(xi.fix) && missing(xi.in)){
    ## posterior Mode
    XXK <- crossprod(X) / sig.in + K_inv / tau.in
    Amat <- diag(N)
    z_star.in <- solve.QP(Dmat = XXK, dvec = as.vector(t(X) %*% y) / sig.in,
                          Amat = t(Amat), bvec = rep(0, N), meq = 0)$solution
    xi.in <- z_star.in 
  }
  tau <- tau.in
  sig <- sig.in
  xi_in <- xi.in
  
  xi_sam <- matrix(NA, nrow = N, ncol = ef)
  tau_sam <- rep(NA, ef)
  sig_sam <- rep(NA, ef)
  fhat_sam <- matrix(NA, nrow = n, ncol = ef)
  
  if (verbose)
    print("MCMC sample draws:")
  
  ptm <- proc.time()
  for (i in 1 : em) {
    # sampling Xi:
    if (missing(xi.fix)) {
      nu.ess <- #LS.KLE(u = my_knots, M = M, N1 = N1, nu = nu, l = l, tausq = tau, tol = tol, sseedLS = i)
        Fast.LS(u = my_knots, M = M, N1 = N1, nu = nu, l = l, tausq = tau, tol = tol, sseedLS = i)
      xi_out <- ESS(beta = xi_in, nu_ess = nu.ess, y = y, X = X, sigsq = sig, eta = eta, seeds = i)
      xi_out <- sapply(xi_out, function(z) return(pmax(0, z)))
    } else {
      xi_out <- xi_in
    }
    set.seed(2 * i)
    # sampling \sigma^2:
    Xxi <- as.vector(X %*% xi_out)
    y_star <- y - Xxi
    if (missing(sig.fix))
      sig <- 1 / rgamma(1, shape = n/2, rate = sum(y_star^2) / 2)
    
    # sampling \tau^2:
    if (missing(tau.fix))
      tau <- 1 / rgamma(1, shape = N/2, rate = (t(xi_out) %*% K_inv %*% xi_out) / 2)
    
    # storing MCMC samples:
    if (i > brn && i%%thin == 0) {
      xi_sam[,(i - brn) / thin] <- xi_out
      sig_sam[(i - brn) / thin] <- sig
      tau_sam[(i - brn) / thin] <- tau
      fhat_sam[, (i - brn) / thin] <- Xxi
    }
    
    if (i %% 1000 == 0 && verbose){
      print(i)
    }
    
    # renewing the intial value:
    xi_in <- xi_out
  } 
  tm <- proc.time() - ptm
  
  ## posterior Mode
  XXK <- crossprod(X) / mean(sig_sam) + K_inv / mean(tau_sam)
  Amat <- diag(N)
  z_star <- solve.QP(XXK, dvec = as.vector(t(X) %*% y) / mean(sig_sam),
                     Amat = t(Amat), bvec = c(rep(0, N)), meq = 0)$solution
  
  MAP <- X %*% z_star # MAP estimate
  fmean <- rowMeans(fhat_sam) # mAP estimate
  qnt <- apply(fhat_sam, 1, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))
  f_low <- qnt[1, ]
  f_upp <- qnt[2, ]
  ub <- max(f_low, f_upp, fmean, MAP)
  lb <- min(f_low, f_upp, fmean, MAP)
  
  if (return.plot) {
    par(mfrow = c(1, 1))
    par(mar = c(2.1, 2.1, 2.1, 1.1)) # adapt margins
    plot(x, y, pch = '*', lwd = 2, lty = 1, col = 'black',
         ylim = range(ub, lb, y), xlab = '', ylab = '')
    polygon(c(x, rev(x)), y = c(f_low, rev(f_upp)), border = F, col = 'gray')
    lines(x, fmean, type = 'l', lty = 4, lwd = 2, col = 'blue')
    lines(x, MAP, type = 'l', lty = 2, lwd = 2, col = 'red')
    points(x, y, pch = '*')
    abline(h = 0, lty = 2, lwd = 2)
  }
  
  return(list("time" = tm, "xi_sam" = xi_sam, "sig_sam" = sig_sam, "tau_sam" = tau_sam,
              "fhat_sam" = fhat_sam, "fmean" = fmean, "f_low" = f_low, "f_upp" = f_upp, "z_star" = z_star))
}
####################################################


#---------------------------------------------------
### Function for drawing posterior samples using ESS and WC with fixed hyperparameters:
### For nonnegative functions estimation 
#---------------------------------------------------
pos.WC.ESS <- function(y, x, N, nu, l, eta, mcmc, brn, thin, tau.in, sig.in, xi.in,
                       tau.fix, sig.fix, xi.fix, sseed, verbose, return.plot, tol){
  # y: Response variable; x: vector to form design matrix X (n x N)
  # N: the number of knots
  # nu:smoothness parameter of Matern; l:length-scale parameter of Matern
  # eta: parameter of the approximation function of the indicator functions
  # mcmc, brn, thin : mcmc samples, burning and thinning for MCMC
  # tau.in, sig.in, xi.in : initial values (supplied by user or use the default values)
  # tau.fix,sig.fix,xi.fix : if fixed values of the parameters are to use
  # verbose : logical; if TRUE, prints current status; default is TRUE
  # return.plot : logical; if true a plot of estimate with 95% CI is returned; default is TRUE
  
  # OUTPUT: Posterior samples on xi,tau,sig and fhat with posterior mean, 95% CI of fhat, mode of posterior distribution z_star
  
  if (length(y) != length(x))
    stop("y and x should be of same length!")
  y <- y[order(x)]
  x <- sort(x)
  n <- length(y)
  if (missing(N))
    N <- ceiling(n / 2) - 1
  
  delta <- 1 / (N - 1)
  my_knots <- seq(from = 0, to = 1, by = delta)
  X <- fcth(x = x, u = my_knots, N = N)
  
  if (missing(nu))
    stop ("nu needs to be supplied")
  if (nu == 0)
    stop ("nu cannot be zero")
  if (!missing(l) && l == 0)
    stop ("l cannot be zero")
  if (missing(return.plot))
    return.plot <- TRUE
  if (!missing(sseed))
    set.seed(sseed)
  if (missing(sseed))
    set.seed(Sys.Date())
  if (missing(verbose))
    verbose <- TRUE
  if (missing(l))
    l <- l_est(nu = nu, range = range(my_knots), 0.05)
  
  # prior covariance K:
  K <- covmat(my_knots, nu, l)
  # prior precision:
  K_inv <- tinv(K)
  
  if (missing(tol))
    tol <- 1e-6
  if (missing(eta))
    eta <- 50
  if (missing(mcmc))
    mcmc <- 5000
  if (missing(brn))
    brn <- 1000
  if (missing(thin))
    thin <- 1
  em <- mcmc + brn
  ef <- mcmc / thin
  
  if (!missing(tau.fix))
    tau.in <- tau.fix
  if (!missing(sig.fix))
    sig.in <- sig.fix
  if (!missing(xi.fix))
    xi.in <- xi.fix
  if (missing(tau.fix) && missing(tau.in))
    tau.in <- 1
  if (missing(sig.fix) && missing(sig.in))
    sig.in <- 1
  if (missing(xi.fix) && missing(xi.in)){
    ## posterior Mode
    XXK <- crossprod(X) / sig.in + K_inv / tau.in
    Amat <- diag(N)
    z_star.in <- solve.QP(Dmat = XXK, dvec = as.vector(t(X) %*% y) / sig.in,
                          Amat = t(Amat), bvec = c(rep(0, N)), meq = 0)$solution
    xi.in <- z_star.in
  }
  tau <- tau.in
  sig <- sig.in
  xi_in <- xi.in
  
  xi_sam <- matrix(NA, nrow = N, ncol = ef)
  tau_sam <- rep(NA, ef)
  sig_sam <- rep(NA, ef)
  fhat_sam <- matrix(NA, nrow = n, ncol = ef)
  
  if(verbose)
    print("MCMC sample draws:")
  
  ptm <- proc.time()
  for (i in 1 : em) {
    # sampling Xi:
    if (missing(xi.fix)) {
      nu.ess <- as.vector(samp.WC(knot = my_knots, nu = nu, l = l, tausq = tau, sseedWC = i))
      xi_out <- ESS(beta = xi_in, nu_ess = nu.ess, y = y, X = X, sigsq = sig, eta = eta, seeds = i)
      xi_out <- sapply(xi_out, function(z) return(pmax(0, z)))
    } else {
      xi_out <- xi_in
    }
    set.seed(2 * i)
    # sampling \sigma^2:
    Xxi <- as.vector(X %*% xi_out)
    y_star <- y - Xxi
    if (missing(sig.fix))
      sig <- 1 / rgamma(1, shape = n/2, rate = sum(y_star^2)/2)
    
    # sampling \tau^2:
    if (missing(tau.fix))
      tau <- 1 / rgamma(1, shape = N / 2, rate = (t(xi_out) %*% K_inv %*% xi_out) / 2)
    
    # storing MCMC samples:
    if (i > brn && i%%thin == 0) {
      xi_sam[, (i - brn) / thin] <- xi_out
      sig_sam[(i - brn) / thin] <- sig
      tau_sam[(i - brn) / thin] <- tau
      fhat_sam[, (i - brn) / thin] <- Xxi
    }
    
    if (i %% 1000 == 0 && verbose) {
      print(i)
    }
    
    # renewing the intial value:
    xi_in <- xi_out
  } 
  tm <- proc.time() - ptm
  
  ## posterior Mode
  XXK <- crossprod(X) / mean(sig_sam) + K_inv / mean(tau_sam)
  Amat <- diag(N)
  z_star <- solve.QP(Dmat = XXK, dvec = as.vector(t(X) %*% y) / mean(sig_sam),
                     Amat = t(Amat), bvec = c(rep(0, N)), meq = 0)$solution
  
  MAP <- X %*% z_star # MAP estimate
  fmean <- rowMeans(fhat_sam) # mAP estimate
  qnt <- apply(fhat_sam, 1, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))
  f_low <- qnt[1, ]
  f_upp <- qnt[2, ]
  ub <- max(f_low, f_upp, fmean, MAP)
  lb <- min(f_low, f_upp, fmean, MAP)
  
  if (return.plot) {
    par(mfrow = c(1, 1))
    par(mar = c(2.1, 2.1, 2.1, 1.1)) # adapt margins
    plot(x, y, pch = '*', lwd = 2, lty = 1, col = 'black',
         ylim = range(ub, lb, y), xlab = '', ylab = '')
    polygon(c(x, rev(x)), y = c(f_low, rev(f_upp)), border = F, col = 'gray')
    lines(x, fmean, type = 'l', lty = 4, lwd = 2, col = 'blue')
    lines(x, MAP, type = 'l', lty = 2, lwd = 2, col = 'red')
    points(x, y, pch = '*')
    abline(h = 0, lty = 2, lwd = 2)
  }
  
  return(list("time" = tm, "xi_sam" = xi_sam, "sig_sam" = sig_sam, "tau_sam" = tau_sam,
              "fhat_sam" = fhat_sam, "fmean" = fmean, "f_low" = f_low, "f_upp" = f_upp, "z_star" = z_star))
}
####################################################







#---------------------------------------------------
### Function for drawing posterior samples using ESS and LS with hyperparameter updates:
### For nonnegative functions estimation 
#---------------------------------------------------
pos.LS.ESS_hyp <- function(y, x, N1, M, nu.in, l.in, eta, mcmc, brn, thin, tau.in, sig.in, xi.in,
                           tau.fix, sig.fix, xi.fix, sseed, verbose, return.plot, tol) {
  # y: Response variable; x: vector to form design matrix X (n x N)
  # N1: number of knots first subdomain; M: nb of subdomain
  # eta:parameter of the approximation function of the indicator functions
  # mcmc, brn, thin : mcmc samples, burning and thinning for MCMC
  # nu.in, l.in, tau.in, sig.in, xi.in : initial values (supplied by user or use the default values)
  # tau.fix,sig.fix,xi.fix : if fixed values of the parameters are to use
  # verbose : logical; if TRUE, prints currnet status; default is TRUE
  # return.plot : logical; if true a plot of estimate with 95% CI is returned; default is TRUE
  # tol: tolerance for numerical stability
  
  # OUTPUT: Posterior samples on xi,tau,sig and fhat with posterior mean, 95% CI of fhat, z_star: mode of posterior distribution
  
  if (length(y) != length(x))
    stop("y and x should be of same length!")
  y <- y[order(x)]
  x <- sort(x)
  n <- length(y)
  N <- N1 * M
  delta <- 1 / (N - 1)
  my_knots <- seq(from = 0, to = 1, by = delta)
  X <- fcth(x = x, u = my_knots, N = N)
  
  if (nu.in == 0)
    stop("nu cannot be zero")
  if (!missing(l.in) && l.in == 0)
    stop("l cannot be zero")
  if (missing(return.plot))
    return.plot <- TRUE
  if (!missing(sseed))
    set.seed(sseed)
  if (missing(sseed))
    set.seed(Sys.Date())
  if (missing(verbose))
    verbose <- TRUE
  if (missing(nu.in))
    nu.in <- 0.75
  if (missing(l.in))
    l.in <- l_est(nu = nu.in, range = range(my_knots), val = 0.05)
  if (missing(tol))
    tol <- 1e-6
  if (missing(eta))
    eta <- 50
  if (missing(mcmc))
    mcmc <- 5000
  if (missing(brn))
    brn <- 1000
  if (missing(thin))
    thin <- 1
  em <- mcmc + brn
  ef <- mcmc / thin
  
  if (!missing(tau.fix))
    tau.in <- tau.fix
  if (!missing(sig.fix))
    sig.in <- sig.fix
  if (!missing(xi.fix))
    xi.in <- xi.fix
  if (missing(tau.fix) && missing(tau.in))
    tau.in <- 1
  if (missing(sig.fix) && missing(sig.in))
    sig.in <- 1
  if (missing(xi.fix) && missing(xi.in)) {
    ## posterior Mode
    K <- covmat(knot = my_knots, n = nu.in, l = l.in)
    K_inv <- tinv(K)
    XXK <- crossprod(X) / sig.in + K_inv / tau.in
    Amat <- diag(N)
    z_star.in <- solve.QP(Dmat = XXK, dvec = as.vector(t(X) %*% y) / sig.in,
                          Amat = t(Amat),bvec = c(rep(0, N)), meq = 0)$solution
    xi.in <- z_star.in # starting point McMC
  }
  tau <- tau.in
  sig <- sig.in
  xi_in <- xi.in
  nu_in <- nu.in
  l_in <- l.in
  
  xi_sam <- matrix(NA, nrow = N, ncol = ef)
  tau_sam <- rep(NA, ef)
  sig_sam <- rep(NA, ef)
  nu_sam <- rep(NA, ef)
  ell_sam <- rep(NA, ef)
  fhat_sam <- matrix(NA, nrow = n, ncol = ef)
  
  if (verbose)
    print("MCMC sample draws:")
  
  ptm <- proc.time()
  for (i in 1 : em) {
    # sampling from \nu and \ell
    MH.out <- nu.MH2(nu_in, l_in, tau, xi_in, my_knots, range.nu = c(0.5, 2.5), range.l = c(0.1, 1), seed = i)
    nu_out <- MH.out$nu
    l_out <- MH.out$l
    L_inv <- MH.out$L_inv
    
    # sampling Xi:
    if (missing(xi.fix)) {
      nu.ess <- Fast.LS(u = my_knots, M = M, N1 = N1, nu = nu_out, l = l_out, tausq = tau, tol = tol, sseedLS = i)
      xi_out <- ESS(beta = xi_in, nu_ess <- nu.ess, y = y, X = X, sigsq = sig, eta = eta, seeds = i)
      xi_out <- sapply(xi_out, function(z) return(pmax(0, z)))
    } else {
      xi_out <- xi_in
    }
    set.seed(2 * i)
    
    # sampling \sigma^2:
    Xxi <- as.vector(X %*% xi_out)
    y_star <- y - Xxi
    if (missing(sig.fix))
      sig <- 1/rgamma(1,shape = n / 2, rate = sum(y_star^2) / 2)
    
    # sampling \tau^2:
    if (missing(tau.fix))
      tau <- 1/rgamma(1, shape = N / 2, rate = (sum((t(L_inv) %*% xi_out)^2)) / 2)
    
    # storing MCMC samples:
    if (i > brn && i%%thin == 0) {
      xi_sam[, (i - brn) / thin] <- xi_out
      sig_sam[(i - brn) / thin] <- sig
      tau_sam[(i - brn) / thin] <- tau
      nu_sam[(i - brn) / thin] <- nu_out
      ell_sam[(i - brn) / thin] <- l_out
      fhat_sam[, (i - brn) / thin] <- Xxi
    }
    
    if (i %% 1000 == 0 && verbose) {
      print(i)
    }
    # renewing the intial value:
    xi_in <- xi_out
    nu_in <- nu_out
    l_in <- l_out
  } 
  tm  <- proc.time() - ptm
  ## Posterior Mode
  K <- covmat(knot = my_knots, nu = mean(nu_sam), l = mean(ell_sam))
  K_inv <- tinv(K)
  XXK <- crossprod(X) / mean(sig_sam) + K_inv / mean(tau_sam)
  Amat <- diag(N)
  z_star <- solve.QP(Dmat = XXK, dvec = as.vector(t(X) %*% y) / mean(sig_sam),
                     Amat = t(Amat), bvec = c(rep(0, N)), meq = 0)$solution
  
  MAP <- X %*% z_star # MAP estimate
  fmean <- rowMeans(fhat_sam) # mAP estimate
  qnt <- apply(fhat_sam, 1, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))
  f_low <- qnt[1, ]
  f_upp <- qnt[2, ]
  ub <- max(f_low, f_upp, fmean, MAP)
  lb <- min(f_low, f_upp, fmean, MAP)
  
  if (return.plot) {
    par(mfrow = c(1, 1))
    par(mar = c(2.1, 2.1, 2.1, 1.1)) # adapt margins
    plot(x, y, pch = '*', lwd = 2, lty = 1, col = 'black',
         ylim = range(ub, lb, y), xlab = '', ylab = '')
    polygon(c(x, rev(x)), y = c(f_low, rev(f_upp)), border = F, col = 'gray')
    lines(x, fmean, type = 'l', lty = 4, lwd = 2, col = 'blue')
    lines(x, MAP, type = 'l', lty = 2, lwd = 2, col = 'red')
    points(x, y, pch = '*')
    abline(h = 0, lty = 2, lwd = 2)
  }
  return(list("time" = tm, "xi_sam" = xi_sam, "sig_sam" = sig_sam, "tau_sam" = tau_sam, "nu_sam" = nu_sam,
              "ell_sam" = ell_sam, "fhat_sam" = fhat_sam, "fmean" = fmean, "f_low" = f_low, "f_upp" = f_upp, "z_star" = z_star))
}
####################################################



#---------------------------------------------------
### For nonnegative functions estimation HMC sampler
#---------------------------------------------------
pos.HMC <- function(y, x, N, nu, l, mcmc, brn, thin, tau.in, sig.in, xi.in,
                    tau.fix, sig.fix, xi.fix, sseed, verbose, return.plot, tol){
  # y: Response variable; x: vector to form design matrix X (n x N)
  # N: the number of knots
  # nu:smoothness parameter of Matern; l:length-scale parameter of Matern
  # mcmc, brn, thin : mcmc samples, burning and thinning for MCMC
  # tau.in, sig.in, xi.in : initial values (supplied by user or use the default values)
  # tau.fix,sig.fix,xi.fix : if fixed values of the parameters are to use
  # verbose : logical; if TRUE, prints current status; default is TRUE
  # return.plot : logical; if true a plot of estimate with 95% CI is returned; default is TRUE
  
  # OUTPUT: Posterior samples on xi,tau,sig and fhat with posterior mean, 95% CI of fhat, mode of posterior distribution z_star
  
  if (length(y) != length(x))
    stop("y and x should be of same length!")
  y <- y[order(x)]
  x <- sort(x)
  n <- length(y)
  if (missing(N))
    N <- ceiling(n / 2) - 1
  
  delta <- 1 / (N - 1)
  my_knots <- seq(from = 0, to = 1, by = delta)
  X <- fcth(x = x, u = my_knots, N = N)
  
  if (missing(nu))
    stop("nu needs to be supplied")
  if (nu == 0)
    stop("nu cannot be zero")
  if (!missing(l) && l == 0)
    stop("l cannot be zero")
  if (missing(return.plot))
    return.plot <- TRUE
  if (!missing(sseed))
    set.seed(sseed)
  if (missing(sseed))
    set.seed(Sys.Date())
  if (missing(verbose))
    verbose <- TRUE
  if (missing(l))
    l <- l_est(nu = nu, range = range(my_knots), val = 0.05)
  
  # prior covariance K:
  K <- covmat(knot = my_knots, nu = nu, l = l)
  # prior precision:
  K_inv <- tinv(K)
  
  if (missing(tol))
    tol <- 1e-6
  if (missing(mcmc))
    mcmc <- 5000
  if (missing(brn))
    brn <- 1000
  if (missing(thin))
    thin <- 1
  em <- mcmc + brn
  ef <- mcmc / thin
  
  if (!missing(tau.fix))
    tau.in <- tau.fix
  if (!missing(sig.fix))
    sig.in <- sig.fix
  if (!missing(xi.fix))
    xi.in <- xi.fix
  if (missing(tau.fix) && missing(tau.in))
    tau.in <- 1
  if (missing(sig.fix) && missing(sig.in))
    sig.in <- 1
  if (missing(xi.fix) && missing(xi.in)){
    ## posterior Mode
    XXK <- crossprod(X) / sig.in + K_inv / tau.in
    Amat <- diag(N)
    z_star.in <- solve.QP(Dmat = XXK, dvec = as.vector(t(X) %*% y) / sig.in,
                          Amat = t(Amat), bvec = c(rep(0, N)), meq = 0)$solution
    xi.in <- z_star.in # starting point McMC
  }
  tau <- tau.in
  sig <- sig.in
  xi_in <- xi.in
  
  xi_sam <- matrix(NA, nrow = N, ncol = ef)
  tau_sam <- rep(NA, ef)
  sig_sam <- rep(NA, ef)
  fhat_sam <- matrix(NA, nrow = n, ncol = ef)
  f <- diag(N)
  g <- rep(tol, N)
  
  if (verbose)
    print("MCMC sample draws:")
  
  ptm <- proc.time()
  for (i in 1 : em) {
    # sampling Xi:
    if (missing(xi.fix)) {
      M <- crossprod(X) / sig + K_inv / tau
      r <- as.vector(t(X) %*% y) / sig
      xi_out <- as.vector(rtmg(n = 1, M = M, r = r, initial = xi_in, f = f, g = g, burn.in = 0))
    } else {
      xi_out <- xi_in
    }
    set.seed(2 * i)
    # sampling \sigma^2:
    Xxi <- as.vector(X %*% xi_out)
    y_star <- y - Xxi
    if (missing(sig.fix))
      sig <- 1 / rgamma(1, shape = n / 2, rate = sum(y_star^2) / 2)
    
    # sampling \tau^2:
    if (missing(tau.fix))
      tau <- 1 / rgamma(1, shape = N / 2, rate = (t(xi_out) %*% K_inv %*% xi_out) / 2)
    
    # storing MCMC samples:
    if (i > brn && i %% thin == 0) {
      xi_sam[, (i - brn) / thin] <- xi_out
      sig_sam[(i - brn) / thin] <- sig
      tau_sam[(i - brn) / thin] <- tau
      fhat_sam[, (i - brn) / thin] <- Xxi
    }
    if (i %% 1000 == 0 && verbose) {
      print(i)
    }
    
    # renewing the intial value:
    xi_in <- xi_out
  }
  tm <- proc.time() - ptm
  
  ## posterior Mode
  XXK <- crossprod(X) / mean(sig_sam) + K_inv / mean(tau_sam)
  Amat <- diag(N)
  z_star <- solve.QP(Dmat = XXK, dvec = as.vector(t(X) %*% y) / mean(sig_sam),
                     Amat = t(Amat), bvec = c(rep(0, N)), meq = 0)$solution
  
  MAP <- X %*% z_star # MAP estimate
  fmean <- rowMeans(fhat_sam) # mAP estimate
  qnt <- apply(fhat_sam, 1, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))
  f_low <- qnt[1, ]
  f_upp <- qnt[2, ]
  ub <- max(f_low, f_upp, fmean, MAP)
  lb <- min(f_low, f_upp, fmean, MAP)
  
  if (return.plot) {
    par(mfrow = c(1, 1))
    par(mar = c(2.1, 2.1, 2.1, 1.1)) # adapt margins
    plot(x, y, pch = '*', lwd = 2, lty = 1, col = 'black',
         ylim = range(ub, lb, y), xlab = '', ylab = '')
    polygon(c(x, rev(x)), y = c(f_low, rev(f_upp)), border = F, col = 'gray')
    lines(x, fmean, type = 'l', lty = 4, lwd = 2, col = 'blue')
    lines(x, MAP, type = 'l', lty = 2, lwd = 2, col = 'red')
    points(x, y, pch = '*')
    abline(h = 0, lty = 2, lwd = 2)
  }
  return(list("time" = tm, "xi_sam" = xi_sam, "sig_sam" = sig_sam, "tau_sam" = tau_sam,
              "fhat_sam" = fhat_sam, "fmean" = fmean, "f_low" = f_low, "f_upp" = f_upp, "z_star" = z_star))
}
####################################################





#---------------------------------------------------
######## monotonicity constraints models ###########
#---------------------------------------------------
### Function for drawing posterior samples using ESS and LS with fixed hyperparameters:
mon.inc.LS.ESS <- function(y, x, N1, M, nu, l, eta, mcmc, brn, thin, tau.in, sig.in, xi0.in, xi.in, xi0.fix,
                           tau.fix, sig.fix, xi.fix, verbose, return.plot, tol, sseed){
  # y: Response variable; x: vector to form design matrix X (n x N)
  # N1: the number of knots of the first subdomain
  # M: nb of subdomain
  # nu:smoothness parameter of Matern; l:length-scale parameter of Matern
  # tau.in, sig.in, xi0.in, xi.in : initial values (supplied by user or the defaul values)
  # eta:parameter of the approximation function of the indicator functions
  # mcmc, brn, thin : mcmc samples, burning and thinning for MCMC
  # xi0.fix,tau.fix,sig.fix,xi.fix : if fixed values of the parameters are to use
  # verbose : logical; if TRUE, prints currnet status; default is TRUE
  # return.plot : logical; if true a plot of estimate with 95% CI is returned; default is TRUE
  # tol: tolerance for numerical stability
  
  # OUTPUT: Posterior samples on xi,xi0,tau,sig and fhat with posterior mean, 95% CI of fhat, z_star: mode posterior distribution
  
  if (length(y) != length(x))
    stop("y and x should be of same length!")
  y <- y[order(x)]
  x <- sort(x)
  n <- length(y)
  N <- N1 * M
  delta <- 1 / (N - 1)
  my_knots <- seq(from = 0, to = 1, by = delta)
  X <- fctphi(x, u = my_knots, N = N)
  
  if (missing(nu))
    stop("nu needs to be supplied")
  if (nu == 0)
    stop("nu cannot be zero")
  if (!missing(l) && l == 0)
    stop("l cannot be zero")
  if (missing(return.plot))
    return.plot <- TRUE
  if (!missing(sseed))
    set.seed(sseed)
  if (missing(sseed))
    set.seed(Sys.Date())
  if (missing(verbose))
    verbose <- TRUE
  if (missing(l))
    l <- l_est(nu = nu, range = range(my_knots), val = 0.05)
  
  # prior covariance K:
  K <- covmat(knot = my_knots, nu = nu, l = l)
  # prior precision:
  K_inv <- tinv(K)
  
  if (missing(tol))
    tol <- 1e-6
  if (missing(eta))
    eta <- 50
  if (missing(mcmc))
    mcmc <- 5000
  if (missing(brn))
    brn <- 1000
  if (missing(thin))
    thin <- 1
  em <- mcmc + brn
  ef <- mcmc / thin
  
  if (!missing(tau.fix))
    tau.in <- tau.fix
  if (!missing(sig.fix))
    sig.in <- sig.fix
  if (!missing(xi0.fix))
    xi0.in <- xi0.fix
  if (!missing(xi.fix))
    xi.in <- xi.fix
  if (missing(tau.fix) && missing(tau.in))
    tau.in <- 1
  if (missing(sig.fix) && missing(sig.in))
    sig.in <- 1
  if (missing(xi0.fix) && missing(xi0.in))
    xi0.in <- 0
  if (missing(xi.fix) && missing(xi.in)) {
    xi.in <- pmax(0, mvrnorm(1, rep(0, N), K))
  }
  
  tau <- tau.in
  sig <- sig.in
  xi0 <- xi0.in
  xi_in <- xi.in
  
  xi_sam <- matrix(NA, nrow = N, ncol = ef)
  xi0_sam <- rep(NA, ef)
  tau_sam <- rep(NA, ef)
  sig_sam <- rep(NA, ef)
  fhat_sam <- matrix(NA, nrow = n, ncol = ef)
  
  if (verbose)
    print("MCMC sample draws:")
  
  ptm  <- proc.time()
  for (i in 1 : em) { 
    # sampling Xi:
    y_tilde <- y - xi0
    if (missing(xi.fix)) {
      nu.ess <- Fast.LS(u = my_knots, M = M, N1 = N1, nu = nu, l = l, tausq = tau, tol = tol, sseedLS = i)
      xi_out <- ESS(beta = xi_in, nu_ess = nu.ess, y = y_tilde, X = X, sigsq = sig, eta = eta, seeds = i)
      xi_out <- sapply(xi_out, function(z) return(pmax(0, z)))
    } else {
      xi_out <- xi_in
    }
    set.seed(2 * i)
    # sampling xi_0:
    Xxi <- as.vector(X %*% xi_out)
    y_star <- y - Xxi
    if (missing(xi0.fix))
      xi0 <- rnorm(n = 1, mean = mean(y_star), sd = sqrt(sig/n))
    
    # sampling \sigma^2:
    y0 <- y_star - xi0
    if (missing(sig.fix))
      sig <- 1 / rgamma(1, shape = n / 2, rate = sum(y0^2) / 2)
    
    # sampling \tau^2:
    if (missing(tau.fix))
      tau <- 1 / rgamma(1, shape = N/2, rate = (t(xi_out) %*% K_inv %*% xi_out) / 2)
    
    # storing MCMC samples:
    if (i > brn && i %% thin == 0) {
      xi_sam[, (i - brn) / thin] <- xi_out
      xi0_sam[(i - brn) / thin] <- xi0
      sig_sam[(i - brn) / thin] <- sig
      tau_sam[(i - brn) / thin] <- tau
      fhat_sam[, (i - brn) / thin] <- xi0 + Xxi
    }
    if (i %% 1000 == 0 && verbose) {
      print(i)
    }
    
    # renewing the intial value:
    xi_in <- xi_out
  }
  tm <- proc.time() - ptm
  
  ## Posterior Mode
  Xp <- cbind(rep(1, n), X)
  K <- matrix(data = NA, nrow = N + 1, ncol = N + 1)
  K[1, 1] <- k(0, nu = nu, l = l)
  K[1, 2 : (N + 1)] <- rep(0, N)
  K[2 : (N + 1), 1] <- rep(0, N)
  K[2 : (N + 1), 2 : (N + 1)] <- covmat(knot = my_knots, nu = nu, l = l)
  K_invp <- bdiag(1 / K[1, 1], K_inv)
  XXK <- crossprod(Xp) / mean(sig_sam) + K_invp / mean(tau_sam)
  Amat <- diag(N + 1)
  Amat <- Amat[-1, ]
  z_star <- solve.QP(Dmat = XXK, dvec = as.vector(t(Xp) %*% y) / mean(sig_sam),
                     Amat = t(Amat), bvec = c(rep(0, N)), meq = 0)$solution
  MAP <- Xp %*% z_star # MAP estimate
  fmean <- rowMeans(fhat_sam) # mAP estimate
  qnt <- apply(fhat_sam, 1, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))
  f_low <- qnt[1, ]
  f_upp <- qnt[2, ]
  ub <- max(f_low, f_upp, fmean, MAP)
  lb <- min(f_low, f_upp, fmean, MAP)
  
  if (return.plot) {
    par(mfrow = c(1, 1))
    par(mar = c(2.1, 2.1, 2.1, 1.1)) # adapt margins
    plot(x, y, pch = '*', lwd = 2, lty = 1, col = 'black',
         ylim = range(ub, lb, y), xlab = '', ylab = '')
    polygon(c(x, rev(x)), y = c(f_low, rev(f_upp)), border = F, col = 'gray')
    lines(x, fmean, type = 'l', lty = 4, lwd = 2, col = 'blue')
    lines(x, MAP, type = 'l', lty = 2, lwd = 2, col = 'red')
    points(x, y, pch = '*')
  }
  
  return(list("time" = tm, "xi_sam" = xi_sam, "xi0_sam" = xi0_sam, "sig_sam" = sig_sam, "tau_sam" = tau_sam,
              "fhat_sam" = fhat_sam, "fmean" = fmean, "f_low" = f_low, "f_upp" = f_upp, "z_star" = z_star))
}
####################################################


#---------------------------------------------------
### Function for drawing posterior samples using ESS and WC with fixed hyperparameters
#---------------------------------------------------
mon.inc.WC.ESS <- function(y, x, N, nu, l, eta, mcmc, brn, thin, tau.in, sig.in, xi0.in, xi.in, xi0.fix,
                           tau.fix, sig.fix, xi.fix, sseed, verbose, return.plot){
  # y: Response variable; x: vector to form design matrix X (n x N)
  # N: the number of knots
  # nu:smoothness parameter of Matern; l:length-scale parameter of Matern
  # eta:parameter of the approximation function of the indicator functions
  # mcmc, brn, thin : mcmc samples, burning and thinning for MCMC
  # tau.in, sig.in, xi0.in, xi.in : initial values (supplied by user or use the default values)
  # xi0.fix,tau.fix,sig.fix,xi.fix : if fixed values of the parameters are to use
  # verbose : logical; if TRUE, prints currnet status; default is TRUE
  # return.plot : logical; if true a plot of estimate with 95% CI is returned; default is TRUE
  
  # OUTPUT: Posterior samples on xi,xi0,tau,sig and fhat with posterior mean, 95% CI of fhat, z_star: mode posterior distribution
  
  if (length(y) != length(x))
    stop("y and x should be of same length!")
  n <- length(y)
  y <- y[order(x)]
  x <- sort(x)
  if (missing(N))
    N <- ceiling(n / 2) - 1
  
  delta <- 1 / (N - 1)
  my_knots  <- seq(from = 0, to = 1, by = delta)
  X <- fctphi(x = x, u = my_knots, N = N)
  
  if (missing(nu))
    stop("nu needs to be supplied")
  if (nu == 0)
    stop("nu cannot be zero")
  if (!missing(l) && l == 0)
    stop("l cannot be zero")
  if (missing(return.plot))
    return.plot <- TRUE
  if (!missing(sseed))
    set.seed(sseed)
  if (missing(sseed))
    set.seed(Sys.Date())
  if (missing(verbose))
    verbose <- TRUE
  if (missing(l))
    l <- l_est(nu = nu, range = range(my_knots), val = 0.05)
  
  # prior covariance K:
  K <- covmat(knot = my_knots, nu = nu, l = l)
  # prior precision:
  K_inv <- tinv(K)
  
  if (missing(eta))
    eta <- 50
  if (missing(mcmc))
    mcmc <- 5000
  if (missing(brn))
    brn <- 1000
  if (missing(thin))
    thin <- 1
  em <- mcmc + brn
  ef <- mcmc / thin
  
  if (!missing(tau.fix))
    tau.in <- tau.fix
  if (!missing(sig.fix))
    sig.in <- sig.fix
  if (!missing(xi0.fix))
    xi0.in <- xi0.fix
  if (!missing(xi.fix))
    xi.in <- xi.fix
  if (missing(tau.fix) && missing(tau.in))
    tau.in <- 1
  if (missing(sig.fix) && missing(sig.in))
    sig.in <- 1
  if (missing(xi0.fix) && missing(xi0.in))
    xi0.in <- 0
  if (missing(xi.fix) && missing(xi.in)) {
    xi.in <- pmax(0, mvrnorm(1, rep(0, N), K))
  }
  
  tau <- tau.in
  sig <- sig.in
  xi0 <- xi0.in
  xi_in <- xi.in
  
  xi_sam <- matrix(NA, nrow = N, ncol = ef)
  xi0_sam <- rep(NA, ef)
  tau_sam <- rep(NA, ef)
  sig_sam <- rep(NA, ef)
  fhat_sam <- matrix(NA, nrow = n, ncol = ef)
  
  if (verbose)
    print("MCMC sample draws:")
  
  ptm <- proc.time()
  for (i in 1 : em) {
    # sampling Xi:
    y_tilde <- y - xi0
    if (missing(xi.fix)) {
      nu.ess <- as.vector(samp.WC(knot = my_knots, nu = nu, l = l, tausq = tau, sseedWC = i))
      xi_out <- ESS(beta = xi_in, nu_ess = nu.ess, y = y_tilde, X = X, sigsq = sig, eta = eta, seeds = i)
      xi_out <- sapply(xi_out, function(z) return(pmax(0, z)))
    } else {
      xi_out <- xi_in
    }
    set.seed(2 * i)
    # sampling xi_0:
    Xxi <- as.vector(X %*% xi_out)
    y_star <- y - Xxi
    if (missing(xi0.fix))
      xi0 <- rnorm(n = 1, mean = mean(y_star), sd = sqrt(sig / n))
    
    # sampling \sigma^2:
    y0 <- y_star - xi0
    if (missing(sig.fix))
      sig <- 1 / rgamma(1, shape = n / 2, rate = sum(y0^2) / 2)
    
    # sampling \tau^2:
    if (missing(tau.fix))
      tau <- 1 / rgamma(1, shape = N / 2, rate = (t(xi_out) %*% K_inv %*% xi_out) / 2)
    
    # storing MCMC samples:
    if (i > brn && i %% thin == 0) {
      xi_sam[, (i - brn) / thin] <- xi_out
      xi0_sam[(i - brn) / thin] <- xi0
      sig_sam[(i - brn) / thin] <- sig
      tau_sam[(i - brn) / thin] <- tau
      fhat_sam[, (i - brn) / thin] <- xi0 + Xxi
    }
    
    if (i %% 1000 == 0 && verbose) {
      print(i)
    }
    
    # renewing the intial value:
    xi_in <- xi_out
  }
  tm <- proc.time() - ptm
  
  ## Posterior Mode
  X <- cbind(rep(1, n), X)
  K <- matrix(data = NA, nrow = N + 1, ncol = N + 1)
  K[1, 1] <- k(0, nu = nu, l = l)
  K[1, 2 : (N + 1)] <- rep(0, N)
  K[2 : (N + 1), 1] <- rep(0, N)
  K[2 : (N + 1), 2 : (N + 1)] <- covmat(knot = my_knots, nu = nu, l= l)
  K_inv <- solve(K)
  XXK <- crossprod(X) / mean(sig_sam) + K_inv / mean(tau_sam)
  Amat <- diag(N + 1)
  Amat <- Amat[-1, ]
  z_star <- solve.QP(Dmat = XXK, dvec = as.vector(t(X) %*% y) / mean(sig_sam),
                     Amat = t(Amat), bvec = c(rep(0, N)), meq = 0)$solution
  MAP <- X %*% z_star # MAP estimate
  fmean <- rowMeans(fhat_sam) # mAP estimate
  qnt <- apply(fhat_sam, 1, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))
  f_low <- qnt[1, ]
  f_upp <- qnt[2, ]
  ub <- max(f_low, f_upp, fmean, MAP)
  lb <- min(f_low, f_upp, fmean, MAP)
  
  if (return.plot) {
    par(mfrow = c(1, 1))
    par(mar = c(2.1, 2.1, 2.1, 1.1)) # adapt margins
    plot(x, y, pch = '*', lwd = 2, lty = 1, col = 'black',
         ylim = range(ub, lb, y), xlab = '', ylab = '')
    polygon(c(x, rev(x)), y = c(f_low, rev(f_upp)), border = F, col = 'gray')
    lines(x, fmean, type = 'l', lty = 4, lwd = 2, col = 'blue')
    lines(x, MAP, type = 'l', lty = 2, lwd = 2, col = 'red')
    points(x, y, pch = '*')
  }
  
  return(list("time" = tm, "xi_sam" = xi_sam, "xi0_sam" = xi0_sam, "sig_sam" = sig_sam, "tau_sam" = tau_sam,
              "fhat_sam" = fhat_sam, "fmean" = fmean, "f_low" = f_low, "f_upp" = f_upp, "z_star" = z_star))
}
####################################################




#---------------------------------------------------
### Function for drawing posterior samples using HMC with fixed hyperparameters:
#---------------------------------------------------
mon.HMC  <- function(y, x, N, nu, l, mcmc, brn, thin, tau.in, sig.in, xi.in, xi0.in,
                     tau.fix, sig.fix, xi.fix, xi0.fix, sseed, verbose, return.plot, tol){
  # y: Response variable; x: vector to form design matrix X (n x N)
  # N: number of knots 
  # nu: smoothness parameter of Matern; l:length-scale parameter of Matern
  # mcmc, brn, thin : mcmc samples, burning and thinning for MCMC
  # tau.in, sig.in, xi.in : initial values (supplied by user or use the default values)
  # tau.fix,sig.fix,xi.fix : if fixed values of the parameters are to use
  # verbose : logical; if TRUE, prints current status; default is TRUE
  # return.plot : logical; if true a plot of estimate with 95% CI is returned; default is TRUE
  
  # OUTPUT: Posterior samples on xi,tau,sig and fhat with posterior mean, 95% CI of fhat, mode of posterior distribution z_star
  
  if (length(y) != length(x))
    stop("y and x should be of same length!")
  n <- length(y)
  y <- y[order(x)]
  x <- sort(x)
  if (missing(N))
    N <- ceiling(n / 2) - 1
  
  delta <- 1 / (N - 1)
  my_knots <- seq(from = 0, to = 1, by = delta)
  X <- fctphi(x = x, u = my_knots, N = N)
  
  if (missing(nu))
    stop("nu needs to be supplied")
  if (nu == 0)
    stop("nu cannot be zero")
  if (!missing(l) && l == 0)
    stop("l cannot be zero")
  if (missing(return.plot))
    return.plot <- TRUE
  if (!missing(sseed))
    set.seed(sseed)
  if (missing(sseed))
    set.seed(Sys.Date())
  if (missing(verbose))
    verbose <- TRUE
  if (missing(l))
    l <- l_est(nu = nu, range = range(my_knots), val = 0.05)
  
  # prior covariance K:
  K <- covmat(knot = my_knots, nu = nu, l = l)
  # prior precision:
  K_inv <- tinv(K)
  
  if (missing(tol))
    tol <- 1e-6
  if (missing(mcmc))
    mcmc <- 5000
  if (missing(brn))
    brn <- 1000
  if (missing(thin))
    thin <- 1
  em <- mcmc + brn
  ef <- mcmc / thin
  
  if (!missing(tau.fix))
    tau.in <- tau.fix
  if (!missing(sig.fix))
    sig.in <- sig.fix
  if (!missing(xi0.fix))
    xi0.in <- xi0.fix
  if (!missing(xi.fix))
    xi.in <- xi.fix
  
  if (missing(tau.fix) && missing(tau.in))
    tau.in <- 1
  if (missing(sig.fix) && missing(sig.in))
    sig.in <- 1
  if (missing(xi0.fix) && missing(xi0.in))
    xi0.in <- 0
  if (missing(xi.fix) && missing(xi.in))
    xi.in <- rep(1, N)
  
  tau <- tau.in
  sig <- sig.in
  xi0 <- xi0.in
  xi_in <- xi.in
  
  xi_sam <- matrix(NA, nrow = N, ncol = ef)
  xi0_sam <- rep(NA, ef)
  tau_sam <- rep(NA, ef)
  sig_sam <- rep(NA, ef)
  fhat_sam <- matrix(NA,nrow = n, ncol = ef)
  ## matrix and vector of constr, respectively
  f <- diag(N)
  g <- rep(tol, N)
  
  if (verbose)
    print("MCMC sample draws:")
  
  ptm <- proc.time()
  for (i in 1 : em) {
    y_tilde <- y - xi0
    # sampling Xi:
    if (missing(xi.fix)) {
      M <- crossprod(X) / sig + K_inv / tau
      r <- as.vector(t(X) %*% y_tilde) / sig
      xi_out <- as.vector(rtmg(n = 1, M = M, r = r, initial = xi_in, f = f, g = g, burn.in = 0))
    } else {
      xi_out <- xi_in
    }
    set.seed(2 * i)
    Xxi <- as.vector(X %*% xi_out)
    y_star <- y - Xxi
    if (missing(xi0.fix))
      xi0 <- rnorm(n = 1, mean = mean(y_star), sd = sqrt(sig / n))
    
    # sampling \sigma^2:
    if (missing(sig.fix))
      y0 <- y_star - xi0
    sig <- 1 / rgamma(1, shape = n / 2, rate = sum(y0^2) / 2)
    
    # sampling \tau^2:
    if (missing(tau.fix))
      tau <- 1 / rgamma(1, shape = N / 2, rate = (t(xi_out) %*% K_inv %*% xi_out) / 2)
    
    # storing MCMC samples:
    if (i > brn && i %% thin == 0) {
      xi_sam[, (i - brn) / thin] <- xi_out
      xi0_sam[(i - brn) / thin] <- xi0
      sig_sam[(i - brn) / thin] <- sig
      tau_sam[(i - brn) / thin] <- tau
      fhat_sam[, (i - brn) / thin] <- xi0 + Xxi
    }
    
    if (i %% 1000 == 0 && verbose) {
      print(i)
    }
    
    # renewing the intial value:
    xi_in <- xi_out
  }
  tm <- proc.time() - ptm
  
  ## Posterior Mode
  X <- cbind(rep(1, n), X)
  K <- matrix(data = NA, nrow = N + 1, ncol = N + 1)
  K[1, 1] <- k(h = 0, nu = nu, l = l)
  K[1, 2 : (N + 1)] <- rep(0, N)
  K[2 : (N + 1), 1] <- rep(0, N)
  K[2 : (N + 1), 2 : (N + 1)] <- covmat(knot = my_knots, nu = nu, l = l)
  K_inv <- solve(K)
  XXK <- crossprod(X) / mean(sig_sam) + K_inv / mean(tau_sam)
  Amat <- diag(N + 1)
  Amat <- Amat[-1, ]
  z_star <- solve.QP(Dmat = XXK, dvec = as.vector(t(X) %*% y) / mean(sig_sam),
                     Amat = t(Amat), bvec = c(rep(0, N)), meq = 0)$solution
  MAP <- X %*% z_star # MAP estimate
  fmean <- rowMeans(fhat_sam) # mAP estimate
  qnt <- apply(fhat_sam, 1, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))
  f_low <- qnt[1, ]
  f_upp <- qnt[2, ]
  ub <- max(f_low, f_upp, fmean, MAP)
  lb <- min(f_low, f_upp, fmean, MAP)
  
  if (return.plot) {
    par(mfrow = c(1, 1))
    par(mar = c(2.1, 2.1, 2.1, 1.1)) # adapt margins
    plot(x, y, pch = '*', lwd = 2, lty = 1, col = 'black',
         ylim = range(ub, lb, y), xlab = '', ylab = '')
    polygon(c(x, rev(x)), y = c(f_low, rev(f_upp)), border = F, col = 'gray')
    lines(x, fmean, type = 'l', lty = 4, lwd = 2, col = 'blue')
    lines(x, MAP, type = 'l', lty = 2, lwd = 2, col = 'red')
    points(x, y, pch = '*')
  }
  
  return(list("time" = tm, "xi_sam" = xi_sam, "sig_sam" = sig_sam, "tau_sam" = tau_sam,
              "fhat_sam" = fhat_sam, "fmean" = fmean, "f_low" = f_low, "f_upp" = f_upp))
}
####################################################







#---------------------------------------------------
### Function for drawing posterior samples using ESS and LS with hyperparameter updates:
#---------------------------------------------------

mon.inc.LS.ESS_hyp <- function(y, x, N1, M, nu.in, l.in, eta, mcmc, brn, thin, tau.in, sig.in, xi0.in, xi.in, xi0.fix,
                               tau.fix, sig.fix, xi.fix, verbose, return.plot, tol, sseed){
  # y: Response variable; x: vector to form design matrix \phi (n X N)
  # N1: number of knots first subdomain; M: nb of subdomain
  # eta: parameter of the spproximation function of the indicator functions
  # mcmc, brn, thin : mcmc samples, burning and thinning for MCMC
  # nu.in, l.in, tau.in, sig.in, xi0.in, xi.in : initial values (supplied by user or the defaul values)
  # xi0.fix, tau.fix, sig.fix, xi.fix : if fixed values of the parameters are to use
  # verbose : logical; if TRUE, prints currnet status; default is TRUE
  # return.plot : logical; if true a plot of estimate with 95% CI is returned; default is TRUE
  # tol: tolerance for numerical stability
  
  # OUTPUT: Posterior samples on xi,xi0,tau,sig and fhat with posterior mean, 95% CI of fhat, z_star: mode posterior distribution
  
  if (length(y) != length(x))
    stop("y and x should be of same length!")
  y <- y[order(x)]
  x <- sort(x)
  n <- length(y)
  N <- N1 * M
  delta <- 1 / (N - 1)
  my_knots <- seq(from = 0, to = 1, by = delta)
  X <- fctphi(x = x, u = my_knots, N = N)
  
  if (missing(return.plot))
    return.plot=TRUE
  if (!missing(sseed))
    set.seed(sseed)
  if (missing(sseed))
    set.seed(Sys.Date())
  if (missing(verbose))
    verbose <- TRUE
  if (missing(eta))
    eta <- 50
  if (missing(mcmc))
    mcmc <- 5000
  if (missing(brn))
    brn <- 1000
  if (missing(thin))
    thin <- 1
  em <- mcmc + brn
  ef <- mcmc / thin
  
  if (!missing(tau.fix))
    tau.in <- tau.fix
  if (!missing(sig.fix))
    sig.in <- sig.fix
  if (!missing(xi0.fix))
    xi0.in <- xi0.fix
  if (!missing(xi.fix))
    xi.in <- xi.fix
  if (missing(nu.in))
    nu.in <- 0.75
  if (missing(l.in))
    l.in <- l_est(nu = nu.in, range = range(my_knots), val = 0.05)
  if (missing(tau.fix) && missing(tau.in))
    tau.in <- 1
  if (missing(sig.fix) && missing(sig.in))
    sig.in <- 1
  if (missing(xi0.fix) && missing(xi0.in))
    xi0.in <- 0
  if (missing(xi.fix) && missing(xi.in)){
    # prior covariance K:
    K <- covmat(knot = my_knots, nu = nu.in, l = l.in)
    xi.in <- pmax(0, mvrnorm(1, rep(0, N), K))
  }
  
  tau <- tau.in
  sig <- sig.in
  xi0 <- xi0.in
  xi_in <- xi.in
  nu_in <- nu.in
  l_in <- l.in
  
  xi_sam <- matrix(NA, nrow = N, ncol = ef)
  xi0_sam <- rep(NA, ef)
  tau_sam <- rep(NA, ef)
  sig_sam <- rep(NA, ef)
  nu_sam <- rep(NA, ef)
  ell_sam <- rep(NA, ef)
  fhat_sam <- matrix(NA, nrow = n, ncol = ef)
  
  if (verbose)
    print("MCMC sample draws:")
  
  ptm <- proc.time()
  for (i in 1 : em) {
    # sampling from \nu and \ell
    MH.out <- nu.MH2(nu_in, l_in, tau, xi_in, my_knots)
    nu_out <- MH.out$nu
    l_out <- MH.out$l
    L_inv <- MH.out$L_inv
    
    # sampling Xi:
    y_tilde <- y - xi0
    if (missing(xi.fix)) {
      nu.ess <- Fast.LS(u = my_knots, M = M, N1 = N1, nu = nu_out, l = l_out, tausq = tau, tol = tol, sseedLS = i)
      xi_out <- ESS(beta = xi_in, nu_ess = nu.ess, y = y_tilde, X= X, sigsq = sig, eta = eta, seeds = i)
      xi_out <- sapply(xi_out, function(z) return(pmax(0, z)))
    } else {
      xi_out <- xi_in
    }
    set.seed(2 * i)
    # sampling xi_0:
    Xxi <- as.vector(X %*% xi_out)
    y_star <- y - Xxi
    if (missing(xi0.fix))
      xi0 <- rnorm(n = 1, mean = mean(y_star), sd = sqrt(sig / n))
    
    # sampling \sigma^2:
    y0 <- y_star - xi0
    if (missing(sig.fix))
      sig <- 1 / rgamma(1, shape = n / 2, rate = sum(y0^2) / 2)
    
    # sampling \tau^2:
    if (missing(tau.fix))
      tau <- 1 / rgamma(1, shape = N / 2, rate = (sum((t(L_inv) %*% xi_out)^2)) / 2)
    
    # storing MCMC samples:
    if (i > brn && i %% thin == 0) {
      xi_sam[, (i - brn) / thin] <- xi_out
      xi0_sam[(i - brn) / thin] <- xi0
      sig_sam[(i - brn) / thin] <- sig
      tau_sam[(i - brn) / thin] <- tau
      nu_sam[(i - brn) / thin] <- nu_out
      ell_sam[(i - brn) / thin] <- l_out
      fhat_sam[, (i - brn) / thin] <- xi0 + Xxi
    }
    
    if (i %% 1000 == 0 && verbose) {
      print(i)
    }
    
    # renewing the intial value:
    xi_in <- xi_out
    nu_in <- nu_out
    l_in <- l_out
  }
  tm <- proc.time() - ptm
  
  ## Posterior Mode
  X <- cbind(rep(1, n), X)
  K <- matrix(data = NA, nrow = N + 1, ncol = N + 1)
  K[1, 1] <- k(h = 0,nu = mean(nu_sam), l = mean(ell_sam))
  K[1, 2 : (N + 1)] <- rep(0, N)
  K[2 : (N + 1), 1] <- rep(0, N)
  K[2 : (N + 1), 2 : (N + 1)] <- covmat(knot = my_knots, nu = mean(nu_sam), l = mean(ell_sam))
  K_inv <- solve(K)
  XXK <- crossprod(X) / mean(sig_sam) + K_inv / mean(tau_sam)
  Amat <- diag(N + 1)
  Amat <- Amat[-1, ]
  z_star <- solve.QP(Dmat = XXK, dvec = as.vector(t(X) %*% y) / mean(sig_sam),
                     Amat = t(Amat), bvec = c(rep(0, N)), meq = 0)$solution
  MAP <- X %*% z_star # MAP estimate
  fmean <- rowMeans(fhat_sam) # mAP estimate
  qnt <- apply(fhat_sam, 1, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))
  f_low <- qnt[1, ]
  f_upp <- qnt[2, ]
  ub <- max(f_low, f_upp, fmean, MAP)
  lb <- min(f_low, f_upp, fmean, MAP)
  
  if (return.plot) {
    par(mfrow = c(1, 1))
    par(mar = c(2.1, 2.1, 2.1, 1.1)) # adapt margins
    plot(x, y,pch = '*', lwd = 2, lty = 1, col = 'black',
         ylim = range(ub, lb, y), xlab = '', ylab = '')
    polygon(c(x, rev(x)), y = c(f_low, rev(f_upp)), border = F, col = 'gray')
    lines(x, fmean, type = 'l', lty = 4, lwd = 2, col = 'blue')
    lines(x, MAP, type = 'l', lty = 2, lwd = 2, col = 'red')
    points(x, y, pch = '*')
  }
  
  return(list("time" = tm, "xi_sam" = xi_sam, "xi0_sam" = xi0_sam, "sig_sam" = sig_sam, "tau_sam" = tau_sam,
              "fhat_sam" = fhat_sam, "fmean" = fmean, "f_low" = f_low, "f_upp" = f_upp, "z_star" = z_star))
}






#---------------------------------------------------
### Function for drawing posterior samples using ESS and WC with hyperparameter updates:
#---------------------------------------------------
mon.inc.WC.ESS_hyp <- function(y, x, N, nu.in, l.in, eta, mcmc, brn, thin, tau.in, sig.in, xi0.in, xi.in, xi0.fix,
                               tau.fix, sig.fix, xi.fix, verbose, return.plot, sseed){
  # y: Response variable; x: vector to form design matrix \phi (n X N)
  # N: number of knots 
  # eta: parameter of the spproximation function of the indicator functions
  # mcmc, brn, thin : mcmc samples, burning and thinning for MCMC
  # nu.in, l.in, tau.in, sig.in, xi0.in, xi.in : initial values (supplied by user or the defaul values)
  # xi0.fix, tau.fix, sig.fix, xi.fix : if fixed values of the parameters are to use
  # verbose : logical; if TRUE, prints currnet status; default is TRUE
  # return.plot : logical; if true a plot of estimate with 95% CI is returned; default is TRUE
  
  # OUTPUT: Posterior samples on xi,xi0,tau,sig and fhat with posterior mean, 95% CI of fhat, z_star: mode posterior distribution
  
  if (length(y) != length(x))
    stop("y and x should be of same length!")
  y <- y[order(x)]
  x <- sort(x)
  n <- length(y)
  delta <- 1 / (N - 1)
  my_knots <- seq(from = 0, to = 1, by = delta)
  X <- fctphi(x = x, u = my_knots, N = N)
  
  if (missing(N))
    N <- ceiling(n / 2) - 1
  if (missing(return.plot))
    return.plot <- TRUE
  if (!missing(sseed))
    set.seed(sseed)
  if (missing(sseed))
    set.seed(Sys.Date())
  if (missing(verbose))
    verbose <- TRUE
  if (missing(eta))
    eta <- 50
  if (missing(mcmc))
    mcmc <- 5000
  if (missing(brn))
    brn <- 1000
  if (missing(thin))
    thin <- 1
  em <- mcmc + brn
  ef <- mcmc / thin
  
  if (!missing(tau.fix))
    tau.in <- tau.fix
  if (!missing(sig.fix))
    sig.in <- sig.fix
  if (!missing(xi0.fix))
    xi0.in <- xi0.fix
  if (!missing(xi.fix))
    xi.in <- xi.fix
  if (missing(nu.in))
    nu.in <- 0.75
  if (missing(l.in))
    l.in <- l_est(nu = nu.in, range = range(my_knots), val = 0.05)
  if (missing(tau.fix) && missing(tau.in))
    tau.in <- 1
  if (missing(sig.fix) && missing(sig.in))
    sig.in <- 1
  if (missing(xi0.fix) && missing(xi0.in))
    xi0.in <- 0
  if (missing(xi.fix) && missing(xi.in)){
    # prior covariance K:
    K <- covmat(knot = my_knots, nu = nu.in, l = l.in)
    xi.in <- pmax(0, mvrnorm(1, rep(0, N), K))
  }
  
  tau <- tau.in
  sig <- sig.in
  xi0 <- xi0.in
  xi_in <- xi.in
  nu_in <- nu.in
  l_in <- l.in
  
  xi_sam <- matrix(NA, nrow = N, ncol = ef)
  xi0_sam <- rep(NA, ef)
  tau_sam <- rep(NA, ef)
  sig_sam <- rep(NA, ef)
  nu_sam <- rep(NA, ef)
  ell_sam <- rep(NA, ef)
  fhat_sam <- matrix(NA, nrow = n, ncol = ef)
  
  if (verbose)
    print("MCMC sample draws:")
  
  ptm <- proc.time()
  for (i in 1 : em) {
    # sampling from \nu and \ell
    MH.out <- nu.MH2(nu_in, l_in, tau, xi_in, my_knots)
    nu_out <- MH.out$nu
    l_out <-MH.out$l
    L_inv <- MH.out$L_inv
    
    # sampling Xi:
    y_tilde <- y - xi0
    if (missing(xi.fix)) {
      nu.ess <- samp.WC(knot = my_knots, nu = nu_out, l = l_out, tausq = tau, sseedWC = i)
      xi_out <- ESS(beta = xi_in, nu_ess = nu.ess, y = y_tilde, X= X, sigsq = sig, eta = eta, seeds = i)
      xi_out <- sapply(xi_out, function(z) return(pmax(0, z)))
    } else {
      xi_out <- xi_in
    }
    set.seed(2 * i)
    # sampling xi_0:
    Xxi <- as.vector(X %*% xi_out)
    y_star <- y - Xxi
    if (missing(xi0.fix))
      xi0 <- rnorm(n = 1, mean = mean(y_star), sd = sqrt(sig / n))
    
    # sampling \sigma^2:
    y0 <- y_star - xi0
    if (missing(sig.fix))
      sig <- 1 / rgamma(1, shape = n / 2, rate = sum(y0^2) / 2)
    
    # sampling \tau^2:
    if (missing(tau.fix))
      tau <- 1 / rgamma(1, shape = N / 2, rate = (sum((t(L_inv) %*% xi_out)^2)) / 2)
    
    # storing MCMC samples:
    if (i > brn && i %% thin == 0) {
      xi_sam[, (i - brn) / thin] <- xi_out
      xi0_sam[(i - brn) / thin] <- xi0
      sig_sam[(i - brn) / thin] <- sig
      tau_sam[(i - brn) / thin] <- tau
      nu_sam[(i - brn) / thin] <- nu_out
      ell_sam[(i - brn) / thin] <- l_out
      fhat_sam[, (i - brn) / thin] <- xi0 + Xxi
    }
    
    if (i %% 1000 == 0 && verbose) {
      print(i)
    }
    
    # renewing the intial value:
    xi_in <- xi_out
    nu_in <- nu_out
    l_in <- l_out
  }
  tm <- proc.time() - ptm
  
  ## Posterior Mode
  X <- cbind(rep(1, n), X)
  K <- matrix(data = NA, nrow = N + 1, ncol = N + 1)
  K[1, 1] <- k(0,nu = mean(nu_sam), l = mean(ell_sam))
  K[1, 2 : (N + 1)] <- rep(0, N)
  K[2 : (N + 1), 1] <- rep(0, N)
  K[2 : (N + 1), 2 : (N + 1)] <- covmat(knot = my_knots, nu = mean(nu_sam), l = mean(ell_sam))
  K_inv <- solve(K)
  XXK <- crossprod(X) / mean(sig_sam) + K_inv / mean(tau_sam)
  Amat <- diag(N + 1)
  Amat <- Amat[-1, ]
  z_star <- solve.QP(Dmat = XXK, dvec = as.vector(t(X) %*% y) / mean(sig_sam),
                     Amat = t(Amat), bvec = c(rep(0, N)), meq = 0)$solution
  MAP <- X %*% z_star # MAP estimate
  fmean <- rowMeans(fhat_sam) # mAP estimate
  qnt <- apply(fhat_sam, 1, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))
  f_low <- qnt[1, ]
  f_upp <- qnt[2, ]
  ub <- max(f_low, f_upp, fmean, MAP)
  lb <- min(f_low, f_upp, fmean, MAP)
  
  if (return.plot) {
    par(mfrow = c(1, 1))
    par(mar = c(2.1, 2.1, 2.1, 1.1)) # adapt margins
    plot(x, y, pch = '*', lwd = 2, lty = 1, col = 'black',
         ylim = range(ub, lb, y), xlab = '', ylab = '')
    polygon(c(x, rev(x)), y = c(f_low, rev(f_upp)), border = F, col = 'gray')
    lines(x, fmean, type = 'l', lty = 4, lwd = 2, col = 'blue')
    lines(x, MAP, type = 'l', lty = 2, lwd = 2, col = 'red')
    points(x, y, pch = '*')
  }
  
  return(list("time" = tm, "xi_sam" = xi_sam, "xi0_sam" = xi0_sam, "sig_sam" = sig_sam, "tau_sam" = tau_sam,
              "fhat_sam" = fhat_sam, "fmean" = fmean, "f_low" = f_low, "f_upp" = f_upp, "z_star" = z_star))
}


## end               
