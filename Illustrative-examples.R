## Illustrative examples shown in
## the section `Numerical performance' of the article: Large-scale constrained GPs for shape-restricted function estimation
## Comparison between Fast.LS and Cholesky
## Comparison between Fast.LS and FFT Wood and Chan in different situations
## Comparison between Fast.LS and LS.KLE
## Illustration of correlation error due to the approximation


source('all_base_functions.R')
source('all_models.R')
library(ggplot2) ## plot zoom inset
library(grid) ## viewport function

## choice of numerical examples illustrated 
## in the article (put 'yes' between quotations)
illustr5split=''
CholvsLS=''
KLEvsLS_N=''
WCvsLS_N=''
WCvsLS_nu=''
WCvsLS_theta=''
LSvsLS.WC_N=''
LSvsLS.WC_N1=''
LS1M=''
correlation=''
monotoneC0Approx=''
comp_pos_syn_LS_WC=''
comp_monot_syn_LS_WC=''


if(illustr5split=='yes'){
  ###################################################
  ############# Illustration 5 split ################
  ###################################################
  par(mar=c(2.1,2.1,1.1,1.1)) # adapt margins
  N1 <- 50
  M <- 5  # nb of subdomains
  nu <- 1.5 # smoothness parameter Matern Kernel (MK)
  theta <- l_est(nu,c(0,1),0.05)#0.4 # length-scale parameter MK
  u <- seq(from = 0, to = 1, length = M * N1)
  U <- function(i){
    u[((i-1)*N1+1):(i*N1)]
  }
  u1 <- U(1)
  u2 <- U(2)
  u3 <- U(3)
  u4 <- U(4)
  u5 <- U(5)
  
  samp <- Fast.LS(u,M=M,N1=N1,nu=nu,l=theta,tausq=1,tol=1e-10,sseed=12345)
  plot(u,u,type='n',ylim=range(samp))
  # title(ylab='Y(x)',xlab='x',line=2)
  lines(u1,samp[1:N1],type='l',lwd=2,col='black')
  lines(u2,samp[(N1+1):(2*N1)],type='l',lwd=2,lty=2,col='gray')
  lines(u3,samp[(2*N1+1):(3*N1)],type='l',lwd=2,lty=2,col='black')
  lines(u4,samp[(3*N1+1):(4*N1)],type='l',lwd=2,lty=2,col='gray')
  lines(u5,samp[(4*N1+1):(5*N1)],type='l',lwd=2,lty=2,col='black')
  ###################################################
}


if(CholvsLS=='yes'){
  ###################################################
  ###### mvrnorm (Cholesky) versus Fast.LS ##########
  ###################################################
  par(mar=c(3.1,3.1,1.1,1.1)) # adapt margins
  nu <- 1.5 # smoothness parameter Matern Kernel (MK)
  theta <- l_est(nu,c(0,1),0.05) # length-scale parameter MK
  N <- seq(100,1000,length=10) # size of the MVN
  N1 <- 50
  M <- N/N1
  trial <- 25
  timeChol <- matrix(NA, nrow = trial, ncol = length(N))
  timeLS <- matrix(NA, nrow = trial, ncol = length(N))
  for(j in 1 : trial) {
    print(j)
    for(i in 1 : length(N)){
      u=seq(0,1,length=N[i])
      timeLS[j,i] <- system.time(Fast.LS(u=u,M=M[i],N1=N1,
                                         nu=nu,l=theta,tausq=1,tol=1e-7,sseed=j))[3]
      timeChol[j,i] <- system.time(mvtnorm::rmvnorm(1,rep(0,N[i]),
                                                    k(outer(u,u,'-'),
                                                      nu=nu,l=theta),method='chol'))[3]
    }
  } 
  averageLS <- colMeans(timeLS)
  averageChol <- colMeans(timeChol)
  plot(N,averageLS,type='l',lwd=2,ylim=range(averageLS,averageChol),xlab='',ylab='')
  title(xlab='Dimension', ylab='Average Time (s)',line=2)
  lines(N,averageChol,type='l',lty=2,lwd=2)
  legend(200,0.25,c('Cholesky','Fast.LS'),
         lty=c(2,1),cex=1,
         text.font = 2, box.lty=0,lwd=2)
  ##############################################################
}


if(KLEvsLS_N=='yes'){
  ####################################################
  ####### Comparison between samp.WC & Fast.LS #######
  ####################################################
  ## The Matérn smoothness parameter nu is fixed
  par(mar=c(3.1,3.1,1.5,1.1)) # adapt margins
  nu <- 1.5 # smoothness parameter Matern Kernel (MK)
  theta <- l_est(nu,c(0,1),0.05) # length-scale parameter MK
  N <- seq(100,10000,length=10) # sizes of MVN
  N1 <- 100 # size of 1st subdomain
  M <- N/N1 # nb of subdomains
  tol <- 1e-6 # tolerance
  trial <- 25 # nb of replicates
  timeKLE <- matrix(NA, nrow = trial, ncol = length(N))
  timeLS <- matrix(NA, nrow = trial, ncol = length(N))
  for(j in 1 : trial){
    print(j)
    for(i in 1 : length(N)){
      u <- seq(0,1,length=N[i])
      timeLS[j,i] <- system.time(Fast.LS(u=u,M=M[i],N1=N1,nu=nu,l=theta,tausq=1,tol=tol,sseed=12345))[3]
      timeKLE[j,i] <- system.time(LS.KLE(u=u,N1=N1,M=M[i],nu=nu,l=theta,tausq=1,sseed=12345))[3]
    }
  } 
  averageLS <- colMeans(timeLS)
  averageKLE <- colMeans(timeKLE)
  plot(N,averageLS,type='l',lwd=2,ylim=range(averageLS,averageKLE),xlab='',ylab='')
  title(xlab='Dimension', ylab='Average Time (s)',line=2)
  mtext(text =  'LS.KLE versus Fast.LS', side = 3, line = 0.3, cex = 1)
  lines(N,averageKLE,type='l',lty=2,lwd=2)
  legend(4000,0.0055,c('LS.KLE','Fast.LS'),
         lty=c(2,1),cex=1,
         text.font = 2, box.lty=0,lwd=2,bg='transparent')
  ####################################################
}





if(WCvsLS_N=='yes'){
  ####################################################
  ####### Comparison between samp.WC & Fast.LS #######
  ####################################################
  ## The Matérn smoothness parameter nu is fixed
  par(mar=c(3.1,3.1,1.1,1.1)) # adapt margins
  nu <- 0.5 # smoothness parameter Matern Kernel (MK)
  theta <- l_est(nu,c(0,1),0.05) # length-scale parameter MK
  N <- seq(100,10000,length=10) # sizes of MVN
  N1 <- 100 # size of 1st subdomain
  M <- N/N1 # nb of subdomains
  tol <- 1e-6 # tolerance
  trial <- 25 # nb of replicates
  timeWC <- matrix(NA, nrow = trial, ncol = length(N))
  timeLS <- matrix(NA, nrow = trial, ncol = length(N))
  for(j in 1 : trial){
    print(j)
    for(i in 1 : length(N)){
      u <- seq(0,1,length=N[i])
      timeLS[j,i] <- system.time(Fast.LS(u,M=M[i],N1=N1,nu=nu,l=theta,tausq=1,tol=tol,sseed=12345))[3]
      timeWC[j,i] <- system.time(samp.WC(u,nu=nu,l=theta,tausq=1,sseed=12345))[3]
    }
  } 
  averageLS <- colMeans(timeLS)
  averageWC <- colMeans(timeWC)
  plot(N,averageLS,type='l',lwd=2,ylim=range(averageLS,averageWC),xlab='',ylab='')
  title(xlab='Dimension', ylab='Average Time (s)',line=2)
  lines(N,averageWC,type='l',lty=2,lwd=2)
  legend(200,0.07,c('samp.WC','Fast.LS'),
         lty=c(2,1),cex=1,
         text.font = 2, box.lty=0,lwd=2)
  ####################################################
}


if(WCvsLS_nu=='yes'){
  ### running time as a fct of the smoothness parameter \nu 
  par(mar=c(3.1,3.1,1.1,1.1)) # adapt margins
  nu <- seq(0.5,1.5,length=10) # smoothness parameter Matern Kernel (MK)
  N <- 10000 # sizes of MVN
  N1 <- 100 # size of 1st subdomain
  M <- N/N1 # nb of subdomains
  trial <- 50 # nb of replicates
  timeWC <- matrix(NA, nrow = length(nu), ncol = trial)
  timeLS <- matrix(NA, nrow = length(nu), ncol = trial)
  u <- seq(0, 1, length = N)
  for(i in 1 : length(nu)){
    theta <- l_est(nu[i],c(0,1),0.05) # length-scale parameter MK
    print(i)
    for(j in 1 : trial){
      timeLS[i,j] <- system.time(Fast.LS(u,M=M,N1=N1,nu=nu[i],l=theta,tausq=1,tol=1e-6,sseed=j))[3]
      timeWC[i,j] <- system.time(samp.WC(u,nu=nu[i],l=theta,tausq=1,sseed=j))[3]
    }
  } 
  averageLS <- rowMeans(timeLS)
  averageWC <- rowMeans(timeWC)
  plot(nu,averageLS,type='l',lwd=2,ylim=range(averageLS,averageWC),xlab='',ylab='')
  title(ylab='Average Time (s)',line=2)
  mtext(text=expression("smoothness parameter" ~ nu),
        side=1,line=2)
  mtext(text =  paste("Dimension N = ", N), side = 3, line = 0.2, cex = 1)
  lines(nu,averageWC,type='l',lty=2,lwd=2)
  legend(0.5,0.4,c('samp.WC','Fast.LS'),
         lty=c(2,1),cex=1,text.font = 2, 
         box.lty=0,lwd=2,bg='transparent')
  #####################################################
}


if(WCvsLS_theta=='yes'){
  ### running time as a fct of the length-scale
  par(mar=c(3.1,3.1,1.1,1.1)) # adapt margins
  nu <- 0.75 # smoothness parameter Matern Kernel (MK)
  theta <- seq(0.1,1,length=10) # length-scale parameter MK
  N <- 10000 # sizes of MVN
  N1 <- 100 # size of 1st subdomain
  M <- N/N1 # nb of subdomains
  trial <- 50 # nb of replicates
  timeWC <- matrix(NA,length(theta),trial)
  timeLS <- matrix(NA,length(theta),trial)
  u <- seq(0,1,length=N)
  for(i in 1 : length(theta)){
    print(i)
    for(j in 1 : trial){
      timeLS[i,j] <- system.time(Fast.LS(u,M=M,N1=N1,nu=nu,l=theta[i],tausq=1,tol=1e-6,sseed=j))[3]
      timeWC[i,j] <- system.time(samp.WC(u,nu=nu,l=theta[i],tausq=1,sseed=j))[3]
    }
  } 
  averageLS <- rowMeans(timeLS)
  averageWC <- rowMeans(timeWC)
  plot(theta,averageLS,type='l',lwd=2,ylim=range(averageLS,averageWC),xlab='',ylab='')
  title(ylab='Average Time (s)',line=2)
  mtext(text=expression("length-scale parameter"),
        side=1,line=2)
  # mtext(text=expression("length-scale parameter" ~ paste('\u2113')),
  #       side=1,line=2)
  mtext(text =  paste("Dimension N = ", N), side = 3, line = 0.2, cex = 1)
  lines(theta,averageWC,type='l',lty=2,lwd=2)
  legend(0.1,0.4,c('samp.WC','Fast.LS'),
         lty=c(2,1),cex=1,text.font = 2, 
         box.lty=0,lwd=2,bg='transparent')
  #####################################################
}



if(LSvsLS.WC_N=='yes'){
  #####################################################
  ### Fast.LS & Fast.LS.WC as function of dimension ###
  #####################################################
  par(mar=c(3.1,3.1,1.1,1.1)) # adapt margins
  nu <- 0.5 # smoothness parameter Matern Kernel (MK)
  theta <- l_est(nu,c(0,1),0.05) # length-scale parameter MK
  N <- seq(1000,5000,length=5) # sizes of MVN
  N1 <- 1000 # size of 1st subdomain
  M <- N/N1 # nb of subdomains
  trial <- 25 # nb of replicates
  timeLS.WC <- matrix(NA,trial,length(N))
  timeLS <- matrix(NA,trial,length(N))
  for(j in 1 : trial){
    print(j)
    for(i in 1 : length(N)){
      u <- seq(0,1,length=N[i])
      timeLS[j,i] <- system.time(Fast.LS(u,M=M[i],N1=N1,nu=nu,l=theta,tausq=1,tol=1e-7,sseed=12345))[3]
      timeLS.WC[j,i] <- system.time(Fast.LS.WC(u,M=M[i],N1=N1,nu=nu,l=theta,tausq=1,tol=1e-7,sseed=12345))[3]
    }
  } 
  averageLS <- colMeans(timeLS)
  averageLS.WC <- colMeans(timeLS.WC)
  plot(N,averageLS,type='l',lwd=2,ylim=range(averageLS,averageLS.WC),xlab='',ylab='')
  title(xlab='Dimension', ylab='Average Time (s)',line=2)
  lines(N,averageLS.WC,type='l',lty=2,lwd=2)
  legend(2000,2,c('Fast.LS','Fast.LS.WC'),
         lty=c(1,2),cex=1,
         text.font = 2, box.lty=0,lwd=2,bg='transparent')
  #####################################################
}


if(LSvsLS.WC_N1=='yes'){
  ####################################################
  ##### Fast.LS & Fast.LS.WC as function of N_1 ######
  ####################################################
  par(mar=c(3.1,3.1,1.1,1.1)) # adapt margins
  nu <- 0.5 # smoothness parameter Matern Kernel (MK)
  theta <- l_est(nu,c(0,1),0.05) # length-scale parameter MK
  N1 <- seq(200,1000,length=5) # sizes of 1st subdomain
  M <- 5 # nb of blocks
  N <- M*N1
  trial <- 25
  timeLS.WC <- matrix(NA,trial,length(N))
  timeLS <- matrix(NA,trial,length(N))
  for(j in 1 : trial){
    print(j)
    for(i in 1 : length(N)){
      u <- seq(0,1,length=N[i])
      timeLS[j,i] <- system.time(Fast.LS(u,M=M,N1=N1[i],nu=nu,l=theta,tausq=1,tol=1e-7,sseed=12345))[3]
      timeLS.WC[j,i] <- system.time(Fast.LS.WC(u,M=M,N1=N1[i],nu=nu,l=theta,tausq=1,tol=1e-7,sseed=12345))[3]
    }
  } 
  averageLS <- colMeans(timeLS)
  averageLS.WC <- colMeans(timeLS.WC)
  plot(N,averageLS,type='l',lwd=2,ylim=range(averageLS,averageLS.WC),xlab='',ylab='')
  title(xlab='Dimension', ylab='Average Time (s)',line=2)
  lines(N,averageLS.WC,type='l',lty=2,lwd=2)
  legend(2000,2,c('Fast.LS','Fast.LS.WC'),
         lty=c(1,2),cex=1,
         text.font = 2, box.lty=0,lwd=2,bg='transparent')
  #######################################################
}



if(LS1M=='yes'){
  ####################################################
  ##### Performance Fast.LS with N=1,000,000 #########
  ####################################################
  ## second test N1=100
  par(mar=c(3.1,3.1,1.5,1.1)) # adapt margins
  nu <- 0.5 # smoothness parameter Matern Kernel (MK)
  theta <- l_est(nu,c(0,1),0.05) # length-scale parameter MK
  N1 <- 500 # size of 1st subdomain (left panel)
  # N1 <- 100 # size of 1st subdomain (right panel)
  # M <- seq(1000,10000,length=10) # nbs of blocks (right panel)
  M <- seq(500,2000,length=11) # nbs of blocks (left panel)
  N <- M*N1 # sizes of the MVN
  trial <- 25
  timeLS <- matrix(NA,trial,length(M))
  for(i in 1 : trial){
    print(i)
    for(j in 1 : length(M)){
      u <- seq(0,1,length=N[j])
      timeLS[i,j] <- system.time(Fast.LS(u,M=M[j],N1=N1,nu=nu,l=theta,tausq=1,tol=1e-7,sseed=12345))[3]
    }
  }
  averageLS <- colMeans(timeLS)
  plot(N,averageLS,type='l',lwd=2,ylim=range(averageLS),
       xlab='',ylab='')
  title(xlab='Dimension N', ylab='Average Time (s)',line=2)
  mtext(bquote(N[1] == .(N1)), side = 3, line = 0.1, cex = 1)
  # mtext(text =  paste("Dimension N1 =", N1), side = 3, line = 0.1, cex = 1)
  # legend(100000,0.24,c('Fast.LS'),
  #        lty=1,cex=1,lwd=2,
  #        text.font = 2, box.lty=0,bg='transparent')
  ########################################################
}


if(correlation=='yes'){
  #####################################################
  ###### correlation structure verification ###########
  #####################################################
  par(mar=c(3.1,3.1,1.1,1.1)) # adapt margins
  nu <- 0.75 # smoothness parameter Matern Kernel (MK)
  theta <- l_est(nu,c(0,1),0.05) # length-scale parameter MK
  N1 <- 50 # size of 1st subdomain
  M <- 5 # nbs of blocks
  N <- M*N1 # sizes of the MVN
  nbsim <- 15000
  trial <- 25 ## nb of replicates
  u <- seq(0,1,length=N)
  FastLS <- Fast.LS_v(nbsim,u,M,N1,nu=nu,l=theta,tausq=1,tol=0,sseed=12345)
  sampWC <- samp.WC_v(nbsim,u,nu=nu,l=theta,tausq=1,sseed=12345)
  corLS <- matrix(NA,N/2,trial)
  corWC <- matrix(NA,N/2,trial)
  for(j in 1 : trial){
    print(j)
    for(i in 1 : dim(corLS)[1]){
      corLS[i,j] <- cor(FastLS[1,],FastLS[N/2+i,])
      corWC[i,j] <- cor(sampWC[1,],sampWC[N/2+i,])
    }
  }
  averageCorLS <- rowMeans(corLS)
  averageCorWC <- rowMeans(corWC)
  t <- u[-(1:N/2)]
  plot(t,averageCorLS,type='l',lty=2,col='blue',
       lwd=2,ylim=range(corLS,k(t,nu=nu,l=theta),
                        corWC),xlab='',ylab='')
  title(xlab='distance',ylab='correlation',line=2)
  lines(t,k(t,nu=nu,l=theta),lty=1,col='black',lwd=2)
  lines(t,averageCorWC,type='l',lwd=2,col='red',lty=4)
  legend(0.7,0.25,c('True','Fast.LS','samp.WC'),
         lty=c(1,2,4),col=c('black','blue','red'),
         cex=1,text.font = 2, box.lty=0,lwd=2,bg='transparent')
  ## Mean square error (MSE):
  MSE_LS <- mean(sum((averageCorLS-k(t,nu=nu,l=theta))^2))
  MSE_WC <- mean(sum((averageCorWC-k(t,nu=nu,l=theta))^2))
  MSE_LS;MSE_WC
}




if(monotoneC0Approx=='yes'){
  #################################################
  ########## monotone C0 approximation h_j ########
  #################################################
  ############################################################
  ##### C0 approximation using the hat basis fct #######
  ############################################################
  
  par(mar=c(2.1,2.1, 1.1, 1.1)) # adapt margins
  N <- 5
  delta <- 1/(N-1)
  l <- 500 # les pas de x
  u <- seq(0, 1, by = delta) # discretized input set
  f <- function(x){
    x^3#exp(2*x)/6
  }
  x <- seq(0, 1, length = l)
  ## Matrix of the basis functions hi_j
  v <- fcth(x,u,N)
  ## illustration
  matplot(x,v,type='l',ylim=range(f(x),v),
          lty=4,col='gray')
  lines(x,f(x),type='l',col='red',lwd=2)
  lines(x,v%*%f(u),type='l',col='black',lwd=2,lty=2)
  points(u,f(u),pch='*',cex=1.5)
  legend(0,0.92,expression(paste('true function ',f),paste('approximated ',tilde(f)[N])),
         lwd=2,lty=c(1,2),col=c('red','black'),
         box.lty=0,bg='transparent')
}
####################################################





if(comp_pos_syn_LS_WC=='yes'){
  ####################################################
  ########## Synthetic positive data #################
  ####################################################
  ntot <- 100 # nb of total data
  ntr <- floor(ntot*0.8) # nb training data
  nte <- ntot-ntr # nb test data
  N1 <- 15#30
  M <- 10#5
  N <- M*N1
  nbsim <- 5000 # mcmc iterations
  brn <- 1000 # burn in
  trial <- 5#25# # nb of replicates
  timeLS <- rep(NA,trial)
  timeWC <- rep(NA,trial)
  rmse_mAP_LS <- rep(NA,trial)
  rmse_MAP_LS <- rep(NA,trial)
  rmse_mAP_WC <- rep(NA,trial)
  rmse_MAP_WC <- rep(NA,trial)
  nu <- 1.5 # smoothness parameter
  l <- l_est(nu,c(0,1),0.05) # length-scale
  sigN <- 0.1 # sd noise
  tol <- 1e-10
  f <- function(x){
    ## Hassan postive function
    #4*x+(cos(10*x))-0.1
    ## Andrew Pensoneault (Nonnegativity enforced GPR 2020)
    1/(1+(10*x)^4)+0.5*exp(-100*(x-0.5)^2) 
  }
  ## split samples
  set.seed(123)
  xtot <- runif(ntot,0,1)
  ytot <- f(xtot) + rnorm(ntot,0,sd=sigN)
  for(Q in 1 : trial){
    set.seed(2*Q)
    print(Q)
    ind <- sort(sample.int(ntot,ntr))
    xtr <- xtot[ind]
    ytr <- ytot[ind]
    yte <- ytot[-ind]
    xte <- xtot[-ind]
    ytr.te <- f(xte)
    
    post_syn_LS <- pos.LS.ESS(y=ytr,x=xtr,N1,M,
                              mcmc=nbsim,brn=brn,thin=1,
                              nu=nu,l=l,tau.in=1,
                              sig.in=sigN^2,prior='Fast.LS',
                              sseed=Q,return.plot=F,tol=tol)
    post_syn_WC <- pos.WC.ESS(y=ytr,x=xtr,N,
                              mcmc=nbsim,brn=brn,thin=1,
                              nu=nu,l=l,tau.in=1,
                              sig.in=sigN^2,
                              sseed=Q,return.plot=F,tol=tol)
    
    timeLS[Q] <- post_syn_LS$time[3]
    timeWC[Q] <- post_syn_WC$time[3]
    
    ## RMSE LS
    delta=(max(xtr)-min(xtr))/(N-1)
    my_knots=seq(min(xtr),max(xtr),by=delta)
    f_meanLS <- fcth(xte,u=my_knots,N)%*%rowMeans(post_syn_LS$xi_sam)
    f_starLS <- fcth(xte,u=my_knots,N)%*%as.vector(post_syn_LS$z_star)
    rmse_mAP_LS[Q]=sqrt(mean((ytr.te-f_meanLS)^2))
    rmse_MAP_LS[Q]=sqrt(mean((ytr.te-f_starLS)^2))
    
    ## RMSE WC
    f_meanWC <- fcth(xte,u=my_knots,N)%*%rowMeans(post_syn_WC$xi_sam)
    f_starWC <- fcth(xte,u=my_knots,N)%*%as.vector(post_syn_WC$z_star)
    rmse_mAP_WC[Q]=sqrt(mean((ytr.te-f_meanWC)^2))
    rmse_MAP_WC[Q]=sqrt(mean((ytr.te-f_starWC)^2))
  }
  
  ## Illustration Fast LS 
  par(mfrow=c(1,1))
  par(mar=c(2.1,2.1,2.1,1.1)) # adapt margins
  t <- seq(min(xtr),max(xtr),length=500)
  my_knots <- seq(min(xtr),max(xtr),length=N)
  Y_LS <- fcth(t,u=my_knots,N)%*%post_syn_LS$xi_sam
  MAP_LS <- fcth(t,u=my_knots,N)%*%post_syn_LS$z_star
  tmp_LS <- apply(Y_LS, 1, quantile, probs=c(0.025, 0.5, 0.975))
  f_low_LS <- tmp_LS[1,]
  fmean_LS <- tmp_LS[2,]
  f_upp_LS <- tmp_LS[3,]
  plot(xtr,ytr,pch='*',lwd=2,lty=1,col='black',
       ylim=range(f_low_LS,f_upp_LS,ytr),xlab='',ylab='')
  polygon(c(t,rev(t)),y=c(f_low_LS, rev(f_upp_LS)),border=F,col='gray')
  lines(t,f(t),lwd=2,lty=1,col='black')
  lines(t,fmean_LS,type='l',lty=4,lwd=2,col='blue')
  lines(t,MAP_LS,type='l',lty=2,lwd=2,col='red')
  mtext(text =  paste("Average Time (s) = ", round(mean(timeLS),2)), side = 3, line = 0.8, cex = 0.8)
  mtext(text =  paste("RMSE mAP = ", round(mean(rmse_mAP_LS),2), "and RMSE MAP = ", round(mean(rmse_MAP_LS),2)), side = 3, line = 0.1, cex = 0.8)
  # mtext(text =  paste("Dimension N = ", M*N1), side = 3, line = 0.1, cex = 1)
  points(xtr,ytr,pch='*')
  abline(h=0,lty=2,lwd=2)
  legend(0.2,1.15,
         c("true function","MAP","mAP"),
         col = c('black','red','blue'), 
         text.col = "black",lty=c(1,2,4),
         lwd = c(2,2,2), text.font=1,box.lty=0,cex=0.8,
         bg='transparent')
  ## Illustration Wood and Chan FFT
  Y_WC <- fcth(t,u=my_knots,N)%*%post_syn_WC$xi_sam
  MAP_WC <- fcth(t,u=my_knots,N)%*%post_syn_WC$z_star
  tmp_WC <- apply(Y_WC, 1, quantile, probs=c(0.025, 0.5, 0.975))
  f_low_WC <- tmp_WC[1,]
  fmean_WC <- tmp_WC[2,]
  f_upp_WC <- tmp_WC[3,]
  plot(xtr,ytr,pch='*',lwd=2,lty=1,col='black',
       ylim=range(f_low_WC,f_upp_WC,ytr),xlab='',ylab='')
  polygon(c(t,rev(t)),y=c(f_low_WC, rev(f_upp_WC)),border=F,col='gray')
  lines(t,f(t),lwd=2,lty=1,col='black')
  lines(t,fmean_WC,type='l',lty=4,lwd=2,col='blue')
  lines(t,MAP_WC,type='l',lty=2,lwd=2,col='red')
  mtext(text =  paste("Average Time (s) = ", round(mean(timeWC),2)), side = 3, line = 0.8, cex = 0.8)
  mtext(text =  paste("RMSE mAP = ", round(mean(rmse_mAP_WC),2), "and RMSE MAP = ", round(mean(rmse_MAP_WC),2)), side = 3, line = 0.1, cex = 0.8)
  # mtext(text =  paste("Dimension N = ", M*N1), side = 3, line = 0.1, cex = 1)
  points(xtr,ytr,pch='*')
  abline(h=0,lty=2,lwd=2)
  legend(0.2,1.15,
         c("true function","MAP","mAP"),
         col = c('black','red','blue'), 
         text.col = "black",lty=c(1,2,4),
         lwd = c(2,2,2), text.font=1,box.lty=0,cex=0.8,
         bg='transparent')
  ###################################################
}




if(comp_monot_syn_LS_WC=='yes'){
  ##################################################
  #### monotone synthethic (LS versus WC) ##########
  ##################################################
  ntot <- 100 # nb of total data
  ntr <- floor(0.8*ntot) # nb training data
  nte <- ntot-ntr # nb test data
  N1 <- 15
  M <- 10
  N <- N1*M
  trial <- 2#5
  sigN <- 0.5 # sd noise
  nbsim <- 5000 # mcmc iterations
  brn <- 1000 # burn in
  tol <- 1e-10
  nu <- 1.5
  l <- l_est(nu,c(0,1),0.05) # length-scale
  f <- function(x){
    # log(20*x+1) # Maatouk & Bay monotone function
    # ifelse(x>=0.6 & x<=1,(5*x-3)^3,0) # fm1
    3/(1+exp(-10*x+2.1)) # fm2
    # sum <- 0 # fm3
    # for(l in 1 : 100){
    #   sum <- sum+(l^(-1.7)*sin(l)*cos(pi*(l-0.5)*(1-x)))
    # }
    # return(sqrt(2)*sum)
    # 5*x^2 # fm4
    # sin(x*10)+10*x # Dunson 2004 (Biometrics)
  }
  ## split samples
  set.seed(123)
  xtot <- runif(ntot,0,1)
  ytot <- f(xtot) + rnorm(ntot,0,sd=sigN)
  timeLS <- rep(NA,trial)
  timeWC <- rep(NA,trial)
  rmse_MAP_LS <- rep(NA,trial)
  rmse_mAP_LS <- rep(NA,trial)
  rmse_MAP_WC <- rep(NA,trial)
  rmse_mAP_WC <- rep(NA,trial)
  
  
  for(Q in 1 : trial){
    print(Q)
    set.seed(2*Q)
    ind <- sample.int(ntot,ntr)
    xtr <- xtot[ind]
    ytr <- ytot[ind]
    xte <- xtot[-ind]
    yte <- ytot[-ind]
    ytr.te <- f(xte)
    
    post_mon_syn_LS <- mon.inc.LS.ESS(ytr,xtr,N1=N1,M=M,nu=nu,l=l,eta=50,mcmc=nbsim,
                                      brn=brn,thin=1,sig.in=sigN^2,tau.in=1,
                                      sseed=Q,return.plot=F,tol=tol)
    post_mon_syn_WC <- mon.inc.WC.ESS(ytr,xtr,N,nu=nu,l=l,eta=50,mcmc=nbsim,
                                      brn=brn,thin=1,sig.in=sigN^2,tau.in=1,
                                      sseed=Q,return.plot=F)
    
    timeLS[Q] <- post_mon_syn_LS$time[3]
    timeWC[Q] <- post_mon_syn_WC$time[3]
    
    ## RMSE LS
    delta=1/(N-1)
    my_knots=seq(0,1,by=delta)
    f_meanLS <- mean(post_mon_syn_LS$xi0_sam)+fctphi(xte,u=my_knots,N)%*%rowMeans(post_mon_syn_LS$xi_sam)
    f_starLS <- cbind(rep(1,nte),fctphi(xte,u=my_knots,N))%*%as.vector(post_mon_syn_LS$z_star)
    rmse_mAP_LS[Q]=sqrt(mean((ytr.te-f_meanLS)^2))
    rmse_MAP_LS[Q]=sqrt(mean((ytr.te-f_starLS)^2))
    
    ## RMSE WC
    f_meanWC <- mean(post_mon_syn_WC$xi0_sam)+fctphi(xte,u=my_knots,N)%*%rowMeans(post_mon_syn_WC$xi_sam)
    f_starWC <- cbind(rep(1,nte),fctphi(xte,u=my_knots,N))%*%as.vector(post_mon_syn_WC$z_star)
    rmse_mAP_WC[Q]=sqrt(mean((ytr.te-f_meanWC)^2))
    rmse_MAP_WC[Q]=sqrt(mean((ytr.te-f_starWC)^2))
    
  }
  
  ## Illustration Fast LS 
  par(mfrow=c(1,1))
  par(mar=c(2.1,2.1,2.1,1.1)) # adapt margins
  t <- seq(min(xtr),max(xtr),length=500)
  my_knots <- seq(min(xtr),max(xtr),length=N)
  Y_LS <- t(post_mon_syn_LS$xi0_sam+t(fctphi(t,u=my_knots,N)%*%post_mon_syn_LS$xi_sam))
  tmp_LS <- apply(Y_LS, 1, quantile, probs=c(0.025, 0.5, 0.975))
  f_low_LS <- tmp_LS[1,]
  fmean_LS <- tmp_LS[2,]
  f_upp_LS <- tmp_LS[3,]
  MAP_LS <- cbind(rep(1,length(t)),fctphi(t,u=my_knots,N))%*%post_mon_syn_LS$z_star
  plot(xtr,ytr,pch='*',lwd=2,lty=1,col='black',
       ylim=range(f_low_LS,f_upp_LS,ytr,MAP_LS),xlab='',ylab='')
  polygon(c(t,rev(t)),y=c(f_low_LS, rev(f_upp_LS)),border=F,col='gray')
  lines(t,f(t),lwd=2,lty=1,col='black')
  lines(t,fmean_LS,type='l',lty=4,lwd=2,col='blue')
  lines(t,MAP_LS,type='l',lwd=2,lty=2,col='red')
  mtext(text =  paste("Average Time (s) = ", round(mean(timeLS),2)), side = 3, line = 0.8, cex = 0.8)
  mtext(text =  paste("RMSE mAP = ", round(mean(rmse_mAP_LS),2), "and RMSE MAP = ", round(mean(rmse_MAP_LS),2)), side = 3, line = 0.1, cex = 0.8)
  # mtext(text =  paste("Dimension N = ", M*N1), side = 3, line = 0.1, cex = 0.8)
  points(xtr,ytr,pch='*')
  legend(0.4,1.5,
         c("true function","MAP","mAP"),
         col = c('black','red','blue'), 
         text.col = "black",lty=c(1,2,4),
         lwd = c(2,2,2), text.font=1,box.lty=0,cex=0.8,
         bg='transparent')
  ## Illustration Wood and Chan FFT
  Y_WC <- t(post_mon_syn_WC$xi0_sam+t(fctphi(t,u=my_knots,N)%*%post_mon_syn_WC$xi_sam))
  tmp_WC <- apply(Y_WC, 1, quantile, probs=c(0.025, 0.5, 0.975))
  f_low_WC <- tmp_WC[1,]
  fmean_WC <- tmp_WC[2,]
  f_upp_WC <- tmp_WC[3,]
  MAP_WC <- cbind(rep(1,length(t)),fctphi(t,u=my_knots,N))%*%post_mon_syn_WC$z_star
  plot(xtr,ytr,pch='*',lwd=2,lty=1,col='black',
       ylim=range(f_low_WC,f_upp_WC,ytr),xlab='',ylab='')
  polygon(c(t,rev(t)),y=c(f_low_WC, rev(f_upp_WC)),border=F,col='gray')
  lines(t,f(t),lwd=2,lty=1,col='black')
  lines(t,fmean_WC,type='l',lty=4,lwd=2,col='blue')
  lines(t,MAP_WC,type='l',lty=2,lwd=2,col='red')
  mtext(text =  paste("Average Time (s) = ", round(mean(timeWC),2)), side = 3, line = 0.8, cex = 0.8)
  mtext(text =  paste("RMSE mAP = ", round(mean(rmse_mAP_WC),2), "and RMSE MAP = ", round(mean(rmse_MAP_WC),2)), side = 3, line = 0.1, cex = 0.8)
  # mtext(text =  paste("Dimension N = ", M*N1), side = 3, line = 0.1, cex = 0.8)
  points(xtr,ytr,pch='*')
  legend(0.4,1.5,
         c("true function","MAP","mAP"),
         col = c('black','red','blue'), 
         text.col = "black",lty=c(1,2,4),
         lwd = c(2,2,2), text.font=1,box.lty=0,cex=0.8,
         bg='transparent')
  ##################################################
}


## end
