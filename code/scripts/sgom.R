library(dplyr)
library(tidyverse)
library(data.table)
library(bedr)
library(readr)
#library(ebpmf)

library(NNLM)
#library(Matrix)
library(CountClust)
library(fastTopics)
library(FNN)

library(robustbase)

library(smashr)
library(ebpmf)
library(Matrix)

library(robustbase)

library(smashr)
library(ashr)
library(wavethresh)

#'@title Smooth Poisson sequence, accounting for nugget effect
#'@param x observed Poisson sequence
#'@param nugget nugget effect
#'@param s Scale factor for Poisson observations: y~Pois(scale*lambda), can be a vector.
#'@param transformation transformation of Poisson data, either 'vst' or 'lik_expansion'; 'vst' for varaince stablizing transformation; 'lik_expansion' for likelihood expansion
#'@param robust whether perform robust wavelet regression
#'@param robust.q quantile to determine outliers
#'@param method smoothing method for Gaussian sequence, either 'smash' or 'ti.thresh'. When n is large, ti.thresh is much faster
#'@param nug.init init value of nugget effect, either a scalar or NULL
#'@param ash.pm If choose lik_expansion, whehter use ash posterior mean approxiamtion if x=0. If not x = x+eps.
#'@param eps If choose lik_expansion, if x=0, set x = x + eps. Either input a numerical value or 'estimate'. If estimate, eps = sum(x==1)/sum(x<=1)
#'@param filter.number,family wavelet basis, see wavethresh pakcage for more details
#'@param maxiter max iterations for estimating nugget effect
#'@param tol tolerance to stop iterations.
#'@return estimated smoothed lambda, estimated nugget effect.
#'@import smashr
#'@import ashr
#'@export

smash.gen.poiss = function(x,nugget=NULL,s=1,
                          transformation = 'vst',
                          method='ti.thresh',
                          robust = T,
                          robust.q = 0.99,
                          nug.init = NULL,
                          ash.pm=FALSE,
                          eps='estimate',
                          filter.number = 1,
                          family = "DaubExPhase",
                          maxiter=10,
                          tol=1e-2){

  if(!ispowerof2(length(x))){
    reflect.x = reflect(x)
    x = reflect.x$x
    idx = reflect.x$idx
  }else{
    idx = 1:length(x)
  }

  n = length(x)

  if(transformation == 'vst'){
    y = sqrt(x+3/8)/sqrt(s)
    st = rep(sqrt(0.25/s),n)
  }

  if(transformation == 'lik_expansion'){

      lambda_tilde = x/s
      if(min(x)<1){
        if(ash.pm){
          x_pm = ash(rep(0,n),1,lik=lik_pois(x,scale=s),
                     optmethod='mixSQP',pointmass=T)$result$PosteriorMean
          lambda_tilde[x<1] = x_pm[x<1]
        }else{
          if(eps == 'estimate'){
            eps = sum(round(x)==1)/sum(round(x)<=1)+0.1
          }
          lambda_tilde[x<1] = (x[x<1]+eps)/s
        }
      }
    # working data
    st=sqrt(1/(s*lambda_tilde))
    y=log(lambda_tilde)+(x-s*lambda_tilde)/(s*lambda_tilde)

  }


  # estimate nugget effect and estimate mu

  if(robust){
    win.size = round(sqrt(n)/2)*2+1
    #win.size = round(log(n,2)/2)*2+1
    #browser()
    y.wd = wd(y,filter.number,family,'station')
    y.wd.coefJ = accessD(y.wd,level = log(n,2)-1)
    y.rmed = runmed(y,win.size)

    robust.z = qnorm(0.5+robust.q/2)

    if(is.null(nugget)){
      nug.init = uniroot(normaleqn,c(-1e6,1e6),y=y,mu=y.rmed,st=st)$root
      nug.init = max(c(0,nug.init))
      outlier.idx = which(abs(y-y.rmed)>=(robust.z*sqrt(st^2+nug.init)))
    }else{
      outlier.idx = which(abs(y-y.rmed)>=robust.z*sqrt(st^2+nugget))
    }
    st[outlier.idx] = abs((y.wd.coefJ)[outlier.idx] - median(y.wd.coefJ))
  }



  if(is.null(nugget)){
    fit = NuggetEst(y,st,nug.init,method,filter.number = filter.number,family = family,maxiter,tol)
    nug.est = fit$nug.est
  }else{
    nug.est = nugget
    if(method=='smash'){
      fit = smash.gaus(y,sigma=sqrt(st^2+nug.est),
                       filter.number = filter.number,family = family,
                       post.var = TRUE)
    }
    if(method == 'ti.thresh'){
      fit = ti.thresh(y,sigma=sqrt(st^2+nug.est),filter.number = filter.number,family = family)
      fit = list(mu.est = fit, mu.est.var = rep(0,length(fit)))
    }
  }

  mu.est = (fit$mu.est)[idx]
  mu.est.var = (fit$mu.est.var)[idx]

  if(transformation == 'vst'){
    lambda.est = mu.est^2-3/(8*s)
  }
  if(transformation == 'lik_expansion'){
    lambda.est = exp(mu.est+mu.est.var/2)
  }

  return(list(lambda.est=lambda.est,mu.est=mu.est,nugget.est=nug.est))
}



normaleqn=function(nug,y,mu,st){
  return(sum((y-mu)^2/(nug+st^2)^2)-sum(1/(nug+st^2)))
}

NuggetEst=function(y,st,nug.init=NULL,method,filter.number,family,maxiter,tol){
  #initialize nugget effect sigma^2
  n=length(y)
  if(is.null(nug.init)){
    x.m=c(y[n],y,y[1])
    st.m=c(st[n],st,st[1])
    nug.init = ((x.m[2:n]-x.m[3:(n+1)])^2+(x.m[2:n]-x.m[1:(n-1)])^2-2*st.m[2:n]^2-st.m[1:(n-1)]^2-st.m[3:(n+1)]^2)/4
    nug.init = nug.init[nug.init>0&nug.init<var(y)]
    nug.init = median(nug.init)
  }
  #given st and nug to estimate mean


  nug.est = nug.init
  for(iter in 1:maxiter){

    #print(nug.est)

    # update mean
    if(method == 'smash'){
      est = smash.gaus(y,sigma=sqrt(st^2+nug.est),filter.number = filter.number,family = family,post.var = TRUE)
      mu.est = est$mu.est
      mu.est.var = est$mu.est.var
    }
    if(method == 'ti.thresh'){
      mu.est = ti.thresh(y,sigma=sqrt(st^2+nug.est),filter.number = filter.number,family = family)
      mu.est.var = rep(0,n)
    }

    # update nugget effect
    nug.est.new=uniroot(normaleqn,c(-1e6,1e6),y=y,mu=mu.est,st=st)$root
    nug.est.new = max(c(0,nug.est.new))


    if(abs(nug.est.new - nug.est)<=tol){
      break
    }else{
      nug.est = nug.est.new
    }
  }

  return(list(mu.est = mu.est,mu.est.var=mu.est.var,nug.est=nug.est))

}

#'@export
ispowerof2 <- function (x){
  x >= 1 & 2^ceiling(log2(x)) == x
}



#'@title Smoothed Poisson Matrix Factorization
#'@description Smoothed Poisson Matrix Factorization/ smoothed topic model
#'@param X count matrix
#'@param K number of factors/ranks
#'@param init initialization methods, 'lee','scd' from package NNLM, or 'uniform' randomly initialize; or provide init as a list with L_init and F_init.
#'@param init_loss loss function of the initialization method, either mkl or mse.
#'@param maxiter maximum iterations
#'@param tol stop criteria
#'@param fix_F if TRUE, F will no be updated.
#'@param bmsm_control control parameters of BMSM, see bmsm_control_default()
#'@param ebpm_method point_gamma or two_gamma
#'@param ebpm_control control parameters of ebpm, see ebpm_control_default()
#'@param nug_control control parameters of smashgen, see nug_control_default()
#'@param smooth_f,smooth_l whether to get smooth estimate of loadings or factors.
#'@param nugget whether to assume nugget effects
#'@param return_all whether return all outputs or simplified ones
#'@return EL,EF: posterior of loadings and factors
#'@import mixsqp
#'@import NNLM
#'@import ebpm
#'@import Matrix
#'@export

stm = function(X,K,
               init = 'scd',init_loss = 'mkl',maxiter=100,tol=1e-3,
               fix_F = FALSE,
                  bmsm_control_l=list(), bmsm_control_f=list(),
                  nug_control_l=list(), nug_control_f=list(),
                  #filter.number = 1,family = "DaubExPhase",
                  ebpm_method='point_gamma',
                  ebpm_control_l=list(), ebpm_control_f=list(),
                  smooth_f=TRUE,smooth_l=FALSE,
                  nugget=FALSE,
               printevery=10){


  n = dim(X)[1]
  p = dim(X)[2]

  res = init_stm(X,K,init,init_loss)
  #plot(res$ql$El[,1])
  #plot(res$ql$El[,2])
  #plot(res$ql$El[,3])
  #inited = list(L_init = res$ql$El,F_init = res$qf$Ef)

  #EZ = array(dim = c(n,p,K))

  #EZ = Calc_EZ(X,K,EZ,ql_hat,qf_hat)

  KL = c()
 # loglik = c()

  #browser()
  KL[1] = mKL(X,tcrossprod(res$ql$El,res$qf$Ef))

  X = Matrix::Matrix(X,sparse = TRUE)
  X_idx = summary(X)

  for(iter in 1:maxiter){

    b_k_max = 0

    for(k in 1:K){

      # get row and col sums of EZ_k
      b_k = res$ql$Elogl[X_idx$i,k]+res$qf$Elogf[X_idx$j,k] - res$a

      EZ_k = sparseMatrix(i=X_idx$i,j=X_idx$j,x = X_idx$x*exp(b_k)/res$b,dims = c(n,p))


      l_seq = rowSums(EZ_k)
      l_scale = sum(res$qf$Ef[,k])


      # adj.ratio = sqrt(l_scale/f_scale)
      #
      # l_scale = l_scale/adj.ratio
      # l_seq = l_seq * adj.ratio
      # f_scale = f_scale*adj.ratio
      # f_seq = f_seq / adj.ratio


      #print(l_scale)
      #print(f_scale)

      # Update L
      if(smooth_l){
        lk_hat = update_smooth(l_seq, l_scale, nugget,bmsm_control_l,nug_control_l)
        res$ql$El[,k] = lk_hat$E
        res$ql$Elogl[,k] = lk_hat$Elog
        #loglikL = loglikL + lk_hat$loglik
        res$nugget_l[k] = lk_hat$nugget

        res$gl[[k]] = lk_hat$pi_weights

      }else{
        lk_hat = update_nsmooth(l_seq,l_scale,ebpm_control_l,ebpm_method)
        res$ql$El[,k] = lk_hat$posterior$mean
        res$ql$Elogl[,k] = lk_hat$posterior$mean_log
        #loglikL = loglikL + lk_hat$log_likelihood
        res$gl[[k]] = lk_hat$fitted_g
      }


      # Update F

      if(!fix_F){
        f_seq = colSums(EZ_k)
        f_scale = sum(res$ql$El[,k])

        if(smooth_f){
          fk_hat = update_smooth(f_seq, f_scale,nugget,bmsm_control_f,nug_control_f)
          res$qf$Ef[,k] = fk_hat$E
          res$qf$Elogf[,k] = fk_hat$Elog
          #loglikR = loglikR + fk_hat$loglik
          res$nugget_f[k] = fk_hat$nugget
          res$gf[[k]] = fk_hat$pi_weight

        }else{
          fk_hat = update_nsmooth(f_seq,f_scale,ebpm_control_f,ebpm_method)
          res$qf$Ef[,k] = fk_hat$posterior$mean
          res$qf$Elogf[,k] = fk_hat$posterior$mean_log
          #loglikR = loglikR + fk_hat$log_likelihood
          res$gf[[k]] = fk_hat$fitted_g
        }
      }

      b_k_new = res$ql$Elogl[X_idx$i,k] + res$qf$Elogf[X_idx$j, k] - res$a
      res$b = res$b - exp(b_k) + exp(b_k_new)
      b_k_max = pmax(b_k_new, b_k_max)


    }
    # # Update Z
    # EZ = Calc_EZ(X,K,EZ,ql_hat,qf_hat)

    res$b = res$b/exp(b_k_max)
    res$a = b_k_max + res$a


    KL[iter+1] = mKL(X,tcrossprod(res$ql$El,res$qf$Ef))
    #loglik[iter] = loglikR+loglikL

    ########
    if(iter%%printevery==0){
      print(sprintf('At iter %d, mean KL: %f',iter,KL[iter+1]))
    }
    ########

    if(abs(KL[iter+1]-KL[iter])<=tol){
      break
    }

  }
  if(iter==maxiter){
    warning('Reached maximum iterations')
  }

  lambda_hat = tcrossprod(res$ql$El,res$qf$Ef)
  #lambda_init = L_init%*%F_init

  # loglik = sum(dpois(X,lambda_hat,log = TRUE))

  ldf = poisson2multinom(res$qf$Ef,res$ql$El)
  fit = list(res = res,EL = ldf$L,EF = ldf$FF,d=ldf$s)
  return(fit)
  # if(return_all){
  #   return(list(ql=ql_hat,qf=qf_hat,gf=gf_hat,gl=gl_hat,KL=KL,Lambda_hat=lambda_hat,
  #               init = inited,
  #               input = list(X=X,K=K),nugget=list(nugget_l=nugget_l,nugget_f=nugget_f)))
  # }else{
  #   return(list(ql=ql_hat$El,qf=qf_hat$Ef,nugget=list(nugget_l=nugget_l,nugget_f=nugget_f),KL=KL))
  # }

}

#'@title initialize the stm model
#'@param X input data matrix
#'@param K number of topics
#'@param init init methods, or a list of init L and F
#'@param init_loss mkl or mse
#'@export
init_stm = function(X,K,init,init_loss){

  if(is.list(init)){
    L_init = init$L_init
    F_init = init$F_init

    if(is.null(L_init)){
      X_init_fit = NNLM::nnmf(as.matrix(X),K,method='lee',
                              loss='mse',show.warning = F,
                              init = list(H=t(F_init)),
                              verbose = F,max.iter = 50)
      L_init = X_init_fit$W
    }

  }else{

    if(init%in%c('scd','lee')){
      X_init_fit = NNLM::nnmf(as.matrix(X),K,method=init,loss=init_loss,show.warning = F,verbose = F,max.iter = 50)
      L_init = X_init_fit$W
      F_init = t(X_init_fit$H)
    }

    # if(init == 'scd'){
    #   X_init_fit = NNLM::nnmf(as.matrix(X),K,method='scd',loss=init_loss,show.warning = F,verbose = F,max.iter = 50)
    #   L_init = X_init_fit$W
    #   F_init = t(X_init_fit$H)
    # }
    # if(init == 'lee'){
    #   X_init_fit = NNLM::nnmf(as.matrix(X),K,method='lee',loss=init_loss,show.warning = F,verbose = F,max.iter = 100)
    #   L_init = X_init_fit$W
    #   F_init = t(X_init_fit$H)
    # }
    if(init == 'uniform'){
      L_init = matrix(runif(n*K),nrow=n,ncol=K)
      F_init = matrix(runif(K*p),nrow=p,ncol=K)
      ratio = median(X)/(median(L_init)*median(F_init))
      L_init = L_init*sqrt(ratio)
      F_init = F_init*sqrt(ratio)
    }
    if(init == 'kmeans'){
      kmeans.init=kmeans(as.matrix(X),K,nstart=5)
      L_init = rep(1,n)%o%normalize(as.vector(table(kmeans.init$cluster)))
      F_init = t(kmeans.init$centers)
      row.names(F_init)=NULL
    }
  }

  # adjust scale of L and F, mainly for stability.
  ratio = adjLF(L_init,F_init)
  L_init = ratio$L_init
  F_init = ratio$F_init

  Elogl = log(L_init+1e-10)
  Elogf = log(F_init+1e-10)

  ql = list(El = L_init, Elogl = Elogl)
  qf = list(Ef = F_init, Elogf = Elogf)

  a = 0
  b = 0

  X = Matrix(X,sparse = TRUE)
  d = summary(X)

  temp = Elogl[d$i,] + Elogf[d$j,]
  a = rowMax(temp)
  b = rowSums(exp(temp-a))


  gl = list()
  gf = list()

  return(list(ql=ql,qf=qf,gl=gl,gf=gf,
              a=a,b=b,
              nugget_l = rep(0,K),
              nugget_f = rep(0,K)))

}

rowMax = function(X){
  do.call(pmax.int, c(na.rm = TRUE, as.data.frame(X)))
}

Calc_EZ = function(X,K,EZ,ql_hat,qf_hat){
  n = nrow(X)
  p = ncol(X)
  for(k in 1:K){
    EZ[,,k] = outer(ql_hat$Elogl[,k], qf_hat$Elogf[k,], "+")
  }
  EZ = softmax3d(EZ)
  EZ = as.vector(EZ)*as.vector(X)
  dim(EZ) = c(n,p,K)
  EZ
}

softmax3d <- function(x){
  score.exp <- exp(x)
  probs <-as.vector(score.exp)/as.vector(rowSums(score.exp,dims=2))
  probs[is.na(probs)] = 0
  dim(probs) <- dim(x)
  return(probs)
}

update_smooth = function(x,sf,nugget,bmsm_control=list(),nug_control=list()){
  if(min(x) < 0){stop ("negative values in x not permitted")}
  if(nugget){

    control0 = nug_control_default()
    if (any(!is.element(names(nug_control),names(control0))))
      stop("Argument \"nug_control\" contains unknown parameter names")

    control1 = modifyList(control0,nug_control,keep.null = TRUE)

    fit = smash.gen.poiss(x,s=sf,filter.number=control1$filter.number,
                          family=control1$family,nugget=control1$nugget,
                          robust=control1$robust,
                          robust.q = control1$robust.q,
                          transformation = control1$transformation,
                          method = control1$method,
                          nug.init = control1$nug.init,
                          ash.pm = control1$ash.pm,
                          eps = control1$eps,
                          maxiter = control1$maxiter,
                          tol = control1$tol)
    est = fit$lambda.est
    if(control1$transformation=='lik_expansion'){
      est_log = fit$mu.est
    }else{
      est_log = log(est)
    }
    pi_weights = NULL
    nugget.est = fit$nugget.est
  }else{
    fit = BMSM(x,sf,bmsm_control)
    est = fit$E
    est_log = fit$Elog
    pi_weights = fit$pi_weights
    nugget.est = 0
  }

  #loglik = fit$loglik

  results = list("E" = est,
                 'Elog' = est_log,
                 "pi_weights" = pi_weights,
                 #"loglik" = loglik,
                 "nugget" = nugget.est)

  return(results)

}


update_nsmooth = function(x,s,ebpm_control = list(),ebpm_method){

  control0 = ebpm_control_default()
  if (any(!is.element(names(ebpm_control),names(control0))))
    stop("Argument \"ebpm_control\" contains unknown parameter names")

  control1 = modifyList(control0,ebpm_control,keep.null = TRUE)

  #scale = control1$scale
  #point_mass=control1$point_mass
  #nullweight=control1$nullweight
  #shape= control1$shape
  g_init = control1$g_init
  fix_g = control1$fix_g
  #m = control1$m
  control =  control1$control
  #low = control1$low
  #d = control1$d
  pi0 = control1$pi0

  if(ebpm_method=='point_gamma'){
    out = ebpm::ebpm_point_gamma(x,s,g_init,fix_g,pi0,control)
  }
  if(ebpm_method=='two_gamma'){
    out = ebpm::ebpm_two_gamma(x,s,g_init,fix_g,pi0,control)
  }

  out



}

#'@title Default parameters of ebpm
#'@export
ebpm_control_default = function(){
  list(pi0 = 'estimate',
       #point_mass=F,
       #nullweight=100,
       #shape=1,
       g_init = NULL,
       fix_g = FALSE,
       #m = 2,
       control =  NULL)
       #low = NULL,
       #d=NULL
}


#'@title Default parameters of smash gen
#'@param filter.number,family wavelet basis, see wavethresh pakcage for more details.
#'@export
nug_control_default = function(){
  list(nugget = NULL,
       robust=T,
       robust.q = 0.99,
       transformation = 'lik_expansion',
       method='ti.thresh',
       nug.init = NULL,
       ash.pm=FALSE,
       eps='estimate',
       maxiter=10,
       tol=1e-2,
       filter.number = 1,
       family = "DaubExPhase")
}




#'@title Xing's smoothed topic model
#'@export
#'
cluster.mix=function(y,smooth=TRUE,pi0=NULL,phi0=NULL,K,tol,maxit,nugget=FALSE,nug_control=list()){
  n=dim(y)[1]
  B=dim(y)[2]
  #if(is.null(pseudocounts)) pseudocounts=10^(round(log10(1/B/100000)))

  if(is.null(pi0)|is.null(phi0)){
    kmeans.init=kmeans(y,K,nstart=5)
  }

  if(is.null(pi0)) pi0=rep(1,n)%o%normalize(as.vector(table(kmeans.init$cluster)))

  if(is.null(phi0)){
    phi0=kmeans.init$centers
    #phi0[phi0==0]=pseudocounts
    #phi0=t(apply(phi0,1,normalize))
    row.names(phi0)=NULL
  }

  out=EMproc.mix(y,smooth,pi0,phi0,n,K,B,tol,maxit,nugget,nug_control)
  return(list(pi=out$pi,phi=out$phi,lambda=out$lambda,gamma=out$gamma,loglik=out$loglik))
}

normalize=function(x){
  #if(sum(abs(x))!=0){
  return(x/sum(x))
  #}else{
  #  return(rep(0,length(x)))
  #}
}

smooth.lambda = function(lambda,nugget,nug_control){
  if(nugget){
    control0 = nug_control_default()
    control1 = modifyList(control0,nug_control,keep.null = TRUE)
    return(t(apply(lambda,1,function(z){
      fit = smash.gen.poiss(z,filter.number=control1$filter.number,
      family=control1$family,nugget=control1$nugget,
      robust=control1$robust,
      robust.q = control1$robust.q,
      transformation = control1$transformation,
      method = control1$method,
      nug.init = control1$nug.init,
      ash.pm = control1$ash.pm,
      eps = control1$eps,
      maxiter = control1$maxiter,
      tol = control1$tol)
      return(fit$lambda.est)
    })))
  }else{
    return(t(apply(lambda,1,smash.poiss,cxx = FALSE)))
  }

}

###mixed membership
EMupd.mix=function(y,smooth,pi,phi,n,K,B,nugget,nug_control){
  #gamma is nB*K, pi is n*K, phi is K*B, y is n*B
  gamma=pi[rep(1:n,each=B),]*t(phi)[rep(1:B,n),]
  gamma=gamma/rowSums(gamma)
  gamma[is.na(gamma)]=1/K
  gammab=(as.vector(t(y))%o%rep(1,K))*gamma
  pi.num=t(apply(array(gammab,dim=c(B,n,K)),2,colSums))
  #pi.num=(diag(1,n)[,rep(1:n,each=B)])%*%gammab
  pi=pi.num/(rowSums(y)%o%rep(1,K))
  ybt=t(apply(array(gammab,dim=c(B,n,K)),1,colSums))
  #ybt=(diag(1,B)[,rep(1:B,n)])%*%gammab
  #ybw=(diag(1,B)[,rep(1:B,n)])%*%gamma
  phi=t(ybt/(rep(1,B)%o%colSums(gammab)))
  #phi[phi==0]=pseudocounts
  #phi=t(apply(phi,1,normalize))
  #ykb=ybt/ybw
  #ykb[is.na(ykb)]=0
  #ykt=colSums(ykb)
  lscale=((colSums(ybt)/colSums(pi))%o%rep(1,B))
  lambda=phi*lscale
  #save(lambda,file="D:/Grad School/projects/sequence_clustering/results/analysis_k562ctcf/debug_lambda.Robj")
  phi.unsmoothed=NULL
  if(smooth==TRUE){
    phi.unsmoothed=phi
    lambda.unsmoothed=lambda
    lambda=smooth.lambda(lambda,nugget,nug_control)
    lambda[is.na(lambda)]=lambda.unsmoothed[is.na(lambda)]
    phi=lambda/lscale
  }


  return(list(pi=pi,phi=phi,phi.unsmoothed=phi.unsmoothed,lambda=lambda,gamma=gamma))
}


negloglik.mix=function(y,pi,phi,n,K,B){
  loglik.ini=log(pi%*%phi)
  yloglik=y*loglik.ini
  yloglik[is.na(yloglik)]=0
  loglik.tot=-sum(yloglik)
  return(loglik.tot)
}



# EMproc.mix=function(y,smooth,pi,phi,n,K,B,tol,maxit){
#   loglik.old=Inf
#   loglik=negloglik.mix(y,pi,phi,n,K,B)
#   cyc=0
#   while(abs(loglik-loglik.old)>tol&cyc<maxit){
#     loglik.old=loglik
#     res=EMupd.mix(y,smooth,pi,phi,n,K,B)
#     pi=res$pi
#     phi=res$phi
#     phi.unsmoothed=res$phi.unsmoothed
#     gamma=res$gamma
#     lambda=res$lambda
#     if(smooth==TRUE){
#       loglik=negloglik.mix(y,pi,phi.unsmoothed,n,K,B)
#     }else{
#       loglik=negloglik.mix(y,pi,phi,n,K,B)
#     }
#     cyc=cyc+1
# #print(cyc)
# #print(pi)
# print(loglik)
#   }
#   return(list(pi=pi,phi=phi,phi.unsmoothed,lambda=lambda,gamma=gamma,loglik=loglik))
# }


rowquantiles = function(x, q) apply(x, 1, quantile, probs = q)

normalized.norm = function(x, y){
  return(rowquantiles(abs(x - y), 0.5)*dim(y)[2])
  #  return(rowquantiles(abs(x - y), 0.95))
}

tol.criterion = function(phi.old, phi, pi.old, pi){
  return(max(max(normalized.norm(phi.old, phi)), mean(normalized.norm(pi.old, pi))))
}

EMproc.mix=function(y,smooth,pi,phi,n,K,B,tol,maxit,nugget,nug_control){
  pi.old=matrix(Inf,nrow=n,ncol=K)
  phi.old=matrix(Inf,nrow=K,ncol=B)
  cyc=0
  while(tol.criterion(phi.old,phi,pi.old,pi)>tol&cyc<maxit){
    pi.old=pi
    phi.old=phi
    res=EMupd.mix(y,smooth,pi,phi,n,K,B,nugget,nug_control)
    pi=res$pi
    phi=res$phi
    phi.unsmoothed=res$phi.unsmoothed
    gamma=res$gamma
    lambda=res$lambda
    if(smooth==TRUE){
      loglik=negloglik.mix(y,pi,phi.unsmoothed,n,K,B)
    }else{
      loglik=negloglik.mix(y,pi,phi,n,K,B)
    }
    cyc=cyc+1
    print("iteration")
    print(cyc)
    #print(pi)
    print("phi difference")
    print(max(normalized.norm(phi.old, phi)))
    print("pi difference")
    print(mean(normalized.norm(pi.old, pi)))
    print("negative loglikelihood")
    print(loglik)
  }
  return(list(pi=pi,phi=phi,phi.unsmoothed,lambda=lambda,gamma=gamma,loglik=loglik))
}

