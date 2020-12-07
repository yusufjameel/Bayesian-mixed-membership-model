SBM=function(dat,size,ngroup.loc,ngroup.spp,gamma.v,gamma.u,ngibbs){
  nloc=nrow(dat)
  nspp=ncol(dat)
  
  #get initial values
  theta=rep(1/ngroup.loc,ngroup.loc)
  phi=rep(1/ngroup.spp,ngroup.spp)
  z=sample(1:ngroup.loc,size=nloc,replace=T)
  w=sample(1:ngroup.spp,size=nspp,replace=T)
  tmp=runif(ngroup.loc*ngroup.spp)
  psi=matrix(tmp,ngroup.loc,ngroup.spp)

  #useful stuff
  datnmy=size-dat
  vec.llk.prior=rep(NA,ngibbs)
  vec.theta=matrix(NA,ngibbs,ngroup.loc)
  vec.z=matrix(NA,ngibbs,nloc)
  vec.w=matrix(NA,ngibbs,nspp)
  vec.phi=matrix(NA,ngibbs,ngroup.spp)
  vec.psi=matrix(NA,ngibbs,ngroup.spp*ngroup.loc)
  param=list(z=z,w=w,psi=psi,theta=theta,phi=phi)
  
  #start gibbs sampler
  options(warn=2)
  for (i in 1:ngibbs){
    print(i)
    param$psi=sample.psi(param=param,dat=dat,datnmy=datnmy,
                         ngroup.loc=ngroup.loc,ngroup.spp=ngroup.spp)
    # param$psi=rbind(cbind(psi.true,0.01,0.01,0.01,0.01,0.01,0.01,0.01),0.01,0.01,0.01,0.01,0.01)  
    lpsi=log(param$psi)
    l1mpsi=log(1-param$psi)
    
    tmp=sample.theta(param=param,ngroup.loc=ngroup.loc,gamma.v=gamma.v)
    param$theta=tmp$theta
    param$vk=tmp$vk
    # param$theta=c(theta.true,rep(0,ngroup.stude-length(theta.true)))
    
    tmp=sample.phi(param=param,ngroup.spp=ngroup.spp,gamma.u=gamma.u)
    param$phi=tmp$phi
    param$uk=tmp$uk
    # param$phi=c(phi.true,rep(0,ngroup.quest-length(phi.true)))
    
    param$z=samplez(ltheta=log(param$theta),dat=dat,datnmy=datnmy,lpsi=lpsi,l1mpsi=l1mpsi,
                   w=param$w-1,runi=runif(nrow(dat)))
    # param$z=z.true
    
    param$w=samplew(lphi=log(param$phi),dat=dat,datnmy=datnmy,lpsi=lpsi,l1mpsi=l1mpsi,
                   z=param$z-1,runi=runif(ncol(dat)))
    # param$w=w.true
    
    #get logl
    psi1=param$psi[param$z,]
    psi2=psi1[,param$w]
    loglikel=dat*log(psi2)+datnmy*log(1-psi2)
    # v=param$vk
    # u=param$uk
    # v[v==1]=NA
    # u[u==1]=NA
    # tmp=sum(loglikel)+
    #     sum(dbeta(v,1,gamma.v,log=T),na.rm=T)+
    #     sum(dbeta(u,1,gamma.u,log=T),na.rm=T)
    
    #store general results
    vec.llk.prior[i]=sum(loglikel)    
    vec.theta[i,]=param$theta
    vec.phi[i,]=param$phi
    vec.psi[i,]=param$psi
    vec.w[i,]=param$w
    vec.z[i,]=param$z
  }
  list(theta=vec.theta,phi=vec.phi,llk=vec.llk.prior,psi=vec.psi,z=vec.z,w=vec.w)
}
