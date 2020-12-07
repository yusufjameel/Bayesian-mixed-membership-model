sample.psi=function(param,dat,ngroup.loc,ngroup.spp,datnmy){
  tmp=getql(z=param$z-1,w=param$w-1,dat=dat,
            ngrloc=ngroup.loc,ngrspp = ngroup.spp,datnmy=datnmy)
  tmp1=rbeta(ngroup.loc*ngroup.spp,tmp$nql1+1,tmp$nql0+1)
  # qqq=tmp$nql1/(tmp$nql1+tmp$nql0)
  matrix(tmp1,ngroup.loc,ngroup.spp)
}
sample.theta=function(param,ngroup.loc,gamma.v){
  nk=rep(0,ngroup.loc)
  tmp=table(param$z)
  nk[as.numeric(names(tmp))]=tmp
  ind=ngroup.loc:1
  invcumsum=cumsum(nk[ind])[ind]
  vk=rbeta(ngroup.loc,nk+1,invcumsum-nk+gamma.v)
  vk[ngroup.loc]=1
  theta=convertSBtoNormal(vk)
  # sum(theta)
  list(vk=vk,theta=theta)
}
sample.phi=function(param,ngroup.spp,gamma.u){
  mk=rep(0,ngroup.spp)
  tmp=table(param$w)
  mk[as.numeric(names(tmp))]=tmp
  ind=ngroup.spp:1
  invcumsum=cumsum(mk[ind])[ind]
  uk=rbeta(ngroup.spp,mk+1,invcumsum-mk+gamma.u)
  uk[ngroup.spp]=1
  phi=convertSBtoNormal(uk)
  # sum(phi)
  list(uk=uk,phi=phi)
}