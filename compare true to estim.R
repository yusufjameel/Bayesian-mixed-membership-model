plot(res$llk,type='l')
seq1=100:ngibbs
plot(res$llk[seq1],type='l')

k=data.frame(ztrue=z.true,zestim=res$z)
table(k)
ind.loc=c(5,4,2,1,6)

plot(res$theta[ngibbs,],type='h')
res$theta[ngibbs,ind.loc];theta.true

rango=range(c(theta.true,res$theta[ngibbs,ind.loc]))
plot(theta.true,res$theta[ngibbs,ind.loc],xlim=rango,ylim=rango)
lines(rango,rango)
#---------------------------------------
k=data.frame(wtrue=w.true,westim=res$w)
table(k)
ind.spp=c(8,5,9)

res$phi[ngibbs,ind.spp];phi.true

rango=range(c(res$phi[ngibbs,ind.spp],phi.true))
plot(phi.true,res$phi[ngibbs,ind.spp],xlim=rango,ylim=rango)
lines(rango,rango)
#----------------------------------------
psi0=matrix(res$psi[ngibbs,],10,10)
psi=psi0[,ind.spp]
psi1=psi[ind.loc,]
rango=c(0,1)
plot(psi1,psi.true,xlim=rango,ylim=rango)
lines(rango,rango)

