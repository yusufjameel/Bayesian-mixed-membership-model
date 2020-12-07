rm(list=ls(all=TRUE))
library('Rcpp')
set.seed(4)

setwd('U:\\STOTEM_SBM')
sourceCpp('rcpp_func.cpp')
source('gibbs functions.R')
source('SBM_main.R')
dat= data.matrix(read.csv('sim y.csv',as.is=T))
size=data.matrix(read.csv('sim n.csv',as.is=T))

ngroup.loc=10
ngibbs=2000
ngroup.spp=10
gamma.v=0.1
gamma.u=0.1

res=SBM(dat=dat,size=size,ngroup.loc=ngroup.loc,ngroup.spp=ngroup.spp,
        gamma.v=gamma.v,gamma.u=gamma.u,ngibbs=ngibbs)

