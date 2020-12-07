#include <Rcpp.h>
#include <iostream>
#include <ctime>
#include <fstream>
using namespace Rcpp;

/***************************************************************************************************************************/
/*********************************                      UTILS          *****************************************************/
/***************************************************************************************************************************/

// This function helps with multinomial draws
int whichLessDVPresence(double value, NumericVector prob) {
  int res=-1;
  double probcum = 0;
  
  for (int i = 0; i < prob.length(); i++) {
    probcum = probcum + prob(i);
    if (value < probcum) {
      res = i;
      break;
    }
  }
  return res;
}

// [[Rcpp::export]]
Rcpp::List getql(IntegerVector z, IntegerVector w, IntegerMatrix dat,
                 int ngrloc, int ngrspp, IntegerMatrix datnmy) {
  IntegerMatrix nql1(ngrloc,ngrspp);
  IntegerMatrix nql0(ngrloc,ngrspp);
  
  for(int i=0; i<dat.nrow();i++){
    for (int j=0; j<dat.ncol(); j++){
      nql1(z[i],w[j])=nql1(z[i],w[j])+dat(i,j);
      nql0(z[i],w[j])=nql0(z[i],w[j])+datnmy(i,j);
    }
  }
  
  Rcpp::List resTemp = Rcpp::List::create(Rcpp::Named("nql1") = nql1,
                                          Rcpp::Named("nql0") = nql0);
  return(resTemp);
}

// [[Rcpp::export]]
NumericVector convertSBtoNormal(NumericVector v) {
  NumericVector res(v.size());
  
  res[0]=v[0];
  double prod=1-v[0];
  for(int j=1; j<v.size();j++){
    res(j)=v[j]*prod;    
    prod=prod*(1-v[j]);
  }
  
  return (res);
}

// [[Rcpp::export]]
NumericVector convertNormaltoSB(NumericVector theta) {
  NumericVector res(theta.size());
  
  res[0]=theta[0];
  double prod=1-res[0];
  for(int j=1; j<theta.size();j++){
    res[j]=theta[j]/prod;    
    prod=prod*(1-res[j]);
  }
  res[theta.size()-1]=1;
  return (res);
}

// [[Rcpp::export]]
IntegerVector samplez(NumericVector ltheta,
                      IntegerMatrix dat,
                      IntegerMatrix datnmy,
                      NumericMatrix lpsi,
                      NumericMatrix l1mpsi,
                      IntegerVector w,
                      NumericVector runi) {
  IntegerVector z(dat.nrow());
  NumericVector prob(ltheta.size());
  double tmp=0;
  
  for(int i=0; i<dat.nrow();i++){
    //calculate probabilities
    for (int k=0; k<ltheta.size(); k++){
      tmp=0;
      for (int j=0; j<dat.ncol(); j++){ //sum over all questions
        tmp=tmp+dat(i,j)*lpsi(k,w[j])+datnmy(i,j)*l1mpsi(k,w[j]);      
      }        
      prob[k]=tmp+ltheta[k];
    }
    prob=prob-max(prob);
    prob=exp(prob);
    prob=prob/sum(prob);
    
    //multinomial draw
    z[i]=whichLessDVPresence(runi[i],prob)+1;
  }
  return (z);
}

// [[Rcpp::export]]
IntegerVector samplew(NumericVector lphi,
                      IntegerMatrix dat,
                      IntegerMatrix datnmy,
                      NumericMatrix lpsi,
                      NumericMatrix l1mpsi,
                      IntegerVector z,
                      NumericVector runi) {
  IntegerVector w(dat.ncol());
  NumericVector prob(lphi.size());
  double tmp=0;
  
  for(int j=0; j<dat.ncol();j++){
    //calculate probabilities
    for (int k=0; k<lphi.size(); k++){
      tmp=0;
      for (int i=0; i<dat.nrow(); i++){ //sum over all students
        tmp=tmp+dat(i,j)*lpsi(z[i],k)+datnmy(i,j)*l1mpsi(z[i],k);      
      }        
      prob[k]=tmp+lphi[k];
    }
    prob=prob-max(prob);
    prob=exp(prob);
    prob=prob/sum(prob);
    
    //multinomial draw
    w[j]=whichLessDVPresence(runi[j],prob)+1;
  }
  return (w);
}

