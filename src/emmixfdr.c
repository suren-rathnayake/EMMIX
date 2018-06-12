#include "emmix.h"

/*

Modified on Nov 1, 2012.

*/


/*

functions to do EMMIX-FDR

*/

void F77_SUB(emmixfdr)(double *y, int *pn, int *pp,int *pg,  int *pncov,
double * pro, double *mu, double * sigma, 
double *tau,  double *sumtau, 
double *loglik, double *lk, 
int *clust, int *pinfo,  int *itmax, double *epsilon)
{
int it;
double zero=0.0;

for(it=0;it< *itmax;it++) lk[it]=zero;
//   start of EM algorithm loop
*pinfo =0;

for(it=0;it< *itmax;it++) {//start of EM loops
//  E-Step
F77_SUB(estepmvn)(y,pn,pp,pg,pro,mu,sigma,tau,sumtau,loglik,pinfo);

if(*pinfo) {*pinfo = *pinfo+10; return;}
lk[it]= *loglik;

//  M-step
F77_SUB(mstepmvn3)(y,pn,pp,pg,pncov,mu,sigma,tau,sumtau);

//  Check if it converges

if(it == (*itmax-1) ) {*pinfo = 1;break;}
	
if(it < imin2(19,*itmax-1)) continue;

if( (fabs(lk[it-10]-lk[it] ) < fabs(lk[it-10]* (*epsilon)))
&& (fabs(lk[it-1]-lk[it]) < fabs(lk[it-1]*(*epsilon)))) break;   
}      
//   end of EM algorithm loop

if(*pinfo <= 1)
{

F77_SUB(tau2clust)(tau, pn, pg, clust);

}


return;
}

void F77_SUB(mstepmvn3)(double *y,int *pn, int *pp, int *pg,int *pncov, 
double *mu, double *sigma,double *tau, double *sumtau)
{
int p= *pp,n= *pn,g= *pg;
int h,i,j,k;
const double two=2.0,zero= 0.0;
double sum;

for(h=0;h<g;h++) {

//   calculate the means 

//#####

   if(*pncov==10 && h==0) continue;

//#####

   for(j=0;j<p;j++) {

	  sum=zero;
      
	  for(i=0;i<n;i++)
	        sum += y[j*n+i]*tau[h*n+i];
      
	  if(sumtau[h] < two) 
            mu[h*p+j] = zero;
	  else
     	    mu[h*p+j] = sum/sumtau[h];}


//  calculate the covariances 
   for(j=0;j<p;j++)     
	  for(i=0;i<=j;i++)
         sigma[h*p*p+j*p+i] = zero;

   for(i=0;i<n;i++)
	  for(k=0;k<p ;k++)
	     for(j=0;j<=k;j++)
sigma[h*p*p+k*p+j]  += (y[j*n+i]-mu[h*p+j])*(y[k*n+i]-mu[h*p+k])*tau[h*n+i];

   for(j=0;j<p;j++)   { 

      for(i=0;i<=j;i++) { 

	  if(sumtau[h] < two)
        sigma[h*p*p+j*p+i] =zero;
	  else
	    sigma[h*p*p+j*p+i]=sigma[h*p*p+j*p+i]/sumtau[h];
	  
	    sigma[h*p*p+i*p+j]=sigma[h*p*p+j*p+i];}
}
}  
// end of loop (h)




return;
}


/*
end of modification
*/
