#include <math.h>
#include <complex>
#include <boost/math/special_functions/bessel.hpp>
#include "kers.h"
#include "prods.h"
#include "geo.h"


#include "mex.h"

using namespace std;
const double pi= 4.0*atan(1.0);


typedef complex<double> dcomp;
//const dcomp i(0.0,1.0);




double p1(dcomp* ft_cofs, int l){
    if(l==0){
        return 0.5*pi*real(ft_cofs[0]-0.5*ft_cofs[2]);
    }else{
        return 0.25*pi*real(ft_cofs[l]-ft_cofs[l+2]);
    }
}

dcomp pld(int l,dcomp* ft_cofs,int dom_op){
    //PORQUE??????
    int sig= pow(-1,l);
    
    //sig=1; 
//     if(dom_op==0){
//         sig=1; 
//     }
    if(l==0){
        return sig*0.5*pi*(ft_cofs[0]-0.5*ft_cofs[2]);
    }else{
        return sig*0.25*pi*(ft_cofs[l]-ft_cofs[l+2]);
    }      
}




///////////////////otras cosas...//////////////////////////////////////////

double *t_polynomial ( int m, int n, double x[] ){
  int i;
  int j;
  double *v;

  if ( n < 0 )
  {
    return NULL;
  }

  v = new double[m*(n+1)];

  for ( i = 0; i < m; i++ )
  {
    v[i] = 1.0;
  }
  if ( n < 1 )
  {
    return v;
  }

  for ( i = 0; i < m; i++ )
  {
    v[i+1*m] = x[i];
  }
  for ( j = 2; j <= n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      v[i+j*m] = 2.0 * x[i] * v[i+(j-1)*m] - v[i+(j-2)*m];
    }
  }
  return v;
}


double chebT ( int n, double x ){
  int m;
  double *v_vec;
  double value;
  double x_vec[1];

  m = 1;
  x_vec[0] = x;

  v_vec = t_polynomial ( m, n, x_vec );

  value = v_vec[n];

  delete [] v_vec;

  return value;
}

double *u_polynomial ( int m, int n, double x[] ){
  int i;
  int j;
  double *v;

  if ( n < 0 )
  {
    return NULL;
  }

  v = new double[m*(n+1)];

  for ( i = 0; i < m; i++ )
  {
    v[i+0*m] = 1.0;
  }

  if ( n < 1 )
  {
    return v;
  }

  for ( i = 0; i < m; i++ )
  {
    v[i+1*m] = 2.0 * x[i];
  }

  for ( i = 0; i < m; i++ )
  {
    for ( j = 2; j <= n; j++ )
    {
      v[i+j*m] = 2.0 * x[i] * v[i+(j-1)*m] - v[i+(j-2)*m];
    }
  }
  return v;
}

double chebU ( int n, double x ){
  int m;
  double *v_vec;
  double value;
  double x_vec[1];

  m = 1;
  x_vec[0] = x;

  v_vec = u_polynomial ( m, n, x_vec );

  value = v_vec[n];

  delete [] v_vec;

  return value;
}


double  chebU_p ( int n, double x ){
  return ((n+1)*chebT(n+1,x)-x*chebU( n, x))/(pow(x,2)-1);
}

double TrialFunction (int n, double x)
{
  //  return chebT(n,x); 
    
  //  return chebU(n,x);
 //     return sqrt(1-x*x)*chebU(n,x); 
//     
     return  chebU(n,x); 
    
}


double chebTWeight ( int n, double x )
{
    return  (chebT(n,x)/ sqrt(1-x*x)); 
    
}


double TrialFunctionDev(int n, double x)
{
    
    //only for W
    return chebU_p (  n,  x );
    
 //   return double(-n-1)*chebT(n+1,x)/sqrt(1-x*x);
    
//     if (n==0)
//         return 0.0; 
//     
//     else 
//         
//         return (double)n*chebU(n-1,x);
    
//     return 0; 
}


///////////////////////////////////////////////////////////////////////////



double TrialFunctionNs (int n, double x)
{
    //sin peso porque esta en la cuadratura.
    
   return chebT(n,x); 
   
   // return chebU(n,x);
    
}


dcomp pVNs(int m, int l ,dcomp** ft_cofs,int Nq,double* xq,double* wq, double k,int domt,int doms){ 
 dcomp u(0,0); 
  
 if(l==0){
    for(int j=0; j<Nq;j++){   

        u+=wq[j]*pi*(ft_cofs[j][0])*TrialFunctionNs(m,xq[j]);
    } 
 }else{
    for(int j=0; j<Nq;j++){        
          u+=wq[j]*pi*(ft_cofs[j][l])*TrialFunctionNs(m,xq[j]);
    } 
 }
 if(m == l)
 {
  if(domt==doms){
     if(l==0){
        u += 0.5*pi*log(2);
     }else{
        u += 0.25*pi/((double)m); 
     } 
  }
 }

 return u;
}


dcomp pldNs(int l,dcomp* ft_cofs,int dom_op){

    if(l==0){
        return pi*(ft_cofs[0]);
    }else{
        return 0.5*pi*(ft_cofs[l]);
    }      
}
