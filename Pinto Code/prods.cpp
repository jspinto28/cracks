#include <math.h>
#include <complex>
#include <boost/math/special_functions/bessel.hpp>
#include "kers.h"
#include "prods.h"
#include "geo.h"
#include <numeric>


#include "mex.h"

using namespace std;
const double pi= 4.0*atan(1.0);


typedef complex<double> dcomp;
//const dcomp i(0.0,1.0);


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


double *t_polynomialD ( int m, int n, double x[], double delta ){
  int i;
  int j;
  double *v;
  
  double *vaux;

  if ( n < 0 )
  {
    return NULL;
  }

  v = new double[m*(n+1)];
  
  vaux = new double[m*(n+1)];

  for ( i = 0; i < m; i++ )
  {
    v[i] = 1.0;
    
    vaux[i] = 1.0;
  }
  if ( n < 1 )
  {
    delete[] vaux; 
      
    return v;
  }

  for ( i = 0; i < m; i++ )
  {
    vaux[i+1*m] = x[i];   
      
    v[i+1*m] = vaux[i+1*m]*pow(2.0,(-2.0-delta));    
    
  }
  for ( j = 2; j <= n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      vaux[i+j*m] = 2.0 * x[i] * vaux[i+(j-1)*m] - vaux[i+(j-2)*m];
      
      v[i+j*m] = vaux[i+j*m]*pow((1.0+j),-2.0-delta); 
    }
  }
  
  delete[] vaux; 
  
  return v;
}


double *t_polynomialDim0 ( int m, int n, double x[], double delta,int dim0 ){
  int i;
  int j;
  double *v;
  
  double *vaux;

  if ( n < 0 )
  {
    return NULL;
  }

  v = new double[m*(n+1)];
  
  vaux = new double[m*(n+1)];

  for ( i = 0; i < m; i++ )
  {
    v[i] = 1.0;
    
    vaux[i] = 1.0;
  }
  if ( n < 1 )
  {
    delete[] vaux; 
      
    return v;
  }

  for ( i = 0; i < m; i++ )
  {
    vaux[i+1*m] = x[i];   
      
    v[i+1*m] = vaux[i+1*m]*pow(2.0,(-2.0-delta));    
    
  }
  for ( j = 2; j <= n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      vaux[i+j*m] = 2.0 * x[i] * vaux[i+(j-1)*m] - vaux[i+(j-2)*m];
      
      v[i+j*m] = vaux[i+j*m]*pow((1.0+j),-2.0-delta); 
    }
  }
  
  for (j = 0; j < dim0; j++)      
  {
    for ( i = 0; i < m; i++ )
    {
        v[i+j*m] = vaux[i+j*m];      
        
    }
      
  }
  
  
  delete[] vaux; 
  
  return v;
}


double *u_polynomialD ( int m, int n, double x[],double delta ){
  int i;
  int j;
  double *v,*vaux;

  if ( n < 0 )
  {
    return NULL;
  }

  v = new double[m*(n+1)];

  vaux = new double[m*(n+1)];
  
  for ( i = 0; i < m; i++ )
  {
    v[i+0*m] = 0.0;
    
    vaux[i+0*m] = 0.0; 
  }

  if ( n < 1 )
  {
    delete[] vaux;   
      
    return v;
    
  }
  
  for ( i = 0; i < m; i++ )
  {
    v[i+1*m] = pow(2.0,(-2.0-delta));
    
    vaux[i+1*m] = 1.0;
  }

  if ( n < 2 )
  {
    delete[] vaux;   
      
    return v;
    
  }
  
  for ( i = 0; i < m; i++ )
  {
    vaux[i+2*m] = 2.0 * x[i];
    
    v[i+2*m]  =  2.0*vaux[i+2*m]*pow(3.0,(-2.0-delta)); 
  }

  for ( i = 0; i < m; i++ )
  {
    for ( j = 3; j <= n; j++ )
    {
      vaux[i+j*m] = 2.0 * x[i] * vaux[i+(j-1)*m] - vaux[i+(j-2)*m];
      
      v[i+j*m] = ((double) j)*vaux[i+j*m]*pow(j+1.0,(-2.0-delta)); 
    }
  }
  
  return v;
  
  delete[] vaux;  
  
}


double *u_polynomialDim0 ( int m, int n, double x[],double delta, int dim0 ){
  int i;
  int j;
  double *v,*vaux;

  if ( n < 0 )
  {
    return NULL;
  }

  v = new double[m*(n+1)];

  vaux = new double[m*(n+1)];
  
  for ( i = 0; i < m; i++ )
  {
    v[i+0*m] = 0.0;
    
    vaux[i+0*m] = 0.0; 
  }

  if ( n < 1 )
  {
    delete[] vaux;   
      
    return v;
    
  }
  
  for ( i = 0; i < m; i++ )
  {
    v[i+1*m] = pow(2.0,(-2.0-delta));
    
    vaux[i+1*m] = 1.0;
  }

  if ( n < 2 )
  {
    delete[] vaux;   
      
    return v;
    
  }
  
  for ( i = 0; i < m; i++ )
  {
    vaux[i+2*m] = 2.0 * x[i];
    
    v[i+2*m]  =  2.0*vaux[i+2*m]*pow(3.0,(-2.0-delta)); 
  }

  for ( i = 0; i < m; i++ )
  {
    for ( j = 3; j <= n; j++ )
    {
      vaux[i+j*m] = 2.0 * x[i] * vaux[i+(j-1)*m] - vaux[i+(j-2)*m];
      
      v[i+j*m] = ((double) j)*vaux[i+j*m]*pow(j+1.0,(-2.0-delta)); 
    }
  }
  
  for (j = 0; j < dim0; j++)      
  {
    for ( i = 0; i < m; i++ )
    {
        v[i+j*m] = vaux[i+j*m];      
        
    }
      
  }
  
  for (j = 0; j < dim0; j++)      
  {
    for ( i = 0; i < m; i++ )
    {
        v[i+j*m] = ((double) j)*vaux[i+j*m];      
        
    }
      
  }
  
  
  return v;
  
  delete[] vaux;  
  
  
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



dcomp pVNs0(int m ,double* OpdR, double* OpdI, double* TN ,int Nq,double wq,
        int Nq0, double wq0, int N )
{
    
    
    double ur=0;
    double ui=0; 
    
    double auxr=0; 
    double auxi=0; 
//     mexPrintf("%d %f , %f | %f \n", m,OpdR[0],OpdI[0], TN[2]);  
// //     
    
    for(int ll(0); ll< Nq0; ++ll)
    {
//         for(int jj(0); jj < Nq; ++jj)
//         {
//             ur += OpdR[ll*Nq+jj]*TN[jj+m*Nq];
//             
//             ui += OpdI[ll*Nq+jj]*TN[jj+m*Nq];
//         }
        
        
        ur+= inner_product((OpdR+ll*Nq),(OpdR+(ll+1)*Nq),(TN+m*Nq),auxr);
        
        ui+= inner_product((OpdI+ll*Nq),(OpdI+(ll+1)*Nq),(TN+m*Nq),auxi);
    }
    

    
//       mexPrintf("%f , %f \n",wq,wq0);
             
//     dcomp u; 
// //    
//     u.real() = ur; 
//     
//     u.imag() = ui; 
    
//      u*=wq*wq0; 
    
    ur*=wq*wq0;
    
    ui*=wq*wq0;
    
    
    return dcomp(ur,ui); 
    
    
}

dcomp pldNs(int l,dcomp* ft_cofs,int dom_op){

    if(l==0){
        return pi*(ft_cofs[0]);
    }else{
        return 0.5*pi*(ft_cofs[l]);
    }      
}


    
