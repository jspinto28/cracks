#include <math.h>
#include <complex>
#include "geo.h"
#include <boost/math/special_functions/hankel.hpp>
#include "mex.h"
#include "prods.h"
#include <numeric>

using namespace std;
const double pi= 4.0*atan(1.0);

typedef complex<double> dcomp;
const dcomp i(0.0,1.0);

const int  Nd = 1; 

const double Cperiod = 2.0;

const double Amplitud = 0.0;


/// num of interfaces for a given domain, used in W for the periodic case. 


void geom(double* x, double* y, double t, double* b, int Nb, double delta){
  //entrega la geometria como funcion de t en -1,1.  
    
//     int dim0=10; 
    
    *x=t; 
    
    double *T_vec;
    
//     T_vec = t_polynomialD ( 1, Nb-1, &t ,delta);
    
        T_vec = t_polynomialD( 1, Nb-1, &t ,delta);
    
    double aux = 0; 
    
//     for(int ii=0; ii<Nb; ii++)
//     {
//       T_vec[ii] *= pow((double)(1+ii),(-2.0-delta));   
//     }
//         
    *y = inner_product(T_vec,T_vec+Nb,b,aux);
    
    delete [] T_vec;       
}



void geomp(double* x, double* y, double t,  double* b, int Nb, double delta){
  //entrega la derivada de la geometria como funcion de t en -1,1. 
    
    int dim0 = 10; 
    
    *x=1; 
    
    double *Tp_vec;
    
    Tp_vec = u_polynomialD ( 1, Nb-1, &t,delta );    
    
//     Tp_vec = u_polynomialDim0 ( 1, Nb-1, &t,delta );    
   
    double aux = 0; 
    
    *y = inner_product(Tp_vec,Tp_vec+Nb,b,aux);
    
    delete[] Tp_vec;

    
}

double J(double t,int domt, int dom_op, double* b, int Nb, double delta){
 double x,y;
 
 geomp(&x,&y,t, b, Nb,delta); 
 
 return sqrt(pow(x,2)+pow(y,2));    
}





