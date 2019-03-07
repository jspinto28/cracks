#include <math.h>
#include <complex>
#include "kers.h"
#include "geo.h"
#include <boost/math/special_functions/hankel.hpp>
#include "mex.h"


#include "greenfunction.h"

using namespace std;
const double pi= 4.0*atan(1.0);
// 
//  typedef complex<double> dcomp;
  const dcomp Ii(0.0,1.0);


dcomp rhsNs(double k,double angle,int domt, double s,double* bgeo, int Nb, double delta){
    
    
      double x,y; 
// //     
      geom(&x,&y,s,bgeo,Nb,delta);      
      
      return GlobalRhsNs(k,angle,x,y);
    
    
}

dcomp GlobalRhsNs(double k,double angle,double x, double y){
    
    double kr=k; 
    
//     if(k < 0.0000001)
//     {
//      //   return x; 
//          kr=5.0; 
//     }       
  
    return exp(-Ii*kr*(cos(angle)*x+sin(angle)*y));
    

    
    
}

dcomp DevrhsNs(double k,double angle,int domt, double s,double* bgeo,
        int Nb, double delta,int p,int component){
    
    
      double x,y;
      
      dcomp g;
// //     
      geom(&x,&y,s,bgeo,Nb,delta);      
      
      g = GlobalRhsNs(k,angle,x,y)*cos(p*acos(s));
      
      if( component ==0 ) 
      {
          
          return -Ii*k*cos(angle)*g; 
      
      }
      else 
      {    
          
          return -Ii*k*sin(angle)*g; 
      }
    
}

dcomp FarOne(double k,double s,double* bgeo,
        int Nb, double delta,int p,int component,
        double y1, double y2){
    
    
      double x,y,aux;
      
//       aux = y1; 
//       
//       if(component > 0) 
//       {
//           aux = y2;
// // //  
//       }
      
      geom(&x,&y,s,bgeo,Nb,delta);  
      
      return exp(-Ii*k*(x*y1+y*y2))*cos(p*acos(s));       
   
}


dcomp GBNs(double t, double s, int domt, int doms,GreenFunctionBase & gf)
{
    if(domt==doms){
        
        if(abs(t-s)>0)
        {
        
            return gf.evaluateGF(t,s,domt,doms,0)
            - gf.evaluateRegularizator2(t,s,doms); 
        
        }
        else
        {
            return gf.LimitRegularized(t,s,doms );
        }
        
        
	}else{
		return gf.evaluateGF(t,s,domt,doms,0);
	}
    
}
