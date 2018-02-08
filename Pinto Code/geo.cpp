#include <math.h>
#include <complex>
#include "geo.h"
#include "prods.h" 
#include <boost/math/special_functions/hankel.hpp>
#include "mex.h"
// #include "Curves.h"

using namespace std;
const double pi= 4.0*atan(1.0);

typedef complex<double> dcomp;
const dcomp i(0.0,1.0);

const int  Nd = 1; 

const double Cperiod = 2.0;

const double Amplitud = 0.0;


/// num of interfaces for a given domain, used in W for the periodic case. 

int NumOfIntefaces(int dom_op)
{
 return 1;    
    
}

////////////////Geometria dada la parametrizacion/////////////////////////

// void curves(double* x, double* y,double t,int curv){
//     
//     if(curv==1) 
//         
//     {
//         *x = t; 
//         
//         *y = sin(pi*t);            
//         
//     }
//     else if (curv == 2) 
//     {
//         *x = 0.5*t-1; 
//         
//         *y= -2; 
//         
//     }
//     else if( curv  ==3) 
//     {
//         *x = t+2; 
//         
//         *y= 2; 
//         
//     }
        
    
    
// 	*x=(-t+1.0)*Cperiod/(2.0);
//     
// 	*y=Amplitud*cos(2.0*pi*(*x)/Cperiod);
    
//     double per = 0.5*Cperiod; 
//     
//     *x = (t+1.0)*per-0.5*per; 
//     
//     *y = Amplitud*cos(2.0*pi*(*x+0.5*per)/per+pi);
    
//     *x = t; 
//     
//     *y=0; 
// }

// void tang_curves(double* x, double* y,double t,int curv){
//     
//         if(curv==1) 
//         
//     {
//         *x = 1; 
//         
//         *y =pi*cos(pi*t);            
//         
//     }
//     else if (curv == 2) 
//     {
//         *x =0.5; 
//         
//         *y= 0; 
//         
//     }
//     else if( curv  ==3) 
//     {
//         *x = 1; 
//         
//         *y= 0; 
//         
//     }
        
    
    
// 	*x=(-t+1.0)*Cperiod/(2.0);
//     
//     double dx = -Cperiod/(2.0);
//     
// 	*y= -Amplitud*2*pi*(dx)/Cperiod*sin(2.0*pi*(*x)/Cperiod);
//     
//     *x = dx;
    
//     double per = 0.5*Cperiod; 
//     
//     double dx = per; 
//     
//     *x = (t+1.0)*per-0.5*per; 
//     
//     *y = -Amplitud*2.0*pi/per*dx*sin(2.0*pi*(*x+0.5*per)/per+pi);
//     
//     *x=dx; 
    
    
//     *x = 1;
//     
//     *y = 0; 
            
//     
//     
// }

void secondd_curves(double* x, double* y,double t,int curv){
    
    *x= 0; 
    
    *y=0; 
    
//     *x=(-t+1.0)*Cperiod/(2.0);
//     
//     double dx = -1.0*Cperiod/(2.0);
//     
// 	*y= -Amplitud*pow(2.0*pi*dx/Cperiod,2.0)*cos(2.0*pi*(*x)/Cperiod);
//     
//     *x=0.0; 
    
//     double per = 0.5*Cperiod; 
//     
//     double dx = per; 
//     
//     *x = (t+1.0)*per-0.5*per; 
//     
//     *y = -Amplitud*pow(2.0*pi/per*dx,2)*cos(2.0*pi*(*x+0.5*per)/per+pi);
//     
//     *x=0.0; 
    
//     
//     *x=0.0; 
    
//     *x=0; 
//     
//     *y=0; 
}


void geomp(double* x, double* y, double t, int domt, int dom_op){
  //entrega la derivada de la geometria como funcion de t en -1,1.    
//     if(dom_op==0){
//         
//         tang_curves(x,y,-t,domt+1);   
//         
//             *x=(-1.0)*(*x);
//             
//             *y=(-1.0)*(*y);
//             
//     }else{
//         
//         tang_curves(x,  y, t,dom_op);
//         
//     }      
//     tang_curves(x,  y, t,domt+1);
  
    return; 
}

void geom(double* x, double* y, double t, int domt, int dom_op){
  //entrega la geometria como funcion de t en -1,1.  
    
    
//     if(dom_op==0){
//         
//         curves( x,  y, -t, domt+1);     
//     
//     }else{
//         
//         curves( x,  y, t, domt+1);
// 
//     }    
//     curves( x,  y, t, domt+1); 
  
    return; 
          
}




double inv_geom(double x, double y, int domt, int dom_op){
    
    //
    
    return x; 
    
}

  

void geompp(double* x, double* y, double t, int domt, int dom_op){
    
    if(dom_op==0){
        
        secondd_curves( x,  y, -t, domt+1);     
    
    }else{
        
        secondd_curves( x,  y, t, dom_op);

    }    

}

int* trial_test(int dom_op_trial, int dom_op_test,int Dt_trial, int Dt_test){
    int *tt= new int[Dt_trial]; 
    for (int j=0; j< Dt_trial; j++){
        for (int l=0; l<Dt_test; l++){   
            double xt,yt,xs,ys,d; 
            geom(&xt,&yt,0,j,dom_op_trial); 
            geom(&xs,&ys,0,l,dom_op_test); 
            d=sqrt(pow(xt-xs,2)+pow(yt-ys,2));
            if(d<pow(10,-10)){
                tt[j]=l; 
                l=Dt_test; 
            }else{
                tt[j]=-1; 
            }        
       }
   }
    return tt; 
}




int reg( int domt, int doms, int dom_op){
    //funcion que determina si se regulariza o no la integral W3, que corresponde 
    //a la integral de los valores de frontera (t=-1,t=1) de Vu_trial(t)
    // 2 se regulariza completamente 
    // 1/-1 se regulariza solo el -1 o el 1 (entrega cual se regulariza)
    //-2 es disjunta no es necesario regularizar. 
    // 0 si se regularian ambos pero no coinsiden los dominios (caso circulo)        
     if(domt==doms){
         return 2;
     }else if ((dist(1,-1,domt,doms,dom_op)<pow(10,-10))&&(dist(-1,1,domt,doms,dom_op)<pow(10,-10))){
         return 0; 
     }else if(dist(1,-1,domt,doms,dom_op)<pow(10,-10)){
         return 1; 
     }else if(dist(-1,1,domt,doms,dom_op)<pow(10,-10)){
         return -1; 
     }else{ 
         return -2; 
     }    
}
//Dummy function. 
int* trial(int Dt){  
 int *tt =new int[Dt];
 for(int j=0; j< Dt; j++){
     tt[j]=j;
 }  
 return tt;       
}

////////////////////////////////////////////////////////////////////////////


//////Jacobianos///////////////////////////////////////////////////////////

double J(double t,int domt, int dom_op){
 double x,y;
 geomp(&x,&y,t,domt,dom_op); 
 return sqrt(pow(x,2)+pow(y,2));    
}

double Jp(double t,int domt, int dom_op){
       
 double x,y;
 geomp(&x,&y,t,domt,dom_op); 
 
 double J= sqrt(pow(x,2)+pow(y,2)); 
 
 double xp, yp; 
 geompp(&xp,&yp,t,domt,dom_op); 
 
 return (xp*x+yp*y)/J;
 
}


///////Otras funciones/////////////////////////////////////////////////////
    
double  dist(double t, double s, int domt, int doms, int dom_op){ 
    double xt,yt,xs,ys; 
    geom(&xt,&yt,t,domt,dom_op); 
    geom(&xs,&ys,s,doms,dom_op); 
    return sqrt(pow(xt-xs,2)+pow(yt-ys,2));
}



void normal(double* nx, double* ny,double t, int domt, int dom_op){
    double x,y,j; 
       
    geomp(&x,&y,t,domt,dom_op);
    j=sqrt(pow(x,2)+pow(y,2));
//      if(dom_op==0){
//          *nx =-y/j;
//          *ny =x/j; 
//      }else{
//          *nx =y/j;
//          *ny =-x/j;     
//      }
   *nx =y/j;
   *ny =-x/j;      


}

//for the periodic problem, the quasiperiodic do not depeend of the orientation
int Orientation(int domt, int dom_op)
{
    
    int  u=0;
    
    if(dom_op ==0) 
    {
        u= -1;
    }
    else
    {
        u= 1; 
    }       

    return u; 
}    


///////////////////////////////////////////////////////////////////////////////////

//perturbable geo: 


void curvesPer(double* x, double* y,double t,
        double* an, double* bn,int NPer){
    
	*x=0.5*(t+1)*Cperiod; 
      
	*y=Amplitud*cos(pi*(t+1));
    
    int Np = NPer/2;
    
    for(int ii(0); ii<Np; ++ii)
    {
     *x += an[ii]*sin((ii+1)*pi*(t+1)); 
     
     *y += bn[ii]*sin((ii+1)*pi*(t+1));
        
    }
    
    for(int ii(0); ii<Np; ++ii)
    {
     *x += an[ii+Np]*cos((ii+1)*pi*(t+1)); 
     
     *y += bn[ii+Np]*cos((ii+1)*pi*(t+1));
        
    }

}

void DcurvesPer(double* x, double* y,double t,
        double* an, double* bn,int NPer){
    
	*x=0.5*Cperiod; 
      
	*y=-Amplitud*pi*sin(pi*(t+1));
    
    int NP = NPer/2;
    
    for(int ii(0); ii<NP; ++ii)
    {
     *x += an[ii]*(ii+1)*pi*cos((ii+1)*pi*(t+1)); 
     
     *y += bn[ii]*(ii+1)*pi*cos((ii+1)*pi*(t+1)); 
        
    }
    
    for(int ii(0); ii<NP; ++ii)
    {
     *x -= an[ii+NP]*(ii+1)*pi*sin((ii+1)*pi*(t+1)); 
     
     *y -= bn[ii+NP]*(ii+1)*pi*sin((ii+1)*pi*(t+1)); 
        
    }

}

void DDcurvesPer(double* x, double* y,double t,
        double* an, double* bn,int NPer){
    
	*x=0; 
      
	*y=-Amplitud*pow(pi,2)*cos(pi*(t+1));
    
    int NP = NPer/2;
    
    for(int ii(0); ii<NP; ++ii)
    {
     *x -= an[ii]*pow((ii+1)*pi,2)*sin((ii+1)*pi*(t+1)); 
     
     *y -= bn[ii]*pow((ii+1)*pi,2)*sin((ii+1)*pi*(t+1)); 
        
    }
    
    for(int ii(0); ii<NP; ++ii)
    {
     *x -= an[ii+NP]*pow((ii+1)*pi,2)*cos((ii+1)*pi*(t+1)); 
     
     *y -= bn[ii+NP]*pow((ii+1)*pi,2)*cos((ii+1)*pi*(t+1)); 
        
    }

}

void geomPer(double* x, double* y, double t, int domt,
        double* Y1, double* Y2,double delta,  int S){
    
    for(int ii(0); ii<S; ++ii)
    {
        double cn = 1.0/pow((double)ii+1.0,(2.0+delta)); 
        
        *x += cn*Y1[ii]*chebT(ii,t); 
        
        *y += cn*Y2[ii]*chebT(ii,t); 
        
    }
    
    

}

void geompPer(double* x, double* y, double t, int domt,
        double* Y1, double* Y2,double delta,  int S){
    
    

    
    for(int ii(1); ii<S; ++ii)
    {
        double dn = ((double)ii)/pow((double)ii+1.0,(2.0+delta)); 
        
        *x += dn*Y1[ii]*chebU((ii-1),t); 
        
        *y += dn*Y2[ii]*chebU((ii-1),t); 
        
    }
    
    

}

void geomppPer(double* x, double* y, double t,int dom_op,
        double* an, double* bn, int NPer){
    
    if(dom_op==0){
        
         DDcurvesPer( x,  y, t, an,bn,NPer);       
    
    }else{
        
         DDcurvesPer( x,  y, -t, an,bn,NPer);   
    }    

}

double JPer(double t,int dom_op,
        double* an, double* bn, int NPer){
    
 //double x,y;
 
 //geompPer(&x,&y,t,dom_op,an,bn,NPer); 
 
 //return sqrt(pow(x,2)+pow(y,2));    
    
    return 0; 
}


double  distPer(double t, double s,int dom_op,
         double* an, double* bn, int NPer){
    
//     double xt,yt,xs,ys; 
//     
//     geomPer(&xt,&yt,t,dom_op,an,bn,NPer); 
//     
//     geomPer(&xs,&ys,s,dom_op,an,bn,NPer); 
//     
//     return sqrt(pow(xt-xs,2)+pow(yt-ys,2));

return 0; 
}



void normalPer(double* nx, double* ny,double t, int dom_op,
        double* an, double* bn, int NPer){
    
    return; 
//     double x,y,j; 
//        
//     geompPer(&x,&y,t,dom_op,an,bn,NPer); 
//     
//     j=sqrt(pow(x,2)+pow(y,2));
// //      if(dom_op==0){
// //          *nx =-y/j;
// //          *ny =x/j; 
// //      }else{
// //          *nx =y/j;
// //          *ny =-x/j;     
// //      }
//    *nx =y/j;
//    *ny =-x/j;      


}

double coorburePer (double t, int dom_op,
        double* an, double* bn, int NPer){
    
    return 0; 
 
//     double xpp,xp,ypp,yp; 
//     
//     geompPer(&xp,&yp,t,dom_op,an,bn,NPer); 
//     
//     geomppPer(&xpp,&ypp,t,dom_op,an,bn,NPer); 
//     
//     double j= sqrt(pow(xp,2)+pow(yp,2));
//     
//     return abs(xp*ypp-yp*xpp)*pow(j,-3.0); 
    
}