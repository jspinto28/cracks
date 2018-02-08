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



    
///////////funciones auxiliares (multiplicadores de kernells)/////////////////

void pk(double* dk, double* pk, double t, double s,  int domt, int doms, int dom_op){
    double xt,yt,xs,ys,nx,ny;    
    geom(&xt,&yt,t,domt,dom_op);
    geom(&xs,&ys,s,doms,dom_op);
    normal(&nx,&ny,t,domt,dom_op);
    *dk=sqrt(pow(xt-xs,2)+pow(yt-ys,2));
    *pk=((xt-xs)*nx+(yt-ys)*ny)/(*dk);    
} 

void pka(double* dk, double* pk, double t, double s,  int domt, int doms, int dom_op){
    double xt,yt,xs,ys,nx,ny;    
    geom(&xt,&yt,t,domt,dom_op);
    geom(&xs,&ys,s,doms,dom_op);
    normal(&nx,&ny,s,doms,dom_op);
    *dk=sqrt(pow(xt-xs,2)+pow(yt-ys,2));
    *pk=((xs-xt)*nx+(ys-yt)*ny)/(*dk);    
} 

void pw(double* dw, double* pw, double t, double s,  int domt, int doms, int dom_op){
    double nxt,nyt,nxs,nys; 
    normal(&nxt,&nyt,t,domt,dom_op);
    normal(&nxs,&nys,s,doms,dom_op);
    *dw=dist(t,s,domt,doms,dom_op);
    *pw=nxt*nxs+nyt*nys;
} 


////////////////////////////////////////////////////////////////////////////////////////

//////Kernells//////////////////////////////////////////////////////////////////


dcomp G0(double t, double s){
	return dcomp(-1/(2*pi)*log(abs(t-s)),0); 
}

dcomp GB(double t, double s, int domt, int doms, int dom_op,
        GreenFunctionBase & gf){
    
	if(domt==doms){
		return gf.evaluateGF(t,s,domt,doms,dom_op)
        - gf.evaluateRegularizator(t,s); 
	}else{
		return gf.evaluateGF(t,s,domt,doms,dom_op);
	}
}

dcomp KB(double t, double s, int domt, int doms, int dom_op,
        GreenFunctionBase & gf){
    
    return gf.evaluateGFDetNx(t,s,domt,doms,dom_op);
    
}

dcomp KAB(double t, double s, int domt, int doms, int dom_op,  
        GreenFunctionBase & gf){    
    
    return gf.evaluateGFDetNy(t,s,domt,doms,dom_op);
    
}   

dcomp WB(double t, double s, int domt, int doms, int dom_op,
        GreenFunctionBase & gf){
    
   	double d,p;
    
   	pw(&d,&p,t,s,domt,doms,dom_op);
    
    dcomp k = gf.GetWaveNumber(); 
    
	if(domt==doms){
		return pow(k,2)*(-1.0*gf.evaluateGF(t,s,domt,doms,dom_op)*p
                +gf.evaluateRegularizator(t,s));
	}else{
		return -1.0*pow(k,2)*gf.evaluateGF(t,s,domt,doms,dom_op)*p;
	}
}

dcomp GBB(double t, double s, int domt, int doms, int dom_op,
        GreenFunctionBase & gf){

	return gf.evaluateGF(t,s,domt,doms,dom_op)
        - gf.evaluateRegularizator(t,s); 

}

dcomp WA1( double t, double s, int domt, int doms, int dom_op,
        GreenFunctionBase & gf){
    
        dcomp u(0,0);    
    
        double js=J(s,doms,dom_op);
        
        u= gf.evaluateGFDevy(t,s,domt,doms,dom_op)/js; 

     	if(domt==doms){         
            
            double jt=J(t,domt,dom_op);
            
            u-= 0.5/pi/(t-s)/jt;
        }
     
        return u; 

}

dcomp WA2(double t, double s, int domt, int doms, int dom_op,
        GreenFunctionBase & gf){
    
    double js;
       
    js=J(s,doms,dom_op);
    
	if(domt==doms){
        
        double jt = J(t,domt,dom_op);
        
		return gf.evaluateGF(t,s,domt,doms,dom_op)/js-
               gf.evaluateRegularizator(t,s)/jt;
	}else{
		return gf.evaluateGF(t,s,domt,doms,dom_op)/js;
	}    
}

dcomp WA3( double t, double s, int domt, int doms, int dom_op,
        GreenFunctionBase & gf){
    double js,jsp,p; 
    
    js=J(s,doms,dom_op);    
    
    jsp=Jp(s,doms,dom_op);
    
    p=jsp*pow(js,-2);
    
	if(domt==doms){
        
        double jt,jtp,p0;
        
        jt=J(t,domt,dom_op);
    
        jtp=Jp(t,domt,dom_op);
        
        p0 = jtp*pow(jt,-2);
        
		return gf.evaluateGF(t,s,domt,doms,dom_op)*p-
                -gf.evaluateRegularizator(t,s)*p0;
        
	}else{
        
		return gf.evaluateGF(t,s,domt,doms,dom_op)*p;
	}       
}
        
dcomp WB1(double alpha, double period,int sig, double t, 
        double s, int domt, int doms, int dom_op,
        GreenFunctionBase & gf){
    
    //Periodic case!!
    
    int Intefs = NumOfIntefaces( dom_op);  
    
    dcomp u(0,0);
    
    if(Intefs == 1)    
    {
               
        double js=J(s,doms,dom_op);
       
        double jt = J(t,doms,dom_op);
    
        double jnt = J(-t,doms,dom_op);

        double k = gf.GetWaveNumber().real();
        
        int au = sig*Orientation(doms,dom_op); 

        double auxx= (-1.0*(double)au)*k*cos(alpha)*period;

        dcomp quasi (cos(auxx),sin(auxx));    

        u= gf.evaluateGF(t,s,domt,doms,dom_op)/js-
                gf.evaluateRegularizator(t,s)/jt
                 -quasi*gf.evaluateRegularizator(-t,s)/jnt;             
        
        
    }
    else
    {    

        double js; 
        int dmx;
        int sign=1; 
        int bol=reg(domt,doms,dom_op);
        if((abs(bol)==1)&&(sig*bol==-1)){
            bol=-2; 
        }
        if(bol>-2){ 
            dmx=doms;        
        }else{
            dmx=domt;
        }       
        if(abs(bol)<2){
            sign=-1; 
        }

        js=J(s,doms,dom_op);
        if(bol>-2){

            double jt=J(t,dmx,dom_op);

            u = gf.evaluateGF(t,s,domt,doms,dom_op)/js-
                    gf.evaluateRegularizator((double)sign*t,s)/jt;

        }else{
            u =   gf.evaluateGF(t,s,domt,doms,dom_op)/js; 
        }
        
        if(abs(domt-doms) == Intefs) 
       //we have an extra singularity of the periodic in the (0,P) or (P,0)                                    
        {
            double jnt = J(-t,doms,dom_op);

            double k = gf.GetWaveNumber().real();

            int orientation = Orientation(doms,dom_op); 

            double auxx= (-(double) (orientation*sig))*k*cos(alpha)*period;

            dcomp quasi (cos(auxx),sin(auxx));       
            
            u -= quasi*gf.evaluateRegularizator(-t,s)/jnt; 

        }
//  
    }
    
//     mexPrintf("%f , %f \n", u.real(), u.imag()); 
    
    return u;
}


// dcomp WB1(double alpha, double period,int sig, double t, 
//         double s, int domt, int doms, int dom_op,
//         GreenFunctionBase & gf){
//     //THIS IS ONLY VALID FOR THE CASE OF ONE INRERFACE !!!!
//     
//     double js=J(s,doms,dom_op);
//        
//     double jt = J(t,doms,dom_op);
//     
//     double jnt = J(-t,doms,dom_op);
//     
//     double k = gf.GetWaveNumber().real();
//     
//     double auxx= (-(double) sig)*k*cos(alpha)*period;
//     
//     dcomp quasi (cos(auxx),sin(auxx));    
// 
//     return gf.evaluateGF(t,s,domt,doms,dom_op)/js-
//             gf.evaluateRegularizator(t,s)/jt
//             -quasi*gf.evaluateRegularizator(-t,s)/jnt; 
// 
// }


dcomp WB2(int sig, double t, double s, int domt, int doms, int dom_op,
        GreenFunctionBase & gf){
    double js,jsp,p; 
    int dmx;
    int sign=1;
    int bol=reg(domt,doms,dom_op); 
    if((abs(bol)==1)&&(sig*bol==-1)){
        bol=-2; 
    }
    if(bol>-2){ 
        dmx=doms;        
    }else{
        dmx=domt;
    }    
    if(abs(bol)<2){
        sign=-1; 
    }
    
    js=J(s,doms,dom_op);
    
    jsp=Jp(s,doms,dom_op);
    
    p=jsp*pow(js,-2);
    
    if(bol>-2){
        
        double jt,jtp,p0; 
        
        jt=J(t,dmx,dom_op);
    
        jtp=Jp(t,dmx,dom_op);
    
        p0=jtp*pow(jt,-2);
        
        
        return gf.evaluateGF(t,s,domt,doms,dom_op)*p-
                      gf.evaluateRegularizator((double)sign*t,s)*p0;
        
    }else{
        
        return gf.evaluateGF(t,s,domt,doms,dom_op)*p;
        
    }
}

//old code without the class for the GF.

// dcomp xG0(double t, double s){
// 	return dcomp(-1/(2*pi)*log(abs(t-s)),0); 
// }
// 
// dcomp xGB(double k,double t, double s, int domt, int doms, int dom_op){
//     double d; 
//     d=dist(t,s,domt,doms,dom_op);      
// 	if(domt==doms){
// 		return 0.25*i*boost::math::cyl_hankel_1(0,k*d)-xG0(t,s);
// 	}else{
// 		return 0.25*i*boost::math::cyl_hankel_1(0,k*d);
// 	}
// }
// 
// dcomp xKB(double k,double t, double s, int domt, int doms, int dom_op){
//     double d,p; 
//   	pk(&d,&p,t,s,domt,doms,dom_op);
// 	return -k*0.25*i*boost::math::cyl_hankel_1(1,k*d)*p;
// }
// 
// dcomp xKAB(double k,double t, double s, int domt, int doms, int dom_op){
//     double d,p; 
//    	pka(&d,&p,t,s,domt,doms,dom_op);
// 	return -k*0.25*i*boost::math::cyl_hankel_1(1,k*d)*p;
// }
// 
// 
// dcomp xWB( double k,double t, double s, int domt, int doms, int dom_op){
//    	double d,p;
//    	pw(&d,&p,t,s,domt,doms,dom_op);
// 	if(domt==doms){
// 		return -pow(k,2)*0.25*i*boost::math::cyl_hankel_1(0,k*d)*p+pow(k,2)*xG0(t,s);
// 	}else{
// 		return -pow(k,2)*0.25*i*boost::math::cyl_hankel_1(0,k*d)*p;
// 	}
// }
// 
// dcomp xGBB(double k,double t, double s, int domt, int doms, int dom_op){
//     double d; 
//     d=dist(t,s,domt,doms,dom_op);
// 	return 0.25*i*boost::math::cyl_hankel_1(0,k*d)-xG0(t,s);
// 
// }
// 
// dcomp xWA1(double k, double t, double s, int domt, int doms, int dom_op){
//     double d,jt,js,p; 
//     d=dist(t,s,domt,doms,dom_op);
//     jt=J(t,domt,dom_op);
//     js=J(s,doms,dom_op);
//     p=1/(jt*js);
// 	if(domt==doms){
// 		return 0.25*i*boost::math::cyl_hankel_1(0,k*d)*p-xWA01(t,s,domt,dom_op);
// 	}else{
// 		return 0.25*i*boost::math::cyl_hankel_1(0,k*d)*p;
// 	}
// }
// 
// dcomp xWA2(double k, double t, double s, int domt, int doms, int dom_op){
//     double d,jt,js,jtp,p; 
//     d=dist(t,s,domt,doms,dom_op);    
//     js=J(s,doms,dom_op);
// 	if(domt==doms){
// 		return 0.25*i*boost::math::cyl_hankel_1(0,k*d)/js-xWA02(t,s,domt,dom_op);
// 	}else{
// 		return 0.25*i*boost::math::cyl_hankel_1(0,k*d)/js;
// 	}    
// }
// 
// dcomp xWA3(double k, double t, double s, int domt, int doms, int dom_op){
//     double d,jt,js,jsp,p; 
//     d=dist(t,s,domt,doms,dom_op);
//     js=J(s,doms,dom_op);    
//     jsp=Jp(s,doms,dom_op);
//     p=jsp*pow(js,-2);
// 	if(domt==doms){
// 		return 0.25*i*boost::math::cyl_hankel_1(0,k*d)*p-xWA03(t,s,domt,dom_op);
// 	}else{
// 		return 0.25*i*boost::math::cyl_hankel_1(0,k*d)*p;
// 	}       
// }
// 
// dcomp xWA4(double k, double t, double s, int domt, int doms, int dom_op){
//     double d,jt,js,jtp,jsp,p; 
//     d=dist(t,s,domt,doms,dom_op);    
//     js=J(s,doms,dom_op);    
//     jsp=Jp(s,doms,dom_op);
//     p=jsp*pow(js,-2); 
// 	if(domt==doms){
// 		return 0.25*i*boost::math::cyl_hankel_1(0,k*d)*p-xWA04(t,s,domt,dom_op);
// 	}else{
// 		return 0.25*i*boost::math::cyl_hankel_1(0,k*d)*p;
// 	}          
// }
// 
// dcomp xWA01(double t, double s, int domt, int dom_op){
//     double jt; 
//     jt=J(t,domt,dom_op);
//     return xG0(t,s)*pow(jt,-2);
// }
// 
// dcomp xWA02(double t, double s, int domt, int dom_op){  
//     double jt;
//     jt=J(t,domt,dom_op);
//     return xG0(t,s)/jt;
// }
// 
// dcomp xWA03(double t, double s, int domt, int dom_op){  
//     double jt,jtp;
//     jt=J(t,domt,dom_op);
//     jtp=Jp(t,domt,dom_op);
//     return xG0(t,s)*jtp*pow(jt,-2);
// }
// 
// dcomp xWA04(double t, double s, int domt, int dom_op){
//     double jt,jtp;
//     jt=J(t,domt,dom_op);
//     jtp=Jp(t,domt,dom_op);
//     return xG0(t,s)*jtp*pow(jt,-2);
// }
//         
// dcomp xWB1(int sig,double k, double t, double s, int domt, int doms, int dom_op){
//     double d,js; 
//     int dmx;
//     int sign=1; 
//     int bol=reg(domt,doms,dom_op);
//     if((abs(bol)==1)&&(sig*bol==-1)){
//         bol=-2; 
//     }
//     if(bol>-2){ 
//         dmx=doms;        
//     }else{
//         dmx=domt;
//     }       
//     if(abs(bol)<2){
//         sign=-1; 
//     }
//     d=dist((double)sign*t,s,dmx,doms,dom_op);    
//     js=J(s,doms,dom_op);
//     if(bol>-2){
//         return 0.25*i*boost::math::cyl_hankel_1(0,k*d)/js-xWB01((double)sign*t,s,dmx,dom_op);
//     }else{
//         return  0.25*i*boost::math::cyl_hankel_1(0,k*d)/js; 
//     }
// }
// 
// dcomp xWB2(int sig,double k, double t, double s, int domt, int doms, int dom_op){
//     double d,js,jsp,p; 
//     int dmx;
//     int sign=1;
//     int bol=reg(domt,doms,dom_op); 
//     if((abs(bol)==1)&&(sig*bol==-1)){
//         bol=-2; 
//     }
//     if(bol>-2){ 
//         dmx=doms;        
//     }else{
//         dmx=domt;
//     }    
//     if(abs(bol)<2){
//         sign=-1; 
//     }
//     d=dist((double)sign*t,s,dmx,doms,dom_op);    
//     js=J(s,doms,dom_op);
//     jsp=Jp(s,doms,dom_op);
//     p=jsp*pow(js,-2);
//     if(bol>-2){
//         return (0.25*i*boost::math::cyl_hankel_1(0,k*d)*p-xWB02((double)sign*t,s,dmx,dom_op));
//     }else{
//         return (0.25*i*boost::math::cyl_hankel_1(0,k*d)*p);
//     }
// }
// 
// dcomp xWB01(double t, double s, int domt, int dom_op){  
//     double jt; 
//     jt=J(t,domt,dom_op);
// 	return xG0(t,s)/jt;
// }
// 
// dcomp xWB02(double t, double s, int domt, int dom_op){
//     double jt,jtp; 
//     jt=J(t,domt,dom_op);
//     jtp=Jp(t,domt,dom_op);
// 	return xG0(t,s)*jtp*pow(jt,-2);
// }

//////////////////////////////////////////////////////////////////////////


///////////////////////////lado derecho/////////////////////////////////

// dcomp uinc_d(double k0, double s, int doms ,int dom_op, double alpha){
//     double x,y; 
//     geom(&x,&y,s,doms,0);
//     return exp(i*k0*(x*cos(alpha)+y*sin(alpha))); 
// }
// 
// dcomp uinc_n(double k0, double s, int doms ,int dom_op, double alpha){
//     double nx,ny; 
//     normal(&nx,&ny,s,doms,0);
//     return i*k0*uinc_d(k0,s,doms,dom_op,alpha)*(cos(alpha)*nx+sin(alpha)*ny); 
//    }
// //ojo signo.


dcomp uinc_d(double k0, double s, int doms ,int dom_op, int l, double R){
    double x,y,norm,unitx,unity,argx; 
    
    geom(&x,&y,s,doms,0);
      
    norm = sqrt(x*x+y*y); 
    
    unitx=x/norm; 
    
    unity=y/norm; 
        
     argx = Argument(unitx,unity);
    
//     argx = atan(unity/unitx); 
       
    dcomp H,Y;
    
    double J; 
            
//     H=boost::math::cyl_hankel_1(abs(l),k0*R);
    
    H=1.0; 
     
    J= boost::math::cyl_bessel_j(l,k0*norm); 
    
    Y = 1/sqrt(2*pi)*exp(i*(double)l*argx); 
           
    return Y*J; 
    
}


double Argument(double x, double y)
{
 
    if((x>=0)&&(y>=0))
        return atan(y/x); 
    else if((x<0)&&(y>0))
        return pi-atan(y/(-x)); 
    else if((x>0) &&(y<0))
        return 2*pi-atan((-y)/x); 
    else 
        return pi+atan(y/x);       
               
    
    
}
// 
// dcomp pointSourcePer(double k0, double s, int doms,int dom_op ,double alpha,
//          int NPer, double* an, double* bn){
//     
// 
//     double x,y; 
//     geomPer(&x,&y,s,dom_op,an,bn,NPer);
// 
//     return exp(i*k0*(x*cos(alpha)-y*sin(alpha))); 
// 
//     
// }



dcomp rhsNs(int domt, double s,double alpha, double kext,
        double* Y1, double* Y2, double delta, int S){
    
 
      double x,y;  
      
      x=0; 
      y=0; 
    
      geomPer(&x,&y,s,domt,Y1,Y2,delta,S);    

      return exp(Ii*kext*(cos(alpha)*x-sin(alpha)*y));
      
}



