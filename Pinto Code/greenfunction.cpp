#include <boost/math/special_functions/hankel.hpp>
#include "greenfunction.h"
#include <vector>
#include <math.h>  
#include "mex.h"


const double pii= 4.0*atan(1.0);

const double Egamma = 0.577215664901533; 
    

dcomp FreeSpace::IntegralKernel()
{
    if(m_WaveNumber.real() > 0.0000001)
    {
    
        return 0.25*i*boost::math::cyl_hankel_1(0,m_WaveNumber.real()*m_distance);      
    
    }
    else
    {
        return dcomp(-1.0/(2.0*pii)*log(m_distance),0);
    }
}
                    
dcomp FreeSpace::evaluateRegularizator2(double t, double s, int doms) 
{
    
    SetPointsFromGeo(t,s,doms,doms,0); 
    
    ComputeDistance(); 
    
    double d = abs(t-s); 
    
    if(m_WaveNumber.real() < 0.0000001)
    {
    
        return dcomp(-0.5/pii*log(d),0);
    
    }
    else 
    {      
               
       
        return dcomp( -0.5/pii*log(d)*
                    boost::math::cyl_bessel_j
                    (0,m_WaveNumber.real()*m_distance),0);
                  
             
    }
    
}


dcomp FreeSpace::LimitRegularized(double t, double s, int doms )
{
    
    double jac = J(s,doms,0,m_bgeo,m_Nb,m_delta);
    
    dcomp u(0.0,0.0); 
    
    if(m_WaveNumber.real() < 0.0000001)
    {
        u = -0.5/pii*log(jac); 
        
    }
    else 
    {
        double k  = m_WaveNumber.real(); 
             
        u.real() = -0.5/pii*(log(0.5*k*jac)+Egamma);
        
        u.imag() = 0.25; 
    }
    
    
    return u; 
    
    
}

dcomp FreeSpace::Jpart(double t, double s, int doms )
{
    SetPointsFromGeo(t,s,doms,doms,0); 
    
    ComputeDistance(); 
    
    return boost::math::cyl_bessel_j
                    (0,m_WaveNumber.real()*m_distance);

    
}

                    
dcomp FreeSpace::evaluateGF( double t, double s, int domt, int doms,
               int dom_op)
{
    SetPointsFromGeo(t,s,domt,doms,dom_op); 

    ComputeDistance(); 

    return IntegralKernel(); 

}

 
  
  dcomp FreeSpace::evaluatePotential(double t, int domt,
             double y1, double y2,int dom_op)
  {
      SetXFromGeo(t,domt,dom_op); 
      
      setYpoint(y1,y2); 
      
      ComputeDistance(); 
      
      return IntegralKernel(); 
      
      
  }
/////////////////////////////////////////////////////////////////////////
  
   double FreeDerivativeSpace::ChebT(double t)
     {
         if( m_p == 0) 
         {    
             return 1; 
         }
         else
         {       
            return cos(m_p*acos(t));
         }
     }
  
   void FreeDerivativeSpace::ComputeDeltaT(double t, double s)
     {
         m_DeltaTp = ChebT(t) - ChebT(s); 
         
     }
  

dcomp FreeDerivativeSpace::IntegralKernel()
{
    if(m_WaveNumber.real() > 0.0000001)
    {
    
        return  -0.25*i*m_WaveNumber.real()*
                boost::math::cyl_hankel_1(1,m_WaveNumber.real()*m_distance)*
                m_DeltaR[m_component]*m_DeltaTp/m_distance;    
    
    }
    else
    {
        return dcomp(-1.0/(2.0*pii)/(m_distance)*
                 m_DeltaR[m_component]*m_DeltaTp/m_distance,0.0);
    }
}
                    
dcomp FreeDerivativeSpace::evaluateRegularizator2(double t, double s, int doms) 
{
    
    SetPointsFromGeo(t,s,doms,doms,0); 
    
    ComputeDistance(); 
    
    double d = abs(t-s); 
    
    ComputeDeltaT(t, s);
    
    ComputeDeltaR(t,doms,s,doms);
        
    if(m_WaveNumber.real() < 0.0000001)
    {
    
        return 0.0;
    
    }
    else 
    {      
               
       
        return dcomp( 0.5*m_WaveNumber.real()/pii*log(d)*
                    boost::math::cyl_bessel_j
                    (1,m_WaveNumber.real()*m_distance)*
                     m_DeltaR[m_component]*m_DeltaTp/m_distance,0);
                  
             
    }
    
//     return 0.0;
    
}


dcomp FreeDerivativeSpace::LimitRegularized(double t, double s, int doms )
{
    
    double jac = J(s,doms,0,m_bgeo,m_Nb,m_delta);
    
    dcomp u(0.0,0.0); 
    
    double x,y,rp; 
    
    geomp(&x,&y,s,m_bgeo,m_Nb,m_delta); 
    
    rp = x; 
    
    if(m_component ==1)
    {
        rp= y; 
    }
        
    if(m_WaveNumber.real() < 0.0000001)
    {
        u =0.0; 
        
    }
    else 
    {   
        if( m_p >0)
        {
            u.real() = -0.5/pii*rp*pow(jac,-2)*((double)m_p)*
                chebU( m_p-1, s );
        }
        else
        {
            return u; 
        }
        

    }
    
    
    return u; 
    
//     return 0.0; 
    
    
}

dcomp FreeDerivativeSpace::Jpart(double t, double s, int doms )
{
    SetPointsFromGeo(t,s,doms,doms,0); 
    
    ComputeDistance(); 
    
    ComputeDeltaT(t, s);
    
    ComputeDeltaR(t,doms,s,doms);        
    
    return  -m_WaveNumber.real()*
                    boost::math::cyl_bessel_j
                    (1,m_WaveNumber.real()*m_distance)*
                     m_DeltaR[m_component]*m_DeltaTp/m_distance;

    
}

                    
dcomp FreeDerivativeSpace::evaluateGF( double t, double s, int domt, int doms,
               int dom_op)
{
    SetPointsFromGeo(t,s,domt,doms,dom_op); 

    ComputeDistance(); 
    
    ComputeDeltaT(t, s);
    
    ComputeDeltaR(t,doms,s,doms);     

    return IntegralKernel(); 

}
  
  dcomp FreeDerivativeSpace::evaluatePotential(double t, int domt,
             double y1, double y2,int dom_op)
  {
      SetXFromGeo(t,domt,dom_op); 
      
      setYpoint(y1,y2); 
      
      ComputeDistance(); 
      
      return IntegralKernel(); 
      
      
  } 
  
 