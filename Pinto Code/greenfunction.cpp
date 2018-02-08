#include <boost/math/special_functions/hankel.hpp>
#include "greenfunction.h"
#include <vector>
#include <math.h>  
#include "mex.h"


const double pii= 4.0*atan(1.0);

dcomp FreeSpace::IntegralKernel()
{
    if(m_WaveNumber.real() > 0.0000001)
    {
    
        return 0.25*i*boost::math::cyl_hankel_1(0,m_WaveNumber.real()*m_distance);      
    
    }
    else
    {
        return dcomp(-1/(2*pii)*log(m_distance),0);
    }
}
                    
dcomp GreenFunctionBase::evaluateRegularizator(double t, double s) 
{

    return dcomp(-0.5/pii*log(abs(t-s)),0);
    
}

                    
dcomp FreeSpace::evaluateGF( double t, double s, int domt, int doms,
               int dom_op)
{
    SetPointsFromGeo(t,s,domt,doms,dom_op); 

    ComputeDistance(); 

    return IntegralKernel(); 

}


 dcomp FreeSpace::IntegralKernelDev()
 {  
     return -m_WaveNumber.real()*0.25*i*boost::math::cyl_hankel_1
             (1,m_WaveNumber.real()*m_distance);
     
 }

 
 void FreeSpace::ComputeGradDistance()
 {
    m_gradDistance[0] = (m_x[0]-m_y[0])/m_distance; 
    
    m_gradDistance[1] = (m_x[1]-m_y[1])/m_distance; 
     
 }

 dcomp FreeSpace::evaluateGFDetNx(double t, double s, int domt, int doms,
               int dom_op)
 {
    SetPointsFromGeo(t,s,domt,doms,dom_op); 
    
    ComputeDistance(); 
    
    ComputeGradDistance(); 
    
    SetXNormalFromGeo(t,domt,dom_op); 
    
    double NormaldotGrad = m_xNormal[0]*m_gradDistance[0]+
            m_xNormal[1]*m_gradDistance[1];
            
    return IntegralKernelDev()*NormaldotGrad;
     
 }
 
  dcomp FreeSpace::evaluateGFDetNy(double  t, double s, int domt, int doms,
               int dom_op)
 {
    SetPointsFromGeo(t,s,domt,doms,dom_op); 
    
    ComputeDistance(); 
    
    ComputeGradDistance(); 
    
    SetYNormalFromGeo(s,doms,dom_op); 
    
    double NormaldotGrad = m_yNormal[0]*(-m_gradDistance[0])+
            m_yNormal[1]*(-m_gradDistance[1]);
            
    return IntegralKernelDev()*NormaldotGrad;
     
 }
  
  dcomp FreeSpace::evaluatePotential(double t, int domt,
             double y1, double y2,int dom_op)
  {
      SetXFromGeo(t,domt,dom_op); 
      
      setYpoint(y1,y2); 
      
      ComputeDistance(); 
      
      return IntegralKernel(); 
      
      
  }    