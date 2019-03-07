#ifndef MY_HEADER_H
#define MY_HEADER_H
#include "IntegrationRules.h"
#include "mex.h"

#include "geo.h"
#include "prods.h"
#include <vector>

using namespace std;

typedef complex<double> dcomp;

const dcomp i(0.0,1.0);

class GreenFunctionBase {
    
protected: 
    
    dcomp m_WaveNumber; 
    
    vector<double> m_x; 
    
    vector<double> m_y; 
    
    vector<double> m_xNormal;
    
    vector<double> m_yNormal; 
    
    double m_distance;
    
    int m_RegOrder;
    
    double* m_bgeo; 
    
    int m_Nb; 
    
    double m_delta;
       
public:

    GreenFunctionBase(dcomp k, double* bgeo, int Nb,double delta):
    m_WaveNumber(k),
    m_RegOrder(0),
    m_bgeo(bgeo),
    m_Nb(Nb), 
    m_delta(delta)        
    {
        m_x.resize(2,0.0);
        m_y.resize(2,0.0);
        m_xNormal.resize(2,0.0);
        m_yNormal.resize(2,0.0);
    }               
        
    GreenFunctionBase(double k, double* bgeo, int Nb,double delta):
    m_RegOrder(0),             
    m_bgeo(bgeo),
    m_Nb(Nb),
    m_delta(delta)        
    {
        m_WaveNumber.real(k); 
        m_WaveNumber.imag(0); 
        m_x.resize(2,0.0);
        m_y.resize(2,0.0);
        m_xNormal.resize(2,0.0);
        m_yNormal.resize(2,0.0);
    }    
    
    void setTwoDVector(double v1, double v2, vector<double>& vec)
    {     
        vec[0] = v1; 
        
        vec[1] = v2; 
        
    }
    
    void setXpoint  (double x1,double x2)
    {
        setTwoDVector(x1,x2,m_x); 
    }
    
    void setYpoint  (double y1,double y2)
    {
        setTwoDVector(y1,y2,m_y); 
    }
    
    void setXNormal  (double x1,double x2)
    {
        setTwoDVector(x1,x2,m_xNormal); 
    }
    
    void setYNormal  (double x1,double x2)
    {
        setTwoDVector(x1,x2,m_yNormal); 
    }    
        
    void setXpoint  (vector<double> data )
    {
        m_x = data;  
    }
    
    void setYpoint  (vector<double> data)
    {
        m_y = data; 
    }
    
    void setXNormal  (vector<double> data)
    {
        m_xNormal = data; 
    }
    
    void setYNormal  (vector<double> data)
    {
       m_yNormal = data;  
    }
        
    virtual void SetXFromGeo(double t, int domt, int dom_op) 
    {   
        double x1,x2; 

        geom(&x1, &x2, t,m_bgeo,m_Nb,m_delta); 
        
        setXpoint(x1,x2); 
        
    }
    
    virtual void SetYFromGeo(double s, int doms, int dom_op) 
    {   
        double y1,y2; 
        
        geom(&y1, &y2, s,m_bgeo,m_Nb,m_delta); 
        
        setYpoint(y1,y2); 
        
    }
    
//     virtual void SetXNormalFromGeo(double t, int domt, int dom_op) 
//     {
//         double nx1,nx2;
//         
//         normal(&nx1,&nx2,t,domt,dom_op);
//         
//         setXNormal(nx1,nx2);
// 
//     }
//     
//     virtual void SetYNormalFromGeo(double s, int doms, int dom_op) 
//     {
//         double ny1,ny2; 
//         
//         normal(&ny1,&ny2,s,doms,dom_op);
//         
//         setYNormal(ny1,ny2);
// 
//     }
    
    void ComputeDistance()
    {
        m_distance = sqrt(pow(m_x[0]-m_y[0],2)+
                pow(m_x[1]-m_y[1],2));
        
        if(m_distance < 1e-12)
        {
            m_distance = 1e-12;        
        
        } 
           
    }
    
    void SetPoints(double x1,double x2,double y1,double y2)
    {
        setXpoint(x1,x2); 
        
        setYpoint(y1,y2);

    }    
    
    void SetPoints(vector<double> datax,vector<double> datay)
    {
        setXpoint(datax); 
        
        setYpoint(datay);
        
    }    
    
    void SetPointsFromGeo(double t, double s, int domt, int doms,
            int dom_op)
    {
        SetXFromGeo(t, domt, dom_op); 
        
        SetYFromGeo(s, doms, dom_op); 
        
    }    
    
    dcomp GetWaveNumber()
    {
        return m_WaveNumber; 
    }
    
//     virtual dcomp evaluate() = 0; 
//     
//     virtual dcomp evaluate(vector<double> dataX, vector<double> dataY)
//     = 0;
    
    
     virtual dcomp IntegralKernel()=0;   
     
     virtual dcomp evaluateRegularizator2(double t, double s, int doms)=0;   
     
     virtual dcomp LimitRegularized(double t, double s, int doms )=0; 
     
     virtual dcomp Jpart(double t, double s, int doms)=0; 

     virtual dcomp evaluateGF( double t, double s, int domt, int doms,
           int dom_op)=0;
     
     virtual dcomp evaluatePotential(double t, int domt,
             double y1, double y2, int dom_op)=0;
     
     void SetRegOrder(int p)
     {
         m_RegOrder = p; 
     }
     
   
//     virtual dcomp evaluateNderX() = 0; 
//     
//     virtual dcomp evaluateNderY() = 0; 
//     
//     virtual dcomp evaluateNderX( double t, double s, int domt, int doms,
//             int dom_op) = 0; 
//     
//     virtual dcomp evaluateNderY( double t, double s, int domt, int doms,
//             int dom_op) = 0; 
        
}; 

class FreeSpace : public GreenFunctionBase
{
              
    public:
        
     FreeSpace(double k,double* bgeo, int Nb, double delta):
         GreenFunctionBase(k,bgeo,Nb,delta) {} 
         
      dcomp IntegralKernel();   
     
      dcomp evaluateRegularizator2(double t, double s, int doms);   
     
      dcomp LimitRegularized(double t, double s, int doms ); 
     
      dcomp Jpart(double t, double s, int doms); 

      dcomp evaluateGF( double t, double s, int domt, int doms,
           int dom_op);
     
      dcomp evaluatePotential(double t, int domt,
             double y1, double y2, int dom_op);

};

class FreeDerivativeSpace : public GreenFunctionBase
{
    protected: 
        
     int m_p;   
             
     double m_DeltaTp; 
     
      double ChebT(double t);

      void ComputeDeltaT(double t, double s); 
    
     vector<double> m_DeltaR; 
             
     void ComputeDeltaR(double t,int domt,double s, int doms)
     {
        SetXFromGeo( t, domt, 0); 
        
        SetYFromGeo( s, doms, 0); 
        
        m_DeltaR[0] = m_x[0] - m_y[0];
        
        m_DeltaR[1] = m_x[1] - m_y[1];
         
     }
     
     int m_component; //0 o 1!!!
              
    public:
        
     FreeDerivativeSpace(double k,double* bgeo, int Nb, double delta
             , int p, int component):
         GreenFunctionBase(k,bgeo,Nb,delta),
         m_component(component),
         m_p(p)
         {
             m_DeltaR.resize(2,0.0);
         }
         
    dcomp evaluatePotential(double t, int domt,
             double y1, double y2, int dom_op);
    
      dcomp IntegralKernel();   
     
      dcomp evaluateRegularizator2(double t, double s, int doms);   
     
      dcomp LimitRegularized(double t, double s, int doms ); 
      
      dcomp Jpart(double t, double s, int doms); 

      dcomp evaluateGF( double t, double s, int domt, int doms,
           int dom_op);
            
   
    
};


#endif