#ifndef MY_HEADER_H
#define MY_HEADER_H
#include "IntegrationRules.h" 

#include "geo.h"
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
    
public:

    GreenFunctionBase(dcomp k):
    m_WaveNumber(k)
    {
        m_x.resize(2,0.0);
        m_y.resize(2,0.0);
        m_xNormal.resize(2,0.0);
        m_yNormal.resize(2,0.0);
    }               
        
    GreenFunctionBase(double k)  
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
        
        geom(&x1, &x2, t,domt,dom_op); 
        
        setXpoint(x1,x2); 
        
    }
    
    virtual void SetYFromGeo(double s, int doms, int dom_op) 
    {   
        double y1,y2; 
        
        geom(&y1, &y2, s,doms,dom_op); 
        
        setYpoint(y1,y2); 
        
    }
    
    virtual void SetXNormalFromGeo(double t, int domt, int dom_op) 
    {
        double nx1,nx2;
        
        normal(&nx1,&nx2,t,domt,dom_op);
        
        setXNormal(nx1,nx2);

    }
    
    virtual void SetYNormalFromGeo(double s, int doms, int dom_op) 
    {
        double ny1,ny2; 
        
        normal(&ny1,&ny2,s,doms,dom_op);
        
        setYNormal(ny1,ny2);

    }
    
    void ComputeDistance()
    {
        m_distance = sqrt(pow(m_x[0]-m_y[0],2)+
                pow(m_x[1]-m_y[1],2));
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
    
    
     virtual dcomp evaluateRegularizator(double t, double s);    

     virtual dcomp evaluateGF( double t, double s, int domt, int doms,
           int dom_op)=0;

     virtual dcomp evaluateGFDetNx(double t, double s, int domt, int doms,
           int dom_op)=0;

     virtual dcomp evaluateGFDetNy(double t, double s, int domt, int doms,
           int dom_op)=0;   
     
     virtual dcomp evaluateGFDevy(double t, double s, int domt, int doms,
           int dom_op) =0; 
     
     virtual double GetAlpha()=0; 
     
     virtual double GetPeriod()=0; 
     
     virtual dcomp evaluatePotential(double t, int domt,
             double y1, double y2, int dom_op)=0;
    
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
    protected: 
            
        vector<double> m_gradDistance; 
          
        dcomp IntegralKernel();    
            
        dcomp IntegralKernelDev();     
        
        void ComputeGradDistance(); 
 
    public:
        
        FreeSpace(double k)
            :GreenFunctionBase(k)
            {
                m_gradDistance.resize(2,0.0); 
            }     
                
//         double XnormalDevGrad(); 
//         
//         double XnormalDevGrad(); 
        
//         virtual dcomp evaluate(); 
//     
//         virtual dcomp evaluate(vector<double> dataX, vector<double>
//                 dataY);
            
     virtual dcomp evaluateGF( double t, double s, int domt, int doms,
           int dom_op);

     virtual dcomp evaluateGFDetNx(double t, double s, int domt, int doms,
           int dom_op);

     virtual dcomp evaluateGFDetNy(double t, double s, int domt, int doms,
           int dom_op);     
     
      dcomp evaluateGFDevy(double t, double s, int domt, int doms,
           int dom_op){return 0;} 
     
     double GetAlpha()
     {
         return 0; 
     }
     
     double GetPeriod()
     {
         return 0; 
     }
     
     virtual dcomp evaluatePotential(double t, int domt,
             double y1, double y2, int dom_op);
         

//         virtual dcomp evaluateNderX(); 
// 
//         virtual dcomp evaluateNderY(); 

//         virtual dcomp evaluateNderX( double t, double s, int domt,
//                 int doms, int dom_op); 
// 
//         virtual dcomp evaluateNderY( double t, double s, int domt,
//                 int doms, int dom_op); 
    



};

class FreeSpacePer: public FreeSpace
{        
        protected: 
            
            int m_S; 
            
            double* m_Y1; 
            
            double* m_Y2;
            
            double m_delta;
            
        public:                
            
            FreeSpacePer(double k, int S,double* Y1, double* Y2,
                    double delta):
                        FreeSpace(k),
                        m_S(S),
                        m_Y1(Y1),
                        m_Y2(Y2),
                        m_delta(delta)
                        {
                        }   
            
           virtual void SetXFromGeo(double t, int domt, int dom_op) 
            {   
                double x1,x2; 

                geomPer(&x1, &x2, t,domt,m_Y1,m_Y2,m_delta,m_S); 

                setXpoint(x1,x2); 


            }

            virtual void SetYFromGeo(double s, int doms, int dom_op) 
            {   
                double y1,y2; 

                geomPer(&y1, &y2, s,doms,m_Y1,m_Y2,m_delta,m_S); 

                setYpoint(y1,y2); 

            }

//             virtual void SetXNormalFromGeo(double t, int domt, int dom_op) 
//             {
//                 double nx1,nx2;
// 
//                 normalPer(&nx1,&nx2,t,dom_op,m_an,m_bn,m_NPer);
// 
//                 setXNormal(nx1,nx2);
// 
//             }
// 
//             virtual void SetYNormalFromGeo(double s, int doms, int dom_op) 
//             {
//                 double ny1,ny2; 
// 
//                 normalPer(&ny1,&ny2,s,dom_op,m_an,m_bn,m_NPer);
// 
//                 setYNormal(ny1,ny2);
// 
//             }
                        
};           

#endif