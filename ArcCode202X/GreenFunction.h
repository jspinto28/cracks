#ifndef GFUN_H
#define GFUN_H

#include <iostream>
#include <complex>
#include "constant.h" //i,pi,eulergamma, vector con la expansion de coefs del log
#include "Arc.h"
#include <boost/math/special_functions/hankel.hpp>
using namespace std;

template <typename T>
class GreenFunctionBase{
	
	
	protected: 
	
		std::complex<T> i = std::complex<T>(0.0,T(1.0));
		ArcBase<T> *m_arcTrial;
		ArcBase<T> *m_arcTest;
		T m_eps; 
					
	
	public: 
	
		GreenFunctionBase(ArcBase<T> *arcTrial,
				ArcBase<T> *arcTest): 
			m_arcTrial(arcTrial),
			m_arcTest(arcTest)
		{
			if(sizeof(T) == 4)
			{
				m_eps = 1e-9; 
			}
			else
			{
				m_eps = 1e-14; 
			}
			
		}
		
		virtual int GetOpDimension()=0; 
		
		virtual complex<T> FundamentalSolution(T t,T s, T dist) = 0; 
		
		virtual complex<T> Jfun(T t, T s, T dist)=0; 
		
		virtual complex<T> JfunLimit(T t)=0; 
		
		virtual complex<T> RegPartLimit(T t)=0; 
		
		virtual T SingPart(T t, T s)=0; 
		
		virtual T SingPartCoefs(int n)=0; 
        
        T ComputeDist(T t, T s)
		{
			T xtrial, ytrial, xtest, ytest; 
			
			m_arcTrial->GetPoint(t,&xtrial,&ytrial); 
			m_arcTest->GetPoint(s,&xtest,&ytest); 
	
			return sqrt((xtrial-xtest)*(xtrial-xtest)+
				(ytrial-ytest)*(ytrial-ytest)); 
		}
								
		void Compute( T t, T s, complex<T>* JPart,
			complex<T>* RegPart)
			{
				
				T dist  = ComputeDist(t,s); 
				
				if(IsSingular())
				{
					if(dist > m_eps)
					{

						*JPart = Jfun(t,s,dist); 
												
						*RegPart =FundamentalSolution(t,s,dist)-
							SingPart(t,s)*(*JPart); 			
						
					}
					else
					{		

						*RegPart = RegPartLimit(t); 
						
						*JPart = JfunLimit(t); 
						
					}
				}
				else
				{
			
					*RegPart = FundamentalSolution(t,s,dist); 
					
					*JPart = complex<T>(0.0,0.0); 
				}
							
			}
		
		complex<T> FundamentalSolution(T t, T s)
		{
			T dist  = ComputeDist(t,s); 
			
			return FundamentalSolution(t,s,dist); 
		}
		
		bool IsSingular()
		{	
			return (m_arcTest->GetIndex() == 
			m_arcTrial->GetIndex());
		}
		
};

template <typename T>
class GFHelmholtz : public GreenFunctionBase<T>
{
	protected: 
	
		T m_WaveNumber; 
		
	public: 
	
		GFHelmholtz(ArcBase<T> *arcTrial,
			ArcBase<T> *arcTest,
			T WaveNumber): 
			GreenFunctionBase<T>(arcTrial,arcTest),
			m_WaveNumber(WaveNumber)
			{}
		
	
		int GetOpDimension()
		{
			return 1; 
		}
		
		complex<T> FundamentalSolution(T t,T s, T dist)
		{
			std::complex<double> H1=
				boost::math::cyl_hankel_1(0,m_WaveNumber*
                  dist);				  		
			
			return T(0.25)*std::complex<T>(-H1.imag(),
				H1.real()); 

		}
		
		complex<T> Jfun(T t, T s, T dist)
		{
			return std::complex<T> 
				( T(boost::math::cyl_bessel_j
                    (0,m_WaveNumber*dist)),0.0);
		}
		
		complex<T> JfunLimit(T t)
		{
			return std::complex<T>(T(1.0),0.0); 
		}			
		
		complex<T> RegPartLimit(T t)
		{
			T jac = this->m_arcTrial->Jacobian(t);
			
			return std::complex<T>(
			-0.5/pi*(log(0.5*jac*m_WaveNumber)+
				eulergamma),0.25); 
		}
		
		T SingPart(T t, T s)
		{
			return -0.5/pi*log(abs(t-s));
		}
		
		T SingPartCoefs(int n)
		{
			return LogCheb[n]; 
		}	
	
};

template <typename T>
class GFElasticityBase : public GreenFunctionBase<T>
{
	protected: 
	
		T m_lambda;
		T m_mu; 
		T m_omega; 
		T m_ks; 
		T m_kp;
		T m_ks2; 
		T m_kp2; 
		T m_omega2; 
		
	public: 
		
		GFElasticityBase(ArcBase<T> *arcTrial,
			ArcBase<T> *arcTest,
			T lambda, T mu, T omega): 
			GreenFunctionBase<T>(arcTrial,arcTest),
			m_lambda(lambda),
			m_mu(mu),
			m_omega(omega)
			{
				m_kp2 = m_omega* m_omega/(m_lambda+2.0*m_mu); 
				m_ks2 = m_omega*m_omega/m_mu; 
				m_ks = sqrt(m_ks2); 
				m_kp = sqrt(m_kp2);
				m_omega2 = m_omega*m_omega; 
				
			}
			
	
		int GetOpDimension()
		{
			return 2; 
		}	
		
		T SingPart(T t, T s)
		{
			return -0.5/pi*log(abs(t-s));
		}
		
		T SingPartCoefs(int n)
		{
			return LogCheb[n]; 
		}	
	
}; 

template <typename T>
class GFElasticity1 : public GFElasticityBase<T> 
{
	public: 
	
	GFElasticity1(ArcBase<T> *arcTrial,
		ArcBase<T> *arcTest,
		T lambda, T mu, T omega): 
		GFElasticityBase<T>(arcTrial,arcTest,
			lambda,mu,omega) 
		{				
		}
		
		complex<T> FundamentalSolution(T t,T s, T dist)
		{
			
			std::complex<double> h0ksd =
				boost::math::cyl_hankel_1(0,this->m_ks*dist); 
			std::complex<double> h1ksd =
				boost::math::cyl_hankel_1(1,this->m_ks*dist); 
			std::complex<double> h1kpd =
				boost::math::cyl_hankel_1(1,this->m_kp*dist); 
				
			std::complex<T> h0ks,h1ks,h1kp;
			
			if(sizeof(T)==4)
			{
				h1ks = std::complex<T>(T(h1ksd.real()),
					T(h1ksd.imag()));
				h0ks = std::complex<T>(T(h0ksd.real()),
					T(h0ksd.imag()));
				h1kp = std::complex<T>(T(h1kpd.real()),
					T(h1kpd.imag()));
			}
			else
			{
				h1ks = h1ksd; 
				h0ks = h0ksd; 
				h1kp = h1kpd;
			}
			
			return T(0.25)*this->i/(this->m_mu)*h0ks
				-T(0.25)*this->i/(this->m_omega2*dist)*
				(this->m_ks*h1ks-this->m_kp*h1kp);
		}
		
		complex<T> Jfun(T t, T s, T dist)
		{
            
            
			return
			complex<T>
			( boost::math::cyl_bessel_j
                    (0,this->m_ks*dist)/(this->m_mu)-
					((this->m_ks)*boost::math::cyl_bessel_j
                    (1,this->m_ks*dist) -(this->m_kp)*
					boost::math::cyl_bessel_j
                    (1,this->m_kp*dist))/(this->m_omega2*dist),0);

		}
		
		complex<T> JfunLimit(T t)
		{
            
			return complex<T>(1.0/(this->m_mu)
				+0.5*(this->m_kp2-this->m_ks2)
				/(this->m_omega2),0.0); 
		}			
		
		complex<T> RegPartLimit(T t)
		{
            T jac = this->m_arcTrial->Jacobian(t);

			return 
			complex<T>(
			-T(0.5/pi)*(
			1.0/(this->m_mu)*(log(0.5*this->m_ks)+eulergamma)
			+0.5/(this->m_omega2)*(
			this->m_kp2*log(0.5*this->m_kp)-
			this->m_ks2*log(0.5*this->m_ks)
			)
			+0.25/(this->m_omega2)*(this->m_ks2-this->m_kp2)*
			(1.0-2.0*eulergamma)
			),
			T(0.25)*(
			1.0/this->m_mu-0.5/this->m_omega2*
				(this->m_ks2-this->m_kp2)
			))
			-T(0.5/pi)*log(jac)*JfunLimit(t); 
		}

};

template <typename T>
class GFElasticity2 : public GFElasticityBase<T> 
{
	
	protected: 
	
	int m_i; 
	int m_j; 
	
	T Ematrix(T t, T s, T dist)
	{
		T x[2]={0.0,0.0}; 
		T y[2] ={0.0,0.0};
		
		this->m_arcTrial->GetPoint(t,&x[0],&x[1]); 
		this->m_arcTest->GetPoint(s,&y[0],&y[1]); 
		
		return (x[m_i]-y[m_i])*(x[m_j]-y[m_j])
				/(dist*dist); 		
				
	}
	
	T EmatrixLimit(T t)
	{
		T x[2] = {0.0,0.0}; 
		
		this->m_arcTrial->GetDevPoint(t,&x[0],&x[1]); 
		
		T jac2 = (x[0]*x[0]+x[1]*x[1]); 
		
		return x[m_i]*x[m_j]/(jac2); 		
		
	}
	
	public: 
	
	GFElasticity2(int indi, int indj,
		ArcBase<T> *arcTrial,
		ArcBase<T> *arcTest,
		T lambda, T mu, T omega): 
		GFElasticityBase<T>(arcTrial,arcTest,
			lambda,mu,omega), 
	    m_i(indi), 
		m_j(indj)
		{}
		
		complex<T> FundamentalSolution(T t,T s, T dist)
		{
			T Eterm = Ematrix(t,s,dist); 	
            
            T dos = T(2.0); 
			
			std::complex<double> h1ksd;
			std::complex<double> h0ksd;
			std::complex<double> h1kpd;
			std::complex<double> h0kpd;
			
			h1ksd =boost::math::cyl_hankel_1(1,this->m_ks*dist);
			h0ksd =boost::math::cyl_hankel_1(0,this->m_ks*dist);
			h1kpd =boost::math::cyl_hankel_1(1,this->m_kp*dist);
			h0kpd =boost::math::cyl_hankel_1(0,this->m_kp*dist);
			
			std::complex<T> h1ks;
			std::complex<T> h0ks;
			std::complex<T> h1kp;
			std::complex<T> h0kp;
			
			if(sizeof(T) == 4)
			{
				h1ks = std::complex<T>(T(h1ksd.real()),
					T(h1ksd.imag()));
				h0ks = std::complex<T>(T(h0ksd.real()),
					T(h0ksd.imag()));
				h1kp = std::complex<T>(T(h1kpd.real()),
					T(h1kpd.imag()));
				h0kp = std::complex<T>(T(h0kpd.real()),
					T(h0kpd.imag()));
			}
			else
			{
				h1ks =h1ksd;
				h0ks = h0ksd;
				h1kp = h1kpd;
			    h0kp = h0kpd;
			}
			
			std::complex<T> z = 
				T(0.25)/(this->m_omega2)*
				(dos*this->m_ks/dist*h1ks-
				this->m_ks2*h0ks-
				dos*this->m_kp/dist*h1kp+
				this->m_kp2*h0kp)*Eterm;
				
// 			std::complex<T> z = 
// 				T(0.25)/(this->m_omega2)*
// 				(dos*this->m_ks/dist*h1ks-
// 				this->m_ks2*h0ks-
// 				dos*this->m_kp/dist*h1kp+
// 				this->m_kp2*h0kp);
			
			return std::complex<T>(-z.imag(),z.real()); 
			
		}
		
		complex<T> Jfun(T t, T s, T dist)
		{
			T Eterm = Ematrix(t,s,dist); 

			complex<T> z = T(1.0)/(this->m_omega2)*				
			complex<T>
			( 
			   this->m_kp2*boost::math::cyl_bessel_j
				(0,this->m_kp*dist)
			  -this->m_ks2*boost::math::cyl_bessel_j
				(0,this->m_ks*dist)+
			  T(2.0)*
			  (
				this->m_ks/dist*
				 boost::math::cyl_bessel_j
				 (1,this->m_ks*dist)
				- 
				this->m_kp/dist*
				 boost::math::cyl_bessel_j
				 (1,this->m_kp*dist)
			  )
			,0.0);

			  return z*Eterm; 
			 
// 			 return z; 

		}
		
		complex<T> JfunLimit(T t)
		{
			return (0.0,0.0); 
		}			
		
		complex<T> RegPartLimit(T t)
		{
            
// 			return 
// 			complex<T>(T(0.25/pi)/(this->m_omega2)*
// 			(this->m_ks2-this->m_kp2),0.0);
//             
			
			T Elim = EmatrixLimit(t);
			
			return 
			T(0.25)/T(pi)/(this->m_omega2)*
			(this->m_ks2-this->m_kp2)*
			Elim; 

		}

};































#endif
