#ifndef RHS_H
#define RHS_H

#include <iostream>
#include <complex>
#include "Arc.h"
#include "ChebTransformer.h"
#include "constant.h" 

const int extraR = 10; 

template< typename T> 
class RhsBase
{
	protected: 
	
		ChebTransformer<T>* m_Transformer; 	
				
		int m_Ndof; 
		
		int m_Nc; 	
		
		virtual std::complex<T> EvalfInside(T t) =0; 
		
	public: 
	
		RhsBase(int Ndof):
			m_Ndof(Ndof)
		{
			m_Nc = m_Ndof+extraR; 
			
			m_Transformer = new ChebTransformer<T>(1,m_Nc); 
		}
		
		~RhsBase()
		{
			delete m_Transformer; 
		}
		
		virtual void Compute(std::complex<T>*  Vector)
		{
			std::complex<T>*  Eval = new std::complex<T>[m_Nc];
			
			T * chebpts = m_Transformer->chebpts(); 
			
			for (int ii(0); ii < m_Nc; ++ii) 
			{
				Eval[ii] = EvalfInside(chebpts[ii]);
				
			}
			
			m_Transformer->VectorCoef(Eval,m_Ndof,Vector); 
			
			delete [] chebpts;
			delete [] Eval;
		}	
};

template< typename T> 
class PlaneWave : public RhsBase<T>
{
	protected:
	
		T m_WaveNumber; 
		
		ArcBase<T> *m_arc; 
		
		T m_alpha; 
		
		std::complex<T> EvalfInside(T t)
		{
			T x[2]={0.0,0.0};
			
			this->m_arc->GetPoint(t,&x[0],&x[1]); 
			
			T prod = x[0]*cos(m_alpha)+x[1]*sin(m_alpha); 
			
			prod *= m_WaveNumber;

			return std::complex<T>(cos(prod),
				sin(prod)); 
		}
		
	public: 
	
		PlaneWave( int Ndof, T WaveNumber, 
			T alpha, ArcBase<T> *arc):
			RhsBase<T>(Ndof),
			m_alpha(alpha),
			m_WaveNumber(WaveNumber),
			m_arc(arc)
			{}			
					
};


template< typename T> 
class PlaneWaveElasticity : public PlaneWave<T>
{
	protected: 
	
		T m_dx; 
		T m_dy; 
		T m_omega, m_lambda, m_mu;
		T m_ks, m_kp; 
	
	public: 
	
		PlaneWaveElasticity( int Ndof, T omega,
				T lambda, T mu, 
				T alpha, ArcBase<T> *arc):
			PlaneWave<T>(Ndof,0,
				alpha,arc),
			m_omega(omega),
			m_lambda(lambda),
			m_mu(mu)
			{
				m_ks = sqrt(m_omega*m_omega/m_mu);
				m_kp = sqrt(m_omega* m_omega/(m_lambda+2.0*m_mu));
			}
		
		void Compute(std::complex<T>*  Vector)
		{			
			RhsBase<T>::Compute(Vector); 
			
			for(int ii(0); ii < this->m_Ndof; ++ii)
			{
				Vector[ii] *= m_dx; 

			}
			
			RhsBase<T>::Compute(Vector+this->m_Ndof); 
			
			for (int ii(0); ii < this->m_Ndof; ++ii)
			{
				Vector[ii+this->m_Ndof] *= m_dy; 
			}						
			
		}
	
};

template< typename T> 
class ElasticityP : public PlaneWaveElasticity<T>
{

	public: 
	
		ElasticityP( int Ndof, T omega,
					T lambda, T mu, 
					T alpha, ArcBase<T> *arc):
			PlaneWaveElasticity<T>(Ndof, omega,
				lambda,mu,alpha,arc)
			{
				this->m_dx = cos(this->m_alpha);
				this->m_dy = sin(this->m_alpha); 
				
				this->m_WaveNumber = this->m_kp; 
				
			}
	
	
}; 

template< typename T> 
class ElasticityS : public PlaneWaveElasticity<T>
{

	public: 
	
		ElasticityS( int Ndof, T omega,
					T lambda, T mu, 
					T alpha, ArcBase<T> *arc):
			PlaneWaveElasticity<T>(Ndof, omega,
				lambda,mu,alpha,arc)
			{
				this->m_dx = -sin(this->m_alpha);
				this->m_dy = cos(this->m_alpha); 
				
				this->m_WaveNumber = this->m_ks; 
				
			}
	
	
}; 


#endif