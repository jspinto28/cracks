#ifndef OPERATOR_H
#define OPERATOR_H
#include <iostream>
#include <omp.h>
#include <complex>
#include "Arc.h"
#include "ChebTransformer.h"
#include "GreenFunction.h"
#include "constant.h" 

const int extra = 28; 

template< typename T> 
class OperatorBase
{
	protected: 
		
		GreenFunctionBase<T>* m_GreenFun;
		
		ChebTransformer<T>* m_Transformer; 
		
		int m_Ndof; 
		
		int m_Nc; 
		
		int m_Nthreads; 
		
		int m_Lmax; 
		
		T m_tol; 
		
		T* m_xc; 
		
		T m_wc;
		
		int m_Nqc;
					
		T Tconstat(int n)
		{
			if(n ==0 )
			{
				return pi;
			}
			
			return 0.5*pi; 
		}			
		
		T Tn(int n,T t)
		{
			return cos(T(n)*acos(t)); 
		}
		
		std::complex<T> DirectQuadrature(int mm,
			std::complex<T> *Gevals)
		{
						
			std::complex<T> result(0.0,0.0); 
			
			for (int n1(0); n1<m_Nqc; ++n1)
			{
				result += 
					Gevals[n1]*
					Tn(mm,m_xc[n1]); 					

			}
						
			
			return result*m_wc*m_wc; 
			
		}
		
		void ComputeFull
			(std::complex<T>* Matrix)
		{
			T* chebpts = m_Transformer->chebpts(); 
			
			std::complex<T>* RegPart = 
				new std::complex<T>[m_Nc*m_Nc]; 
			std::complex<T>* JPart = 
				new std::complex<T>[m_Nc*m_Nc]; 

			 #pragma omp parallel for collapse(2)  
			for(int ii=0; ii < m_Nc; ii++)
			{
				for(int jj=0; jj < m_Nc; jj++)
				{					
					//aqui podemos cambiar para dar la matrix transpuesta
					m_GreenFun->Compute(chebpts[jj],
						chebpts[ii], &JPart[jj+ii*m_Nc],
						&RegPart[jj+ii*m_Nc]); 
																
				}
			}

			m_Transformer->MatrixIntegration(RegPart,m_Ndof,Matrix); 

			if( m_GreenFun->IsSingular())
			{
				std::complex<T>* JPartTrans= 
					new std::complex<T>[4*m_Nc*m_Nc]; 								
				
				for (int ii(0); ii < 4*m_Nc*m_Nc; ++ii)
				{
					JPartTrans[ii] =std::complex<T>(0.0,0.0); 
				}
				
				m_Transformer->
					MatrixCoef(JPart,m_Nc,JPartTrans);
					
			
				// needs to be cheked on paper
				#pragma omp parallel for collapse(2)  
				for (int mm=0; mm < m_Ndof; mm++)
				{
					for( int ll=0; ll < m_Ndof; ll++)
					{
						std::complex<T>aux(0.0,0.0); 
						
						for(int nn(0); nn <m_Nc; ++nn)
						{
							aux += m_GreenFun->SingPartCoefs(nn)*
							(JPartTrans[(ll+nn)+(mm+nn)*m_Nc]*Tconstat((mm+nn))*Tconstat(ll+nn)+
							 JPartTrans[(ll+nn)+abs(mm-nn)*m_Nc]*Tconstat(abs(mm-nn))*Tconstat((ll+nn))+
							 JPartTrans[abs(ll-nn)+(mm+nn)*m_Nc]*Tconstat(mm+nn)*Tconstat(abs(ll-nn))+
							 JPartTrans[abs(ll-nn)+abs(mm-nn)*m_Nc]*Tconstat(abs(mm-nn))*Tconstat(abs(ll-nn)));
							
							// aux += m_GreenFun->SingPartCoefs(nn)*
							// (JPartTrans[abs(ll-nn)+abs(mm-nn)*m_Nc]+
							// JPartTrans[abs(ll-nn)+(mm+nn)*m_Nc]+
							// JPartTrans[(ll+nn)+abs(mm-nn)*m_Nc]+
							// JPartTrans[(ll+nn)+(mm+nn)*m_Nc]); 
						}
						
						aux *=0.25; 
						
						 // aux *= 0.25*Tconstat(mm)*Tconstat(ll); 
						
						Matrix[ll+mm*m_Ndof] += aux; 
					}
				}
				
				delete [] JPartTrans; 
				
			}		

			
			delete [] chebpts; 
			delete [] RegPart; 
			delete[] JPart;
		}
		
		
		void ComputeRegCompresed
			(std::complex<T>* Matrix)
		{
			int Ndofb = m_Ndof; 
			int Ncb = m_Nc; 
            						
			std::complex<T>* RegPart = 
				new std::complex<T>[m_Nqc*m_Nqc]; 
				
			#pragma omp parallel for collapse(2)  
			for(int ii=0; ii < m_Nqc; ii++)
			{
				for(int jj=0; jj < m_Nqc; jj++)
				{					
					RegPart[jj+ii*m_Nqc]=
					m_GreenFun->FundamentalSolution
					(m_xc[ii],
						m_xc[jj]); 
				}
			}
			
			std::complex<T>* RowVector = 
				new std::complex<T>[m_Nqc]; 
				
			std::complex<T>* ColumVector = 
				new std::complex<T>[m_Nqc]; 
			
			for(int ii(0); ii < m_Nqc; ++ii)
			{
				RowVector[ii] =std::complex<T>(0.0,0.0); 
				ColumVector[ii] =std::complex<T>(0.0,0.0);
			}
			
			#pragma omp parallel for collapse(1) 
			for(int ii=0; ii < m_Nqc; ii++)
			{
				for(int jj=0; jj < m_Nqc; jj++)
				{	
					RowVector[ii] += RegPart[jj+ii*m_Nqc]; 
					//not cache efficient...(shold be done bloking)
					ColumVector[ii] +=RegPart[ii+jj*m_Nqc]; 
					
				}
			}
             
			delete [] RegPart; 
			
			int lev = 0; 
			int a = 0; 
			int b = this->m_Ndof-1; 
			
			int ti,tc,td,mid,n1,n2; 
			T vi,vc,vd; 
			
			while(lev < m_Lmax)
			{
				mid = (int)((a+b)/2); 
			
				ti = mid-1; 
				
				tc = mid; 
				
				td = mid+1;
		   
				vi =  abs(DirectQuadrature(ti,RowVector)); 
				
				vc =  abs(DirectQuadrature(tc,RowVector)); 
				
				vd =  abs(DirectQuadrature(td,RowVector));
				
				if( ((vd < 0.5*m_tol)&&(vc < 0.5*m_tol))||
						((vi < 0.5*m_tol)&&(vc < 0.5*m_tol))) 
				{
					
					b = mid; 
					
				}
				else
				{            
					a = mid; 
				}
				
				lev++; 
        
			}   
        
			n1 =b; 
			
			a=0;
			
			b= this->m_Ndof-1; 
			
			lev=0;

			while(lev < m_Lmax)
			{
				mid = (int)((a+b)/2); 
			
				ti = mid-1; 
				
				tc = mid; 
				
				td = mid+1;
						
				vi =  abs(DirectQuadrature(ti,ColumVector)); 
				
				vc =  abs(DirectQuadrature(tc,ColumVector)); 
				
				vd =  abs(DirectQuadrature(td,ColumVector));
				
				if( ((vd < 0.5*m_tol)&&(vc < 0.5*m_tol))||
						((vi < 0.5*m_tol)&&(vc < 0.5*m_tol))) 
				{
					
					b = mid; 
					
				}
				else
				{            
					a = mid; 
				}
				
				lev++; 
				
			}
             
            delete [] ColumVector;
            delete [] RowVector;
			
			n2 =b; 
     
			m_Ndof = max(n1,n2); 
			m_Nc = m_Ndof +extra; 

            m_Transformer->ChangePoints(m_Nc); 
 			 
            std::complex<T>* MatrixC = new std::complex<T>[m_Ndof*m_Ndof]; 
            
            for (int ii(0); ii < m_Ndof*m_Ndof; ++ii)
            {
                MatrixC[ii] = std::complex<T>(0.0,0.0); 
            }
            
 			ComputeFull(MatrixC);
            
            for (int ii(0); ii < m_Ndof; ++ii)
            {
                for(int jj(0); jj < m_Ndof; ++jj)                
                {
                    Matrix[jj+ii*Ndofb] =
                             MatrixC[jj+ii*m_Ndof]; 
                } 
            }
			
			m_Nc = Ncb; 
			m_Ndof = Ndofb;
            
            m_Transformer->ChangePoints(m_Nc); 
			
		}
		
		void ComputeSingCompresed
			(std::complex<T>* Matrix)
		{
			ComputeFull(Matrix); 
            
			for(int mm=0; mm< m_Ndof; ++mm)
			{
				for(int ll=0; ll < m_Ndof; ++ll)
				{
					if(abs(Matrix[ll+mm*m_Ndof]) < m_tol)
					{
						Matrix[ll+mm*m_Ndof] = std::complex<T>(0.0,0.0); 
					}
				}
			}
			
		}
		
	public: 
	
		OperatorBase(int Ndof, int Nthreads):
			m_Ndof(Ndof),
			m_Nthreads(Nthreads),
			m_tol(0),
			m_Nqc(0)
		{
			m_Nc = m_Ndof+extra; 
			
			m_Transformer = new ChebTransformer<T>(m_Nthreads,m_Nc);
			
			omp_set_num_threads(Nthreads);
			
		}
		
		OperatorBase(int Ndof, int Nthreads, T tol, int Lev):
			m_Ndof(Ndof),
			m_Nthreads(Nthreads),
			m_tol(tol),
			m_Lmax(Lev)
		{
            m_Nc = m_Ndof+extra;     
                
			m_Nqc = m_Ndof+extra/2; 
			
			m_xc = new T[m_Nqc]; 
			m_wc = pi/T(m_Nqc);
			
			for(int ii(0); ii < m_Nqc; ++ii)
			{
				m_xc[ii] =cos(T(2*ii+1)*pi/T(2*m_Nqc)); 

			}
			
			m_Transformer = new ChebTransformer<T>(m_Nthreads,m_Nc);
            
            omp_set_num_threads(Nthreads);
			
		}
		
		~OperatorBase()
		{
			delete m_Transformer; 
			
			delete m_GreenFun; 
			
			if(m_tol >0) 
			{
				delete [] m_xc; 

			}
			
		}	
	
		virtual void Compute(std::complex<T>* Matrix)
		{			
          
            
			if(m_tol >0)
			{
				if(m_GreenFun->IsSingular())
				{
 					ComputeSingCompresed(Matrix); 
				}
				else
				{
 					ComputeRegCompresed(Matrix); 
				}
			}						
			else
			{
				ComputeFull(Matrix); 
			}
            

		}
	
};


template< typename T> 
class VHemholtz : public OperatorBase<T>
{
	protected: 
	
	T m_WaveNumber;
	
	ArcBase<T> *m_arcTrial;
	ArcBase<T> *m_arcTest;
	
	public: 
	
		VHemholtz(int Ndof, int Nthreads, 
			T WaveNumber, ArcBase<T>* arcTrial, 
			ArcBase<T> * arcTest):
			OperatorBase<T>(Ndof,Nthreads),
			m_WaveNumber(WaveNumber),
			m_arcTrial(arcTrial),
			m_arcTest(arcTest)
		{
			this->m_GreenFun = new GFHelmholtz<T>(arcTrial,
				arcTest,WaveNumber); 
			
		}
		
		VHemholtz(int Ndof, int Nthreads,
			T tol, int Lev, 
			T WaveNumber, ArcBase<T>* arcTrial, 
			ArcBase<T> * arcTest):
			OperatorBase<T>(Ndof,
			Nthreads,tol,Lev),
			m_WaveNumber(WaveNumber),
			m_arcTrial(arcTrial),
			m_arcTest(arcTest)
		{
			this->m_GreenFun = new GFHelmholtz<T>(arcTrial,
				arcTest,WaveNumber); 
		}
};


template< typename T> 
class VElasticity : public OperatorBase<T>
{
	protected: 
		ArcBase<T> *m_arcTrial;
		ArcBase<T> *m_arcTest;

		T m_omega;
		T m_lambda;
		T m_mu; 
		
		int m_Ndof2; 
		
	public: 

	VElasticity(int Ndof, int Nthreads, 
			T omega, T lambda, T mu
			, ArcBase<T>* arcTrial, 
			ArcBase<T>* arcTest):
			OperatorBase<T>(Ndof,Nthreads),
			m_arcTrial(arcTrial),
			m_arcTest(arcTest),
			m_omega(omega),
			m_lambda(lambda),
			m_mu(mu)
		{			
			m_Ndof2 = 2*Ndof; 
		}
		
	VElasticity(int Ndof, int Nthreads,
			T tol, int Lev, 
			T omega, T lambda, T mu,
			ArcBase<T>* arcTrial, 
			ArcBase<T> * arcTest):
			OperatorBase<T>(Ndof,
			Nthreads,tol,Lev),			
			m_arcTrial(arcTrial),
			m_arcTest(arcTest),
			m_omega(omega),
			m_lambda(lambda),
			m_mu(mu)
		{
			m_Ndof2 = 2*Ndof; 
		}
	
	void Compute(std::complex<T>* Matrix)
	{
			int Ndof = this->m_Ndof; 
		
			std::complex<T>* M11 
				= new std::complex<T>[Ndof*Ndof]; 
			std::complex<T>* M12 
				= new std::complex<T>[Ndof*Ndof]; 
			std::complex<T>* M21 
				= new std::complex<T>[Ndof*Ndof]; 
			std::complex<T>* M22 
				= new std::complex<T>[Ndof*Ndof]; 
			
			for (int ii(0); ii < Ndof*Ndof; ++ii)
			{
				M11[ii] = std::complex<T>(0.0,0.0);
				M12[ii] = std::complex<T>(0.0,0.0);
				M21[ii] = std::complex<T>(0.0,0.0);
				M22[ii] = std::complex<T>(0.0,0.0);
			}
				
			this->m_GreenFun = new GFElasticity1<T>(m_arcTrial,
				m_arcTest,
				m_lambda,m_mu,m_omega); 
				
			OperatorBase<T>::Compute(M11); 
				
			for (int ii(0); ii < Ndof*Ndof; ++ii)
			{
				M22[ii] = M11[ii]; 

			}
			
			delete this->m_GreenFun; 
			
			this->m_GreenFun = new GFElasticity2<T>
				(0,0,m_arcTrial,m_arcTest,
				m_lambda,m_mu,m_omega); 
			
			OperatorBase<T>::Compute(M11); 
			
			delete this->m_GreenFun; 
			
			this->m_GreenFun = new GFElasticity2<T>
				(1,0,m_arcTrial,m_arcTest,
				m_lambda,m_mu,m_omega); 
				
			OperatorBase<T>::Compute(M21);

			if(this->m_GreenFun->IsSingular())
			{
				for (int ii(0); ii < Ndof*Ndof; ++ii)
				{
					M12[ii] = M21[ii]; 
				}
			}
			else
			{
				delete this->m_GreenFun;
				
				this->m_GreenFun = new GFElasticity2<T>
					(0,1,m_arcTrial,m_arcTest,
					m_lambda,m_mu,m_omega); 
				
				OperatorBase<T>::Compute(M12);				
				
			}
			
			delete this->m_GreenFun;
				
			this->m_GreenFun = new GFElasticity2<T>
					(1,1,m_arcTrial,m_arcTest,
					m_lambda,m_mu,m_omega); 
				
			OperatorBase<T>::Compute(M22);	
			
			#pragma omp parallel for collapse(2) 
			for(int ii=0; ii < Ndof; ii++)
			{
				for (int jj=0; jj < Ndof; jj++)
				{
					Matrix[jj+ii*m_Ndof2]
						= M11[jj+ii*Ndof]; 
					Matrix[jj+Ndof+ii*m_Ndof2]
						= M21[jj+ii*Ndof]; 
					Matrix[jj+(ii+Ndof)*m_Ndof2]
						= M12[jj+ii*Ndof]; 
					Matrix[jj+Ndof+(ii+Ndof)*m_Ndof2]
						= M22[jj+ii*Ndof]; 
				}
			}									
		
			delete [] M11; 
			delete [] M12; 
			delete [] M21; 
			delete [] M22;
			
		}

};

#endif