#ifndef CHEBTRANSF_H
#define CHEBTRANSF_H
#include <iostream>
#include <complex>
#include <fftw3.h>
#include <omp.h>
template< typename T> 
class ChebTransformer
{
	protected: 
	
		fftw_plan m_plan; 
		
		int m_N;

		bool m_init; 
		
		int m_Nthreads; 
		
		bool m_getIntegration;
		
		void InitializePlan(std::complex<T>* fxc)
		{
			std::complex<double>* vec=
			new std::complex<double>[2*m_N-2];    
			
			vec[0]=std::complex<double>(
				double(fxc[m_N-1].real()),double(fxc[m_N-1].imag())); 
			vec[m_N-1]=std::complex<double>(
				double(fxc[0].real()),double(fxc[0].imag())); 
			
			for(int jj(1); jj<m_N-1;++jj)
			{
				vec[jj]=std::complex<double>(
				double(fxc[m_N-1-jj].real()),double(fxc[m_N-1-jj].imag())); 
				vec[jj+m_N-1]=std::complex<double>(
				double(fxc[jj].real()),double(fxc[jj].imag())); 
				
				// vec[jj]=fxc[m_N-1-jj];
				// vec[jj+m_N-1]=fxc[jj];
			}	
			
			fftw_complex* vecf =reinterpret_cast<fftw_complex*>(vec);
			
			m_plan=fftw_plan_dft_1d(2*m_N-2, vecf, vecf,
				FFTW_FORWARD, FFTW_ESTIMATE);   
			
			m_init = true; 
			
			delete [] vec; 
		}
		
		void Tchv_cof(std::complex<T>* fxc)
		{ 
			int m; 
			m= m_N-1;    
			
			std::complex<double>* vec=
			new std::complex<double>[2*m_N-2];    
			
			vec[0]=std::complex<double>(
				double(fxc[m_N-1].real()),double(fxc[m_N-1].imag())); 
			vec[m_N-1]=std::complex<double>(
				double(fxc[0].real()),double(fxc[0].imag())); 
			
			for(int jj(1); jj<m_N-1;++jj)
			{
				vec[jj]=std::complex<double>(
				double(fxc[m_N-1-jj].real()),double(fxc[m_N-1-jj].imag())); 
				vec[jj+m_N-1]=std::complex<double>(
				double(fxc[jj].real()),double(fxc[jj].imag())); 
			}	
		
			fftw_execute_dft(m_plan, reinterpret_cast<fftw_complex*>(vec),
				reinterpret_cast<fftw_complex*>(vec)); 

			for (int jj(0); jj<m_N;++jj)
			{				   
				fxc[jj]=std::complex<T>(
					T(vec[jj].real())/T(m),
					T(vec[jj].imag())/T(m))
					;
			}    
			
			fxc[0]=T(0.5)*fxc[0];
			fxc[m_N-1]=T(0.5)*fxc[m_N-1];
			
			if(m_getIntegration)
			{
				fxc[0] *=  T(pi); 
				for (int jj(1); jj<m_N;++jj)
				{
					fxc[jj] *= T(0.5*pi); 
				}
			}
			
			
			delete []  vec;    
		}	
	
		void VectorInside(std::complex<T>* Evals,
			int m_Ndof,
			std::complex<T>* Result)
		{
			if(!m_init)
			{
				InitializePlan(Evals); 
			}
			
			Tchv_cof(Evals); 
			
			for (int ii(0); ii < m_Ndof; ++ii)
			{
				Result[ii] += Evals[ii]; 
			}
		}
	
		void MatrixInside(std::complex<T>* Evals,
			int m_Ndof,
			std::complex<T>* Result)
		{
			
			if(!m_init)
			{
				InitializePlan(Evals); 
			}
            
			//could be done in parallel? (parece que si) 
			#pragma omp parallel for collapse(1) 
			for (int ii=0; ii < m_N; ii++)
			{
				Tchv_cof(Evals+ii*m_N); 

			}		
			
			std::complex<T>** Copy =
				new std::complex<T>*[m_Nthreads];
				
			for(int ii (0); ii <m_Nthreads; ++ii)
			{
				Copy[ii] = new std::complex<T>[m_N]; 
			}
			
		    #pragma omp parallel for collapse(1) 
			for (int ii=0; ii < m_Ndof; ii++)
			{
				int tid = omp_get_thread_num();							
				
				for(int jj(0); jj < m_N; ++jj)
				{
					Copy[tid][jj] = Evals[ii+jj*m_N];

				}
				
				Tchv_cof(Copy[tid]);

				//result dado en filas
				for( int jj(0); jj < m_Ndof; ++jj) 
				{					
					Result[jj+ii*m_Ndof] += Copy[tid][jj]; 
				}
				
			}	
			
			for(int ii(0); ii< m_Nthreads; ++ii)
			{
				delete [] Copy[ii];
			}
			
			delete [] Copy; 
		}
	
	
	public: 
	
		ChebTransformer(int Nthreads, int N):
		m_init(false),
		m_Nthreads(Nthreads),
		m_N(N)
		{
			omp_set_num_threads(Nthreads);			
		}
		
		 ~ChebTransformer()
        {
			if(m_init)
            {
				fftw_destroy_plan(m_plan);   
				fftw_cleanup();	 				
			}
        }
			
		T* chebpts()
		{
			T* pts = new T[m_N]; 
			
			for(int jj(0); jj < m_N; ++jj)
			{
				pts[jj] = cos((m_N-1-jj)*pi/(m_N-1));
			}
			
			return pts; 
		}
				
		void VectorIntegration(std::complex<T>* Evals,
			int m_Ndof,
			std::complex<T>* Result)
		{
			m_getIntegration = true; 
			
			VectorInside(Evals,m_Ndof,Result); 
			
		}
		
		void VectorCoef(std::complex<T>* Evals,
			int m_Ndof,
			std::complex<T>* Result)
		{
			m_getIntegration = false; 
			
			VectorInside(Evals,m_Ndof,Result); 
			
		}			
			
		void MatrixIntegration(std::complex<T>* Evals,
			int m_Ndof,
			std::complex<T>* Result)
		{
			m_getIntegration = true; 
			
			MatrixInside(Evals,m_Ndof,Result); 
				
		}
		
		void MatrixCoef(std::complex<T>* Evals,
			int m_Ndof,
			std::complex<T>* Result)
		{
			m_getIntegration = false; 
			
			MatrixInside(Evals,m_Ndof,Result); 
				
		}
        
        void ChangePoints(int N)
        {
            m_N = N; 
            
            m_init = false; 
        }
}; 

#endif