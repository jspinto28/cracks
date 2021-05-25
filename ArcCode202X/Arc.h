#ifndef ARC_H
#define ARC_H
#include <iostream>
#include "constant.h"

template <typename T> 
class ArcBase
{
	protected: 
	
		int m_index;
		
	public: 
	
		ArcBase(int index): 
			m_index(index)
		{}
		
		int GetIndex()
        {
            return m_index; 
        }
		
		virtual void GetPoint(T t, T *x, T *y) =0;

		virtual void GetDevPoint( T t, T *x, T *y) =0; 
		
		T Jacobian(T t)
		{
			T x,y; 
			
			GetDevPoint(t,&x,&y); 
			
			return sqrt(x*x+y*y); 
		}	
}; 

template <typename T> 
class ArcList : public ArcBase<T>
{
	public:

	ArcList(int index): 
		ArcBase<T>(index)
		{
		}
	
	void GetPoint(T t, T *x, T *y)
	{
		if(this->m_index== 1) 
		{ 
			*x = t; 
			
			*y = 0.2+0.28*t+0.3*(2.0*t*t-1.0)-0.5*(4.0*t*t*t-3.0*t); 
		
			// *x = 0.289413*0.5*(t+1)+ 4.000003; 

			// *y = 1.776822*sin(0.199908*t+2.057481)+16.000003; 
	
			return;
		} 

		if(this->m_index== 2) 
		{ 
			*x = t; 
			
			*y = 2+0.28*t; 
		
			// *x = 0.336164*0.5*(t+1)+ 1.500003; 

			// *y = 1.572857*sin(0.092395*t+3.312582)+16.000003; 
	
			return;
		} 

		if(this->m_index== 3) 
		{ 
			*x = 0.460626*0.5*(t+1)+ 2.500003; 

			*y = 1.754944*sin(0.161934*t+4.450605)+16.000003; 
			
			return;

		} 

		if(this->m_index== 4) 
		{ 
			*x = 0.487897*0.5*(t+1)+ 3.500003; 

			*y = 1.488565*sin(0.066611*t+6.163855)+16.000003; 
			
			return;

		} 

		if(this->m_index== 5) 
		{ 
			*x = 0.298058*0.5*(t+1)+ 3.000003; 

			*y = 1.994647*sin(0.030504*t+3.939598)+16.000003; 
			
			return;

		} 

		if(this->m_index== 6) 
		{ 
			*x = 0.339920*0.5*(t+1)+ 2.000003; 

			*y = 1.306060*sin(0.089160*t+2.819228)+16.000003; 
			
			return;

		} 

		if(this->m_index== 7) 
		{ 
			*x = 0.471085*0.5*(t+1)+ 2.000003; 

			*y = 1.062521*sin(0.190657*t+5.703751)+12.000003; 
			
			return;

		} 

		if(this->m_index== 8) 
		{ 
			*x = 0.381239*0.5*(t+1)+ 1.500003; 

			*y = 1.240953*sin(0.002210*t+5.843232)+12.000003; 
			
			return;

		} 

		if(this->m_index== 9) 
		{ 
			*x = 0.363771*0.5*(t+1)+ 2.500003; 

			*y = 1.071533*sin(0.089550*t+3.159880)+12.000003; 
			
			return;

		} 

		if(this->m_index== 10) 
		{ 
			*x = 0.493564*0.5*(t+1)+ 3.000003; 

			*y = 1.006572*sin(0.069831*t+4.836404)+12.000003; 
			
			return;

		} 

		if(this->m_index== 11) 
		{ 
			*x = 0.299936*0.5*(t+1)+ 3.500003; 

			*y = 1.434835*sin(0.073112*t+5.471992)+12.000003; 
			
			return;

		} 

		if(this->m_index== 12) 
		{ 
			*x = 0.388021*0.5*(t+1)+ 4.000003; 

			*y = 1.523866*sin(0.191748*t+2.059149)+12.000003; 
			
			return;

		} 

		if(this->m_index== 13) 
		{ 
			*x = 0.331134*0.5*(t+1)+ 4.000003; 

			*y = 1.209944*sin(0.190516*t+5.506406)+20.000003; 
			
			return;

		} 

		if(this->m_index== 14) 
		{ 
			*x = 0.288760*0.5*(t+1)+ 2.000003; 

			*y = 1.904723*sin(0.091736*t+3.684866)+20.000003; 
			
			return;

		} 

		if(this->m_index== 15) 
		{ 
			*x = 0.428749*0.5*(t+1)+ 3.000003; 

			*y = 1.251373*sin(0.055478*t+1.349997)+20.000003; 
			
			return;

		} 

		if(this->m_index== 16) 
		{ 
			*x = 0.286259*0.5*(t+1)+ 3.500003; 

			*y = 1.202787*sin(0.110663*t+6.195915)+20.000003; 
			
			return;

		} 

		if(this->m_index== 17) 
		{ 
			*x = 0.351353*0.5*(t+1)+ 1.500003; 

			*y = 1.097227*sin(0.130117*t+2.404943)+20.000003; 
			
			return;

		} 

		if(this->m_index== 18) 
		{ 
			*x = 0.428947*0.5*(t+1)+ 2.500003; 

			*y = 1.270708*sin(0.017917*t+5.641250)+20.000003; 
			
			return;

		} 

		if(this->m_index== 19) 
		{ 
			*x = 0.296572*0.5*(t+1)+ 2.000003; 

			*y = 1.585895*sin(0.023030*t+2.289206)+28.000003; 
			
			return;

		} 

		if(this->m_index== 20) 
		{ 
			*x = 0.451886*0.5*(t+1)+ 4.000003; 

			*y = 1.594124*sin(0.066120*t+3.153163)+28.000003; 
			
			return;

		} 

		if(this->m_index== 21) 
		{ 
			*x = 0.288799*0.5*(t+1)+ 3.000003; 

			*y = 1.897216*sin(0.164679*t+3.578781)+28.000003;

			return;

		} 

		if(this->m_index== 22) 
		{ 
			*x = 0.426558*0.5*(t+1)+ 1.500003; 

			*y = 1.899389*sin(0.191497*t+0.957194)+28.000003; 
			
			return;

		} 

		if(this->m_index== 23) 
		{ 
			*x = 0.348412*0.5*(t+1)+ 2.500003; 

			*y = 1.868069*sin(0.011534*t+3.570715)+28.000003; 
			
			return;

		} 

		if(this->m_index== 24) 
		{ 
			*x = 0.300720*0.5*(t+1)+ 3.500003; 

			*y = 1.053193*sin(0.020675*t+6.095131)+28.000003; 
			
			return;

		} 

		if(this->m_index== 25) 
		{ 
			*x = 0.472982*0.5*(t+1)+ 1.500003; 

			*y = 1.571249*sin(0.144650*t+2.419822)+24.000003; 
			
			return;

		} 

		if(this->m_index== 26) 
		{ 
			*x = 0.467364*0.5*(t+1)+ 2.500003; 

			*y = 1.067616*sin(0.149310*t+0.259148)+24.000003; 
			
			return;

		} 

		if(this->m_index== 27) 
		{ 
			*x = 0.336677*0.5*(t+1)+ 3.000003; 

			*y = 1.661009*sin(0.117396*t+3.694377)+24.000003; 
			
			return;

		} 

		if(this->m_index== 28) 
		{ 
			*x = 0.349155*0.5*(t+1)+ 2.000003; 

			*y = 1.895247*sin(0.136098*t+2.394704)+24.000003; 
			
			return;

		} 
	}

	void GetDevPoint(T t, T *x, T *y)
	{
		if(this->m_index== 1) 
		{ 
			*x = 1.0; 
			
			*y = 0.28+0.3*(4.0*t)-0.5*(12.0*t*t-3.0); 
		
			// *x = 0.289413*0.5; 

			// *y = 1.776822*0.199908*cos(0.199908*t+2.057481); 
			
			return;

		} 

		if(this->m_index== 2) 
		{ 
			*x =1.0;
			
			*y = 0.28; 
		
			// *x = 0.336164*0.5; 

			// *y = 1.572857*0.092395*cos(0.092395*t+3.312582); 
			
			return;

		} 

		if(this->m_index== 3) 
		{ 
			*x = 0.460626*0.5; 

			*y = 1.754944*0.161934*cos(0.161934*t+4.450605); 
			
			return;

		} 

		if(this->m_index== 4) 
		{ 
			*x = 0.487897*0.5; 

			*y = 1.488565*0.066611*cos(0.066611*t+6.163855); 
			
			return;

		} 

		if(this->m_index== 5) 
		{ 
			*x = 0.298058*0.5; 

			*y = 1.994647*0.030504*cos(0.030504*t+3.939598); 
			
			return;

		} 

		if(this->m_index== 6) 
		{ 
			*x = 0.339920*0.5; 

			*y = 1.306060*0.089160*cos(0.089160*t+2.819228); 
			
			return;

		} 

		if(this->m_index== 7) 
		{ 
			*x = 0.471085*0.5; 

			*y = 1.062521*0.190657*cos(0.190657*t+5.703751); 
			
			return;

		} 

		if(this->m_index== 8) 
		{ 
			*x = 0.381239*0.5; 

			*y = 1.240953*0.002210*cos(0.002210*t+5.843232); 
			
			return;

		} 

		if(this->m_index== 9) 
		{ 
			*x = 0.363771*0.5; 

			*y = 1.071533*0.089550*cos(0.089550*t+3.159880); 
			
			return;

		} 

		if(this->m_index== 10) 
		{ 
			*x = 0.493564*0.5; 

			*y = 1.006572*0.069831*cos(0.069831*t+4.836404); 
			
			return;

		} 

		if(this->m_index== 11) 
		{ 
			*x = 0.299936*0.5; 

			*y = 1.434835*0.073112*cos(0.073112*t+5.471992); 
			
			return;

		} 

		if(this->m_index== 12) 
		{ 
			*x = 0.388021*0.5; 

			*y = 1.523866*0.191748*cos(0.191748*t+2.059149); 
			
			return;

		} 

		if(this->m_index== 13) 
		{ 
			*x = 0.331134*0.5; 

			*y = 1.209944*0.190516*cos(0.190516*t+5.506406); 
			
			return;

		} 

		if(this->m_index== 14) 
		{ 
			*x = 0.288760*0.5; 

			*y = 1.904723*0.091736*cos(0.091736*t+3.684866); 
			
			return;

		} 

		if(this->m_index== 15) 
		{ 
			*x = 0.428749*0.5; 

			*y = 1.251373*0.055478*cos(0.055478*t+1.349997); 
			
			return;

		} 

		if(this->m_index== 16) 
		{ 
			*x = 0.286259*0.5; 

			*y = 1.202787*0.110663*cos(0.110663*t+6.195915); 
			
			return;

		} 

		if(this->m_index== 17) 
		{ 
			*x = 0.351353*0.5; 

			*y = 1.097227*0.130117*cos(0.130117*t+2.404943); 
			
			return;

		} 

		if(this->m_index== 18) 
		{ 
			*x = 0.428947*0.5; 

			*y = 1.270708*0.017917*cos(0.017917*t+5.641250); 
			
			return;

		} 

		if(this->m_index== 19) 
		{ 
			*x = 0.296572*0.5; 

			*y = 1.585895*0.023030*cos(0.023030*t+2.289206); 
			
			return;

		} 

		if(this->m_index== 20) 
		{ 
			*x = 0.451886*0.5; 

			*y = 1.594124*0.066120*cos(0.066120*t+3.153163); 
			
			return;

		} 

		if(this->m_index== 21) 
		{ 
			*x = 0.288799*0.5; 

			*y = 1.897216*0.164679*cos(0.164679*t+3.578781); 
			
			return;

		} 

		if(this->m_index== 22) 
		{ 
			*x = 0.426558*0.5; 

			*y = 1.899389*0.191497*cos(0.191497*t+0.957194); 
			
			return;

		} 

		if(this->m_index== 23) 
		{ 
			*x = 0.348412*0.5; 

			*y = 1.868069*0.011534*cos(0.011534*t+3.570715); 
			
			return;

		} 

		if(this->m_index== 24) 
		{ 
			*x = 0.300720*0.5; 

			*y = 1.053193*0.020675*cos(0.020675*t+6.095131); 
			
			return;

		} 

		if(this->m_index== 25) 
		{ 
			*x = 0.472982*0.5; 

			*y = 1.571249*0.144650*cos(0.144650*t+2.419822); 
			
			return;

		} 

		if(this->m_index== 26) 
		{ 
			*x = 0.467364*0.5; 

			*y = 1.067616*0.149310*cos(0.149310*t+0.259148); 
			
			return;

		} 

		if(this->m_index== 27) 
		{ 
			*x = 0.336677*0.5; 

			*y = 1.661009*0.117396*cos(0.117396*t+3.694377); 
			
			return;

		} 

		if(this->m_index== 28) 
		{ 
			*x = 0.349155*0.5; 

			*y = 1.895247*0.136098*cos(0.136098*t+2.394704); 
			
			return;

		} 
	}

};

template <typename T> 
class ArcParametrized : public ArcBase<T>
{
	protected: 
	
		int m_Nx, m_Ny; 
		
		T* m_coefX, *m_coefY;
			
	
	public: 
	
		ArcParametrized( int index, 
			int Nx, T* coefX, 
			int Ny, T* coefY):
			ArcBase<T>(index),
			m_Nx(Nx),
			m_Ny(Ny), 
			m_coefX(coefX),
			m_coefY(coefY)
			{}			
				
		void GetPoint(T t, T *x, T *y)
		{		
			*x =m_coefX[0] + m_coefX[1]*t; 
			*y =m_coefY[0] + m_coefY[1]*t;
            
            int NXdiv2 = (int)((m_Nx-2)/2);
			int NYdiv2 = (int)((m_Ny-2)/2);
                    
			for(int ii(0); ii < NXdiv2; ++ii)
			{
				*x += m_coefX[2*ii+2]*cos(T(ii)*t)+
                        m_coefX[2*ii+3]*sin(T(ii+1)*t);                
			}
            
            if( ((m_Nx-2)%2))
            {
                *x += m_coefX[2*NXdiv2+2]*cos(T(NXdiv2)*t); 
            }
                
            for(int ii(0); ii < NYdiv2; ++ii)
			{
				*y += m_coefY[2*ii+2]*cos(T(ii)*t)+
                        m_coefY[2*ii+3]*sin(T(ii+1)*t);                
			}
            
            if( ((m_Ny-2)%2))
            {
                *y += m_coefY[2*NYdiv2+2]*cos(T(NYdiv2)*t); 
            }			
			

		}
		
		
		void GetDevPoint(T t, T *x, T *y)
		{
			*x = m_coefX[1]; 
			*y = m_coefY[1];
            
            int NXdiv2 = (int)((m_Nx-2)/2);
			int NYdiv2 = (int)((m_Ny-2)/2);
                    
			for(int ii(0); ii < NXdiv2; ++ii)
			{
				*x += -T(ii)*m_coefX[2*ii+2]*sin(T(ii)*t)+
                        T(ii+1)*m_coefX[2*ii+3]*cos(T(ii+1)*t);                
			}
            
            if( ((m_Nx-2)%2))
            {
                *x -= T(NXdiv2)*m_coefX[2*NXdiv2+2]*sin(T(NXdiv2)*t); 
            }
                
            for(int ii(0); ii < NYdiv2; ++ii)
			{
				*y += -T(ii)*m_coefY[2*ii+2]*sin(T(ii)*t)+
                        T(ii+1)*m_coefY[2*ii+3]*cos(T(ii+1)*t);                
			}
            
            if( ((m_Ny-2)%2))
            {
                *y -= T(NYdiv2)*m_coefY[2*NYdiv2+2]*sin(T(NYdiv2)*t); 
            }		
		}
		
};


// template <typename T> 
// class ArcParametrized : public ArcBase<T>
// {
// 	protected: 
// 	
// 		int m_Nx, m_Ny; 
// 		
// 		T* m_coefX, *m_coefY;
// 	
// 		T* PolsT(int N,T t) 
// 		{
// 			T* chebpols = new T[N]; 
// 			
// 			T acost = acos(t); 
// 			
// 			for(int ii(0);ii<N; ++ii)
// 			{
// 				chebpols[ii] = cos(T(ii)*acost); 
// 			}
// 			
// 			return chebpols;
// 		}
// 		
// 		T* PolsnU(int N,T t) 
// 		{
// 			T* chebpols = new T[N]; 
// 			
// 			T acost = acos(t); 
// 			
// 			if( acost==0)
// 			{
// 				for(int ii(0);ii<N; ++ii)
// 				{
// 					chebpols[ii] = (T(ii)+1)*T(ii+1); 
// 				}
// 				
// 				if(t < 0)
// 				{
// 					for(int ii(1);ii<N; ii+=2)
// 					{
// 						chebpols[ii] *= T(-1.0); 
// 					}
// 				}
// 			}
// 			else
// 			{
// 				for(int ii(0);ii<N; ++ii)
// 				{
// 					chebpols[ii] = T(ii+1)*
// 						sin(T(ii+1)*acost)/sin(acost); 
// 				}
// 				
// 			}		
// 			
// 			return chebpols;
// 		}
// 	
// 	public: 
// 	
// 		ArcParametrized( int index, 
// 			int Nx, T* coefX, 
// 			int Ny, T* coefY):
// 			ArcBase<T>(index),
// 			m_Nx(Nx),
// 			m_Ny(Ny), 
// 			m_coefX(coefX),
// 			m_coefY(coefY)
// 			{}			
// 				
// 		void GetPoint(T t, T *x, T *y)
// 		{
// 			T* polsX = PolsT(m_Nx,t);
// 			T* polsY = PolsT(m_Ny,t); 
// 			
// 			*x =0; 
// 			*y =0; 
// 			
// 			for(int ii(0); ii < m_Nx; ++ii)
// 			{
// 				*x += m_coefX[ii]*polsX[ii];  
// 			}
// 			
// 			for(int ii(0); ii < m_Ny; ++ii)
// 			{
// 				*y += m_coefY[ii]*polsY[ii];  
// 			}
// 			
// 			
// 			delete [] polsX; 
// 			delete [] polsY; 
// 		}
// 		
// 		
// 		void GetDevPoint(T t, T *x, T *y)
// 		{
// 			T* polsX = PolsnU(m_Nx,t);
// 			T* polsY = PolsnU(m_Ny,t); 
// 			
// 			*x =0; 
// 			*y=0; 
// 			
// 			for(int ii(1); ii < m_Nx; ++ii)
// 			{
// 				
// 				*x += m_coefX[ii]*polsX[ii-1];  
// 			}
// 			
// 			for(int ii(1); ii < m_Ny; ++ii)
// 			{
// 				*y += m_coefY[ii]*polsY[ii-1];  
// 			}
// 			
// 			
// 			
// 			
// 			delete [] polsX; 
// 			delete [] polsY; 
// 		}
// 		
// };

#endif