#include "mex.h"
#include <iostream>
#include <complex>
#include "Operator.h"
#include "Arc.h"
#include "RightHandSide.h"
#include "ChebTransformer.h"
#include "GreenFunction.h"

using namespace std;

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int fun = (int)mxGetScalar(prhs[0]);
	int usefloat = (int)mxGetScalar(prhs[1]);
	
	if(fun==0)
	{

		ChebTransformer<float>* m_Transformer = new ChebTransformer<float>(1,32);
 
		std::complex<float>* evals = new std::complex<float>[1024]; 
		std::complex<float>* vec = new std::complex<float>[64]; 
		
		float* pts =  m_Transformer->chebpts();
		
		int nt =3; 
		
		int ns = 2; 
		
		for (int ii(0); ii < 32; ++ii)
		{			
		
			for (int jj(0); jj < 32; ++jj)
			{
				evals[jj+ii*32] =
				std::complex<float>(
				cos(acos(pts[ii])*float(nt))*
				cos(acos(pts[jj])*float(ns)),
				0.0);
				
			}
			
		}
		
		m_Transformer->MatrixCoef(evals,8,vec); 
		
		for(int ii(0); ii < 8; ++ii)
		{
			for(int jj(0); jj < 8; ++jj)
			{
				std::cout<<vec[ii+jj*8]<<" "; 
					
			}
			std::cout<<std::endl; 
		}


		double *pR,*pI; 
		mxArray *p_m; 
		
		p_m=plhs[0] =mxCreateDoubleMatrix(1,1,mxCOMPLEX);
		pR=  mxGetPr(p_m);
        pI=  mxGetPi(p_m);
		
		for(int ii(0); ii < 1; ++ii)
		{
			pR[ii] = 2.8; 
			pI[ii] = 0.28; 
		}			


		delete m_Transformer; 
		delete [] pts; 
		delete [] evals;
		delete [] vec; 
		
		
	}
	
    //V helmholtz float list
	else if((fun ==1)&&(usefloat==1))
	{
		int indTrial, indTest, Ndofs,nthreads; 
		float wavenumber; 
		double *pR,*pI; 
		
		indTrial = (int)mxGetScalar(prhs[2]);
		indTest = (int)mxGetScalar(prhs[3]);
		wavenumber = mxGetScalar(prhs[4]);
		Ndofs = (int)mxGetScalar(prhs[5]);
		nthreads = (int)mxGetScalar(prhs[6]);
		
		ArcBase<float>* arcTrial =
			new ArcList<float>(indTrial); 
			
		ArcBase<float>* arcTest =
			new ArcList<float>(indTest); 

		std::complex<float>* Matrix = new 
			std::complex<float>[Ndofs*Ndofs]; 
		
		for(int ii(0); ii < Ndofs*Ndofs; ++ii)
		{
			Matrix[ii] = std::complex<float>(0.0,0.0); 
		}
		
		OperatorBase<float>* V = new VHemholtz<float>
			(Ndofs,nthreads,
			wavenumber,arcTrial,arcTest); 
			
		V->Compute(Matrix); 
		
		mxArray *p_m; 
		
		p_m=plhs[0] =mxCreateDoubleMatrix(Ndofs,Ndofs,mxCOMPLEX);
		pR=  mxGetPr(p_m);
        pI=  mxGetPi(p_m);
		
		for(int ii(0); ii < Ndofs*Ndofs; ++ii)
		{
			pR[ii] = Matrix[ii].real(); 
			pI[ii] = Matrix[ii].imag(); 
		}			
		
		delete V;
		delete arcTrial; 
		delete arcTest; 
		delete [] Matrix; 			
		
		
	}
	
    //V helmholtz double list
	else if((fun ==1)&&(usefloat==0))
	{
		int indTrial, indTest, Ndofs,nthreads; 
		double wavenumber; 
		double *pR,*pI; 
		
		indTrial = (int)mxGetScalar(prhs[2]);
		indTest = (int)mxGetScalar(prhs[3]);
		wavenumber = mxGetScalar(prhs[4]);
		Ndofs = (int)mxGetScalar(prhs[5]);
		nthreads = (int)mxGetScalar(prhs[6]);
		
		ArcBase<double>* arcTrial =
			new ArcList<double>(indTrial); 
			
		ArcBase<double>* arcTest =
			new ArcList<double>(indTest); 

		std::complex<double>* Matrix = new 
			std::complex<double>[Ndofs*Ndofs]; 
		
		for(int ii(0); ii < Ndofs*Ndofs; ++ii)
		{
			Matrix[ii] = std::complex<double>(0.0,0.0); 
		}
		
		OperatorBase<double>* V = new VHemholtz<double>
			(Ndofs,nthreads,
			wavenumber,arcTrial,arcTest); 
			
		V->Compute(Matrix); 
		
		mxArray *p_m; 
		
		p_m=plhs[0] =mxCreateDoubleMatrix(Ndofs,Ndofs,mxCOMPLEX);
		pR=  mxGetPr(p_m);
        pI=  mxGetPi(p_m);
		
		for(int ii(0); ii < Ndofs*Ndofs; ++ii)
		{
			pR[ii] = Matrix[ii].real(); 
			pI[ii] = Matrix[ii].imag(); 
		}			
		
		delete V;
		delete arcTrial; 
		delete arcTest; 
		delete [] Matrix; 		
	}
	
    //V helmholtz double coefs
	else if((fun ==11)&&(usefloat==0))
	{
		int indTrial, indTest, Ndofs,nthreads; 
		int NxTrial,NyTrial, NxTest, NyTest; 
		double *coefxTrial, *coefyTrial, *coefxTest,*coefyTest; 
		double wavenumber; 
		double *pR,*pI; 
		
		indTrial = (int)mxGetScalar(prhs[2]);
		
		NxTrial = (int) mxGetScalar(prhs[3]); 
		NyTrial = (int) mxGetScalar(prhs[4]); 
		
		coefxTrial = mxGetPr(prhs[5]); 
		coefyTrial = mxGetPr(prhs[6]); 
		
		indTest = (int)mxGetScalar(prhs[7]);
		
		NxTest = (int) mxGetScalar(prhs[8]);
		NyTest = (int) mxGetScalar(prhs[9]);
		
		coefxTest = mxGetPr(prhs[10]); 
		coefyTest = mxGetPr(prhs[11]); 

		wavenumber = mxGetScalar(prhs[12]);
		Ndofs = (int)mxGetScalar(prhs[13]);
		nthreads = (int)mxGetScalar(prhs[14]);
		
		ArcBase<double>* arcTrial =
			new ArcParametrized<double>(indTrial,
				NxTrial, coefxTrial,
				NyTrial,coefyTrial); 
			
		ArcBase<double>* arcTest =
			new ArcParametrized<double>(indTest,
				NxTest,coefxTest,
				NyTest,coefyTest); 

		std::complex<double>* Matrix = new 
			std::complex<double>[Ndofs*Ndofs]; 
		
		for(int ii(0); ii < Ndofs*Ndofs; ++ii)
		{
			Matrix[ii] = std::complex<double>(0.0,0.0); 
		}
		
		OperatorBase<double>* V = new VHemholtz<double>
			(Ndofs,nthreads,
			wavenumber,arcTrial,arcTest); 
			
		V->Compute(Matrix); 
		
		mxArray *p_m; 
		
		p_m=plhs[0] =mxCreateDoubleMatrix(Ndofs,Ndofs,mxCOMPLEX);
		pR=  mxGetPr(p_m);
        pI=  mxGetPi(p_m);
		
		for(int ii(0); ii < Ndofs*Ndofs; ++ii)
		{
			pR[ii] = Matrix[ii].real(); 
			pI[ii] = Matrix[ii].imag(); 
		}			
		
		delete V;
		delete arcTrial; 
		delete arcTest; 
		delete [] Matrix; 		
	}
        
	//V elasticity float list
	else if((fun==2)&&(usefloat==1))
	{
		int indTrial, indTest, Ndofs,nthreads; 
		float omega, lambda,mu; 
		double *pR,*pI; 
		
		indTrial = (int)mxGetScalar(prhs[2]);
		indTest = (int)mxGetScalar(prhs[3]);
		omega = mxGetScalar(prhs[4]);
		lambda = mxGetScalar(prhs[5]);
		mu = mxGetScalar(prhs[6]);
		Ndofs = (int)mxGetScalar(prhs[7]);
		nthreads = (int)mxGetScalar(prhs[8]);
				
		ArcBase<float>* arcTrial =
			new ArcList<float>(indTrial); 
			
		ArcBase<float>* arcTest =
			new ArcList<float>(indTest); 

		std::complex<float>* Matrix = new 
			std::complex<float>[2*Ndofs*2*Ndofs]; 
		
		for(int ii(0); ii < 4*Ndofs*Ndofs; ++ii)
		{
			Matrix[ii] = std::complex<float>(0.0,0.0); 
		}
		
		OperatorBase<float>* V = new VElasticity<float>
			(Ndofs,nthreads,
			omega,lambda,mu,arcTrial,arcTest); 
			
		V->Compute(Matrix); 
		
		mxArray *p_m; 
		
		p_m=plhs[0] =mxCreateDoubleMatrix(2*Ndofs,2*Ndofs,mxCOMPLEX);
		pR=  mxGetPr(p_m);
        pI=  mxGetPi(p_m);
		
		for(int ii(0); ii < 4*Ndofs*Ndofs; ++ii)
		{
			pR[ii] = Matrix[ii].real(); 
			pI[ii] = Matrix[ii].imag(); 
		}			
		
		delete V;
		delete arcTrial; 
		delete arcTest; 
		delete [] Matrix; 		
	}
		
    //V elasticity double list
	else if((fun==2)&&(usefloat==0))
	{
		int indTrial, indTest, Ndofs,nthreads; 
		double omega, lambda,mu; 
		double *pR,*pI; 
		
		indTrial = (int)mxGetScalar(prhs[2]);
		indTest = (int)mxGetScalar(prhs[3]);
		omega = mxGetScalar(prhs[4]);
		lambda = mxGetScalar(prhs[5]);
		mu = mxGetScalar(prhs[6]);
		Ndofs = (int)mxGetScalar(prhs[7]);
		nthreads = (int)mxGetScalar(prhs[8]);
				
		ArcBase<double>* arcTrial =
			new ArcList<double>(indTrial); 
			
		ArcBase<double>* arcTest =
			new ArcList<double>(indTest); 

		std::complex<double>* Matrix = new 
			std::complex<double>[2*Ndofs*2*Ndofs]; 
		
		for(int ii(0); ii < 4*Ndofs*Ndofs; ++ii)
		{
			Matrix[ii] = std::complex<double>(0.0,0.0); 
		}
		
		OperatorBase<double>* V = new VElasticity<double>
			(Ndofs,nthreads,
			omega,lambda,mu,arcTrial,arcTest); 
			
		V->Compute(Matrix); 
		
		mxArray *p_m; 
		
		p_m=plhs[0] =mxCreateDoubleMatrix(2*Ndofs,2*Ndofs,mxCOMPLEX);
		pR=  mxGetPr(p_m);
        pI=  mxGetPi(p_m);
		
		for(int ii(0); ii < 4*Ndofs*Ndofs; ++ii)
		{
			pR[ii] = Matrix[ii].real(); 
			pI[ii] = Matrix[ii].imag(); 
		}			
		
		delete V;
		delete arcTrial; 
		delete arcTest; 
		delete [] Matrix; 		
	}
	
    //V elasticity double coefs
	else if((fun==21)&&(usefloat==0))
	{		
		int indTrial, indTest, Ndofs,nthreads; 
		int NxTrial,NyTrial, NxTest, NyTest; 
		double *coefxTrial, *coefyTrial, *coefxTest,*coefyTest; 
		double omega, lambda,mu; 
		double *pR,*pI; 
		
		indTrial = (int)mxGetScalar(prhs[2]);
		
		NxTrial = (int) mxGetScalar(prhs[3]); 
		NyTrial = (int) mxGetScalar(prhs[4]); 
		
		coefxTrial = mxGetPr(prhs[5]); 
		coefyTrial = mxGetPr(prhs[6]); 
		
		indTest = (int)mxGetScalar(prhs[7]);
		
		NxTest = (int) mxGetScalar(prhs[8]);
		NyTest = (int) mxGetScalar(prhs[9]);
		
		coefxTest = mxGetPr(prhs[10]); 
		coefyTest = mxGetPr(prhs[11]); 
		
		omega = mxGetScalar(prhs[12]);
		lambda = mxGetScalar(prhs[13]);
		mu = mxGetScalar(prhs[14]);
		Ndofs = (int)mxGetScalar(prhs[15]);
		nthreads = (int)mxGetScalar(prhs[16]);
				
		ArcBase<double>* arcTrial =
			new ArcParametrized<double>(indTrial,
				NxTrial, coefxTrial,
				NyTrial,coefyTrial);         
							
		ArcBase<double>* arcTest =
			new ArcParametrized<double>(indTest,
				NxTest,coefxTest,
				NyTest,coefyTest); 

		std::complex<double>* Matrix = new 
			std::complex<double>[2*Ndofs*2*Ndofs]; 
		
		for(int ii(0); ii < 4*Ndofs*Ndofs; ++ii)
		{
			Matrix[ii] = std::complex<double>(0.0,0.0); 
		}
		
        double tol = mxGetScalar(prhs[17]);
       
        int Nlevs = (int)mxGetScalar(prhs[18]);      
        
		OperatorBase<double>* V = new VElasticity<double>
			(Ndofs,nthreads,tol,Nlevs,
			omega,lambda,mu,arcTrial,arcTest); 
 		
 		V->Compute(Matrix); 
				
		mxArray *p_m; 
		
		p_m=plhs[0] =mxCreateDoubleMatrix(2*Ndofs,2*Ndofs,mxCOMPLEX);
		pR=  mxGetPr(p_m);
        pI=  mxGetPi(p_m);
		
		for(int ii(0); ii < 4*Ndofs*Ndofs; ++ii)
		{
			pR[ii] = Matrix[ii].real(); 
			pI[ii] = Matrix[ii].imag(); 
		}			
					
		// delete V;
		delete arcTrial; 
		delete arcTest; 
		delete [] Matrix; 		
	}
	
    // plane wave helmholtz;
	else if(fun == 3)
	{
		
		int  indTest, Ndofs;
		double wavenumber,alpha; 
		double *pR,*pI; 

		indTest = (int)mxGetScalar(prhs[2]);
		wavenumber = mxGetScalar(prhs[3]);
		alpha = mxGetScalar(prhs[4]);
		Ndofs = (int)mxGetScalar(prhs[5]);

		ArcBase<double>* arcTest =
			new ArcList<double>(indTest); 
			
		std::complex<double>* Vector = new 
			std::complex<double>[Ndofs]; 
			
		for(int ii(0); ii < Ndofs; ++ii)
		{
			Vector[ii] = std::complex<double>(0.0,0.0); 
		}
		
		RhsBase<double>* V = new PlaneWave<double>
			(Ndofs,wavenumber, alpha,arcTest); 
		
		V->Compute(Vector); 
		
		mxArray *p_m; 
		
		p_m=plhs[0] =mxCreateDoubleMatrix(Ndofs,1,mxCOMPLEX);
		pR=  mxGetPr(p_m);
        pI=  mxGetPi(p_m);
		
		for(int ii(0); ii < Ndofs; ++ii)
		{
			pR[ii] = Vector[ii].real(); 
			pI[ii] = Vector[ii].imag(); 
		}			
		
		delete V;
		delete arcTest; 
		delete [] Vector; 		
		
	}
	
    // plane wave p;
	else if(fun == 4)
	{
		
		int  indTest, Ndofs;
		double omega, lambda , mu ,alpha; 
		double *pR,*pI; 

		indTest = (int)mxGetScalar(prhs[2]);
		omega = mxGetScalar(prhs[3]);
		lambda = mxGetScalar(prhs[4]);
		mu = mxGetScalar(prhs[5]);
		alpha = mxGetScalar(prhs[6]);
		Ndofs = (int)mxGetScalar(prhs[7]);

		ArcBase<double>* arcTest =
			new ArcList<double>(indTest); 
			
		std::complex<double>* Vector = new 
			std::complex<double>[2*Ndofs]; 
			
		for(int ii(0); ii < 2*Ndofs; ++ii)
		{
			Vector[ii] = std::complex<double>(0.0,0.0); 
		}
		
		RhsBase<double>* V = new ElasticityP<double>
			(Ndofs,omega,lambda,mu, alpha,arcTest); 
		
		V->Compute(Vector); 
		
		mxArray *p_m; 
		
		p_m=plhs[0] =mxCreateDoubleMatrix(2*Ndofs,1,mxCOMPLEX);
		pR=  mxGetPr(p_m);
        pI=  mxGetPi(p_m);
		
		for(int ii(0); ii < 2*Ndofs; ++ii)
		{
			pR[ii] = Vector[ii].real(); 
			pI[ii] = Vector[ii].imag(); 
		}			
		
		delete V;
		delete arcTest; 
		delete [] Vector; 		
		
	}
	
    // plane wave p;
	else if(fun == 41)
	{
		
		int  indTest, Ndofs;
		int Nx, Ny; 
		double *coefx, *coefy; 
		double omega, lambda , mu ,alpha; 
		double *pR,*pI; 

		indTest = (int)mxGetScalar(prhs[2]);
		
		Nx = (int)mxGetScalar(prhs[3]);
		Ny = (int)mxGetScalar(prhs[4]);
		
		coefx = mxGetPr(prhs[5]); 
		coefy = mxGetPr(prhs[6]); 
		
		omega = mxGetScalar(prhs[7]);
		lambda = mxGetScalar(prhs[8]);
		mu = mxGetScalar(prhs[9]);
		alpha = mxGetScalar(prhs[10]);
		Ndofs = (int)mxGetScalar(prhs[11]);

		ArcBase<double>* arcTest =
			new ArcParametrized<double>(indTest,Nx,coefx,Ny,coefy); 
			
		std::complex<double>* Vector = new 
			std::complex<double>[2*Ndofs]; 
			
		for(int ii(0); ii < 2*Ndofs; ++ii)
		{
			Vector[ii] = std::complex<double>(0.0,0.0); 
		}
		
		RhsBase<double>* V = new ElasticityP<double>
			(Ndofs,omega,lambda,mu, alpha,arcTest); 
		
		V->Compute(Vector); 
		
		mxArray *p_m; 
		
		p_m=plhs[0] =mxCreateDoubleMatrix(2*Ndofs,1,mxCOMPLEX);
		pR=  mxGetPr(p_m);
        pI=  mxGetPi(p_m);
		
		for(int ii(0); ii < 2*Ndofs; ++ii)
		{
			pR[ii] = Vector[ii].real(); 
			pI[ii] = Vector[ii].imag(); 
		}			
		
		delete V;
		delete arcTest; 
		delete [] Vector; 		
		
	}
	
    // plane wave s;
	else if(fun == 5)
	{
		
		int  indTest, Ndofs;
		double omega, lambda , mu ,alpha; 
		double *pR,*pI; 

		indTest = (int)mxGetScalar(prhs[2]);
		omega = mxGetScalar(prhs[3]);
		lambda = mxGetScalar(prhs[4]);
		mu = mxGetScalar(prhs[5]);
		alpha = mxGetScalar(prhs[6]);
		Ndofs = (int)mxGetScalar(prhs[7]);

		ArcBase<double>* arcTest =
			new ArcList<double>(indTest); 
			
		std::complex<double>* Vector = new 
			std::complex<double>[2*Ndofs]; 
			
		for(int ii(0); ii < 2*Ndofs; ++ii)
		{
			Vector[ii] = std::complex<double>(0.0,0.0); 
		}
		
		RhsBase<double>* V = new ElasticityS<double>
			(Ndofs,omega,lambda,mu, alpha,arcTest); 
		
		V->Compute(Vector); 
		
		mxArray *p_m; 
		
		p_m=plhs[0] =mxCreateDoubleMatrix(2*Ndofs,1,mxCOMPLEX);
		pR=  mxGetPr(p_m);
        pI=  mxGetPi(p_m);
		
		for(int ii(0); ii < 2*Ndofs; ++ii)
		{
			pR[ii] = Vector[ii].real(); 
			pI[ii] = Vector[ii].imag(); 
		}			
		
		delete V;
		delete arcTest; 
		delete [] Vector; 		
		
	}
	
    // plane wave s;
	else if(fun == 51)
	{
		
		int  indTest, Ndofs;
		int Nx, Ny; 
		double *coefx, *coefy; 
		double omega, lambda , mu ,alpha; 
		double *pR,*pI; 

		indTest = (int)mxGetScalar(prhs[2]);
		
		Nx = (int)mxGetScalar(prhs[3]);
		Ny = (int)mxGetScalar(prhs[4]);
		
		coefx = mxGetPr(prhs[5]); 
		coefy = mxGetPr(prhs[6]); 
		
		omega = mxGetScalar(prhs[7]);
		lambda = mxGetScalar(prhs[8]);
		mu = mxGetScalar(prhs[9]);
		alpha = mxGetScalar(prhs[10]);
		Ndofs = (int)mxGetScalar(prhs[11]);

		ArcBase<double>* arcTest =
			new ArcParametrized<double>(indTest,Nx,coefx,Ny,coefy); 
			
		std::complex<double>* Vector = new 
			std::complex<double>[2*Ndofs]; 
			
		for(int ii(0); ii < 2*Ndofs; ++ii)
		{
			Vector[ii] = std::complex<double>(0.0,0.0); 
		}
		
		RhsBase<double>* V = new ElasticityS<double>
			(Ndofs,omega,lambda,mu, alpha,arcTest); 
		
		V->Compute(Vector); 
		
		mxArray *p_m; 
		
		p_m=plhs[0] =mxCreateDoubleMatrix(2*Ndofs,1,mxCOMPLEX);
		pR=  mxGetPr(p_m);
        pI=  mxGetPi(p_m);
		
		for(int ii(0); ii < 2*Ndofs; ++ii)
		{
			pR[ii] = Vector[ii].real(); 
			pI[ii] = Vector[ii].imag(); 
		}			
		
		delete V;
		delete arcTest; 
		delete [] Vector; 		
		
	}
	
	
	//falta farfield helmholtz...
	
    //Farfield P list
	else if( fun == 6) 
	{
		
		int  ind, Ndofs;
		double omega, lambda , mu ,beta; 
		double *pR,*pI; 
		double *LambdaR, *LambdaI; 
		
		ind = (int)mxGetScalar(prhs[2]);
		omega = mxGetScalar(prhs[3]);
		lambda = mxGetScalar(prhs[4]);
		mu = mxGetScalar(prhs[5]);
		beta = mxGetScalar(prhs[6]);
		Ndofs = (int)mxGetScalar(prhs[7]);
		LambdaR = mxGetPr(prhs[8]); 
		LambdaI = mxGetPr(prhs[9]);

		ArcBase<double>* arc =
			new ArcList<double>(ind); 
			
		std::complex<double>* Vector = new 
			std::complex<double>[2*Ndofs]; 
			
		std::complex<double>* Lambda = new 
			std::complex<double>[2*Ndofs]; 
			
		for(int ii(0); ii < 2*Ndofs; ++ii)
		{
			Vector[ii] = std::complex<double>(0.0,0.0); 
			Lambda[ii] = std::complex<double>(LambdaR[ii],
				LambdaI[ii]); 
		}
		
		RhsBase<double>* V = new ElasticityP<double>
			(Ndofs,omega,lambda,mu, beta,arc); 
		
		V->Compute(Vector); 
		
		std::complex<double> z (0.0,0.0); 
		
		for(int ii(0); ii < Ndofs; ++ii)
		{
			Vector[ii] = std::conj(Vector[ii]); 
			Vector[ii+Ndofs] = std::conj(Vector[ii+Ndofs]); 
			z += Vector[ii]*Lambda[ii] + Vector[ii+Ndofs]*Lambda[ii+Ndofs];
		}
		
		double kp = sqrt(omega*omega/(lambda+2*mu)); 
		
		std::complex<double> aux (1,1); 
		
		aux *= 0.5*z/mu/sqrt(2)/sqrt(8*pi*kp); 			
		
		mxArray *p_m; 
		
		p_m=plhs[0] =mxCreateDoubleMatrix(2,1,mxCOMPLEX);
		pR=  mxGetPr(p_m);
        pI=  mxGetPi(p_m);
		
		pR[0] = aux.real()*cos(beta); 
		pI[0] = aux.imag()*cos(beta); 
		pR[1] = aux.real()*sin(beta); 
		pI[1] = aux.imag()*sin(beta); 
		
		delete V; 
		delete [] Lambda; 
		delete [] Vector; 	
		
	}
	
    //Farfield P coef
	else if( fun == 61) 
	{
		
		int  ind, Ndofs;
		int Nx, Ny; 
		double *coefx, *coefy; 
		double omega, lambda , mu ,beta; 
		double *pR,*pI; 
		double *LambdaR, *LambdaI; 
		
		ind = (int)mxGetScalar(prhs[2]);
		
		Nx = (int)mxGetScalar(prhs[3]);
		Ny = (int)mxGetScalar(prhs[4]);
		
		coefx = mxGetPr(prhs[5]); 
		coefy = mxGetPr(prhs[6]); 
		
		omega = mxGetScalar(prhs[7]);
		lambda = mxGetScalar(prhs[8]);
		mu = mxGetScalar(prhs[9]);
		beta = mxGetScalar(prhs[10]);
		Ndofs = (int)mxGetScalar(prhs[11]);
		LambdaR = mxGetPr(prhs[12]); 
		LambdaI = mxGetPr(prhs[13]);


		ArcBase<double>* arc =
			new ArcParametrized<double>(ind,Nx,coefx,Ny,coefy); 
			
		std::complex<double>* Vector = new 
			std::complex<double>[2*Ndofs]; 
			
		std::complex<double>* Lambda = new 
			std::complex<double>[2*Ndofs]; 
			
		for(int ii(0); ii < 2*Ndofs; ++ii)
		{
			Vector[ii] = std::complex<double>(0.0,0.0); 
			Lambda[ii] = std::complex<double>(LambdaR[ii],
				LambdaI[ii]); 
		}
		
		RhsBase<double>* V = new ElasticityP<double>
			(Ndofs,omega,lambda,mu, beta,arc); 
		
		V->Compute(Vector); 
		
		std::complex<double> z (0.0,0.0); 
		
		for(int ii(0); ii < Ndofs; ++ii)
		{
			Vector[ii] = std::conj(Vector[ii]); 
			Vector[ii+Ndofs] = std::conj(Vector[ii+Ndofs]); 
			z += Vector[ii]*Lambda[ii] + Vector[ii+Ndofs]*Lambda[ii+Ndofs];
		}
		
		double kp = sqrt(omega*omega/(lambda+2*mu)); 
		
		std::complex<double> aux (1,1); 
		
		aux *= 0.5*z/mu/sqrt(2)/sqrt(8*pi*kp); 
		
		
		
		mxArray *p_m; 
		
		p_m=plhs[0] =mxCreateDoubleMatrix(2,1,mxCOMPLEX);
		pR=  mxGetPr(p_m);
        pI=  mxGetPi(p_m);
		
		pR[0] = aux.real()*cos(beta); 
		pI[0] = aux.imag()*cos(beta); 
		pR[1] = aux.real()*sin(beta); 
		pI[1] = aux.imag()*sin(beta); 
		
		delete V; 
		delete arc; 
		delete [] Lambda; 
		delete [] Vector; 	
		
	}
    
    //Farfield S list
	else if( fun == 7) 
	{
		
		int  ind, Ndofs;
		double omega, lambda , mu ,beta; 
		double *pR,*pI; 
		double *LambdaR, *LambdaI; 
		
		ind = (int)mxGetScalar(prhs[2]);
		omega = mxGetScalar(prhs[3]);
		lambda = mxGetScalar(prhs[4]);
		mu = mxGetScalar(prhs[5]);
		beta = mxGetScalar(prhs[6]);
		Ndofs = (int)mxGetScalar(prhs[7]);
		LambdaR = mxGetPr(prhs[8]); 
		LambdaI = mxGetPr(prhs[9]);

		ArcBase<double>* arc =
			new ArcList<double>(ind); 
			
		std::complex<double>* Vector = new 
			std::complex<double>[2*Ndofs]; 
			
		std::complex<double>* Lambda = new 
			std::complex<double>[2*Ndofs]; 
			
		for(int ii(0); ii < 2*Ndofs; ++ii)
		{
			Vector[ii] = std::complex<double>(0.0,0.0); 
			
			Lambda[ii] = std::complex<double>(LambdaR[ii],
				LambdaI[ii]); 
		}
		
		RhsBase<double>* Vs = new ElasticityS<double>
			(Ndofs,omega,lambda,mu, beta,arc); 
		
		Vs->Compute(Vector); 

		std::complex<double>* PPS1 = new 
			std::complex<double>[2*Ndofs]; 
		std::complex<double>* PPS2 = new 
			std::complex<double>[2*Ndofs]; 
			
		for(int ii(0); ii < 2*Ndofs; ++ii)
		{
			PPS1[ii] = std::complex<double>(0.0,0.0); 
			PPS2[ii] = std::complex<double>(0.0,0.0); 
		}		
			
		double sbeta = sin(beta); 
		double cbeta = cos(beta); 
		
		if(abs(sbeta) > 1e-14) 
		{
			for(int ii(0); ii < Ndofs; ++ii)
			{
				PPS1[ii] = std::conj(Vector[ii]/(-sbeta)); 
				PPS2[ii] = std::conj(Vector[ii] *(-cbeta/sbeta)); 
			}
		}
		
		if( abs(cbeta) > 1e-14)
		{
			for(int ii(Ndofs); ii < 2*Ndofs; ++ii)
			{
				PPS1[ii] = std::conj(Vector[ii]/(cbeta)); 
				PPS2[ii] = std::conj(Vector[ii] *(sbeta/cbeta)); 
			}
		}
		
		std::complex<double> z2 (0.0,0.0); 
		
		std::complex<double> z1[2];
		
		z1[0] = std::complex<double>(0.0,0.0); 
		z1[1] = std::complex<double>(0.0,0.0); 
		
		
		for(int ii(0); ii < Ndofs; ++ii)
		{
			z1[0] += PPS1[ii]*Lambda[ii];
			z1[1] += PPS1[ii+Ndofs]*Lambda[ii+Ndofs];			
			z2 += PPS2[ii]*Lambda[ii] + PPS2[ii+Ndofs]*Lambda[ii+Ndofs];
		}
		
		double ks = sqrt(omega*omega/mu); 
		
		std::complex<double> aux (1,1); 
		
		aux *= 1.0/mu/sqrt(2)/sqrt(8*pi*ks); 	
		z1[0] *= aux; 
		z1[1] *= aux; 
		
		aux *= z2; 
		
		mxArray *p_m; 
		
		p_m=plhs[0] =mxCreateDoubleMatrix(2,1,mxCOMPLEX);
		pR=  mxGetPr(p_m);
        pI=  mxGetPi(p_m);
		
		pR[0] = z1[0].real()-aux.real()*cos(beta); 
		pI[0] = z1[0].imag()-aux.imag()*cos(beta); 
		pR[1] = z1[1].real()-aux.real()*sin(beta); 
		pI[1] = z1[1].imag()-aux.imag()*sin(beta); 
		
		delete Vs; 
		delete arc; 
		delete [] Lambda; 
		delete [] Vector; 	
		delete [] PPS1; 
		delete [] PPS2; 
		
	}
	
    //Farfield S coef
	else if( fun == 71) 
	{
		
		int  ind, Ndofs;
		int Nx, Ny; 
		double *coefx, *coefy; 
		double omega, lambda , mu ,beta; 
		double *pR,*pI; 
		double *LambdaR, *LambdaI; 
		
		ind = (int)mxGetScalar(prhs[2]);
		
		Nx = (int)mxGetScalar(prhs[3]);
		Ny = (int)mxGetScalar(prhs[4]);
		
		coefx = mxGetPr(prhs[5]); 
		coefy = mxGetPr(prhs[6]); 
		
		omega = mxGetScalar(prhs[7]);
		lambda = mxGetScalar(prhs[8]);
		mu = mxGetScalar(prhs[9]);
		beta = mxGetScalar(prhs[10]);
		Ndofs = (int)mxGetScalar(prhs[11]);
		LambdaR = mxGetPr(prhs[12]); 
		LambdaI = mxGetPr(prhs[13]);


		ArcBase<double>* arc =
			new ArcParametrized<double>(ind,Nx,coefx,Ny,coefy); 
			
			
		std::complex<double>* Vector = new 
			std::complex<double>[2*Ndofs]; 
			
		std::complex<double>* Lambda = new 
			std::complex<double>[2*Ndofs]; 
			
		for(int ii(0); ii < 2*Ndofs; ++ii)
		{
			Vector[ii] = std::complex<double>(0.0,0.0); 
			
			Lambda[ii] = std::complex<double>(LambdaR[ii],
				LambdaI[ii]); 
		}
		
		RhsBase<double>* Vs = new ElasticityS<double>
			(Ndofs,omega,lambda,mu, beta,arc); 
		
		Vs->Compute(Vector); 

		std::complex<double>* PPS1 = new 
			std::complex<double>[2*Ndofs]; 
		std::complex<double>* PPS2 = new 
			std::complex<double>[2*Ndofs]; 
			
		for(int ii(0); ii < 2*Ndofs; ++ii)
		{
			PPS1[ii] = std::complex<double>(0.0,0.0); 
			PPS2[ii] = std::complex<double>(0.0,0.0); 
		}		
			
		double sbeta = sin(beta); 
		double cbeta = cos(beta); 
		
		if(abs(sbeta) > 1e-14) 
		{
			for(int ii(0); ii < Ndofs; ++ii)
			{
				PPS1[ii] = std::conj(Vector[ii]/(-sbeta)); 
				PPS2[ii] = std::conj(Vector[ii] *(-cbeta/sbeta)); 
			}
		}
		
		if( abs(cbeta) > 1e-14)
		{
			for(int ii(Ndofs); ii < 2*Ndofs; ++ii)
			{
				PPS1[ii] = std::conj(Vector[ii]/(cbeta)); 
				PPS2[ii] = std::conj(Vector[ii] *(sbeta/cbeta)); 
			}
		}
		
		std::complex<double> z2 (0.0,0.0); 
		
		std::complex<double> z1[2];
		
		z1[0] = std::complex<double>(0.0,0.0); 
		z1[1] = std::complex<double>(0.0,0.0); 
		
		
		for(int ii(0); ii < Ndofs; ++ii)
		{
			z1[0] += PPS1[ii]*Lambda[ii];
			z1[1] += PPS1[ii+Ndofs]*Lambda[ii+Ndofs];			
			z2 += PPS2[ii]*Lambda[ii] + PPS2[ii+Ndofs]*Lambda[ii+Ndofs];
		}
		
		double ks = sqrt(omega*omega/mu); 
		
		std::complex<double> aux (1,1); 
		
		aux *= 1.0/mu/sqrt(2)/sqrt(8*pi*ks); 	
		z1[0] *= aux; 
		z1[1] *= aux; 
		
		aux *= z2; 
		
		mxArray *p_m; 
		
		p_m=plhs[0] =mxCreateDoubleMatrix(2,1,mxCOMPLEX);
		pR=  mxGetPr(p_m);
        pI=  mxGetPi(p_m);
		
		pR[0] = z1[0].real()-aux.real()*cos(beta); 
		pI[0] = z1[0].imag()-aux.imag()*cos(beta); 
		pR[1] = z1[1].real()-aux.real()*sin(beta); 
		pI[1] = z1[1].imag()-aux.imag()*sin(beta); 
		
		delete Vs; 
		delete arc; 
		delete [] Lambda; 
		delete [] Vector; 	
		delete [] PPS1; 
		delete [] PPS2; 
		
	}
    
    //get geo
    else if ( fun == 100) 
    {
        int  ind;
		int Nx, Ny,Nt; 
		double *coefx, *coefy; 
		double *pR,*pI,*ts; 
		
		ind = (int)mxGetScalar(prhs[2]);
		
		Nx = (int)mxGetScalar(prhs[3]);
		Ny = (int)mxGetScalar(prhs[4]);
		
		coefx = mxGetPr(prhs[5]); 
		coefy = mxGetPr(prhs[6]); 
        
        Nt = (int)mxGetScalar(prhs[7]);
        ts = mxGetPr(prhs[8]); 
		
        double* Xs = new double[Nt]; 
        
        double* Ys = new double[Nt];       
        
        ArcBase<double>* arc =
			new ArcParametrized<double>(ind,Nx,coefx,Ny,coefy); 
        
        for(int ii(0); ii <Nt; ++ii) 
        {
            arc->GetPoint(ts[ii],Xs+ii,Ys+ii); 
        }
        
        mxArray *p_m; 
		
		p_m=plhs[0] =mxCreateDoubleMatrix(Nt,2,mxREAL);
		pR=  mxGetPr(p_m);
       
        for(int ii(0); ii < Nt; ++ii)
		{
			pR[ii] = Xs[ii]; 
            pR[ii+Nt] = Ys[ii]; 
		}        
        
        delete [] Xs; 
        delete [] Ys; 
        delete arc; 
        
		
    }
    else if(fun ==200)
    {
        int indTrial, indTest, Ndofs,nthreads; 
		int NxTrial,NyTrial, NxTest, NyTest; 
		double *coefxTrial, *coefyTrial, *coefxTest,*coefyTest; 
		double omega, lambda,mu; 
		double *pR,*pI; 
		
		indTrial = (int)mxGetScalar(prhs[2]);
		
		NxTrial = (int) mxGetScalar(prhs[3]); 
		NyTrial = (int) mxGetScalar(prhs[4]); 
		
		coefxTrial = mxGetPr(prhs[5]); 
		coefyTrial = mxGetPr(prhs[6]); 
		
		indTest = (int)mxGetScalar(prhs[7]);
		
		NxTest = (int) mxGetScalar(prhs[8]);
		NyTest = (int) mxGetScalar(prhs[9]);
		
		coefxTest = mxGetPr(prhs[10]); 
		coefyTest = mxGetPr(prhs[11]); 
		
		omega = mxGetScalar(prhs[12]);
		lambda = mxGetScalar(prhs[13]);
		mu = mxGetScalar(prhs[14]);
		Ndofs = (int)mxGetScalar(prhs[15]);
		nthreads = (int)mxGetScalar(prhs[16]);
				
		ArcBase<double>* arcTrial =
			new ArcParametrized<double>(indTrial,
				NxTrial, coefxTrial,
				NyTrial,coefyTrial); 
							
		ArcBase<double>* arcTest =
			new ArcParametrized<double>(indTest,
				NxTest,coefxTest,
				NyTest,coefyTest); 
        
        double t = mxGetScalar(prhs[17]);
        
        int Ns = (int) mxGetScalar(prhs[18]);
        
        double *ss = mxGetPr(prhs[19]); 
        
//         GFElasticity2<double>* m_GreenFun = new GFElasticity2<double>
// 				(0,0,arcTrial,arcTest,
// 				lambda,mu,omega); 
        
        
       GFElasticity1<double>* m_GreenFun = new GFElasticity1<double>
				(arcTrial,arcTest,
				lambda,mu,omega); 
        
        std::complex<double>* evals = new std::complex<double>[Ns]; 
        
        for(int ii(0); ii < Ns; ++ii)
        {
            double dist  = m_GreenFun->ComputeDist(t,ss[ii]); 
            
            if( dist > 1e-14)
            {
                std::complex<double> J = m_GreenFun->Jfun(t,ss[ii],dist); 
                
                evals[ii] = m_GreenFun->
                        FundamentalSolution( t, ss[ii],  dist)+
                        0.5/pi*log(abs(t-ss[ii]))*J;
                
            }
            else
            {
                evals[ii] = m_GreenFun->RegPartLimit(t); 
            
                std::cout<<"here "<<evals[ii]<<std::endl; 
            }           
                        
        }
        
        mxArray *p_m; 
		
		p_m=plhs[0] =mxCreateDoubleMatrix(Ns,1,mxCOMPLEX);
		pR=  mxGetPr(p_m);
        pI = mxGetPi(p_m); 
       
        for(int ii(0); ii < Ns; ++ii)
		{
			pR[ii] = evals[ii].real(); 
            pI[ii] = evals[ii].imag(); 
		}        
               
        
        delete [] evals;        
        
        delete m_GreenFun;
        
    }
    
}
	
	
