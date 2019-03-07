#include "mex.h"
#include "Adom.h"
#include "geo.h"
#include "kers.h"
#include <omp.h>
#include <complex>
#include "IntegrationRules.h" 
#include "greenfunction.h"
#include "in_cells.h"
#include "prods.h"


using namespace std;
typedef complex<double> dcomp;

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //mexPrintf("imprime \n");   
    
    //Wich mex function 
    int fun; 
    //fun=1 ==> A
    //fun=2 ==> X
    //fun=3 ==> b
    //fun=4 ==> j
    fun=(int)mxGetScalar(prhs[0]);
   // mexPrintf("%d",fun);
           


    



        
    if(fun==5){
        //geometry
        //input variables 
        double t; 
        double *bgeo,*pR,delta;
        int Nb;
        mxArray *p_m;

        //get scalar variables.     
        t=mxGetScalar(prhs[1]);
        bgeo=mxGetPr(prhs[2]);
        Nb=(int)mxGetScalar(prhs[3]);
        
        delta = mxGetScalar(prhs[4]); 
        

       //output variable
        p_m=plhs[0] =mxCreateDoubleMatrix(2,1,mxREAL);          
        pR=  mxGetPr(p_m);
        geom(&pR[0], &pR[1], t, bgeo,Nb,delta);
    

    }
    
    else if(fun == 106)         
    {
    
          //right hand side 
       //inpuut variables

        int dom_op,NT,Dt,Dt_0;
        bool trace; 
        double *pR,*pI,*N_m,k,angle,*bgeo, delta;
        mxArray *p_m;
        int *N,Nb; 
        //get scalar variables.  
        

        k = mxGetScalar(prhs[1]);
        angle = mxGetScalar(prhs[2]);
        N_m=mxGetPr(prhs[3]);
        dom_op=(int)mxGetScalar(prhs[4]);
        Dt=(int)mxGetScalar(prhs[5]);
        Dt_0=(int)mxGetScalar(prhs[6]);
        
        bgeo = mxGetPr(prhs[7]); 
        
        Nb = (int)mxGetScalar(prhs[8]);
        
        delta = mxGetScalar(prhs[9]); 
        
        N=new int[Dt]; 
        for (int t=0; t <Dt; t++){
            N[t]=(int)N_m[t]; } 
       

        //generete vector
        NT=0; 
        for (int t=0; t<Dt; t++){
            //NT+=2*N[t]+1; 
            NT += N[t];
        }       
        dcomp* b=new dcomp[NT];
        for (int j=0; j<NT; j++){
            b[j]=dcomp(0,0);
        }
        bdomPlaneNs(b,k,angle,N,dom_op,Dt,Dt_0,NT,bgeo,Nb,delta);    

        //output variable
        p_m=plhs[0] =mxCreateDoubleMatrix(NT,1,mxCOMPLEX);
        pR=  mxGetPr(p_m);
        pI=  mxGetPi(p_m);
        for (int j=0; j<NT;j++){          
           pR[j]=real(b[j]);
           pI[j]=imag(b[j]);
        }
    
        //free memory and end 
        delete b;  
        delete N;
        
    }  
    
    else if(fun == 206)         
    {
    
          //right hand side 
       //inpuut variables

        int dom_op,NT,Dt,Dt_0;
        bool trace; 
        double *pR,*pI,*N_m,k,angle,*bgeo, delta;
        mxArray *p_m;
        int *N,Nb,p,component; 
        //get scalar variables.  
        

        k = mxGetScalar(prhs[1]);
        angle = mxGetScalar(prhs[2]);
        N_m=mxGetPr(prhs[3]);
        dom_op=(int)mxGetScalar(prhs[4]);
        Dt=(int)mxGetScalar(prhs[5]);
        Dt_0=(int)mxGetScalar(prhs[6]);
        
        bgeo = mxGetPr(prhs[7]); 
        
        Nb = (int)mxGetScalar(prhs[8]);
        
        delta = mxGetScalar(prhs[9]); 
        
        p = (int)mxGetScalar(prhs[10]);
        
        component = (int)mxGetScalar(prhs[11]);
        
        N=new int[Dt]; 
        for (int t=0; t <Dt; t++){
            N[t]=(int)N_m[t]; } 
       

        //generete vector
        NT=0; 
        for (int t=0; t<Dt; t++){
            //NT+=2*N[t]+1; 
            NT += N[t];
        }       
        dcomp* b=new dcomp[NT];
        for (int j=0; j<NT; j++){
            b[j]=dcomp(0,0);
        }
        
        bdomDevPlaneNs( b, k, angle,  N, dom_op,
         Dt, Dt_0, NT, bgeo, Nb,  delta, p, component);     

        //output variable
        p_m=plhs[0] =mxCreateDoubleMatrix(NT,1,mxCOMPLEX);
        pR=  mxGetPr(p_m);
        pI=  mxGetPi(p_m);
        for (int j=0; j<NT;j++){          
           pR[j]=real(b[j]);
           pI[j]=imag(b[j]);
        }
    
        //free memory and end 
        delete b;  
        delete N;
        
    }  
    else if(fun ==207) 
    {
        
        double k, *bgeo, delta, *y,*pR,*pI;
        int p, component,N,Nb,Nobs; 
        
        mxArray *p_m;
        
        k = mxGetScalar(prhs[1]);
        
        N = (int)mxGetScalar(prhs[2]);
        
        bgeo = mxGetPr(prhs[3]); 
        
        Nb = (int)mxGetScalar(prhs[4]);
        
        delta = mxGetScalar(prhs[5]);
        
        p = (int)mxGetScalar(prhs[6]);
        
        component = (int)mxGetScalar(prhs[7]);
        
        y = mxGetPr(prhs[8]); 
              
        Nobs  = (int)mxGetScalar(prhs[9]);
                
        int Nc=2*(2*N+1)+28;  
        
        dcomp* A = new dcomp[Nc*Nobs];
        
        FarOneComputation(A,k,N,
        bgeo,Nb,delta,p,component,y,Nobs); 
        
        p_m=plhs[0] =mxCreateDoubleMatrix(Nobs,N,mxCOMPLEX);
        
        pR=  mxGetPr(p_m);
        
        pI=  mxGetPi(p_m);    
        
        for( int s(0); s< Nobs ; s++)
        {
            for (int n(0); n< N; n++)
            {
               pR[s+n*Nobs]=real(A[s*Nc+n]);
               
               pI[s+n*Nobs]=imag(A[s*Nc+n]);    
                               
            }
            
        }
        
        
 
        delete[] A;  

                
    }
    else if(fun == 109)         
    {
        
        double k,tol;
        
        int Nc,Nt,segmentTrial,segmentTest,N,Nb;         

        //poiters for arrays coming from matlab. 
        double *pR,*pI,*N_m,*bgeo, delta;
        
        mxArray *p_m;
          
        k= mxGetScalar(prhs[1]);        
             
        N= (int)mxGetScalar(prhs[2]);     
        
        Nc=(int)mxGetScalar(prhs[3]);
        
        Nt=(int)mxGetScalar(prhs[4]);
        
        segmentTrial=(int)mxGetScalar(prhs[5]);
        
        segmentTest=(int)mxGetScalar(prhs[6]);
        
        tol = mxGetScalar(prhs[7]);
        
        bgeo = mxGetPr(prhs[8]); 
        
        Nb = (int)mxGetScalar(prhs[9]);
        
        delta = mxGetScalar(prhs[10]); 
               
        dcomp** A=new dcomp*[N];
        
        for(int j=0; j<N;j++){ 
            
            A[j]= new dcomp[N]; 
            
            for (int i(0); i <N; ++i)
            {
                A[j][i]=0.0; 
            }
        }   
        
 //       printf("%d \n",NT); 
        
        
        if(tol > 0)
        {

            VdomNsPerSegment3(A,k,N,Nc,Nt,
                    segmentTrial,segmentTest,tol,bgeo,Nb,delta);
        }
        else 
        {
            
            VdomNsPerSegment2(A,k,N,Nc,Nt,
                    segmentTrial,segmentTest,bgeo,Nb,delta);
        }   
        

        int Mt = N; 
        
        if( k < pow(1.0,-6.0))
        {     
            Mt = N-1;             
        }    
        
        //output varaibles, and whrite the matrix.       
        p_m=plhs[0] =mxCreateDoubleMatrix(Mt,Mt,mxCOMPLEX);
        
        pR=  mxGetPr(p_m);
        
        pI=  mxGetPi(p_m);      
        
        if( k < pow(1.0,-6.0))
        {    
          for (int j=1; j<N;j++){    
            
            for (int l=1; l<N;l++){   
                
               pR[(l-1)*(N-1)+(j-1)]=real(A[j][l]);
               
               pI[(l-1)*(N-1)+(j-1)]=imag(A[j][l]);          
              

             }
            } 
        }
        else 
        {
           for (int j=0; j<N;j++){    
            
            for (int l=0; l<N;l++){   
                
               pR[(l)*(N)+(j)]=real(A[j][l]);
               
               pI[(l)*(N)+(j)]=imag(A[j][l]);               
              
             }
            } 
        }
            
        for(int j=0; j<N;j++){ 
            
             delete[] A[j]; 
        }              
        
        delete[] A;      
    }
     else if(fun == 209)         
    {
        
        double k,tol;
        
        int Nc,Nt,segmentTrial,segmentTest,N,Nb,p,component;         

        //poiters for arrays coming from matlab. 
        double *pR,*pI,*N_m,*bgeo, delta;
        
        mxArray *p_m;
          
        k= mxGetScalar(prhs[1]);        
             
        N= (int)mxGetScalar(prhs[2]);     
        
        Nc=(int)mxGetScalar(prhs[3]);
        
        Nt=(int)mxGetScalar(prhs[4]);
        
        segmentTrial=(int)mxGetScalar(prhs[5]);
        
        segmentTest=(int)mxGetScalar(prhs[6]);
        
        tol = mxGetScalar(prhs[7]);
        
        bgeo = mxGetPr(prhs[8]); 
        
        Nb = (int)mxGetScalar(prhs[9]);
        
        delta = mxGetScalar(prhs[10]); 
        
        p = (int)mxGetScalar(prhs[11]);
        
        component = (int)mxGetScalar(prhs[12]);
               
        dcomp** A=new dcomp*[N];
        
        for(int j=0; j<N;j++){ 
            
            A[j]= new dcomp[N]; 
            
            for (int i(0); i <N; ++i)
            {
                A[j][i]=0.0; 
            }
        }   
        
 //       printf("%d \n",NT); 
                 
        
        VDerivedNsPerSegment2(A,k,N,Nc,Nt,
                    segmentTrial,segmentTest,bgeo,Nb,delta,p,component);
        
        int Mt = N; 
        
        if( k < pow(1.0,-6.0))
        {     
            Mt = N-1;             
        }    
        
        //output varaibles, and whrite the matrix.       
        p_m=plhs[0] =mxCreateDoubleMatrix(Mt,Mt,mxCOMPLEX);
        
        pR=  mxGetPr(p_m);
        
        pI=  mxGetPi(p_m);      
        
        if( k < pow(1.0,-6.0))
        {    
          for (int j=1; j<N;j++){    
            
            for (int l=1; l<N;l++){   
                
               pR[(l-1)*(N-1)+(j-1)]=real(A[j][l]);
               
               pI[(l-1)*(N-1)+(j-1)]=imag(A[j][l]);          
              

             }
            } 
        }
        else 
        {
           for (int j=0; j<N;j++){    
            
            for (int l=0; l<N;l++){   
                
               pR[(l)*(N)+(j)]=real(A[j][l]);
               
               pI[(l)*(N)+(j)]=imag(A[j][l]);               
              
             }
            } 
        }
            
        for(int j=0; j<N;j++){ 
            
             delete[] A[j]; 
        }              
        
        delete[] A;      
    }
  
    else if(fun == 112)         
    {
        
        double k,tol,wq,wq0,Nb,delta;
        
        int Nq,Nq0,segmentTrial,segmentTest,N,Lmax;    
        
        

        //poiters for arrays coming from matlab. 
        double *pR,*pI,*N_m,*xq,*xq0,*bgeo;
        
        mxArray *p_m;
          
        k= mxGetScalar(prhs[1]);        
             
        N= (int)mxGetScalar(prhs[2]);     
        
        segmentTrial=(int)mxGetScalar(prhs[3]);
        
        segmentTest=(int)mxGetScalar(prhs[4]);
        
        Nq = (int)mxGetScalar(prhs[5]);
        
        xq = mxGetPr(prhs[6]); 
        
        wq = mxGetScalar(prhs[7]); 
        
        Nq0 = (int)mxGetScalar(prhs[8]);
        
        xq0 = mxGetPr(prhs[9]); 
        
        wq0 = mxGetScalar(prhs[10]); 
        
        tol = mxGetScalar(prhs[11]);
        
        Lmax = (int)mxGetScalar(prhs[12]);
        
        bgeo = mxGetPr(prhs[13]); 
        
        Nb = (int)mxGetScalar(prhs[14]);
        
        delta = mxGetScalar(prhs[15]); 
               
        dcomp** A=new dcomp*[N];
        
        for(int j=0; j<N;j++){ 
            
            A[j]= new dcomp[N]; 
            
            for (int i(0); i <N; ++i)
            {
                A[j][i]=0.0; 
            }
        }   
        
        VdomNsPerSegment5(A, k,N, Nq, 
         xq, wq, Nq0,  xq0,  wq0,
         segmentTrial, segmentTest,
         tol,  Lmax, bgeo,Nb,delta);

        int Mt = N; 
        
        if( k < pow(1.0,-6.0))
        {     
            Mt = N-1;             
        }    
        
        //output varaibles, and whrite the matrix.       
        p_m=plhs[0] =mxCreateDoubleMatrix(Mt,Mt,mxCOMPLEX);
        
        pR=  mxGetPr(p_m);
        
        pI=  mxGetPi(p_m);      
        
        if( k < pow(1.0,-6.0))
        {    
          for (int j=1; j<N;j++){    
            
            for (int l=1; l<N;l++){   
                
               pR[(l-1)*(N-1)+(j-1)]=real(A[j][l]);
               
               pI[(l-1)*(N-1)+(j-1)]=imag(A[j][l]);          
              

             }
            } 
        }
        else 
        {
           for (int j=0; j<N;j++){    
            
            for (int l=0; l<N;l++){   
                
               pR[(l)*(N)+(j)]=real(A[j][l]);
               
               pI[(l)*(N)+(j)]=imag(A[j][l]);               
              
             }
            } 
        }
            
        for(int j=0; j<N;j++){ 
            
             delete[] A[j]; 
        }              
        
        delete[] A;  
    
    }
    
    return; 
}

