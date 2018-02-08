#include "mex.h"
#include "Adom.h"
#include "geo.h"
#include "kers.h"
#include <omp.h>
#include <complex>
#include "IntegrationRules.h" 
#include "greenfunction.h"
#include "prods.h"


using namespace std;
typedef complex<double> dcomp;    
void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int fun; 
    fun=(int)mxGetScalar(prhs[0]);
         
    if (fun ==11)
    {
        //trial function 
        
        int n = (int) mxGetScalar(prhs[1]);
        
        double *pR; 
        
        mxArray *p_m;
        
        double t=mxGetScalar(prhs[2]);
        
        p_m=plhs[0] =mxCreateDoubleMatrix(1,1,mxREAL);
        
        pR=  mxGetPr(p_m);
        
        pR[0] =  TrialFunction(n,t);
        
    }
   
        
    else if(fun ==21)  
    {
        //perturbed geo
        
//         double t; 
//         int dom_op,NPer;
//         double *an,*bn; 
//         double *pR;
//         mxArray *p_m;
// 
//         //get scalar variables.     
//         t=mxGetScalar(prhs[1]);
//         dom_op=(int)mxGetScalar(prhs[2]);
//         NPer=(int)mxGetScalar(prhs[3]);
//         
//         //get vector varibales.
//         an = mxGetPr(prhs[4]);
//         bn = mxGetPr(prhs[5]);
// 
//        //output variable
//         p_m=plhs[0] =mxCreateDoubleMatrix(2,1,mxREAL);          
//         pR=  mxGetPr(p_m);
//         geomPer(&pR[0], &pR[1], t, dom_op,an,bn,NPer);      
//         
    }
    
//     else if(fun ==22)  
//     {
//         //perturbed normal
//         
//         double t; 
//         int dom_op,NPer;
//         double *an,*bn; 
//         double *pR;
//         mxArray *p_m;
// 
//         //get scalar variables.     
//         t=mxGetScalar(prhs[1]);
//         dom_op=(int)mxGetScalar(prhs[2]);
//         NPer=(int)mxGetScalar(prhs[3]);
//         
//         //get vector varibales.
//         an = mxGetPr(prhs[4]);
//         bn = mxGetPr(prhs[5]);
// 
//        //output variable
//         p_m=plhs[0] =mxCreateDoubleMatrix(2,1,mxREAL);          
//         pR=  mxGetPr(p_m);
//         normalPer(&pR[0], &pR[1], t, dom_op,an,bn,NPer);      
//         
//     }
//     
//      else if(fun ==23)  
//     {
//         //perturbed jacobian
//         
//         double t; 
//         int dom_op,NPer;
//         double *an,*bn; 
//         double *pR;
//         mxArray *p_m;
// 
//         //get scalar variables.     
//         t=mxGetScalar(prhs[1]);
//         dom_op=(int)mxGetScalar(prhs[2]);
//         NPer=(int)mxGetScalar(prhs[3]);
//         
//         //get vector varibales.
//         an = mxGetPr(prhs[4]);
//         bn = mxGetPr(prhs[5]);
// 
//        //output variable
//         p_m=plhs[0] =mxCreateDoubleMatrix(1,1,mxREAL);          
//         pR=  mxGetPr(p_m);
//         pR[0]=JPer(t, dom_op,an,bn,NPer);      
//         
//     } 
   

    else if(fun == 105)         
    {
        
        double k,delta;
        
        int dom_op,Nq,Nc,Dt,NT,*N,S;         

        //poiters for arrays coming from matlab. 
        double *xq,*wq,*pR,*pI,*N_m,*Y1,*Y2;
        
        mxArray *xg_m,*p_m,*wg_m;
          
        k= mxGetScalar(prhs[1]);        
             
        N_m= mxGetPr(prhs[2]);     
        
        Nc=(int)mxGetScalar(prhs[3]);
        
        Nq=(int)mxGetScalar(prhs[4]);
        
        dom_op=(int)mxGetScalar(prhs[5]);
        
        Dt=(int)mxGetScalar(prhs[6]);
        
        xq = mxGetPr(prhs[7]); 
        
        wq = mxGetPr(prhs[8]); 
        
        Y1 = mxGetPr(prhs[9]); 
        
        Y2 = mxGetPr(prhs[10]); 
        
        delta = mxGetScalar(prhs[11]); 
        
        S = (int)mxGetScalar(prhs[12]); 
        
        N=new int[Dt];         
        
        for (int t=0; t <Dt; t++)
            N[t]=(int)N_m[t];         

        //generate matrix
        NT=0;        
        for (int t=0; t<Dt; t++){
            NT+=2*N[t]+1; 
        }
                        
        dcomp** A=new dcomp*[NT];
        
        for(int j=0; j<NT;j++){ 
            
            A[j]= new dcomp[NT]; 
        }   
       
        VdomNs(A,k,N,Nc,Nq,xq,wq,dom_op,Dt,NT,Y1,Y2,delta,S);           
        
        //output varaibles, and whrite the matrix.       
        p_m=plhs[0] =mxCreateDoubleMatrix(NT,NT,mxCOMPLEX);
        
        pR=  mxGetPr(p_m);
        
        pI=  mxGetPi(p_m);
        
        for (int j=0; j<NT;j++){    
            
            for (int l=0; l<NT;l++){   
                
               pR[l*NT+j]=real(A[j][l]);
               
               pI[l*NT+j]=imag(A[j][l]);
             }
        } 

       //free memory and end     
        
        delete A,N;      
    
    }
    
    else if(fun == 106)         
    {
    
          //right hand side 
       //inpuut variables

        int NT,Dt,S;
        bool trace; 
        double *pR,*pI,*N_m,alpha,kext,*Y1,*Y2,delta;
        mxArray *p_m;
        int *N; 
        //get scalar variables.  
        

        
        N_m=mxGetPr(prhs[1]);
        alpha =mxGetScalar(prhs[2]);
        kext=mxGetScalar(prhs[3]);
        Dt=mxGetScalar(prhs[4]);
        
        Y1 = mxGetPr(prhs[5]); 
        
        Y2 = mxGetPr(prhs[6]); 
        
        delta = mxGetScalar(prhs[7]); 
        
        S = (int)mxGetScalar(prhs[8]); 
        
        N=new int[Dt]; 
        for (int t=0; t <Dt; t++){
            N[t]=(int)N_m[t]; } 
       

        //generete vector
        NT=0; 
        for (int t=0; t<Dt; t++){
            NT+=2*N[t]+1; 
        }       
        dcomp* b=new dcomp[NT];
        for (int j=0; j<NT; j++){
            b[j]=dcomp(0,0);
        }
        bdomPlaneNs(b,N,alpha,kext,Dt,NT,Y1,Y2,delta,S);   

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
    
      else if(fun == 107)         
    {
        
        double k;
        
        int Nq,Nc,segmentTrial,segmentTest,S;         

        //poiters for arrays coming from matlab. 
        double *xq,*wq,*pR,*pI,*N_m,*Y1,*Y2,delta;
        
        mxArray *xg_m,*p_m,*wg_m;
          
        k= mxGetScalar(prhs[1]);        
             
        N_m= mxGetPr(prhs[2]);     
        
        Nc=(int)mxGetScalar(prhs[3]);
        
        Nq=(int)mxGetScalar(prhs[4]);
        
        segmentTrial=(int)mxGetScalar(prhs[5]);
        
        segmentTest=(int)mxGetScalar(prhs[6]);
        
        xq = mxGetPr(prhs[7]); 
        
        wq = mxGetPr(prhs[8]); 
        
        Y1 = mxGetPr(prhs[9]); 
        
        Y2 = mxGetPr(prhs[10]); 
        
        delta = mxGetScalar(prhs[11]); 
        
        S = (int)mxGetScalar(prhs[12]); 
                
        int NT = 2*N_m[segmentTrial]+1; 
        
        dcomp** A=new dcomp*[NT];
        
        for(int j=0; j<NT;j++){ 
            
            A[j]= new dcomp[NT]; 
            
            for (int i(0); i <NT; ++i)
            {
                A[j][i]=0.0; 
            }
        }   
        
 //       printf("%d \n",NT); 
       
        VdomNsPerSegment(A,k,N_m[segmentTrial],Nc,Nq,xq,wq,
                segmentTrial,segmentTest,Y1,Y2,delta,S);     
        
//          printf("%f,%f !\n", A[0][0].real(),A[0][0].imag()); 
        
        int Mt = NT; 
        
        if( k < pow(1.0,-6.0))
        {     
            Mt = NT-1;             
        }    
        
        //output varaibles, and whrite the matrix.       
        p_m=plhs[0] =mxCreateDoubleMatrix(Mt,Mt,mxCOMPLEX);
        
        pR=  mxGetPr(p_m);
        
        pI=  mxGetPi(p_m);      
        
        if( k < pow(1.0,-6.0))
        {    
          for (int j=1; j<NT;j++){    
            
            for (int l=1; l<NT;l++){   
                
               pR[(l-1)*(NT-1)+(j-1)]=real(A[j][l]);
               
               pI[(l-1)*(NT-1)+(j-1)]=imag(A[j][l]);
               
               //mexPrintf("%d,%d | %f , %f\n",j,l, A[j][l].real(),A[j][l].imag()); 

             }
            } 
        }
        else 
        {
           for (int j=0; j<NT;j++){    
            
            for (int l=0; l<NT;l++){   
                
               pR[(l)*(NT)+(j)]=real(A[j][l]);
               
               pI[(l)*(NT)+(j)]=imag(A[j][l]);
               
               //mexPrintf("%d,%d | %f , %f\n",j,l, A[j][l].real(),A[j][l].imag()); 

             }
            } 
        }
            
        


       //free memory and end     
        
        delete A;      
    
    }   
    
    return; 
}


