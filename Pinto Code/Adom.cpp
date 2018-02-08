#include <omp.h>
#include "in_cells.h"
#include "prods.h"
#include "geo.h"
#include "tchv_cof.h"

#include "mex.h"
using namespace std;
typedef complex<double> dcomp;
const double pi= 4.0*atan(1.0);

void VdomNs(dcomp** A, double k,int* N,int Nc, int Nq,
        double*xq, double* wq,int dom_op,int Dt, int NT,double* Y1,
        double* Y2,double delta, int S){

    
    ///////////////guardamos los coeficientes de fourier para este domonio
     
    dcomp**** Opd =new dcomp***[Dt];       
    for(int t=0;t< Dt; t++){
        
        Opd[t]=new dcomp**[Dt]; 
        
        for(int s=0; s<Dt; s++){            

             Opd[t][s]=new dcomp*[Nq];
                 
             for(int j=0; j<Nq;j++)                 
                 Opd[t][s][j]=new dcomp[Nc];
         }
     }    
      
   double* xc;
    
   xc = chebpts(Nc); 

   in_cellsVNs(k,Opd,Dt,dom_op,Nq, xq, Nc,xc,Y1,Y2,delta,S);
    
   int nthread = omp_get_max_threads();
   omp_set_num_threads(nthread);
    
    int ns=0;
    int nt=0;      
    
    for(int t=0;t<Dt;t++){
        
        if(t>0) 
        {
            nt+=2*N[t-1]+1;
        }
        
        ns=0; 
        
        for(int s=0;s<Dt;s++){ 
            
            if(s>0)
            {
                ns+=2*N[s-1]+1;
            }
            
            #pragma omp parallel for collapse(2)  
                for(int l=0;l<2*N[s]+1;l++){
                 for(int m=-N[t];m<N[t]+1;m++){     
                              
                     A[l+ns][m+N[t]+nt]=pVNs(m+N[t],l ,Opd[t][s],Nq,xq,wq,k,t,s);                     

                 }}       

            

        }
    }
    

     delete [] xc;   
     xc=0;

    
    
    for(int t=0;t< Dt; t++){       
        for(int s=0; s<Dt; s++){                             
             for(int j=0; j<Nq;j++){
                 delete [] Opd[t][s][j];}                 
            delete [] Opd[t][s];}
        delete [] Opd[t];}
    delete [] Opd; 
    Opd= 0;
    

    return; 
    
    
}

void bdomPlaneNs(dcomp* b, int* N,double alpha, double k,int Dt,int NT,
        double* Y1,double* Y2,double delta,int S){
 // b inicializado en 0. 
     int Nc=2*(2*N[0]+1)+28;  
//      int* tt;
//      int sign_dom=1;
//      if(dom_op==0){ 
//          sign_dom=-1; 
//          tt=trial(Dt);
//      }else{
//          tt=trial(Dt);
//      }
     double* xc;
     xc = chebpts(Nc); 
     dcomp** Opd =new dcomp*[Dt];       
     for(int t=0;t< Dt; t++){
            Opd[t]=new dcomp[Nc];    
 
      }   

     in_cells_Ns(Opd,Dt,alpha,k,Nc,xc,Y1,Y2,delta,S);
     int ns=0;
     for (int s=0; s<Dt;s++){
         if( s>0) ns+=2*N[s-1]+1;
                
        for(int l=0;l<2*N[s]+1;l++){
            b[l+ns]=pldNs(l,Opd[s],0);               
           
         }
     }

          for(int t=0;t< Dt; t++){
            delete [] Opd[t];     
   
      }   
      delete [] Opd;

      return;    
}

void VdomNsPerSegment(dcomp** A, double k,int N,int Nc, int Nq,
        double*xq, double* wq,int segmentTrial,int segmentTest,
        double* Y1,double* Y2,double delta,int S){

    //
        
    dcomp** Opd =new dcomp*[Nq];       

    for(int j=0; j<Nq;j++)   
    {
         Opd[j]=new dcomp[Nc];

    }       
      
   double* xc;
    
   xc = chebpts(Nc); 
 
   in_cellsVNsPerSegment(k,Opd,segmentTrial,segmentTest,Nq, xq, Nc,xc,
           Y1,Y2,delta,S);
   
   int nthread = omp_get_max_threads();
   omp_set_num_threads(nthread);
   
   
   //not compresed code. 
//         #pragma omp parallel for collapse(2)  
//      for(int l=0;l<2*N+1;l++){
//       for(int m=-N;m<N+1;m++){     
// 
//          A[l][m+N]=pVNs(m+N,l ,Opd,Nq,xq,wq,k,
//                  segmentTrial,segmentTest);                     
// 
//      }}


        int l =0; 
        
         double tol = 1e-10; 
         int Ntol = 2*N+1;
         int tolFound =0; 
              
         for(int m=-N;m<N+1;m++)
         {     

             A[0][m+N]=pVNs(m+N,0 ,Opd,Nq,xq,wq,k,
                      segmentTrial,segmentTest);                    

             if ( abs(A[0][m+N]) < tol) 
             {
                 ++tolFound ;

                 if(tolFound ==2) 
                 {
                  Ntol = m+N; 
                  break; 
                 }
             }                          
             else 
             {
                 tolFound =0; 
             }

         }
         
         tolFound =0; 
         int Ntol2 =2*N+1;
         
         for(int l=0;l<2*N+1;l++)
         {     

             A[l][0]=pVNs(0,l ,Opd,Nq,xq,wq,k,
                      segmentTrial,segmentTest);                    

             if ( abs(A[l][0]) < tol) 
             {
                 ++tolFound ;

                 if(tolFound ==2) 
                 {
                  Ntol2 = l; 
                  break; 
                 }
             }                          
             else 
             {
                 tolFound =0; 
             }

         }
               
         int Nmax1 = min(Ntol,2*N+1);
         
         int Nmax2 = min(Ntol2,2*N+1);
                  
         //terms arround the diagonal. 
         if (segmentTrial == segmentTest) 
         {
             int count = 0; 
             
             int l0 = 0; 
             
             int l1 = 0; 
             
             int n = N;
             
             if( Nmax1  <2*N+1)
             {
                 for (int l(0); l<2*N+1;++l) 
                 {
                     A[l][n]=
                            pVNs(n,l ,Opd,Nq,xq,wq,k,
                            segmentTrial,segmentTest); 
                     
                     if((count == 0)&&(abs(A[l][n])>tol))
                     {
                         count = 1; 
                         l0 = l; 
                     }
                     else if ((count > 0))
                     {
                        if((count == 1)&&(abs(A[l][n])<tol)) 
                        {
                            count = 2;
                            
                            continue; 
                        }
                        else if((count == 2)&&(abs(A[l][n])<tol))
                        {
                            l1 = l; 
                            
                            break; 
                        }
                        else 
                        {

                            count =1; 
                        }
                     }
                     else if(abs(A[l][n])<tol)
                     {
                         A[l][n] =0.0; 
                     }                  
                 
                 }
                 
             }
                 
             if(l1 == 0)
             {
                 l1 = 2*N+1; 
             }

             int nc = max(n-l0,l1+1-n); 

             #pragma omp parallel for  
             for (int m=1; m<2*N+1; ++m )
             {
                 int lmin = max(0,m-nc); 

                 if(m < Nmax1) 
                 {
                     lmin=0; 

                 }                 

                 int lmax = max(Nmax2,nc+m);

                 lmax = min(2*N+1,lmax);

                // printf("%d , %d | %d \n",m,lmin,lmax);

                 for( int l=lmin; l<lmax;++l)
                 {
                        A[l][m]=
                        pVNs(m,l ,Opd,Nq,xq,wq,k,
                        segmentTrial,segmentTest); 
                 }

             }

             
             
                     
                     
         
         }
         else if((Nmax2>0)&&(Nmax1>0))
         {
            #pragma omp parallel for collapse(2)  
            for(int l=0;l<Nmax2;l++){
                for(int m=0;m<Nmax1;m++){    

                    A[l][m]=
                            pVNs(m,l ,Opd,Nq,xq,wq,k,
                            segmentTrial,segmentTest);                     

                   }}
             
            
         }    


     delete [] xc;   
     xc=0;

    
    
                         
    for(int j=0; j<Nq;j++){
        delete [] Opd[j];} 
    delete [] Opd; 
    Opd= 0;
    

    return; 
    
    
}