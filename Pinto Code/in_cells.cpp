#include "kers.h"
#include "prods.h"
#include "geo.h"
#include "tchv_cof.h"
#include "in_cells.h"
#include <omp.h>
#include <fftw3.h>
#include "greenfunction.h"

#include <iostream>
#include <fstream>

#include "mex.h"

using namespace std;
typedef complex<double> dcomp;



void in_cells_Ns(dcomp** Opd,double k,double angle,int Dt,int dom_op,
        int Nc,double* xc, double* bgeo, int Nb, double delta){    
    fftw_plan p;
    int nthread = omp_get_max_threads();
    omp_set_num_threads(nthread);       


         #pragma omp parallel for collapse(2)  
        for(int t=0; t<Dt;t++){         
            for(int l=0;l<Nc;l++){          
                        Opd[t][l]=rhsNs(k,angle,t,xc[l],bgeo,Nb,delta);   
    
            }
        }
 

          
    
   p=plan_0(Nc, Opd[0]);   
   #pragma omp parallel for collapse(1)   
   for(int t=0; t<Dt;t++){                      
            Tchv_cof(Nc,Opd[t],p);                            

   }
   

   fftw_destroy_plan(p);
   fftw_cleanup();	  
   return;

}


void in_cells_NsDev(dcomp** Opd,double k,double angle,int Dt,int dom_op,
        int Nc,double* xc, double* bgeo, int Nb, double delta,
        int p, int component){    
    fftw_plan pp;
    int nthread = omp_get_max_threads();
    omp_set_num_threads(nthread);       


         #pragma omp parallel for collapse(2)  
        for(int t=0; t<Dt;t++){         
            for(int l=0;l<Nc;l++){          
                        Opd[t][l]=DevrhsNs(k,angle,t,xc[l],bgeo,Nb,delta
                                ,p,component);   
    
            }
        }
 

          
    
   pp=plan_0(Nc, Opd[0]);   
   #pragma omp parallel for collapse(1)   
   for(int t=0; t<Dt;t++){                      
            Tchv_cof(Nc,Opd[t],pp);                            

   }
   

   fftw_destroy_plan(pp);
   fftw_cleanup();	  
   return;

}

// void in_cellsBcofs(int N,dcomp*** Opd,int segmentTrial,
//          int Nc,double* xc,int order){        
//    
//     fftw_plan p;
//     int nthread = omp_get_max_threads();
//     omp_set_num_threads(nthread);
// 
//     #pragma omp parallel
//         {    
// 
//          #pragma omp for collapse(3)    
//                    for (int ii=0; ii<=order; ii++)
//                    {
//                     for(int m=0; m<N;m++){                    
//                         for(int l=0;l<Nc;l++){ 
//                             
// 
//                             //falta el orden 
//                             Opd[ii][m][l]=chebT(m,xc[l])*
//                                     pow(J(xc[l],segmentTrial,0),2*ii);          
//                             
// 
//                         }
//                     }
//                 }
//         }
//     
// 
//       p=plan_0(Nc, Opd[0][0]  );             
//      #pragma omp parallel for collapse(2)   
//       for (int ii=0; ii<=order; ii++){
//              for(int m=0; m<N;m++){
//                  Tchv_cof(Nc,Opd[ii][m],p);                  
//                 
//     }
//       }
//       
//     fftw_destroy_plan(p);
//     fftw_cleanup();	  
//     
//     
//     return;
// 
// }

void FarOneCofs(double k,dcomp* Opd,int Nc,double* xc, 
        double* bgeo, int Nb, double delta, 
        int p, int component,
        double *y, int Nobs)
        
{
    
    fftw_plan pp;
    int nthread = omp_get_max_threads();
    omp_set_num_threads(nthread);
    
    //Opd de tamño Nc*Nc

    #pragma omp parallel
        {
                 
         #pragma omp for collapse(2)       
                    for(int s=0; s<Nobs;s++){                    
                        for(int n=0;n<Nc;n++){ 
                            
                            Opd[s*Nc+n]=FarOne(k,xc[n],bgeo,Nb,delta,p,
                                    component,y[s], y[s+Nobs]);  

                             }
                        }
                    
                    
          }
    
    
     pp=plan_0(Nc, Opd  ); 

     #pragma omp parallel for collapse(1)   
            for(int s=0; s<Nobs;s++){
                 Tchv_cof(Nc,Opd+s*Nc,pp); 
                                    
                 Opd[s*Nc+0]*= 3.141592653589793*y[s+component*Nobs];
                 
                 for(int n=1;n<Nc;n++){ 
                  
                     Opd[s*Nc+n]*= 0.5*3.141592653589793*y[s+component*Nobs];
                     
                 }
                
    }
          
   
    
    fftw_destroy_plan(pp);
    fftw_cleanup();	  
    
    
    return;          
 
    
    
}





void VNsRegularPart(double k,dcomp* Opd, int domt,
        int doms, int Nc,double* xc, double* bgeo, int Nb, double delta)
        
{
    
    fftw_plan p;
    int nthread = omp_get_max_threads();
    omp_set_num_threads(nthread);
    
    //Opd de tamño Nc*Nc

    #pragma omp parallel
        {
        
         FreeSpace gf(k,bgeo,Nb,delta); 
                 
         #pragma omp for collapse(2)       
                    for(int j=0; j<Nc;j++){                    
                        for(int l=0;l<Nc;l++){ 
                            
                            Opd[j+l*Nc]=GBNs(xc[l],xc[j],domt,
                                    doms,gf);  

                             }
                        }
          }
    
    
     p=plan_0(Nc, Opd  ); 

     #pragma omp parallel for collapse(1)   
            for(int l=0; l<Nc;l++){
                 Tchv_cof(Nc,Opd+l*Nc,p);                  
                
    }
          
    
     #pragma omp parallel
     {
        dcomp *v=new dcomp[Nc]; 

        #pragma omp for collapse(1)  
          for(int j=0; j<Nc;j++)
            {
             for(int i=0; i<Nc;i++){

                 v[i] = Opd[j+i*Nc];

             }

              Tchv_cof(Nc,v,p); 

              for(int i=0; i<Nc;i++)
             {

                 Opd[j+i*Nc] = v[i];

             }

         }
        
        delete[] v; 
     }
    
    fftw_destroy_plan(p);
    fftw_cleanup();	  
    
    
    return;          
 
    
    
}


void JpartKernel(double k,dcomp* Opd, int domt,
       int Nc,double* xc, double* bgeo, int Nb, double delta)
        
{
    
    fftw_plan p;
    int nthread = omp_get_max_threads();
    omp_set_num_threads(nthread);
    
    //Opd de tamño Nc*Nc

    #pragma omp parallel
        {
    
        
         FreeSpace gf(k,bgeo,Nb,delta); 

         #pragma omp for collapse(2)       
                    for(int j=0; j<Nc;j++){                    
                        for(int l=0;l<Nc;l++){ 
                            
                             Opd[j+l*Nc]=gf.Jpart(xc[j], xc[l],domt ); 

                             }
                        }
          }

     p=plan_0(Nc, Opd  ); 

     #pragma omp parallel for collapse(1)   
            for(int l=0; l<Nc;l++){
                 Tchv_cof(Nc,Opd+l*Nc,p);                  
                
    }

     
     #pragma omp parallel
     {
        dcomp *v=new dcomp[Nc]; 

        #pragma omp for collapse(1)  
          for(int j=0; j<Nc;j++)
            {
             for(int i=0; i<Nc;i++){

                 v[i] = Opd[j+i*Nc];

             }

              Tchv_cof(Nc,v,p); 

              for(int i=0; i<Nc;i++)
             {

                 Opd[+j+i*Nc] = v[i];

             }

         }
        
        delete[] v; 
     }
    
    fftw_destroy_plan(p);
    fftw_cleanup();	  
    
    
    return;          
 
    
    
}

// void SlTcofsComp (double k, int Nc, int Ncurves, 
//         int Nxy, double* x, double* y, dcomp* Cofs)
// {
//     fftw_plan p;
//     int nthread = omp_get_max_threads();
//     omp_set_num_threads(nthread);
//  
//     double* xc;
//     
//     xc = chebpts(Nc);
//     
//     FreeSpace gf(k); 
//     
//     for(int m=0; m<Nc; m++)
//     {
//         for(int ii=0; ii<Nxy; ii++)
//         {
//             for(int jj=0; jj<Ncurves; jj++)
//             {
//                 Cofs[m+ii*Nc+jj*Nc*Nxy] = 
//                 gf.evaluatePotential(xc[m],jj,
//               x[ii], y[ii], 0); 
//             }
//         }
//     }
//     
//      p=plan_0(Nc, Cofs  ); 
// 
//      #pragma omp parallel for collapse(2)   
//         for(int ii=0; ii<Nxy; ii++)
//         {
//             for(int jj=0; jj<Ncurves; jj++)
//             {
//                  Tchv_cof(Nc,Cofs+ii*Nc+jj*Nc*Nxy,p);  
//             }
//                 
//         }
//      
//     fftw_destroy_plan(p);
//     fftw_cleanup();	  
//     
//     
//     return;   
//     
// }
        

void in_cellsVNsPerSegment(double k,double* OpdmR, double* OpdmI,
        double* OpdlR, double* OpdlI,double* TN,
        int segmentTrial,int segmentTest,int N,
        int Nq, double*xq,
        int Nq0, double*xq0, double* bgeo, int Nb, double delta){        
   
    fftw_plan p;
    int nthread = omp_get_max_threads();
    omp_set_num_threads(nthread);


       #pragma omp parallel
        {
        
         dcomp u(0.0,0.0);
    
        
         FreeSpace gf(k,bgeo,Nb,delta); 

         #pragma omp for collapse(2)    
            for(int ll=0; ll<Nq; ++ll)
            {                                
                for (int jj=0; jj < Nq0 ; ++jj)
                {
                    
                    u=GBNs(xq[ll],xq0[jj],segmentTrial,
                                    segmentTest, gf);
                    
                    OpdmR[ll+jj*Nq] = u.real(); 
                    
                    OpdmI[ll+jj*Nq] = u.imag(); 
                    
                    u = GBNs(xq0[jj],xq[ll],segmentTrial,
                                    segmentTest,gf);
                    
                    OpdlR[ll+jj*Nq] = u.real(); 
                    
                    OpdlI[ll+jj*Nq] = u.imag(); 
                    
                 //   mexPrintf("%f \n", Opd[jj][ll].real()); 
                    
                }
            }
        }
         
//         TN = t_polynomial( Nq, N-1, xq );
        
//         for (int m=0;m<Nq;m++)
//         {
//             mexPrintf(" %f \n",TN[m]);
//         }
// // 
//          mexPrintf("\n");


}       


void VNsDerivedRegularPart(double k,dcomp* Opd, int domt,
        int doms, int Nc,double* xc, double* bgeo, int Nb, double delta, 
        int p, int component)
        
{
    
    fftw_plan pp;
    int nthread = 1;
    omp_set_num_threads(nthread);
    
    //Opd de tamño Nc*Nc

    #pragma omp parallel
        {
        
         FreeDerivativeSpace gf(k,bgeo,Nb,delta,p,component); 
                 
         #pragma omp for collapse(2)       
                    for(int j=0; j<Nc;j++){                    
                        for(int l=0;l<Nc;l++){ 
                            
                            Opd[j+l*Nc]=GBNs(xc[l],xc[j],domt,
                                    doms,gf);  

                             }
                        }
          }
    
    
     pp=plan_0(Nc, Opd  ); 

     #pragma omp parallel for collapse(1)   
            for(int l=0; l<Nc;l++){
                 Tchv_cof(Nc,Opd+l*Nc,pp);                  
                
    }
          
    
     #pragma omp parallel
     {
        dcomp *v=new dcomp[Nc]; 

        #pragma omp for collapse(1)  
          for(int j=0; j<Nc;j++)
            {
             for(int i=0; i<Nc;i++){

                 v[i] = Opd[j+i*Nc];

             }

              Tchv_cof(Nc,v,pp); 

              for(int i=0; i<Nc;i++)
             {

                 Opd[j+i*Nc] = v[i];

             }

         }
        
        delete[] v; 
     }
    
    fftw_destroy_plan(pp);
    fftw_cleanup();	  
    
    
    return;          
 
    
    
}

void JpartDevKernel(double k,dcomp* Opd, int domt,
       int Nc,double* xc, double* bgeo, int Nb, double delta, int p, int component)
        
{
    
    fftw_plan pp;
    int nthread = omp_get_max_threads();
    omp_set_num_threads(1);
    
    //Opd de tamño Nc*Nc

    #pragma omp parallel
        {
    
        
         FreeDerivativeSpace gf(k,bgeo,Nb,delta,p,component);

         #pragma omp for collapse(2)       
                    for(int j=0; j<Nc;j++){                    
                        for(int l=0;l<Nc;l++){ 
                            
                             Opd[j+l*Nc]=gf.Jpart(xc[j], xc[l],domt ); 

                             }
                        }
          }

     pp=plan_0(Nc, Opd  ); 

     #pragma omp parallel for collapse(1)   
            for(int l=0; l<Nc;l++){
                 Tchv_cof(Nc,Opd+l*Nc,pp);                  
                
    }

     
     #pragma omp parallel
     {
        dcomp *v=new dcomp[Nc]; 

        #pragma omp for collapse(1)  
          for(int j=0; j<Nc;j++)
            {
             for(int i=0; i<Nc;i++){

                 v[i] = Opd[j+i*Nc];

             }

              Tchv_cof(Nc,v,pp); 

              for(int i=0; i<Nc;i++)
             {

                 Opd[+j+i*Nc] = v[i];

             }

         }
        
        delete[] v; 
     }
    
    fftw_destroy_plan(pp);
    fftw_cleanup();	  
    
    
    return;          
 
    
    
}