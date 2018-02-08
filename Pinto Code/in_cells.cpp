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




void in_cellsVNs(double k,dcomp**** Opd,int Dt, int dom_op,
        int Nq, double*xq, int Nc,double* xc,double* Y1, double* Y2,
        double delta,int S){        
   
    fftw_plan p;
    int nthread = omp_get_max_threads();
    omp_set_num_threads(nthread);


    #pragma omp parallel
        {
        
         FreeSpacePer gf(k,S,Y1,Y2,delta); 


         #pragma omp for collapse(4)          
            for(int t=0; t<Dt;t++){             
                for(int s=0; s<Dt;s++){                
                    for(int j=0; j<Nq;j++){                    
                        for(int l=0;l<Nc;l++){                    

                            Opd[t][s][j][l]=GB(xc[l],xq[j],t,s, dom_op,gf);      
                            
//                             if(t ==1) 
//                             {
//                                 mexPrintf("%f %f",Opd[t][s][j][l].real(), 
//                                         Opd[t][s][j][l].imag()); 
//                             }
                        }
                    }
                }
            }
        }   

      p=plan_0(Nc, Opd[0][0][0]  );             
     #pragma omp parallel for collapse(3)   
       for(int t=0; t<Dt;t++){
         for(int s=0; s<Dt;s++){
            for(int j=0; j<Nq;j++){
                 Tchv_cof(Nc,Opd[t][s][j],p);                  
                
    }}}  
    fftw_destroy_plan(p);
    fftw_cleanup();	  
    
    
    return;
    
    //ESTA FUNCION SE PODRIA MEJORAR, como se necesitan guardar arrays 
    // de tamaño Nc, bastan los primeros 2N+1 coeficientes. 
}


void in_cells_Ns(dcomp** Opd,int Dt,double alpha,
        double kext, int Nc,double* xc,double* Y1, double* Y2,
        double delta,int S){    
    fftw_plan p;
    int nthread = omp_get_max_threads();
    omp_set_num_threads(nthread);       


         #pragma omp parallel for collapse(2)  
        for(int t=0; t<Dt;t++){         
            for(int l=0;l<Nc;l++){          
                        Opd[t][l]=rhsNs(t,xc[l],alpha,kext,Y1,Y2,delta,S);   
    
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

void in_cellsVNsPerSegment(double k,dcomp** Opd,int segmentTrial,
        int segmentTest,int Nq, double*xq, int Nc,double* xc,
        double* Y1, double* Y2,double delta,int S){        
   
    fftw_plan p;
    int nthread = omp_get_max_threads();
    omp_set_num_threads(nthread);


    #pragma omp parallel
        {
    
        
        FreeSpacePer gf(k,S,Y1,Y2,delta); 


         #pragma omp for collapse(2)       
                    for(int j=0; j<Nq;j++){                    
                        for(int l=0;l<Nc;l++){                    

                            Opd[j][l]=GB(xc[l],xq[j],segmentTrial,
                                    segmentTest, 0,gf);      
                            
//                             if(t ==1) 
//                             {
//                                 mexPrintf("%f %f",Opd[t][s][j][l].real(), 
//                                         Opd[t][s][j][l].imag()); 
//                             }
                        }
                    }
                }


      p=plan_0(Nc, Opd[0]  );             
     #pragma omp parallel for collapse(1)   
            for(int j=0; j<Nq;j++){
                 Tchv_cof(Nc,Opd[j],p);                  
                
    } 
    fftw_destroy_plan(p);
    fftw_cleanup();	  
    
    
    return;
    
    //ESTA FUNCION SE PODRIA MEJORAR, como se necesitan guardar arrays 
    // de tamaño Nc, bastan los primeros 2N+1 coeficientes. 
}
