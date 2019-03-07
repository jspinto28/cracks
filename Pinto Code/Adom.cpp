#include <omp.h>
#include "in_cells.h"
#include "prods.h"
#include "geo.h"
#include "tchv_cof.h"
#include "greenfunction.h"

#include "mex.h"
using namespace std;
typedef complex<double> dcomp;
const double pi= 4.0*atan(1.0);




void bdomPlaneNs(dcomp* b,double k,double angle, int* N,int dom_op,
        int Dt,int Dt_0,int NT,double* bgeo, int Nb, double delta ){
 // b inicializado en 0. 
     int Nc=2*(2*N[0]+1)+28;  
     int sign_dom=-1;
     double* xc;
     xc = chebpts(Nc); 
     dcomp** Opd =new dcomp*[Dt];       
     for(int t=0;t< Dt_0; t++){
            Opd[t]=new dcomp[Nc];    
 
      }   
     
     in_cells_Ns(Opd,k,angle,Dt,dom_op,Nc,xc,bgeo,Nb,delta);
     int ns=0;
     for (int s=0; s<Dt;s++){
         if( s>0) ns+=N[s-1];
           
            for(int l=0;l<N[s];l++){
                b[l+ns]=pldNs(l,Opd[s],dom_op);               
            }
       
     }
     
       

          for(int t=0;t< Dt_0; t++){
            delete [] Opd[t];     
   
      }   
      delete [] Opd;

      return;    
}

void bdomDevPlaneNs(dcomp* b,double k,double angle, int* N,int dom_op,
        int Dt,int Dt_0,int NT,double* bgeo, int Nb, double delta,
        int p, int component){
 // b inicializado en 0. 
     int Nc=2*(2*N[0]+1)+28;  
     int sign_dom=-1;
     double* xc;
     xc = chebpts(Nc); 
     dcomp** Opd =new dcomp*[Dt];       
     for(int t=0;t< Dt_0; t++){
            Opd[t]=new dcomp[Nc];    
 
      }   
     
     in_cells_NsDev(Opd,k,angle,Dt,dom_op,Nc,xc,bgeo,Nb,delta,p,component);
               
     int ns=0;
     for (int s=0; s<Dt;s++){
         if( s>0) ns+=N[s-1];
           
            for(int l=0;l<N[s];l++){
                b[l+ns]=pldNs(l,Opd[s],dom_op);               
            }
       
     }
     
       

          for(int t=0;t< Dt_0; t++){
            delete [] Opd[t];     
   
      }   
      delete [] Opd;

      return;    
}


void FarOneComputation(dcomp* A,double k,int N,
        double* bgeo, int Nb, double delta,
        int p, int component, double *y, int Nobs){

    
     int Nc=2*(2*N+1)+28;  
     double* xc;
     xc = chebpts(Nc); 
     
     FarOneCofs(k,A,Nc,xc,bgeo,Nb,delta,p,component,y,Nobs);
     
     return;    
}



//////////////////////////////////////////////////////////////////////

void VdomNsPerSegment2(dcomp** A, double k,int N,int Nc, int Nt,
        int segmentTrial,int segmentTest, double* bgeo, int Nb, double delta){
  
   double* xc;
    
   xc = chebpts(Nc); 
   
   dcomp* Opd =new dcomp[Nc*Nc];     


   VNsRegularPart(k,Opd,segmentTrial,segmentTest,Nc,xc,bgeo,Nb,delta);      
   
   int nthread = omp_get_max_threads();
   omp_set_num_threads(nthread);

	
   double* cm = new double[N];

   cm[0] =pi; 
   
   for(int ii(1); ii<N; ++ii)
   {
     cm[ii] = 0.5*pi;
   }
   
	
//    Full matrix from the regular part. 	
   #pragma omp parallel for collapse(2)  	
   for(int l=0;l<N;l++){ 

	 for(int m=0;m<N;m++){  

		A[l][m]=Opd[m+l*Nc]*cm[m]*cm[l]; 

	}

   }

   if( segmentTrial == segmentTest)
   {
	if(k>=0.0000001)	
	{
//              
	   int Nc2 = 2*(Nc+Nt);
       
       double* xc2;

	   xc2 = chebpts(Nc2);    
		
	   dcomp* Opd2 =new dcomp[Nc2*Nc2]; 

	   JpartKernel(k,Opd2,segmentTrial,Nc2,xc2,bgeo,Nb,delta);

 	   double* an = new double[Nt];

       double* cn = new double[N+Nt];

       an[0] = 0.5/pi*log(2); 

	   for(int ii(1); ii<Nt; ++ii)
	   {
	     an[ii] = 1.0/(pi*((double)ii));
	   }

       cn[0] = pi; 

	   for(int ii(1); ii<Nt+N; ++ii)
	   {
	     cn[ii] = 0.5*pi;
	   }

       #pragma omp parallel for collapse(2)  
	   for(int m=0; m<N; m++)
	   {
             for(int l=0; l<N; l++)
             {
               for(int n=0; n<Nt; n++)
               {
                 A[l][m] += 0.25*an[n]*(
                    Opd2[(m+n)+(l+n)*Nc2]*cn[m+n]*cn[l+n]+
                    Opd2[(m+n)+abs(l-n)*Nc2]*cn[m+n]*cn[abs(l-n)]+
                    Opd2[abs(m-n)+(l+n)*Nc2]*cn[abs(m-n)]*cn[l+n]+
                    Opd2[abs(m-n)+abs(l-n)*Nc2]*cn[abs(m-n)]*cn[abs(l-n)]
                    );

               }
             }
           }
	  

	   delete[] Opd2;
	   
	   delete[] an;

           delete[] cn; 

    }
	else
	{

	  A[0][0] += 0.5*pi*log(2);

       //   #pragma omp parallel for collapse(1) 
	   for(int ii=1; ii<N; ii++)
	   {
	     A[ii][ii] += 1.0*pi/(4.0*((double)ii));
	   }

	}
   }	       
 
   delete[] Opd; 
   
   delete[] cm; 

   return; 
    

}

void VdomNsPerSegment3(dcomp** A, double k,int N,int Nc, int Nt,
        int segmentTrial,int segmentTest,
        double tol, double* bgeo, int Nb, double delta){
  
   double* xc;
    
   xc = chebpts(Nc); 
   
   dcomp* Opd =new dcomp[Nc*Nc];     


   VNsRegularPart(k,Opd,segmentTrial,segmentTest,Nc,xc,bgeo,Nb,delta);      
   
   int nthread = omp_get_max_threads();
   omp_set_num_threads(nthread);
   
   double* cm = new double[N];

   cm[0] =pi; 
   
   for(int ii(1); ii<N; ++ii)
   {
     cm[ii] = 0.5*pi;
   }
   
	
   //Full matrix from the regular part. 	
   //#pragma omp parallel for collapse(2)  	
   for(int l=0;l<N;l++){ 

	 for(int m=0;m<N;m++){  

		A[l][m]=Opd[m+l*Nc]*cm[m]*cm[l]; 

	}

   }

   if( segmentTrial == segmentTest)
   {       
	if(k>=0.0000001)	
	{
             
	   int Nc2 = 2*(Nc+Nt);
       
       double* xc2;

	   xc2 = chebpts(Nc2);    
		
	   dcomp* Opd2 =new dcomp[Nc2*Nc2]; 

	   JpartKernel(k,Opd2,segmentTrial,Nc2,xc2,bgeo,Nb,delta);

 	   double* an = new double[Nt];

       double* cn = new double[N+Nt];

       an[0] = 0.5/pi*log(2); 

	   for(int ii(1); ii<Nt; ++ii)
	   {
	     an[ii] = 1.0/(pi*((double)ii));
	   }

       cn[0] = pi; 

	   for(int ii(1); ii<Nt+N; ++ii)
	   {
	     cn[ii] = 0.5*pi;
	   }

       #pragma omp parallel for collapse(2)  
	   for(int m=0; m<N; m++)
	   {
             for(int l=0; l<N; l++)
             {
               for(int n=0; n<Nt; n++)
               {
                 A[l][m] += 0.25*an[n]*(
                    Opd2[(m+n)+(l+n)*Nc2]*cn[m+n]*cn[l+n]+
                    Opd2[(m+n)+abs(l-n)*Nc2]*cn[m+n]*cn[abs(l-n)]+
                    Opd2[abs(m-n)+(l+n)*Nc2]*cn[abs(m-n)]*cn[l+n]+
                    Opd2[abs(m-n)+abs(l-n)*Nc2]*cn[abs(m-n)]*cn[abs(l-n)]
                    );

               }
             }
           }
	  

	   delete[] Opd2;
	   
	   delete[] an;

           delete[] cn; 

    }
	else
	{

	  A[0][0] += 0.5*pi*log(2);

       //#pragma omp parallel for collapse(1) 
	   for(int ii=1; ii<N; ii++)
	   {
	     A[ii][ii] += 1.0*pi/(4.0*((double)ii));
	   }

	}
   }	   
   
   delete [] cm; 

   delete [] Opd; 
   
   //#pragma omp parallel for collapse(2)  
   for(int m=0; m<N; m++)
   {
         for(int l=0; l<N; l++)
         {
             if(abs(A[l][m])<tol)
             {
                 A[l][m]=0.0; 
             }
         }
   }
   

   return; 
   

    

}

void VdomNsPerSegment5(dcomp** A, double k,int N,int Nq, 
        double* xq,double wq,int Nq0, double* xq0, double wq0,
        int segmentTrial,int segmentTest,
        double tol, int Lmax, double* bgeo, int Nb, double delta){
  
  
    double* OpdmR =new double[Nq*Nq0];       
    
    double* OpdlR = new double [Nq*Nq0]; 
    
    double* OpdmI =new double[Nq*Nq0];       
    
    double* OpdlI = new double [Nq*Nq0]; 
    
    double* TN = new double [Nq*N]; 
    
    int n1=N; 
    
    int n2=N; 
    
    in_cellsVNsPerSegment( k, OpdmR,OpdmI,  OpdlR,OpdlI, TN,
         segmentTrial, segmentTest, N,
         Nq, xq,
         Nq0, xq0,bgeo,Nb, delta);
    
 
    
    TN = t_polynomial( Nq, N-1, xq );   
    
//     for(int jj =0; jj< N ; ++jj)       
//     {
//         for( int ii = 0; ii<Nq; ++ii)
//         {
//             mexPrintf("%f , %f | %f \n",xq[ii],xq0[jj],TN[ii+Nq*jj]);   
//         }
//     }   
    
    int lev =0; 
    
    int a = 0; 
    
    int b= N-1; 
   
    int ti,tc,td; 
    
    double vi,vc,vd;
        
    int mid;
        
    while(lev < Lmax)
    {
        mid = (int)((a+b)/2); 
    
        ti = mid-1; 
        
        tc = mid; 
        
        td = mid+1;
        
        vi =  abs(pVNs0(ti ,OpdmR,OpdmI, TN ,Nq, wq,Nq0,wq0,N ));         
          
        vc =  abs(pVNs0(tc ,OpdmR,OpdmI, TN ,Nq, wq,Nq0,wq0,N ));   
        
        vd =  abs(pVNs0(td ,OpdmR,OpdmI, TN ,Nq, wq,Nq0,wq0,N ));    
        
//         if((segmentTrial == 20)&&  (segmentTest==0))
//         {
//             mexPrintf("Row 0:  %d , %f; %f ; %f \n", mid,vi,vc,vd); 
//             
//         }

        if( ((vd < 0.5*tol)&&(vc < 0.5*tol))||
                ((vi < 0.5*tol)&&(vc < 0.5*tol))) 
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
    
    b= N-1; 
    
    lev=0;

    while(lev < Lmax)
    {
        mid = (int)((a+b)/2); 
    
        ti = mid-1; 
        
        tc = mid; 
        
        td = mid+1;
                
        vi =  abs(pVNs0(ti ,OpdlR,OpdlI, TN ,Nq, wq,Nq0,wq0,N ));         
          
        vc =  abs(pVNs0(tc ,OpdlR,OpdlI, TN ,Nq, wq,Nq0,wq0,N ));   
        
        vd =  abs(pVNs0(td ,OpdlR,OpdlI, TN ,Nq, wq,Nq0,wq0,N ));   
        
//                 if((segmentTrial == 20)&&  (segmentTest==0))
//         {
//             mexPrintf("Col 0:  %d , %f; %f ; %f \n", mid,vi,vc,vd); 
//             
//         }
        
        if( ((vd < 0.5*tol)&&(vc < 0.5*tol))||
                ((vi < 0.5*tol)&&(vc < 0.5*tol))) 
        {
            
            b = mid; 
            
        }
        else
        {            
            a = mid; 
        }
        
        lev++; 
                
        
    }
    
     n2 =b; 
    
    
     int n = max(n1,n2); 
    
//     mexPrintf("%d,",n); 
    
    int Nc = 2*n+28; 
    
    VdomNsPerSegment2(A, k,n,Nc, 0,
         segmentTrial, segmentTest, bgeo,Nb,delta);
    
  
        

    delete [] OpdmR; 
    
    delete [] OpdmI; 
    
    delete [] OpdlR;
    
    delete [] OpdlI;
    
    delete [] TN; 


}

void VDerivedNsPerSegment2(dcomp** A, double k,int N,int Nc, int Nt,
        int segmentTrial,int segmentTest, double* bgeo, int Nb, double delta,
        int p, int component){
  
   double* xc;
    
   xc = chebpts(Nc); 
   
   dcomp* Opd =new dcomp[Nc*Nc];     

   VNsDerivedRegularPart
           (k,Opd,segmentTrial,segmentTest,Nc,xc,bgeo,Nb,delta,p,component);      
   
   int nthread = omp_get_max_threads();
   omp_set_num_threads(1);

	
   double* cm = new double[N];

   cm[0] =pi; 
   
   for(int ii(1); ii<N; ++ii)
   {
     cm[ii] = 0.5*pi;
   }
   
	
//    Full matrix from the regular part. 	
   #pragma omp parallel for collapse(2)  	
   for(int l=0;l<N;l++){ 

	 for(int m=0;m<N;m++){  

		A[l][m]=Opd[m+l*Nc]*cm[m]*cm[l]; 

	}

   }

   if( segmentTrial == segmentTest)
   {
	if(k>=0.0000001)	
	{
//              
	   int Nc2 = 2*(Nc+Nt);
       
       double* xc2;

	   xc2 = chebpts(Nc2);    
		
	   dcomp* Opd2 =new dcomp[Nc2*Nc2]; 

	   JpartDevKernel(k,Opd2,segmentTrial,Nc2,xc2,bgeo,Nb,delta,p,component);

 	   double* an = new double[Nt];

       double* cn = new double[N+Nt];

       an[0] = 0.5/pi*log(2); 

	   for(int ii(1); ii<Nt; ++ii)
	   {
	     an[ii] = 1.0/(pi*((double)ii));
	   }

       cn[0] = pi; 

	   for(int ii(1); ii<Nt+N; ++ii)
	   {
	     cn[ii] = 0.5*pi;
	   }

       #pragma omp parallel for collapse(2)  
	   for(int m=0; m<N; m++)
	   {
             for(int l=0; l<N; l++)
             {
               for(int n=0; n<Nt; n++)
               {
                 A[l][m] += 0.25*an[n]*(
                    Opd2[(m+n)+(l+n)*Nc2]*cn[m+n]*cn[l+n]+
                    Opd2[(m+n)+abs(l-n)*Nc2]*cn[m+n]*cn[abs(l-n)]+
                    Opd2[abs(m-n)+(l+n)*Nc2]*cn[abs(m-n)]*cn[l+n]+
                    Opd2[abs(m-n)+abs(l-n)*Nc2]*cn[abs(m-n)]*cn[abs(l-n)]
                    );

               }
             }
           }
	  

	   delete[] Opd2;
	   
	   delete[] an;

           delete[] cn; 

    }
   }	       
 
   delete[] Opd; 
   
   delete[] cm; 

   return; 
    

}