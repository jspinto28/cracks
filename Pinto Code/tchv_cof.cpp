#include <complex>
#include <fftw3.h>
#include "tchv_cof.h"
#include "mex.h"
using namespace std;
typedef complex<double> dcomp;
const dcomp i(0.0,1.0);
const double pi= 4.0*atan(1.0);


void *Tchv_cof(int n, dcomp* fxc,fftw_plan p){
    // n length of fxc 
    // fxc vector of the function valued at cheby points    
    int m; 
    m= n-1;     
    dcomp* vec =new dcomp[2*n-2];     
    vec[0]=fxc[n-1]; 
    vec[n-1]=fxc[0];
    for(int j=1; j<n-1;j++){
        vec[j]=fxc[n-1-j];
        vec[j+n-1]=fxc[j];
    }          
//     
//     for(int ii(0); ii<n;++ii)
//     {
//         mexPrintf("%f %f\n",fxc[ii].real(), fxc[ii].imag()); 
//     }
//     
//     mexPrintf("------------------------\n");
    

    fftw_execute_dft(p, reinterpret_cast<fftw_complex*>(vec),reinterpret_cast<fftw_complex*>(vec));   
    for (int j=0; j<n;j++){
           
        fxc[j]=vec[j]/(double)m;
    }    
    fxc[0]=0.5*fxc[0];
    fxc[n-1]=0.5*fxc[n-1];
    delete vec;           
   
}

fftw_plan plan_0(int n, dcomp* fxc){
 //se ejecuta una sola vez para obtener el plan.  
    fftw_plan p;
    dcomp* vec= new dcomp[2*n-2];    
    vec[0]=fxc[n-1]; 
    vec[n-1]=fxc[0];
    for(int j=1; j<n-1;j++){
        vec[j]=fxc[n-1-j];
        vec[j+n-1]=fxc[j];
    }
    fftw_complex* vecf =reinterpret_cast<fftw_complex*>(vec);
    p=fftw_plan_dft_1d(2*n-2, vecf, vecf, FFTW_FORWARD, FFTW_ESTIMATE);   
    
    delete [] vec ;
    return p; 
    
}

double *chebpts(int n){
    double* x= new double[n]; 
    for(int j=0; j<n;j++){
        x[j]=cos((n-1-j)*pi/(n-1));
    }
    return x;
}

