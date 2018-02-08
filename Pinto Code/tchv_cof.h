#include <complex>
#include <fftw3.h>
using namespace std;
typedef complex<double> dcomp;
extern void *Tchv_cof(int n, dcomp* fxc,fftw_plan p); 
extern fftw_plan plan_0(int n, dcomp* fxc); 
extern  double *chebpts(int n);
