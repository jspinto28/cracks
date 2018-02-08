#include <math.h> 
#include "Curves.h" 
#include "prods.h" 
using namespace std; 

void curves(double* x, double* y,double t,int curv){ 

    *x = 
    0.000000*1.000000*chebT(0,t)+ 
    1.000000*0.176777*chebT(1,t)+ 
    0; 
    *y = 
    1.000000*1.000000*chebT(0,t)+ 
    0.000000*0.176777*chebT(1,t)+ 
    0; 
}

void tang_curves(double* x, double* y,double t,int curv){ 

    *x = 
    1.000000*0.176777*chebU(0,t)+ 
    0; 
    *y = 
    0.000000*0.176777*chebU(0,t)+ 
    0; 
}
