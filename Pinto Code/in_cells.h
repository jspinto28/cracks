#include <complex>
using namespace std;
typedef complex<double> dcomp;

extern void in_cells_Ns(dcomp** Opd,double k,double angle,int Dt,
        int dom_op, int Nc,double* xc, double* bgeo, int Nb, double delta ); 

// extern void in_cellsBcofs(int N,dcomp*** Opd,int segmentTrial,
//          int Nc,double* xc, int order); 

extern void VNsRegularPart(double k,dcomp* Opd, int domt,
        int doms, int Nc,double* xc, double* bgeo, int Nb, double delta);

extern void JpartKernel(double k, dcomp* Opd, int domt,
       int Nc,double* xc, double* bgeo, int Nb, double delta); 

// extern void SlTcofsComp (double k, int Nc, int Ncurves, 
//         int Nxy, double* x, double* y, dcomp* Cofs); 

extern void in_cellsVNsPerSegment(double k,double* OpdmR, double* OpdmI,
        double* OpdlR, double* OpdlI,double* TN,
        int segmentTrial,int segmentTest,int N,
        int Nq, double*xq,
        int Nq0, double*xq0, double* bgeo, int Nb, double delta);


extern void VNsDerivedRegularPart(double k,dcomp* Opd, int domt,
        int doms, int Nc,double* xc, double* bgeo, int Nb, double delta, 
        int p, int component); 

extern void JpartDevKernel(double k,dcomp* Opd, int domt,
       int Nc,double* xc, double* bgeo, int Nb, double delta, int p, int component); 

extern void in_cells_NsDev(dcomp** Opd,double k,double angle,int Dt,int dom_op,
        int Nc,double* xc, double* bgeo, int Nb, double delta,
        int p, int component); 

extern void FarOneCofs(double k,dcomp* Opd,int Nc,double* xc, 
        double* bgeo, int Nb, double delta, 
        int p, int component,
        double *y, int Nobs); 