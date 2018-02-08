#include <complex>
using namespace std;
typedef complex<double> dcomp;

extern void in_cellsVNs(double k,dcomp**** Opd,int Dt, int dom_op,
        int Nq, double*xq, int Nc,double* xc,double* Y1, double* Y2,
        double delta,int S); 

extern void in_cells_Ns(dcomp** Opd,int Dt,double alpha,
        double kext, int Nc,double* xc,double* Y1, double* Y2,
        double delta,int S); 

extern void in_cellsVNsPerSegment(double k,dcomp** Opd,int segmentTrial,
        int segmentTest,int Nq, double*xq, int Nc,double* xc,
        double* Y1, double* Y2,double delta,int S); 