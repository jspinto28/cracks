#include <complex>
using namespace std;
typedef complex<double> dcomp;

void bdomPlaneNs(dcomp* b,double k,double angle, int* N,int dom_op,
        int Dt,int Dt_0,int NT,double* bgeo, int Nb, double delta );

extern void VdomNsPerSegment2(dcomp** A, double k,int N,int Nc, int Nt,
        int segmentTrial,int segmentTest, double* bgeo, int Nb, double delta);

extern void VdomNsPerSegment3(dcomp** A, double k,int N,int Nc, int Nt,
        int segmentTrial,int segmentTest,double tol, double* bgeo, int Nb, double delta); 

extern void VdomNsPerSegment5(dcomp** A, double k,int N,int Nq, 
        double* xq,double wq,int Nq0, double* xq0, double wq0,
        int segmentTrial,int segmentTest,
        double tol, int Lmax, double* bgeo, int Nb, double delta);

extern void VDerivedNsPerSegment2(dcomp** A, double k,int N,int Nc, int Nt,
        int segmentTrial,int segmentTest, double* bgeo, int Nb, double delta,
        int p, int component); 

extern void bdomDevPlaneNs(dcomp* b,double k,double angle, int* N,int dom_op,
        int Dt,int Dt_0,int NT,double* bgeo, int Nb, double delta,
        int p, int component); 

extern void FarOneComputation(dcomp* A,double k,int N,
        double* bgeo, int Nb, double delta,
        int p, int component, double *y, int Nobs); 