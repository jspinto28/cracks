#include <complex>
using namespace std;
typedef complex<double> dcomp;

extern void VdomNs(dcomp** A, double k,int* N,int Nc, int Nq,
        double*xq, double* wq,int dom_op,int Dt, int NT,double* Y1,
        double* Y2,double delta, int S);

extern void bdomPlaneNs(dcomp* b, int* N,double alpha, double k,int Dt,int NT,
        double* Y1,double* Y2,double delta,int S);
        
extern void VdomNsPerSegment(dcomp** A, double k,int N,int Nc, int Nq,
        double*xq, double* wq,int segmentTrial,int segmentTest,
        double* Y1,double* Y2,double delta,int S);