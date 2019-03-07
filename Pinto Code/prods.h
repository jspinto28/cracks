#include <complex>
using namespace std;
typedef complex<double> dcomp;

extern double *t_polynomial ( int m, int n, double x[] ); 
extern double chebT ( int n, double x );
extern double *u_polynomial ( int m, int n, double x[] );
extern double chebU ( int n, double x );
extern double chebU_p ( int n, double x );

extern double *t_polynomialD ( int m, int n, double x[], double delta ); 

extern double *u_polynomialD ( int m, int n, double x[],double delta ); 

// extern double *t_polynomial ( int m, int n, double* x ); 


extern dcomp pWonly(double alpha, double period,
        int m, int l ,dcomp*** ft_cofs,dcomp** wb,int Ng,double* xg,double* wg, double k,int domt,int doms,int dom_op);


extern dcomp pldNs(int l,dcomp* ft_cofs,int dom_op); 

extern dcomp pVNs0(int m ,double* OpdR, double* OpdI, double* TN ,int Nq,double wq,
        int Nq0, double wq0, int N );

extern double *t_polynomialDim0 ( int m, int n, double x[], double delta,int dim0 ); 

extern double *u_polynomialDim0 ( int m, int n, double x[],double delta, int dim0 );