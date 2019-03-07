#include <complex>
#include "greenfunction.h"

using namespace std;
typedef complex<double> dcomp;

extern dcomp rhsNs(double k,double angle,int domt, double s,double* bgeo, int Nb, double delta); 

extern dcomp GlobalRhsNs(double k,double angle, double x, double y); 

extern dcomp GBNs(double t, double s, int domt, int doms,GreenFunctionBase & gf);

extern dcomp DevrhsNs(double k,double angle,int domt, double s,double* bgeo,
        int Nb, double delta,int p,int component);

extern dcomp FarOne(double k,double s,double* bgeo,
        int Nb, double delta,int p,int component,
        double y1, double y2); 