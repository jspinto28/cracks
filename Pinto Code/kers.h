#include <complex>
#include "greenfunction.h"

using namespace std;
typedef complex<double> dcomp;
extern void pk(double* dk, double* pk, double t, double s,  int domt, int doms, int dom_op); 
extern void pka(double* dk, double* pk, double t, double s,  int domt, int doms, int dom_op); 
extern void pw(double* dw, double* pw, double t, double s,  int domt, int doms, int dom_op); 





extern dcomp G0(double t, double s);

extern dcomp GB(double t, double s, int domt, int doms, int dom_op,
        GreenFunctionBase & gf);
extern dcomp KB(double t, double s, int domt, int doms, int dom_op,
         GreenFunctionBase & gf);
extern dcomp KAB(double t, double s, int domt, int doms, int dom_op,
         GreenFunctionBase & gf);
extern dcomp WB(double t, double s, int domt, int doms, int dom_op,
         GreenFunctionBase & gf);
extern dcomp GBB(double t, double s, int domt, int doms, int dom_op,
         GreenFunctionBase & gf);
extern dcomp WA1(double t, double s, int domt, int doms, int dom_op,
         GreenFunctionBase & gf);
extern dcomp WA2(double t, double s, int domt, int doms, int dom_op,
         GreenFunctionBase & gf);
extern dcomp WA3(double t, double s, int domt, int doms, int dom_op,
         GreenFunctionBase & gf);
extern dcomp WA4(double t, double s, int domt, int doms, int dom_op,
         GreenFunctionBase & gf);
extern dcomp WB1(double alpha, double period, int sig,double t, double s,
        int domt, int doms, int dom_op,GreenFunctionBase & gf);
extern dcomp WB2(int sig,double t, double s, int domt, int doms, int dom_op,
         GreenFunctionBase & gf);

//olde code (no class for the GF) 
// extern dcomp xG0(double t, double s);
// extern dcomp xGB(double k,double t, double s, int domt, int doms, int dom_op);
// extern dcomp xKB(double k,double t, double s, int domt, int doms, int dom_op);
// extern dcomp xKAB(double k,double t, double s, int domt, int doms, int dom_op);
// extern dcomp xWB( double k,double t, double s, int domt, int doms, int dom_op);
// extern dcomp xGBB(double k,double t, double s, int domt, int doms, int dom_op);
// extern dcomp xWA1(double k, double t, double s, int domt, int doms, int dom_op);
// extern dcomp xWA2(double k, double t, double s, int domt, int doms, int dom_op);
// extern dcomp xWA3(double k, double t, double s, int domt, int doms, int dom_op);
// extern dcomp xWA4(double k, double t, double s, int domt, int doms, int dom_op);
// extern dcomp xWA01( double t, double s, int domt,  int dom_op);
// extern dcomp xWA02( double t, double s, int domt,  int dom_op);
// extern dcomp xWA03( double t, double s, int domt,  int dom_op);
// extern dcomp xWA04( double t, double s, int domt,  int dom_op);
// extern dcomp xWB1(int sig,double k, double t, double s, int domt, int doms, int dom_op);
// extern dcomp xWB2(int sig,double k, double t, double s, int domt, int doms, int dom_op);
// extern dcomp xWB01(double t, double s, int domt,  int dom_op);
// extern dcomp xWB02(double t, double s, int domt,  int dom_op);


extern dcomp uinc_d(double k0, double s, int doms ,int dom_op,int l, double R);
extern dcomp uinc_n(double k0, double s, int doms ,int dom_op, int l, double R);
extern dcomp pointSource(double k0, double s, int doms,int dom_op ,double alpha); 
extern dcomp pointSourceN(double k0, double s, int doms,int dom_op ,double alpha); 

extern double Argument(double x, double y); 

extern double TrialFunction (int n, double x); 

extern double TrialFunctionDev (int n, double x); 


extern dcomp pointSourcePer(double k0, double s, int doms,int dom_op ,double alpha,
        int NPer, double* an, double* bn); 

extern dcomp pointSourceNPer(double k0, double s, int doms,int dom_op ,double alpha,
         int NPer, double* an, double* bn); 

extern dcomp rhsNs(int domt, double s,double alpha, double kext,
        double* Y1, double* Y2, double delta, int S); 

extern dcomp GlobalRhsNs(double x, double y); 
        