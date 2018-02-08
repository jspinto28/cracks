/////////////////////////////////////////////////////////////////
// Interface for quadratures: 

//Gauss-Legendre:
extern void lgwt(int NumPoints, double a, double b, double  x[],
        double  w[]); 

//Laguerre: 
extern void laguerre(int NumPoints, double a, double b, double  x[],
        double  w[]); 

//////////////////////////////////////////////////////////////////
//Extern functions: 

extern void cdgqf ( int nt, int kind, double alpha, double beta, 
        double t[], double wts[] );

extern void cgqf ( int nt, int kind, double alpha, double beta, double a,
        double b,  double t[], double wts[] );

extern double class_matrix ( int kind, int m, double alpha, double beta, 
        double aj[],   double bj[] );

extern void imtqlx ( int n, double d[], double e[], double z[] );

extern void parchk ( int kind, int m, double alpha, double beta );

extern double r8_epsilon ( );

extern double r8_sign ( double x );

extern void scqf ( int nt, double t[], int mlt[], double wts[], int nwts,
        int ndx[],   double swts[], double st[], int kind, double alpha, 
        double beta, double a,   double b );

extern void sgqf ( int nt, double aj[], double bj[], double zemu,
        double t[],   double wts[] );

