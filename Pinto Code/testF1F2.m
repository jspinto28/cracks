%  mex mex_funs.cpp IntegrationRules.cpp Adom.cpp kers.cpp prods.cpp geo.cpp in_cells.cpp tchv_cof.cpp greenfunction.cpp CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" -I/home/matlab/boost_1_55_0 -lfftw3 -lm ;

 Nb = 10; 
 
 k=10; 
 
 alpha = pi/3; 
 
 delta = -2;
 
 Nsol = 128;
 
 Nobs = 128; 

 Thetas = linspace(0,2*pi,Nobs+1); 

 Thetas = Thetas(1:Nobs); 

 y = [cos(Thetas'), sin(Thetas')]; 
 
 bgeo = [1,1,0.5,0.5,0.28,0.28,0,0,0,0]; 
 
 eps = 1e-6; 
 
 hgeo = [eps,0,0,0,0,0,0,0,0,0]; 
 
 [tu1,L1] = Solve(k,Nsol,alpha,bgeo+hgeo,delta);
 
 B1 = Farfield(tu1,k,bgeo+hgeo,delta,y); 
 
 [tu0,L0] = Solve(k,Nsol,alpha,bgeo,delta);
 
 B0 = Farfield(tu0,k,bgeo,delta,y); 
 
 [f1,F1] = GetF1(k,tu0,bgeo,delta,y); 
 
 [f2,F2] = GetF2(k,alpha,tu0,L0,bgeo,delta,y); 
 
 norm( (B1-B0)/eps - (F1+F2)*(hgeo'/eps))
 