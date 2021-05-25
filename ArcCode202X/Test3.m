clear all; 

mex ArcSolver.cpp CXXFLAGS="\$CXXFLAGS -std=c++0x -O3 -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" CXXOPTIMFLAGS="-O3"...
		   LDOPTIMFLAGS=" -O3" -I/home/matlab/boost_1_55_0 -lfftw3 -lm ; 

Ndof = 32;
omega = 2.8; 
lambda = -0.8;
mu = 1; 
alpha = pi*0.25; 

Nobs = 6;
xtheta = linspace(0,2*pi,Nobs+1); 
xtheta = xtheta(1:end-1); 
wtheta = zeros(Nobs,1)+2*pi/Nobs; 

bdecayexp = 1.5; 
% bdcayscale = 0.125/2; 
bdcayscale = zeros(129,1)+1;
bdcayscale(1:4) = [0.1;0.1;0.2;0.2]; 

Narcs = 2; 
coefx0 = cell(Narcs,1); 
coefy0 = cell(Narcs,1); 
Nx = zeros(Narcs,1); 
Ny = zeros(Narcs,1); 

coefx0{1}= [0;1];  
coefy0{1}= [0;0];
coefx0{2}= [0.5;1];  
coefy0{2}= [1.5;0];

Nx(1) = 1; 
Ny(1) = 2; 
Nx(2) = 3; 
Ny(2) = 4; 

totalN = sum(Nx)+sum(Ny); 


dcoefs = 0;
if(length(bdcayscale) ==1)
    dcoefs= bdcayscale*[1;1./(1:128)'].^bdecayexp; 
else
     dcoefs= bdcayscale.*([1;1./(1:128)'].^bdecayexp); 
end 
     

checkGeoConfig; 

% random points
% Y = 2*rand(totalN,1)-1; 
Y= [-0.4023    0.0742    0.2930   -0.1523    0.4414    0.3086   -0.3711    0.1680   -0.1992   -0.3477]'; 


[coefx,coefy] =GetGeoCofs(coefx0,coefy0,Nx,Ny,dcoefs,Y);

[u,V]= DirectSolverFull('s',false,coefx,coefy,omega,lambda,mu,alpha,Ndof,1e-10);

fobs = ObservedFunc('both',xtheta,wtheta,u,...
        coefx,coefy,omega,lambda,mu,Ndof); 