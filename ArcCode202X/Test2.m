% clear all; 

mex ArcSolver.cpp CXXFLAGS="\$CXXFLAGS -std=c++0x -O3 -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" CXXOPTIMFLAGS="-O3"...
		   LDOPTIMFLAGS=" -O3" -I/home/matlab/boost_1_55_0 -lfftw3 -lm ; 

Ndof = 64;
omega = 2.8; 
lambda = -0.8;
mu = 1; 
alpha = pi*0.25; 

bdecayexp = 3; 
% bdcayscale = 0.125/2; 
bdcayscale = zeros(129,1)+1;
bdcayscale(1:4) = [0.1;0.1;0.2;0.2]; 

Narcs = 1; 
coefx0 = cell(Narcs,1); 
coefy0 = cell(Narcs,1); 
Nx = cell(Narcs,1); 
Ny = cell(Narcs,1); 

coefx0{1}= [0;1];  
coefy0{1}= [0;0];
coefx0{2}= [0.5;1];  
coefy0{2}= [1.5;0];

Nx{1} = 8; 
Ny{1} = 8; 
Nx{2} = 8; 
Ny{2} = 8; 

dcoefs = 0;
if(length(bdcayscale) ==1)
    dcoefs= bdcayscale*[1;1./(1:128)'].^bdecayexp; 
else
     dcoefs= bdcayscale.*([1;1./(1:128)'].^bdecayexp); 
end 
     
t0 = 0; 
Nt =128+1; 
ts = linspace(-1,1,Nt); 


G200 = ArcSolver(200,false,1,length(coefx{1}),length(coefy{1}),coefx{1},coefy{1},...
1,length(coefx{1}),length(coefy{1}),coefx{1},coefy{1},...
omega,lambda,mu,Ndof,8,t0,Nt,ts);

plot(ts,real(G200),'x');