% clear all; 

% mex ArcSolver.cpp CXXFLAGS="\$CXXFLAGS -std=c++0x -O3 -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" CXXOPTIMFLAGS="-O3"...
% 		   LDOPTIMFLAGS=" -O3" -I/home/matlab/boost_1_55_0 -lfftw3 -lm ; 

Ndof = 64;
omega = 2.8; 
lambda = -0.8;
mu = 1; 
alpha = pi*0.25; 

bdecayexp = 3; 
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
Y = 2*rand(totalN,1)-1; 

% for ii=1:Narcs
% 
%     %points for the x part
%     Y{ii}{1} = 2*rand(Nx{ii},1)-1; 
%     %points for the y part
%     Y{ii}{2} = 2*rand(Ny{ii},1)-1; 
%     
%     
% end 

%think about how produce the points... only one large vector or
%tenzosization???

[coefx,coefy] =GetGeoCofs(coefx0,coefy0,Nx,Ny,dcoefs,Y);


%can be used with compresion parameters!
% [u,V]= DirectSolverFull('s',false,coefx,coefy,omega,lambda,mu,alpha,Ndof);
% plot(1:length(u),log10(abs(u(:,1))))

hold on; 
Nplot = 64;
ts = linspace(-1,1,Nplot); 
cc =['b','r','g','y','k']; 

for ii=1:Narcs
    
    XY = ArcSolver(100,0,ii,length(coefx{ii}),length(coefy{ii}),coefx{ii},coefy{ii},Nplot,ts);
    plot(XY(:,1),XY(:,2),cc(ii)); 
end 

 axis equal;