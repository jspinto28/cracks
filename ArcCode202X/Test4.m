clear all; 

mex ArcSolver.cpp CXXFLAGS="\$CXXFLAGS -std=c++0x -O3 -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" CXXOPTIMFLAGS="-O3"...
		   LDOPTIMFLAGS=" -O3" -I/home/matlab/boost_1_55_0 -lfftw3 -lm ; 

Ndof = 32;
omega = 2.8; 
lambda = -0.8;
mu = 1; 
alpha = pi*0.25; 

Nobs = 3;
xtheta = linspace(0,2*pi,Nobs+1); 
xtheta = xtheta(1:end-1); 
wtheta = zeros(Nobs,1)+2*pi/Nobs; 

bdecayexp = 3; 
% bdcayscale = 0.125/2; 
bdcayscale = zeros(129,1)+1;
bdcayscale(1:4) = [0.2;0.2;0.3;0.3]; 

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

%%generate quadrature: 
Npoints = 512; 
latticeseq_b2('init0');
Niters = 8;
Quads = zeros(4,Niters);
sums = zeros(4,Niters+1);

for ii=1:Niters

    x1 = latticeseq_b2(totalN, Npoints);
    
    sums(:,ii+1)= sums(:,ii); 

    for jj=1:Npoints 
    
        % random points
        Y = x1(:,jj)-0.5;  
        
%         fileID = fopen('debug.txt','w'); 
%         fprintf(fileID,'%f \n',Y);
%         fprintf(fileID,'%d %d ',ii,jj);
%         fclose(fileID); 

        [coefx,coefy] =GetGeoCofs(coefx0,coefy0,Nx,Ny,dcoefs,Y);
        
        [u,V]= DirectSolverFull('s',false,coefx,coefy,omega,lambda,mu,alpha,Ndof,1e-10);

        sums(:,ii+1) =  sums(:,ii+1)+ObservedFunc('both',xtheta,wtheta,u,...
                coefx,coefy,omega,lambda,mu,Ndof);           

            
    end 
    
    Quads(:,ii) = sums(:,ii+1)/(Npoints*ii); 
    
end 

errors = zeros(4,Niters-1); 

for ii=1:(Niters-1)

    errors(:,ii) = abs(Quads(:,ii)-Quads(:,end));
    
end 

plot(log10(Npoints*(1:(Niters-1))),log10(errors(1,:)),'-x',...
log10(Npoints*(1:(Niters-1))),log10(errors(2,:)),'-x',...
log10(Npoints*(1:(Niters-1))),log10(errors(3,:)),'-x',...
log10(Npoints*(1:(Niters-1))),log10(errors(4,:)),'-x'); 