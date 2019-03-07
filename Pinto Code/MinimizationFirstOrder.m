clear all; 

mex mex_funs.cpp IntegrationRules.cpp Adom.cpp kers.cpp prods.cpp geo.cpp in_cells.cpp tchv_cof.cpp greenfunction.cpp CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" -I/home/matlab/boost_1_55_0 -lfftw3 -lm ;

% fixed: 

k= 1; 

alpha = pi\4; 

Nb = 10; 

delta = -2; 

Nobs = 28; 

Thetas = linspace(0,2*pi,Nobs+1); 

Thetas = Thetas(1:Nobs); 

y = [cos(Thetas'), sin(Thetas')]; 

MaxIter = 1024;

bgeo = cell(MaxIter,1); 

hgeo = cell(MaxIter,1);

Nsol = 28; 

step = 0.01; 

%initialization 

br = zeros(Nb,1); 

br(1) = 1; 

br(2) = 0.5; 
% 
br(3) = 0.28; 
% 
br(4) = 0.28; 
% 
% br(5) = 0.0; 

% bgeo{1} = zeros(Nb,1); 

% bgeo{1} =[       0.8624;   -0.0156;   -0.0633;    0.0667;    0.0121;   -0.1046;   -0.1709;    0.0791;   -0.3551;   -0.2144];


bgeo{1} = [        0.9544; 0.5584; 0.2577;-0.0643;-0.1487;-0.1044;0.0137;0.0164;0.0003;-0.0648];


[tur,Lr] = Solve(k,Nsol,alpha,br,delta);

f = Farfield(tur,k,br,delta,y); 

ts = linspace(-1,1,1000); 

zr = GetGeo(ts',br,delta); 

%first stage: 

[tu,L] = Solve(k,Nsol,alpha,bgeo{1},delta); 

f0 = Farfield(tu,k,bgeo{1},delta,y); 

[f1,F1] = GetF1(k,tu,bgeo{1},delta,y); 

[f2,F2] = GetF2(k,alpha,tu,L,bgeo{1},delta,y); 

F3 = F1+F2; 

hgeo{1} = -2*real(F3'*(f0-f));

hgeo{1} = step*hgeo{1}/norm(hgeo{1});

iters = 2; 

error = norm(f-f0)/norm(f0); 

while( (iters <= MaxIter) )
    
    bgeo{iters} = bgeo{iters-1} + hgeo{iters-1};
    
%     error = norm( hgeo{iters-1});
    
    if( error < 1e-12) 
        
        break; 
        
    end 
    
    [tu,L] = Solve(k,Nsol,alpha,bgeo{iters},delta);
    
    f0 = Farfield(tu,k,bgeo{iters},delta,y); 
    
    [f1,F1] = GetF1(k,tu,bgeo{iters},delta,y); 

    [f2,F2] = GetF2(k,alpha,tu,L,bgeo{iters},delta,y); 
    
    error = norm(f-f0)/norm(f0);
    
    [iters,error]

    F3 = F1+F2; 

    hgeo{iters} = -2*real(F3'*(f0-f));

    hgeo{iters} = step*hgeo{iters}*norm(bgeo{iters})/norm(hgeo{iters});
    
    zg = GetGeo(ts',bgeo{iters},delta);     
     
    iters = iters +1 ; 
    
end 

