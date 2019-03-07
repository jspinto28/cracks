 Nb = 10; 
 
 k=10; 
 
 alpha = pi/3; 
 
 delta = -2;
 
 Nsol = 128;
 
 tu = zeros(Nsol,1) +1; 
 
 Nobs = 128; 

 Thetas = linspace(0,2*pi,Nobs+1); 

 Thetas = Thetas(1:Nobs); 

 y = [cos(Thetas'), sin(Thetas')]; 
 
 bgeo = [1,1,0.5,0.5,0.28,0.28,0,0,0,0]; 
 
 eps = 1e-4; 
 
 hgeo = [eps,0,0,0,0,0,0,0,0,0]; 
 
 B1 = Farfield(tu,k,bgeo+hgeo,delta,y); 
 
 B0 = Farfield(tu,k,bgeo,delta,y); 
 
 [f1,F1] = GetF1(k,tu,bgeo,delta,y); 
 
[norm((B1-B0)/eps ) ,  norm( (B1-B0)/eps - (F1)*(hgeo'/eps))]
 