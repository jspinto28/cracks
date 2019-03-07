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
 
 eps = 1e-8; 
 
 hgeo = [eps,0,0,0,0,0,0,0,0,0]; 
 
 [tu1,L1] = Solve(k,Nsol,alpha,bgeo+hgeo,delta); 
 
 B1 = Farfield(tu1,k,bgeo,delta,y);
 
 [tu0, L0] = Solve(k,Nsol,alpha,bgeo,delta);
 
 B0 = Farfield(tu0,k,bgeo,delta,y); 
 
 [f2,F2] = GetF2(k,alpha,tu0,L0,bgeo,delta,y);  
 
[norm((B1-B0)/eps ) ,  norm( (B1-B0)/eps - (F2)*(hgeo'/eps))]
 