%%
%Generate a geometry contaning a open arc, given by the expansion parameter
% y1, y2
%
%           Smax-1
%           --
% R(t) =    \    (Y1(k))Tn(t)Cn 
%           /    (y2(k)) 
%           ----
%           k=0
%
% where Cn =  1/(n+1)**(2+delta)

function GenerateLine(Y1,Y2,delta)


Smax = length(Y1); 
cn = 1./(1:(Smax)).^(2+delta); 
dn = (0.000:(Smax-1))./(1:(Smax)).^(2+delta); 

fileID = fopen('Curves.cpp','w');

fprintf(fileID,'#include <math.h> \n');
fprintf(fileID,'#include "Curves.h" \n');
fprintf(fileID,'#include "prods.h" \n');
fprintf(fileID,'using namespace std; \n');
fprintf(fileID,'\n');
fprintf(fileID,'void curves(double* x, double* y,double t,int curv){ \n');
fprintf(fileID,'\n');    
fprintf(fileID,'    *x = \n'); 

for ii=0:(Smax-1)
   
    fprintf(fileID,'    %f*%f*chebT(%d,t)+ \n',Y1(ii+1),cn(ii+1),ii); 
end
fprintf(fileID,'    0; \n');
fprintf(fileID,'    *y = \n'); 

for ii=0:(Smax-1)
   
    fprintf(fileID,'    %f*%f*chebT(%d,t)+ \n',Y2(ii+1),cn(ii+1),ii); 
end
fprintf(fileID,'    0; \n');
fprintf(fileID,'}\n'); 


fprintf(fileID,'\n');

fprintf(fileID,'void tang_curves(double* x, double* y,double t,int curv){ \n');
fprintf(fileID,'\n');    

fprintf(fileID,'    *x = \n'); 

for ii=1:(Smax-1)
   
    fprintf(fileID,'    %f*%f*chebU(%d,t)+ \n',Y1(ii+1),dn(ii+1),(ii-1)); 
end
fprintf(fileID,'    0; \n');
fprintf(fileID,'    *y = \n'); 

for ii=1:(Smax-1)
   
    fprintf(fileID,'    %f*%f*chebU(%d,t)+ \n',Y2(ii+1),dn(ii+1),(ii-1)); 
end
fprintf(fileID,'    0; \n');
fprintf(fileID,'}\n'); 



fclose(fileID);
    
end