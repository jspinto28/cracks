function q=QoI(Y1,Y2,delta,dim) 
% Quantity of interest, the result is a vector in (R**dim,2) 
% implemented as a evaluation of the geometry (one curve) in some points. 

xc = chebpts(dim); 

q = lineGeo(Y1,Y2,delta,xc); 


end