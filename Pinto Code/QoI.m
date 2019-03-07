function q=QoI(b,dim,delta) 
% Quantity of interest, the result is a vector in (R**dim,2) 
% implemented as a evaluation of the geometry (one curve) in some points. 

xc = chebpts(dim); 

q = GetGeo(xc,b,delta); 

end