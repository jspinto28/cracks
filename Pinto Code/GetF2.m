function [f2,F2] = GetF2(k,alpha,tu,L,b,delta,y)

    Nsol = length(tu); 

    L = inv(L); 
    
    FT = FarChebyshev(y,Nsol,2*Nsol+28,b,delta,k); 
    
    [v1,V1] = GetV1(k,Nsol,alpha,b,delta); 
    
    [v2,V2] = GetV2(k,tu,b,delta); 
    
    f2 = FT*L*(v1-v2); 
    
    F2 = FT*L*(V1-V2); 
    
end