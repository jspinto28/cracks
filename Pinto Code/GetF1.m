function [f1,F1] = GetF1(k,tu,b,delta,y)

    Nsol = length(tu); 
    
    Nb = length(b); 
    
    Nobs= length(y(:,1)); 
    
    const = sqrt(k/(8*pi))*exp(-1i*pi/4); 
    
    f1=0;
    
%     f1 = const*mex_funs(207,k,Nsol,b,Nb,delta,1,0,[y(:,1);y(:,2)],Nobs)*tu;
    
    F1 = zeros(Nobs,Nb); 
    
    for p=0:(Nb-1)
        
        F1(:,p+1) = const*mex_funs(207,k,Nsol,b,Nb,delta,p,1,[y(:,1);y(:,2)],Nobs)*tu;
        
    end 



end