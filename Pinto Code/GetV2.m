function [v2,V2] = GetV2(k,tu,b,delta)

    Nsol = length(tu); 
    
    Nb = length(b); 
    
    v2 = mex_funs(209,k,Nsol,2*Nsol+28,2*Nsol+28,0,0,-1,b,Nb,delta,1,0)*tu;
    
    V2 = zeros(Nsol,Nb);
    
    for p=0:(Nb-1)

        V2(:,p+1) = mex_funs(209,k,Nsol,2*Nsol+28,2*Nsol+28,0,0,-1,b,Nb,delta,p,1)*tu;
                    
    end


end 