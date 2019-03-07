function [v1,V1] = GetV1(k,Nsol,alpha,b,delta)
   
    Nb = length(b); 
    
    v1 = mex_funs(206,k,alpha,[Nsol],0,1,1,b,Nb,delta,1,0); 
    
    V1 = zeros(Nsol,Nb); 
    
    for p=0:(Nb-1)
        
        V1(:,p+1) = mex_funs(206,k,alpha,[Nsol],0,1,1,b,Nb,delta,p,1); 
        
    end 
      
end 