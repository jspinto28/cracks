function [tu,L] = Solve(k,Nsol,alpha,bgeo,delta)

    Nb = length(bgeo); 
        
    b =  mex_funs(106,k,alpha,Nsol,0,1,1,bgeo,Nb,delta); 
    
    Nc = 2*Nsol+28;
    
    L = mex_funs(109,k,Nsol,Nc,Nc,... 
                      1,1,-1,bgeo,Nb,delta);
    
    tu = L\b; 


end