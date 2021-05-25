function f = ObservedFunc(type,xtheta,wtheta,u,...
        Xcoefs,Ycoefs,omega,lambda,mu,Ndof)
        
    Nobs= length(xtheta); 
    
    f = zeros(2,1); 
    
    if(strcmp(type,'both'))
                 
         f = zeros(4,1);
         
    end 
    
    for ii=1:Nobs
        
        f = f + wtheta(ii)*Farfield(type,u,xtheta(ii),...
            Xcoefs,Ycoefs,omega,lambda,mu,Ndof);
        
    end 
    
    




end 