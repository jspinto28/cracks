function f = Farfield(tu,k,bgeo,delta,y)

    Nsol = length(tu); 
    
    Nobs = length(y(:,1)); 
    
    Nb = length(bgeo); 
    
    Fun =@(t,z) exp(-1i*k*GetGeo(t,bgeo,delta)*z);
    
    Nc = 2*Nsol+28; 
    
    A = zeros(Nobs,Nsol); 
    
    xc = chebpts(Nc); 
    
    for ii=1:Nobs
        
        v = Fun(xc,y(ii,:)');
        
        aux = zeros(2*Nc-2,1); 
        
        aux(1) = v(Nc); 
        
        aux(Nc) = v(1); 
        
        aux(2:(Nc-1)) = flip(v(2:(Nc-1))); 
        
        aux((Nc+1):end) = v(2:(Nc-1)); 
        
        aux = fft(aux);
        
        v = aux(1:Nc)/(Nc-1); 
        
        v(1) = 0.5*v(1); 
        
        v(1) = v(1)*pi;
        
        v(end) = 0.5*v(end); 
        
        v(2:end) = 0.5*pi*v(2:end); 
        
        A(ii,:) = v(1:Nsol); 
        
    end 
    
    f = exp(1i*pi*0.25)/sqrt(8*pi*k)*A*tu; 
    
    f = reshape(f,Nobs,1);


end