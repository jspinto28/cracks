function V = FarChebyshev(y,Nsol,Nc,bgeo,delta,k)

    Nobs = length(y(:,1)); 

    V = zeros(Nsol*Nobs,1); 
    
    xc = chebpts(Nc);
    
    const = exp(1i*pi/4)/sqrt(8*pi*k); 
    
    for s=1:Nobs
       
        r = GetGeo(xc,bgeo,delta); 
        
        r = bsxfun(@times,r,y(s,:));
        
        r = exp(-1i*k*sum(r,2)); 
        
        aux = zeros(2*Nc-2,1);
        
        aux(1) = r(Nc); 
        
        aux(Nc) = r(1); 
        
        aux(2:(Nc-1)) = flip(r(2:(Nc-1))); 
        
        aux((Nc+1):end) = r(2:(Nc-1)); 
        
        aux = fft(aux);
        
        r = aux(1:Nc)/(Nc-1); 
        
        r(1) = 0.5*r(1); 
        
        r(end) = 0.5*r(end); 
        
        r(1) = r(1)*pi; 
        
        r(2:end) = r(2:end)*pi/2; 
        
        V(((s-1)*Nsol+1):s*Nsol) = r(1:Nsol);         
    
    end
    
    V = V*const; 
    
    V = reshape( V, Nsol,Nobs).'; 
end