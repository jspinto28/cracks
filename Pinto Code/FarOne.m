function V = FarOne(y,Nsol,Nc,bgeo,delta,k,p,component)

    Tn = @(t) cos(p*acos(t)); 

    Nobs = length(y(:,1)); 

    V = zeros(Nobs*Nsol,1); 
    
    xc = chebpts(Nc);
    
    for s=1:Nobs
       
        r = GetGeo(xc,bgeo,delta); 
        
        r = bsxfun(@times,r,y(s,:));
        
        r = exp(-1i*k*sum(r,2)).*Tn(xc); 
        
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
        
        r(2:end) = r(2:end)*pi*0.5; 
        
        V(((s-1)*Nsol+1):s*Nsol) = r(1:Nsol)*y(s,component+1);         
    
    end
end