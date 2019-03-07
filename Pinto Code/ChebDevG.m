function V = ChebDevG(Nsol,Nc,bgeo,delta,k,alpha,p,component)

    g =@(z) exp(-1i*k*sum(bsxfun(@times,z,[cos(alpha),sin(alpha)]),2));
    
    dg =@(z) [-1i*k*cos(alpha)*g(z),-1i*k*sin(alpha)*g(z)];
    
    Tn =@(t) cos(p*acos(t));

    xc = chebpts(Nc);
       
    r = GetGeo(xc,bgeo,delta); 
    
    r = dg(r);
    
    r = Tn(xc).*r(:,component+1);        
       
    aux = zeros(2*Nc-2,1);
        
    aux(1) = r(Nc); 
        
    aux(Nc) = r(1); 
        
    aux(2:(Nc-1)) = flip(r(2:(Nc-1))); 
        
    aux((Nc+1):end) = r(2:(Nc-1)); 
        
    aux = fft(aux);
        
    r = aux(1:Nc)/(Nc-1); 
        
    r(1) = 0.5*r(1); 
        
    r(end) = 0.5*r(end); 
        
    V = r(1:Nsol);         
    
    V(1) = V(1)*pi; 
    
    V(2:end) = V(2:end)*pi/2; 
    
end
