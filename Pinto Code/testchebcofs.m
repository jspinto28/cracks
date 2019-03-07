Nc = 40; 

xc = chebpts(Nc); 

Tn =@(t,n) cos(n*acos(t)); 

A=zeros(10,10); 


    for ii=0:9
        
        v = Tn(xc,ii);
        
        aux = zeros(2*Nc-2,1); 
        
        aux(1) = v(Nc); 
        
        aux(Nc) = v(1); 
        
        aux(2:(Nc-1)) = flip(v(2:(Nc-1))); 
        
        aux((Nc+1):end) = v(2:(Nc-1)); 
        
        aux = fft(aux);
        
        v = aux(1:Nc)/(Nc-1); 
        
        v(1) = 0.5*v(1); 
        
        v(end) = 0.5*v(end);
        
        A(ii+1,:) = v(1:10); 
        
        
    end