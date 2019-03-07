function g = Obserbable(bt,delta, kt,angles,thetsObs, Nsol)
    %%Compute the Observable (Far field in this case) for a given geometry,
    % and solution. 

    Nb = length(bt);
    
    Na = length(angles); 
    
    Nc = 2*Nsol+28; 
    
    Nt = 2*Nsol+28; 
    
    Vst = mex_funs(109,kt,Nsol,Nc,Nt,... 
                      1,1,-1,bt,Nb,delta);
                  
    b = zeros(length(Vst(:,1)),Na);    
    
    for ii=1:Na
                  
        b(:,ii) =  mex_funs(106,kt,angles(ii),Nsol,0,1,1,bt,Nb,delta);   

    end
    
    tu = Vst \b; 
    
    %%%%
    %FarField: 
    
    Nobs = length(thetsObs); 
    
    N=length(tu(:,1)); 

    xobs = cos(thetsObs); 
    
    yobs = sin(thetsObs);  
    
    Fun =@(t,z) exp(-1i*kt*GetGeo(t,bt,delta)*z);
    
    Nc = 4*N+2*ceil(kt); 
    
    A = zeros(Nobs,N); 
    
    xc = chebpts(Nc); 
    
    for ii=1:Nobs
        
        v = Fun(xc,[xobs(ii);yobs(ii)]);
        
        aux = zeros(2*Nc-2,1); 
        
        aux(1) = v(Nc); 
        
        aux(Nc) = v(1); 
        
        aux(2:(Nc-1)) = flip(v(2:(Nc-1))); 
        
        aux((Nc+1):end) = v(2:(Nc-1)); 
        
        aux = fft(aux);
        
        v = aux(1:Nc)/(Nc-1); 
        
        v(1) = 0.5*v(1); 
        
        v(end) = 0.5*v(end); 
        
        A(ii,:) = v(1:N); 
        
    end 
    
    g = exp(1i*pi*0.25)/sqrt(8*pi*kt)*A*tu; 
    
    g = reshape(g,Nobs*Na,1);
    
end