function z = GetGeo(t,bn,delta)

    Tn = @(n,t) cos(n*acos(t)); 
    
    x = t; 
  
    ns = (0:(length(bn)-1))'; 
    
    t = reshape(t,length(t),1); 
    
    bn = reshape(bn,length(bn),1); 
    
    v = 1./((1+ns((1:end))').^(2+delta));
    
    M = Tn(ns,t')'; 
    
    M = bsxfun(@times,M,v);
    
    y = M*bn; 
    
    z=[x,y]; 

end 