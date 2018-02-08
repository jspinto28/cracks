function z = lineGeo(Y1,Y2,delta,t)
    % see GenerateLine, this just do the same on matlab.
        
    z=zeros(length(t),2); 
    Smax = length(Y1); 
    cn = 1./(1:(Smax)).^(2+delta); 
    Tn=@(n,t) cos(n*acos(t)); 
    
    for i=1:length(t)
       
        z(i,1) = cn*(Y1.*Tn((0:(Smax-1))',t(i))); 
        
        z(i,2) = cn*(Y2.*Tn((0:(Smax-1))',t(i))); 
        
    end

end 