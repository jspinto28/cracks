function u=Un(n,t)

if(abs(t)==1)
   
    u=(n+1)*t^n;
    
    return;  
    
end

if(n == 0)
    
    u=1;
    
    return; 
    
else if(n==1) 
        
    u=2.*t; 

    return; 
end 

    u = 2*Un(n-1,t)-Un(n-2,t); 
    

end
