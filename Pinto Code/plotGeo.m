function plotGeo(dom_op,ts,an,bn)

    Nplot = length(ts); 
    
    NPer = length(an); 
    
    xs = zeros(Nplot,1); 
    
    ys = zeros(Nplot,1); 
    
    xs0 = zeros(Nplot,1); 
    
    ys0 = zeros(Nplot,1); 
    
    for i=1:Nplot
        
        z = mex_funs(21,ts(i),dom_op,NPer,an,bn); 
        
        xs(i)=z(1); 
        
        ys(i)=z(2);      
        
        z = mex_funs(21,ts(i),dom_op,1,0,0); 
        
        xs0(i)=z(1); 
        
        ys0(i)=z(2);   
        
    end 
    
     figure; 
    
    plot(xs,ys,'r');
    
    hold on; 
    
    plot(xs0,ys0,'-.'); 

end