function [Xcoefs,Ycoefs] = GetGeoCofs(coefx0,coefy0,Nx,Ny,dcoefs,Y)

    Narcs = length(Nx); 
    Xcoefs = cell(Narcs,1); 
    Ycoefs = cell(Narcs,1); 

    
    for ii=1:Narcs
        
        Xcoefs{ii} = zeros(Nx(ii),1); 
        Ycoefs{ii} = zeros(Ny(ii),1);     
        
    end 
    

    shiftD = 0; 
    for ii=1:max(max(Nx),max(Ny))
       
        NarcsX = sum(Nx >= ii); 
        shiftX = shiftD+1; 
        shiftY = NarcsX+shiftD+1; 
        Totalii = 0; 
        Totalx = 0; 
        Totaly = 0; 
        
        for jj=1:Narcs
            
            if(Nx(jj) >= ii)
            
                Xcoefs{jj}(ii) = Y(shiftX+Totalx);
                Totalii = Totalii+1; 
                Totalx = Totalx + 1; 

                
            end 
            
            if(Ny(jj) >= ii)
                
                Ycoefs{jj}(ii) = Y(shiftY+Totaly);
                Totalii = Totalii+1; 
                Totaly = Totaly + 1; 

                
                
            end
            
        end 
        
        shiftD = shiftD+ Totalii; 
        
        
    end    
        
    for ii=1:Narcs             

        for jj=1:floor(length(Xcoefs{ii})/2)
            
            Xcoefs{ii}(2*jj-1) = Xcoefs{ii}(2*jj-1)* dcoefs(jj); 
            Xcoefs{ii}(2*jj) =  Xcoefs{ii}(2*jj)*dcoefs(jj); 
            
        end
        
        if(mod(length( Xcoefs{ii}),2))
           
            Xcoefs{ii}(end)  = Xcoefs{ii}(end) *...
                dcoefs(floor(length(Xcoefs{ii})/2)+1);
            
        end
        
        for jj=1:floor(length(  Ycoefs{ii})/2)
            
            Ycoefs{ii}(2*jj-1) = Ycoefs{ii}(2*jj-1)* dcoefs(jj); 
            Ycoefs{ii}(2*jj) =  Ycoefs{ii}(2*jj)*dcoefs(jj); 
            
        end
        
        if(mod(length( Ycoefs{ii}),2))
           
            Ycoefs{ii}(end)  = Ycoefs{ii}(end) *...
                dcoefs(floor(length( Xcoefs{ii})/2)+2);
            
        end
        
        Xcoefs{ii} = [coefx0{ii}; Xcoefs{ii}]; 
        Ycoefs{ii} = [coefy0{ii}; Ycoefs{ii}]; 
        
    end  




end 