function [u,V] = DirectSolverFull(type,useFloat,Xcoefs,Ycoefs,omega,lambda,mu,alpha,Ndof,tolerance,level)

    if ~exist('tolerance','var')        
        tolerance = 0; 
        level = 3; 
        
    end

    if ~exist('level','var')
        level = 3;         
    end

    Narcs = length(Xcoefs); 
        
    V = sparse(2*Ndof*Narcs,2*Ndof*Narcs); 
    
    g = 0; 
    
    if(strcmp(type,'both'))
        
        g= zeros(2*Ndof*Narcs,2);
        
        g(:,1) = GetRighHandSide('p',useFloat,Xcoefs,Ycoefs,omega,lambda,mu,alpha,Ndof);
        g(:,2) = GetRighHandSide('s',useFloat,Xcoefs,Ycoefs,omega,lambda,mu,alpha,Ndof);
        
    else 
        
        g = GetRighHandSide(type,useFloat,Xcoefs,Ycoefs,omega,lambda,mu,alpha,Ndof);
        
    end        
    
    
    for ii=1:Narcs
        for jj=1:ii
            
%             fileID = fopen('debug.txt','w'); 
%             fprintf(fileID,'%d %d ',ii,jj);
%             fclose(fileID); 
            
        	V((1:2*Ndof)+(ii-1)*2*Ndof,(1:2*Ndof)+(jj-1)*2*Ndof) = ...
               sparse( ArcSolver(21,useFloat,jj,length(Xcoefs{jj}),length(Ycoefs{jj}),Xcoefs{jj},Ycoefs{jj},...
                ii,length(Xcoefs{ii}),length(Ycoefs{ii}),Xcoefs{ii},Ycoefs{ii},...
                omega,lambda,mu,Ndof,8,tolerance,level));
			
            if(ii > jj)

                V((1:2*Ndof)+(jj-1)*2*Ndof,(1:2*Ndof)+(ii-1)*2*Ndof) =...
                V((1:2*Ndof)+(ii-1)*2*Ndof,(1:2*Ndof)+(jj-1)*2*Ndof);
            end    
            
        end         
    end 
    
    u = V\g;      
    

end