function V = GetRighHandSide(type,useFloat,Xcoefs,Ycoefs,omega,lambda,mu,alpha,Ndof)

    Narcs = length(Xcoefs); 

    V = zeros(2*Ndof*Narcs,1); 
    
    indx = 41;    
    
    if(type =='s')
    
        indx = 51; 
        
    end 
    
    for ii=1:Narcs       

        V((1:2*Ndof)+(ii-1)*2*Ndof) = ...
        ArcSolver(indx,useFloat,ii,length(Xcoefs{ii}),length(Ycoefs{ii}),Xcoefs{ii},Ycoefs{ii},...
            omega,lambda,mu,alpha,Ndof);             
        
    end 



end