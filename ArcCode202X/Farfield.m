function f = Farfield(type,u,theta,Xcoefs,Ycoefs,omega,lambda,mu,Ndof)

    f = zeros(2,1); 
    
    indx = 61; 
    
    if (type == 's')
    
        indx = 71; 
        
    end
    
    if(strcmp(type,'both'))
                 
         f = zeros(4,1);
         
    end 
    
    Narcs = length(Xcoefs); 
    
    for ii=1:Narcs
        
        f(1:2) = f(1:2) +   ...
           ArcSolver(indx,false,ii,length(Xcoefs{ii}),length(Ycoefs{ii}),Xcoefs{ii},Ycoefs{ii},...
                   omega,lambda,mu,theta,Ndof,...
                   real(u((1:2*Ndof)+2*Ndof*(ii-1))),...
                   imag(u((1:2*Ndof)+2*Ndof*(ii-1))));
        
    end 
    
    if(strcmp(type,'both'))
         
        indx = 71; 
        
        for ii=1:Narcs
        
            f(3:4) = f(3:4) +   ...
               ArcSolver(indx,false,ii,length(Xcoefs{ii}),length(Ycoefs{ii}),Xcoefs{ii},Ycoefs{ii},...
                       omega,lambda,mu,theta,Ndof,...
                       real(u((1:2*Ndof)+2*Ndof*(ii-1))),...
                       imag(u((1:2*Ndof)+2*Ndof*(ii-1))));
               
        end 
         
    end 



end