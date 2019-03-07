%maximum geometry basis 
Nb = 10; 

%Temporal stages
Tmax = 10; 

kt = [10,10,10,10,20,20,20,20,28,28]

anglet =[ 0, pi/4, pi/2, 3*pi/2,0, pi/4, pi/2, 3*pi/2,0,pi/2]; 

LosesExperts = zeros(Tmax,Nb);

Experts = eye(Nb); 

bt = [0.5,0.28,0,0.1,0,0,0.12,0,0,0];

eta = 1; 

bte = Experts(:,1); 

Ls =  zeros(Nb,1); 

for t =2: Tmax  
    
    for n=1:Nb
        
        Ls(n) = Ls(n) + loss(Experts(:,n),bt,kt(t-1),anglet(t-1)); 
        
    end 
    
    P = sum(exp(-eta*Ls)); 
    
    bte = exp(-eta*Ls)/P;     
    
end 
