function g = Obserbable(Y1,Y2,delta, kArray, angleArray, Nobs)
    %%Compute the Observable (Far field in this case) for a given geometry
    % given by Y1,Y2,delta.
    % The obserbable correspond to point evaluation of the far field,
    % observed in Nobs points (equally spaces in the unitary circle). 
    % We construct a matrix of obserbables by looping in al the freqs in
    % kArray, and all the angles in angleArray. The incident field in each
    % case is consider as a plane wave (notices that we have to assume that
    % k >0). 
    
    %for now fixed.
    Ncurves = 1;
    
    S =length(Y1); 
    
    domt=Ncurves;
    
%     GenerateLine(Y1,Y2,delta);     
       
    %use compression:
    fast = true; 
    
    Tu = cell( length(kArray),1);     
    
    %solve:
    for ii=1:length(kArray) 
        
        %matrix setup
        k_0=kArray(ii);       
        
        N = ceil(k_0/10 +10); 
        
        Nc=128*2+2*N;
        
        Nq=128*2+2*N;
        
        [xq,wq] = chebyshev1_rule ( Nq(1), -1, 1);      
        
        V=sparse((2*N+1)*domt,(2*N+1)*domt);  
        
        Nshif = 2*N+1; 
        
        if(fast)
        
            shiftt =0;      
        
            for tt=1:domt(1)
            
                shifss =0;        
            
        
                for ss=1:domt(1)                   
                                                
                    Vst = mex_funs(107,k_0,N,Nc,Nq,...
                      tt-1,ss-1,xq,wq,Y1,Y2,delta,S); 
                  
                     V(1+shifss:Nshif+shifss,...
                        1+shiftt:Nshif+shiftt) = Vst; 
%                   add transpose???
%                    if(ss ~= tt)                    
%                       V(1+shifss:Ntt+shifss,
%                           ...
%                           1+shiftt:Ntt+shiftt) = Vst; 
%                    end 
                
                shifss = shifss+ Nshif;            
            
                end 
            
            shiftt = shiftt+(Nshif); 
            
            end
            
        else
            
            V=mex_funs(105,k_0,N,Nc(1),Nq(1),0,domt,xq,wq,Y1,Y2,delta,S);
            
        end 
        
        Bs = zeros(2*N+1,length(angleArray)); 

        for jj =1:length(angleArray)
            
            %%rhs:          
            
            Bs(:,jj) =  mex_funs(106,N,angleArray(jj),k_0,domt,Y1,Y2,delta,S);     
                        
            
        end 
        
        Tu{ii}= V\Bs; 
        
    end 
    
 
    
    %%compute far fields: 
        
    thets = linspace(0,2*pi,Nobs+1); 
    
    thets = thets(1:Nobs);
    
    xobs = cos(thets); 
    
    yobs = sin(thets);  
    
    g = zeros(length(kArray) ,length(angleArray), Nobs); 
    
    Tn =@(n,t) cos(n*acos(t)); 
    
    for mm=1:Nobs
        
        zobs = [xobs(mm);yobs(mm)]; 
        
        for kk=1:length(kArray) 
            
            kext = kArray(kk); 
            
            N = (kext/10 +10);
            
            M= zeros(2*N+1,4*N);
     
            [xq,wq] = chebyshev1_rule ( 4*N, -1, 1);              
    
            for ii=1:(2*N+1);

                M(ii,:) = Tn((ii-1),xq)'; 

            end     
            
            ExpVec = zeros(4*N,domt); 
            
            for tt=1:domt

                for jj=1:4*N
                    
                    
                 %check numeration from 0 or 1.
                 z = lineGeo(Y1,Y2,delta,xq(jj));              

                 ExpVec(jj,tt) = ...
                     exp(-1i*kext*z*zobs);

                end 
            end 
            
            A = M*ExpVec;
            
            A = reshape(A,[(2*N+1)*domt,1]); 
                                    
            for aa=1:length(angleArray)                
       
                g(kk,aa,mm) = wq(1)*exp(1i*0.25*pi)...
                    /sqrt(8*pi*kext)*(conj(Tu{kk}(:,aa)')*A);
                
            end 
            
        end 
        
        
    end 
    
    g = reshape(g,[Nobs*length(kArray)*length(angleArray),1]); 
  


end