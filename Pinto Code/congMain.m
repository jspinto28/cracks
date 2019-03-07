clear all

tic;

% BASIC PARAMETERS

%permitivity array
ep=zeros(1,2)+1; 
% ks=[0,10,25,30,35,40,45,50,52,54,56,58,60,70,75,78,80,90,92,95,97,100,110,115,120,128,150,151,155,160,164,170,175,180,190,200,205,210]; 
% ks=[0,10];
%   ks=10; 
%  k_0=[20]; 
%ks=[0,25,50,100];
    ks=100;
% ks=100; 

% ks=25;

% ks=0; 


 Ncurves = 40; 
 fmin=pi/3; 
 fmax=pi/2; 
 ampMin=1;
 ampMax=1.5; 
 lmin=0.28;
 lmax=0.5;

%  [freqs,xmin,xmax,ymin,ymax,anchs,amps,phases,xs,ys] = ...
%      GenerateGeo(Ncurves,fmin,fmax,ampMin,ampMax,lmin,lmax);

% [xmin,xmax,ymin,ymax]...
%     = GenerateSpirals(Ncurves,fmin,fmax,ampMin,ampMax);

 mex mex_funs.cpp IntegrationRules.cpp Adom.cpp kers.cpp prods.cpp geo.cpp Curves.cpp in_cells.cpp tchv_cof.cpp greenfunction.cpp CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" -I/home/matlab/boost_1_55_0 -lfftw3 -lm ;


% [xs2,ys2,zs]=PlotField(Ncurves,xmin,xmax,ymin,ymax,anchs,amps,freqs,phases,...
%                     xs,ys,tu{1},k(1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


DOMS=0; 
N0 = 0;


domt=zeros(1,DOMS+1)+Ncurves;
ca=1; 
ck=1;

%changues in N_0 is aplied to all the domain and interfaces

Ns=[2,5,10,20,40,60,80,100,105,110,115,120,125,130,135,140,160,180,200,220,240,260,280,300,340,400,460,480];

% Ns=[Ns;Ns];

% Ns = [2,5,10,20,40,60,80];

% Ns=[20]; 
    
tol = -1;

errorsOk  = zeros(length(ks),length(Ns)-1); 

normsOk = zeros(length(ks),length(Ns)-1); 

% Ns=[Ns;Ns;Ns;Ns];

for jj=1:length(ks)
    
k_0 = ks(jj);

N_p =Ns(1,:); 



% Ns(jj,:) = N_p;
%  
 
cN=length(N_p);

norms = zeros(cN,1);
% ercomp = zeros(cN,1); 
% 
% patern = zeros(cN,1);

itersFull = zeros(cN,1); 

itersSparse = zeros(cN,1);

itersFullLh = zeros(cN,1); 

itersSparseLh = zeros(cN,1);

itersFullPh = zeros(cN,1); 

itersSparsePh = zeros(cN,1);
 
k_p=[]; 

k_p=[0,k_p];

dom_op = 0; 

indx=0;

fast = true; 

    
    NN = N_p(end); 
    
    Nc=zeros(1,DOMS+1)+2*NN+128;
  
    Nt=2*NN+128; 
    
    Nq = NN; 
    
    [xq,wq] = chebyshev1_rule ( Nq, -1, 1);
    
    Nq0 = ceil(Nq-1); 
    
    [xq0,wq0] = chebyshev1_rule ( Nq0, -1, 1);

    indx2erase = zeros(domt(1),1);
    
    N=cell(1,DOMS+1); 
    
    k=zeros(1,DOMS+1); 

    for d=0:DOMS
        k(d+1)=(k_0+k_p(1))*ep(d+1);
    end
    
    Mk=max(k);
    
    N=cell(1,DOMS+1); 
          
    for d=0:DOMS
        
        N{1}=zeros(1,domt(1))+NN;
 
    end
   
     [N{1},ks(jj)]
    
    s=0; 
    for i=1:domt(1) 
        
        indx2erase(i) = 1+(i-1)*(N{1}(i)); 

    end          
    
    if(k(1) < 1e-6) 
        NT= sum(N{1})-Ncurves;
    else 
        NT= sum(N{1});
    end 
    
%     V=zeros(NT,NT); 
%     
%     Vp =zeros(NT); 


        V=zeros(NT); 
%     
%         Vp =sparse(NT,NT); 

        shiftt =0;      
        
        for tt=1:domt(1)

            shifss =0; 
            
            Ntt = N{1}(tt);
            
            if (k(1) < 1e-6) 
            
                Ntt = N{1}(tt)-1;                              
                
            end             
            
        
            for ss=1:tt
                
                if (k(1) < 1e-6) 
                
                    Nss = N{1}(ss)-1;
                else
                                    
                    Nss = N{1}(ss);
                
                end
                
                if(ss == tt)
% 
                Vst = mex_funs(109,k(1),N{1}(1),Nc(1),Nt,... 
                      tt-1,ss-1,-1);
                  
                V(1+shifss:Nss+shifss,...
                    1+shiftt:Ntt+shiftt) = Vst;    
                
%                 Vp(1+shifss:Nss+shifss,...
%                     1+shiftt:Ntt+shiftt) = inv(Vst); 
                
                else
                    
                    Vst=0;                    
%                     
                    if(tol >0)
                    
%                         Vst = mex_funs(110,k(1),N{1}(1),tt-1,ss-1,Nq(1),xq,wq,... 
%                           tol,3);

                          Vst = mex_funs(112,k(1),N{1}(1),tt-1,ss-1,Nq(1),xq,wq(1),... 
                             Nq0,xq0,wq0,tol,3);

                    else 
                  
                      Vst = mex_funs(109,k(1),N{1}(1),Nc(1),Nt,... 
                          tt-1,ss-1,tol);
                  
                    end
                    
                   V(1+shifss:Nss+shifss,...
                    1+shiftt:Ntt+shiftt) = Vst;  
                    
                    V(1+shiftt:Ntt+shiftt,...
                        1+shifss:Nss+shifss) = Vst.'; 
         
                end 
                
                shifss = shifss+ Nss; 
            
            end 
            
            shiftt = shiftt+(Ntt); 
            
        end

    b =  mex_funs(106,k(1),N{1},dom_op,domt(1),domt(1));
        
    if (k(1) < 1e-6) 
    
         indx = setdiff(1:length(b),indx2erase);

        b=b(indx); 
    
    end
    
    toc
    
    tic;
    
    tu= V\b; 

    toc
%     Vpp= V*Vp; 
%     
%     tu = gmres(Vpp,b,(length(Vp(:,1))+1),1e-15,[]); 
%     
%     tu = Vp*tu; 



for n=1:(length(N_p)-1) 
    
    Np = N_p(n); 
    
    if( k(1) < 1e-6) 
        
        Np = Np-1; 
    end 
    
    Vn = zeros(Ncurves*Np);
    
    Vnp = zeros(Ncurves*Np);
    
    bn = zeros(Ncurves*Np,1);
    
    Nend = length(b)/Ncurves;
    
    for m=1:Ncurves
       
        for l=1:Ncurves
        
            Vn((m-1)*Np+1:m*Np,(l-1)*Np+1:l*Np) = ...
                V((m-1)*Nend+1:(m-1)*Nend+Np,(l-1)*Nend+1:(l-1)*Nend+Np); 

%             Vn((m-1)*Np+1:m*Np,(l-1)*Np+1:l*Np) = ...
%                 Vpp((m-1)*Nend+1:(m-1)*Nend+Np,(l-1)*Nend+1:(l-1)*Nend+Np); 
%             
%             Vnp((m-1)*Np+1:m*Np,(l-1)*Np+1:l*Np) = ...
%                 Vp((m-1)*Nend+1:(m-1)*Nend+Np,(l-1)*Nend+1:(l-1)*Nend+Np); 
%             
        end
        
%         bn((m-1)*Np+1:m*Np) = b((m-1)*Nend+1:(m-1)*Nend+Np);  

          bn((m-1)*Np+1:m*Np) = b((m-1)*Nend+1:(m-1)*Nend+Np);  
        
    end
    
    tic; 

     tun = Vn\bn; 

     toc
%     tun = gmres(Vn,bn,(length(Vn(:,1))+1),1e-15,[]);
%     
%     tun = Vnp*tun; 
%        
    MperCurve = N_p(end); 
    
    NperCurve = N_p(n); 
    
    if(k(1) < 1e-6)
        
        MperCurve = MperCurve-1; 
        
        NperCurve = NperCurve-1; 
        
    end 
    
    err = 0 ; 
    
    nor = 0; 
    
    if(k(1)<1e-6)
    
        for ii=1:Ncurves

            for m=1:NperCurve

                err = err + abs(tun(m+(ii-1)*NperCurve)-...
                    tu(m+(ii-1)*MperCurve))^2/m;
                
                nor = nor + abs(tu(m+(ii-1)*MperCurve))^2/m;

            end
            
            for m=(NperCurve+1):(MperCurve)
                
                err = err + abs(tu(m+(ii-1)*MperCurve))^2/m; 
                
                 nor = nor + abs(tu(m+(ii-1)*MperCurve))^2/m;
                
            end 


        end
    else 
        
        for ii=1:Ncurves
            
            err = err + abs(tun(1+(ii-1)*NperCurve)-...
                    tu(1+(ii-1)*MperCurve))^2;
                
            nor = nor + abs(tu(1+(ii-1)*MperCurve))^2;

            for m=2:NperCurve

                err = err + abs(tun(m+(ii-1)*NperCurve)-...
                    tu(m+(ii-1)*MperCurve))^2/(m-1);
                
                nor = nor + abs(tu(m+(ii-1)*MperCurve))^2/(m-1);

            end
            
            for m=(NperCurve+1):(MperCurve)
                
                err = err + abs(tu(m+(ii-1)*MperCurve))^2/(m-1); 
                
                 nor = nor + abs(tu(m+(ii-1)*MperCurve))^2/(m-1);
                
            end 


        end
        
        
    end
    
    errorsOk(jj,n)  =  err; 
    
    normsOk(jj,n) = nor; 
      
    
%     
 end 

end 

figure; 
hold on;
for jj=1:length(ks)
           
    plot((Ns(jj,1:(end-1))),log10(sqrt(errorsOk(jj,1:end))),'-x'); 
    
end 

hold off;


