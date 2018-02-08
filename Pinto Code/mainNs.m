clear all

% BASIC PARAMETERS

%permitivity array
ep=zeros(1,2)+1; 

k_0=28; 
Ncurves = 50; 
fmin=pi; 
fmax=2*pi; 
ampMin=1; 
ampMax=2; 
lmin=0.28;
lmax=0.5;

% [freqs,xmin,xmax,ymin,ymax,anchs,amps,phases,xs,ys] = ...
%     GenerateGeo(Ncurves,fmin,fmax,ampMin,ampMax,lmin,lmax);

[xmin,xmax,ymin,ymax]...
    = GenerateSpirals(Ncurves,fmin,fmax,ampMin,ampMax);

mex mex_funs.cpp IntegrationRules.cpp Adom.cpp kers.cpp prods.cpp geo.cpp Curves.cpp in_cells.cpp tchv_cof.cpp greenfunction.cpp CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" -I/home/matlab/boost_1_55_0 -lfftw3 -lm ;


% [xs2,ys2,zs]=PlotField(Ncurves,xmin,xmax,ymin,ymax,anchs,amps,freqs,phases,...
%                     xs,ys,tu{1},k(1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%exterior epsilon.
% ep(1)=2; 

% ADVANCED PARAMETERS (dont touch)
DOMS=0; 
N0 = 0;
Nc=zeros(1,DOMS+1)+128*10;
Nq=zeros(1,DOMS+1)+128*10;

[xq,wq] = chebyshev1_rule ( Nq(1), -1, 1);



domt=zeros(1,DOMS+1)+Ncurves;
ca=1; 
ck=1;

%changues in N_0 is aplied to all the domain and interfaces
N_p=[28*2]; 
cN=length(N_p);

% ercomp = zeros(cN,1); 
% 
% patern = zeros(cN,1);

itersFull = zeros(cN,1); 

itersSparse = zeros(cN,1);

itersFullLh = zeros(cN,1); 

itersSparseLh = zeros(cN,1);

itersFullPh = zeros(cN,1); 

itersSparsePh = zeros(cN,1);

tu=cell(cN,1);

tud=cell(cN,1);
%changues in k_0, will also changue the other wave numbers. 
k_p=[]; 

k_p=[0,k_p];

dom_op = 0; 

indx=0;

fast = true; 

for n=1:cN
    
    indx2erase = zeros(domt(1),1);
    
    N=cell(1,DOMS+1); 
    
    k=zeros(1,DOMS+1); 

    for d=0:DOMS
        k(d+1)=(k_0+k_p(1))*ep(d+1);
    end
    
    Mk=max(k);
    
    N=cell(1,DOMS+1); 
    
    for d=0:DOMS
        
        N{1}=zeros(1,domt(1))+N_p(n);
%         N{d+1}=N_p(n); 
    end
    
    s=0; 
    for i=1:domt(1) 
        
        indx2erase(i) = 1+(i-1)*(2*N{1}(i)+1); 

    end   
    
    V=sparse(2*sum(N{1}),2*sum(N{1})); 
    
    if(fast)
        
        shiftt =0;      
        
        for tt=1:domt(1)
            
            shifss =0; 
            
            if (k(1) < 1e-6) 
            
                Ntt = 2*N{1}(tt);
           
            else
                
                Ntt = 2*N{1}(tt)+1;
                
            end             
            
        
            for ss=1:domt(1)
                
                if (k(1) < 1e-6) 
                
                    Nss = 2*N{1}(ss);
                else
                                    
                    Nss = 2*N{1}(ss)+1;
                
                end   
                    
                                                
                Vst = mex_funs(107,k(1),N{1},Nc(1),Nq(1),...
                      tt-1,ss-1,xq,wq);    
                  
                V(1+shifss:Nss+shifss,...
                    1+shiftt:Ntt+shiftt) = Vst; 
%            add transpose???
%                 if(ss ~= tt)                    
%                     V(1+shifss:Ntt+shifss,
...
%                         1+shiftt:Ntt+shiftt) = Vst; 
%                 end 
                
                shifss = shifss+ Nss;            
            
            end 
            
            shiftt = shiftt+(Ntt); 
            
        end
        
    else
        
        V=mex_funs(105,k(1),N{1},Nc(1),Nq(1),dom_op,domt(1),xq,wq);         
        
    end   
    
%     Vd=mex_funs(105,k(1),N{1},Nc(1),Nq(1),dom_op,domt(1),xq,wq);  
    
    b =  mex_funs(106,N{1},dom_op,domt(1),domt(1));
        
    if (k(1) < 1e-6) 
    
        indx = setdiff(1:length(b),indx2erase);
        
%         if (~fast)
% 
%             V = V(indx,indx); 
%             
%         else 
%             
%             Vd = Vd(indx,indx); 
%             
%         end

        b=b(indx); 
    
    end
    
%     Lh =Lhat(k,N);
%     
%     Ph = sparse(diag(diag(V))); 
    
%     [x,flag,relres,iter]=gmres(Vd,b,[],1e-8,length(b));
%     
%     itersFull(n)= iter(2); 
%     
%     [x,flag,relres,iter]=gmres(V,b,[],1e-8,length(b));
% 
%     itersSparse(n) = iter(2); 
    
%     [x,flag,relres,iter]=gmres(Vd,b,[],1e-8,length(b),Lh);
%     
%     itersFullLh(n) = iter(2); 
    
%     [x,flag,relres,iter]=gmres(V,b,[],1e-8,length(b),Lh);
%     
%     itersSparseLh(n) = iter(2); 
 
%     [x,flag,relres,iter]=gmres(Vd,b,[],1e-8,length(b),Ph);
%     
%     itersFullPh(n) = iter(2); 
    
%     [x,flag,relres,iter]=gmres(V,b,[],1e-8,length(b),Ph);
%     
%     itersSparsePh(n) = iter(2); 
%    
    
     tu{n} = V\(b); 
%     
%     tud{d} = Vd\b; 
%     
%     ercomp(n) = norm(tu{n}-tud)/norm(tud); 
%     
%     patern(n)  = 100*nnz(V)/(length(V(:,1))^2);

    
%     tuTest0 = zeros(length(b),1); 
% 
%     tuTest0(1) = 2/log(2)*b(1)/pi; 
% 
%     m = length(b)-1;
% 
%     tuTest0(2:end) = 4*(1:m)'.*b(2:end)/pi; 
%                              printf("count %d  \n",count); 
                            
%     norm(tu{n}-tuTest0);


end 

NpC = length(tu{1})/Ncurves;

[xg,wg] = chebyshev1_rule(NpC+10,-1,1);

tic;
% 
 xmin = xmin-15.0;
% 
 xmax = xmax+15.0; 
% 
 ymin=ymin-15; 
% 
 ymax=ymax+15; 

 Ps=mex_funs(108,Ncurves,xmin,xmax,ymin,ymax,length(xg),xg,wg,tu{1},NpC,k);

 toc
 

 
PlotField2(Ps,Ncurves);


%  PlotTrace(Ps,1,Ncurves,tu{1},NpC,k);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% V0=mex_funs(105,0,N{1},Nc(1),Nq(1),dom_op,domt(1),xq,wq);

% ft =@(n,t) cos(n*acos(t));
% 
% M0 =zeros(length(b)); 
% 
% dtt=length(b)/domt(1);

% shift=1; 

% for ii=1:dtt
%     
%    shift=shift+(dtt-1)*2*N{1}(i)
%    
%    for m=1:dtt
% 
%     for l=1:dtt
%         
%         fun = @(t) ft(m,t).*ft(l,t); 
%         
%         M0(l,m) = wq'*fun(xq);
%         
%     end 
%     
%   end 
%     
%     
% % end
% 
% plot error overkill:
% 
% if (k(1) < 1e-6) 
%     
%     V0 = V0(indx,indx);
%     
% end
% 
% V0=diag(diag(V0)); 
% 
% error = zeros(cN-1,1); 
% 
% errorRel = zeros(cN-1,1); 
% 
% overkillnorm = sqrt(tu{end}'*V0*tu{end});
% 
% errorl2 = zeros(cN-1,1); 
% 
% for n=1:cN-1
%     
%     tuaux = zeros(length(tu{end}),1); 
%     
%     Ndm = length(tu{n})/domt(1); 
%     
%     Nd0 = length(tu{end})/domt(1); 
%     
% %     tuaux(1:length(tu{n})) = tu{n}; 
% 
%     for i=1:domt(1)
%         
%         tuaux(1+(i-1)*Nd0:(i-1)*Nd0+Ndm) = ...
%             tu{n}(1+(i-1)*Ndm:(i-1)*Ndm+Ndm);
%         
%     end 
%     
%     error(n) = sqrt((tu{end}-tuaux)'*V0*(tu{end}-tuaux));
%     
%     errorRel(n) =error(n)/overkillnorm; 
%     
% %     errorl2(n) = sqrt((tu{end}-tuaux)'*M0*(tu{end}-tuaux));
% %     
% end 
% format shortE
% % 
% % 
% % [2*N_p(1:(end-1))'+1,abs(error),abs(errorRel)]
% 
% % coefficients = polyfit(2*N_p((end-2):(end-1))'+1, abs(errorl2((end-1):end)), 1);
% % slope = coefficients(1);
% % yl2= slope*(2*N_p((end-2):(end-1))+1)+coefficients(2);
% 
%  semilogy(2*N_p(1:(end-1))+1,abs(error),'b');
%  leg= legend('$\tilde{H}^{-1/2}$ -error');
%  set(leg,'Interpreter','latex');
%  set(leg,'FontSize',14);
%  xlabel('N');
%  ylabel('Error');

