r1 =@(t) t; 

r1p = @(t) zeros(size(t))+1;

r2 =@(t) 2*t.^2+t; 

r2p =@(t) 4*t+1; 

r=@(t) [r1(t),r2(t)]; 

rp =@(t) [r1p(t),r2p(t)]; 

rr=@(t,s) (bsxfun(@minus,r(t),r(s))).^2; 

j= @(t) sqrt(sum(rp(t).^2,2)); 

d=@(t,s) sqrt(sum(rr(t,s),2)); 

daux =@(t,s) abs(bsxfun(@minus,t,s));

k=10; 

Tn =@(n,t) cos(n*acos(t)); 

Un =@(n,t) sin((n+1)*acos(t))./sqrt(1-t.*t);

p= 4; 

Fsing1 =@(t,s) -1i*k*0.25*besselh(1,k*d(t,s)).*(r1(t)-r1(s)).*(Tn(p,t)-Tn(p,s))./d(t,s);

Freg1 = @(t,s) Fsing1(t,s) -k*0.5/pi*besselj(1,k*d(t,s)).*(r1(t)-r1(s)).*(Tn(p,t)-Tn(p,s))./d(t,s).*log(daux(t,s));

Reg1 =@(t,s) k*besselj(1,k*d(t,s)).*(r1(t)-r1(s)).*(Tn(p,t)-Tn(p,s))./d(t,s);

Limit1 =@(t) -0.5/pi*r1p(t).*p.*Un(p-1,t)./(j(t).^2);

Nq = 512; 

[xq,wq] = chebyshev1_rule(Nq,-1,1);

u = 0; 

n = 0; 

l=1; 

A = zeros(Nq,Nq); 

for i=1:Nq 
    
    v = Fsing1(xq,xq(i)).*Tn(n,xq).*Tn(l,xq(i)); 
    
    v(d(xq,xq(i))==0) = Limit1(xq(i));
    
    A(:,i) = v; 
    
end 

u = wq'*A*wq; 

% bt=[1,1,1];
% 
% mex mex_funs.cpp IntegrationRules.cpp Adom.cpp kers.cpp prods.cpp geo.cpp in_cells.cpp tchv_cof.cpp greenfunction.cpp CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" -I/home/matlab/boost_1_55_0 -lfftw3 -lm ;
% 
% Vst = mex_funs(209,10,10,256,256,1,1,-1,bt,length(bt),-2,4,0);