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

p= 4; 

Fsing1 =@(t,s) -1i*k*0.25*besselh(1,k*d(t,s)).*(r1(t)-r1(s)).*(Tn(p,t)-Tn(p,s))./d(t,s);

Freg1 = @(t,s) Fsing1(t,s) -k*0.5/pi*besselj(1,k*d(t,s)).*(r1(t)-r1(s)).*(Tn(p,t)-Tn(p,s))./d(t,s).*log(daux(t,s));

Reg1 =@(t,s) k*besselj(1,k*d(t,s)).*(r1(t)-r1(s)).*(Tn(p,t)-Tn(p,s))./d(t,s);

Limit1 =@(t) -0.5/pi*r1p(t).*p.*Un(p-1,t)./(j(t).^2);

% ts = linspace(-1,1,10000); 
% 
% s0 = 0.28; 
% 
% plot(ts,Fsing1(ts',s0));

Nc = 256; 
    
xc = chebpts(Nc); 

A= zeros(Nc); 

for i=1:Nc 
    
%     A(i,:) = Reg1(xc,xc(i))';    

%     v = Reg1(xc,xc(i));
%     
%     v( d(xc,xc(i))==0)=0; 


    v = Fsing1(xc,xc(i)); 
    
    v( d(xc,xc(i))==0) = Limit1(xc(i)); 
    
    aux = zeros(2*Nc-2,1); 
    
    aux(1) = v(Nc); 
        
    aux(Nc) = v(1); 
        
    aux(2:(Nc-1)) = flip(v(2:(Nc-1))); 

    aux((Nc+1):end) = v(2:(Nc-1)); 

    aux = fft(aux);

    v = aux(1:Nc)/(Nc-1); 

    v(1) = 0.5*v(1); 

    v(end) = 0.5*v(end); 
    
    A(i,:) = v(1:Nc)';
    
    
end 

A(1,:) = A(1,:)*pi; 

A(2:end,1) = A(2:end,1)*pi; 

A(2:end,2:end) = A(2:end,2:end)*0.5*pi; 