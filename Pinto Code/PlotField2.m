function PlotField2(Ps,Ncurves)

    figure;

    Xs = Ps(:,1); 
    
    Ys = Ps (:,2); 
    
    Zs = Ps(:,3:end);     
    
    hold on; 
    
    contourf(Xs,Ys,Zs,40,'LineColor','none'); 
    
    Muestr = 100;
    
    ts = linspace(-1,1,Muestr); 
    
    C2plotX=zeros(Ncurves,Muestr);
    
    C2plotY=zeros(Ncurves,Muestr);

    for m=1:Ncurves
        
        for tn = 1:Muestr
            
        z = mex_funs(5,ts(tn),m-1,0);     
        
        C2plotX(m,tn) = z(1);
        
        C2plotY(m,tn) = z(2);       
        
        
        end 
        
        plot(C2plotX(m,:) , C2plotY(m,:)); 
        
    end 
    hold off
    




% [xs,ys,zs]=PlotField(Ncurves,xmin,xmax,ymin,ymax,anchs,amps,freqs,phases,...
%                     xdesps,ydesps,tu,k)
    
%     Ngrid =100;
%     
%     xs = linspace(3*xmin,3*xmax,Ngrid); 
%     
%     ys = linspace(3*ymin,3*ymax,Ngrid); 
%     
%     Gk = @(tx,ty,x,y) -0.5/pi*log(sqrt((x-tx).^2+(y-ty).^2));
%     
%     N =length(tu)/Ncurves; 
%     
%     if(k>1e-6)
%      %change green fun
%     end
%     
%     Tn = @(n,t) cos(n*acos(t))./sqrt(1-t.^2); 
%     
%     funY=@(t,amp,freq,phase,ymin) amp*sin(freq*t+phase)+ymin; 
%     
%     funX=@(t,anch,xmin) anch*0.5*(t+1)+xmin; 
%     
%     zs=zeros(Ngrid); 
%     
%     for m=1:Ncurves
%         
%         ['curve %d of %d',m,Ncurves]
%         
%         for j=1:Ngrid
%             
%             for i=1:Ngrid
%                 
%                 for n=1:N
%                 
%                     fun =@(t) Gk(funX(t,anchs(m),xdesps(m)),...
%                         funY(t,amps(m),freqs(m),phases(m),ydesps(m)),...
%                         xs(i),ys(j)).*Tn(n,t)*tu((m-1)*N+n);
%                     
%                     zs(j,i) = zs(j,i)+quadgk(fun,-1,1); 
%                     
%                 end 
%                 
%             end
%         end 
%     end 
%         
% 



end
