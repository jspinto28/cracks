function PlotCurves(Ncurves)
    ts=linspace(-1,1,100);
    x=zeros(100,1);
    y=zeros(100,1);
    hold on
    for dt=0:Ncurves(1)-1
        for t=1:length(ts)
            v=mex_funs(5,ts(t),dt,0);       
            x(t)=v(1); 
            y(t)=v(2);
        end
        plot(x,y,'LineWidth',3);
    end 
  
%     xlim([-2,4]);
%     ylim([-3,3]);
    
%     axis off;
    hold off;
    

end 