x=linspace(-1,1,100); 

uIntVectorD = zeros(length(x), 2*Nt+1); 
uIntVectorN = zeros(length(x), 2*Nt+1);
uExtVectorD = zeros(length(x), 2*Nt+1); 
uExtVectorN = zeros(length(x), 2*Nt+1);
uIncVectorD = zeros(length(x), 2*Nt+1); 
uIncVectorN = zeros(length(x), 2*Nt+1);

interface= 1; 



for ll=-Nt:Nt
    
    lShifh = ll+Nt+1;    
    
    domt_trial = interface;

    dom_test=-1;
    dom_op_test=-1;
    for d=1:DOMS        
        trial_test=mex_funs(8,0,d,domt(1),domt(d+1));
        if(trial_test(domt_trial+1)>-1)
           domt_test=trial_test(domt_trial+1);
           dom_op_test=d;
           break;
        end        
    end
    for i=1:length(x)

        uIntVectorD(i,lShifh)=u_d(-x(i),tu(:,lShifh),domt_test,dom_op_test,N,domt);
        uIntVectorN(i,lShifh)=u_n(-x(i),tu(:,lShifh),domt_test,dom_op_test,N,domt);
        uExtVectorD(i,lShifh)=u_d(x(i),tu(:,lShifh),domt_trial,0,N,domt);
        uIncVectorD(i,lShifh)= mex_funs(9,k(1),x(i),ll,diam,domt_trial); 
        uExtVectorN(i,lShifh)=u_n(x(i),tu(:,lShifh),domt_trial,0,N,domt); 
        uIncVectorN(i,lShifh) = mex_funs(10,k(1),x(i),ll,diam,domt_trial); 
    end     

end 

plotL = Nt;

plot(x,uIntVectorD(:,plotL),x,uExtVectorD(:,plotL)+uIncVectorD(:,plotL)); 
figure; 
plot(x,uIntVectorN(:,plotL),x,-uExtVectorN(:,plotL)-uIncVectorN(:,plotL)); 

uIncMat =@(t) besselh(1,1*2)*besselj(1,1)/sqrt(2*pi)*exp(-1i*pi/2*(-t+2*interface+1)); 