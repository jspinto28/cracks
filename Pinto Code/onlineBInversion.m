clear all; 
mex mex_funs.cpp IntegrationRules.cpp Adom.cpp kers.cpp prods.cpp geo.cpp in_cells.cpp tchv_cof.cpp greenfunction.cpp CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" -I/home/matlab/boost_1_55_0 -lfftw3 -lm ;

%geometry dimension
s0 = 10; 
%quadrature points (agrandar) 
n0 = 1024; 

%geometry decay parameter
delta = 2.28; 

%Stages: 
T =  5;

kT= [sqrt(2)*pi,2*sqrt(2)*pi,4*sqrt(2)*pi,8*sqrt(2)*pi,16*sqrt(2)*pi,32*sqrt(2)*pi];

% kT = [8*sqrt(2)*pi,32*sqrt(2)*pi];

angles = [-pi/4,0,pi/4,pi/6,pi/3];

Na = length(angles);

% angles = angles(1:Na);

resultT= cell(T,1); 

bComputed = cell(T,1); 

errorT = zeros(T,1); 

FarerrorT = zeros(T,1); 

%Nsol: 
NsolT =[50,50,50,50,50,150]; 

%points of observation, for farfield 
% ThetaObs = [0, pi/3,pi/4, pi/2,3*pi/4, pi, 3*pi/2];

Nobs =300;

ThetaObs = linspace(0,2*pi,Nobs+1); 

ThetaObs = ThetaObs(1:end-1); 

Nobs= length(ThetaObs); 


%generate quadratures:
 latticeseq_b2('init0');

% xmc = latticeseq_b2(s0,n0)-0.5; 

%test points 
dimtest = 2*s0; 
ts = chebpts(dimtest); 

%%generate the real geometry 
bR = zeros(s0,1);
bR(1) = 0.5; 
bR(2) = -0.5; 
bR(4) = 0.1; 
bR(8) = -0.1; 
% bR(16) = 0.05; 
% bR(32) = -0.05; 
% bR = rand(s0,1)-0.5;

% bR = [   0.1894;   -0.4495;   -0.3156;   -0.4543;    0.3850;    0.3398;   -0.3818;   -0.0896;   -0.3798;   0.0721;    0.4494;   -0.2436;    0.4899;   -0.1502;   -0.2915;    0.1658;    0.4733;    0.1227;   -0.4365;   -0.1265];

%real geomtry
RealGeo = GetGeo(ts,bR,delta); 

%observations
ObservationsQuad = zeros(T*Nobs,n0); 
Observations = zeros(T*Nobs,1); 


for t=1:T
    
   Nsol=NsolT(t);

    %get the observation
    Observations((Nobs*Na*(t-1)+1):(Nobs*Na*(t)))...
        = Obserbable(bR,delta,kT(t),angles,ThetaObs,Nsol); 

    %Add Gausian noise (que pasa si no se le añade?) 
    sigma = 1e-4;
    Observations((Nobs*Na*(t-1)+1):(Nobs*Na*(t))) ...
        = Observations((Nobs*Na*(t-1)+1):(Nobs*Na*(t)))+sigma*randn(Nobs*Na,1);

    %obtain the prediction: 
    DnumA =0; 
    NnumA =zeros(dimtest,2); 

    %quadrature loop. 

%     numerator =zeros(dimtest,2); 
%     denominator =0;
    prevResl = -1; 
    result =0; 
    iters = 6; 
    DnumA =0; 
    NnumA =zeros(dimtest,2); 

    
    while( (norm(prevResl-result)/norm(result) > 1e-1)&&(iters < 7))  
    
    xmc = latticeseq_b2(s0,n0)-0.5; 
    
    prevResul = result;
      
    QoIV =zeros(dimtest*n0,2); 
    
    exponents = zeros(n0,1); 
    
        for ii=1:n0
            
%             tic;

            bq = xmc(:,ii);

            ObservationsQuad((Nobs*Na*(t-1)+1):(Nobs*Na*(t)),ii)...
                = Obserbable(bq,delta,kT(t),angles,ThetaObs',Nsol); 
            
%             toc

            ObservationsQuad((Nobs*Na*(t-1)+1):(Nobs*Na*(t)),ii)...
                =ObservationsQuad((Nobs*Na*(t-1)+1):(Nobs*Na*(t)),ii)-...
                    Observations((Nobs*Na*(t-1)+1):(Nobs*Na*(t))) ; 

            Gy = ObservationsQuad(1:Nobs*Na*(t),ii);

            exponents(ii) = -0.5*Gy'*inv(sigma)*Gy; 
            
            QoIV((1+dimtest*(ii-1)):dimtest*ii,:) = QoI(bq,dimtest,delta); 
            
            
% %             Gy = expinters(-0.5*Gy'*inv(sigma)*Gy);
%             
%             
% 
%             denominator = denominator+Gy; 
% 
%             numerator = numerator+QoI(bq,dimtest,delta)*Gy;    

        end 
    
        QoIV = reshape(QoIV,dimtest,2*n0);
        
        expmin = min(exponents); 
        
        exponents = exponents-expmin; 
        
        exponents = exp(exponents); 
        
        result = result + ...
            1/sum(exponents)*[QoIV(:,1:n0)*exponents,QoIV(:,(n0+1):end)*exponents];
        
        iters = iters +1; 

    end 
    
    resultT{t} = real(result); 
    
    %%compute farerror
    
    aux = zeros(2*dimtest-2,1); 
    
    aux(1) = resultT{t}(end,2);
    
    aux(dimtest) = resultT{t}(1,2);
    
    aux(2:(dimtest-1)) = flip(resultT{t}(2:(end-1),2));
    
    aux((dimtest+1):end) = resultT{t}(2:(end-1),2); 
    
    aux = fft(aux); 
    
    bComputed{t} = aux(1:s0)/(dimtest-1).*(((1:s0).').^(2+delta)); 
    
    bComputed{t}(1) = bComputed{t}(1)*0.5; 
    
    for tt=1:T
        
    FarerrorT(t) =  FarerrorT(t)+ norm(Obserbable(bComputed{t},delta,kT(tt),angles,ThetaObs',Nsol)...
        -  Obserbable(bR,delta,kT(tt),angles,ThetaObs,Nsol));

    end 
    
    errorT(t) = norm(RealGeo(:,2)-resultT{t}(:,2)); 
     
end 