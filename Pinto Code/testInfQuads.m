clear all; 
mex mex_funs.cpp IntegrationRules.cpp Adom.cpp kers.cpp prods.cpp geo.cpp in_cells.cpp tchv_cof.cpp greenfunction.cpp CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" -I/home/matlab/boost_1_55_0 -lfftw3 -lm ;

%geometry dimension
s0 = 40; 
%quadrature points (agrandar) 
% nQs = [64,128,256,512,1024]; 

nQs = [64,128,256,512,1024];

%cantidad de ciclos
nLoops = [4,4,4,4,4]; 

%geometry decay parameter
delta = 0.28; 

%Stages: 
kT= [sqrt(2)*pi];

Nsol = 50; 

angles = [-pi/4,0,pi/4,pi/6,pi/3];

Na = length(angles);

Integrals = cell(length(nQs),1); 

Nobs =300;

ThetaObs = linspace(0,2*pi,Nobs+1); 

ThetaObs = ThetaObs(1:end-1); 

Nobs= length(ThetaObs); 

bR = [   0.1894;   -0.4495;   -0.3156;   -0.4543;    0.3850;    0.3398;   -0.3818;   -0.0896;   -0.3798;   0.0721;    0.4494;   -0.2436;    0.4899;   -0.1502;   -0.2915;    0.1658;    0.4733;    0.1227;   -0.4365;   -0.1265];

dimtest = 2*s0; 

ts = chebpts(dimtest); 

Observations(1:Nobs*Na*(1),1)...
        = Obserbable(bR,delta,kT(1),angles,ThetaObs,Nsol).'; 
    
sigma = 1; 


for nn=1:length(nQs)
    
    Integrals{nn} = zeros(nLoops(nn),1); 

    latticeseq_b2('init0');

    ObservationsQuad = zeros(Nobs,nQs(nn)); 
    
    exponents = zeros(nQs(nn)*nLoops(nn),1); 
    
    QoIV = zeros(nQs(nn)*nLoops(nn),1); 
        
    for mm = 1:nLoops(nn)
       
        xmc = latticeseq_b2(s0,nQs(nn))-0.5;
        
        for ii=1:nQs(nn)
            
            bq = xmc(:,ii);

            ObservationsQuad(1:(Nobs*Na*(1)),ii)...
                = Obserbable(bq,delta,kT(1),angles,ThetaObs',Nsol); 
            
            ObservationsQuad(1:(Nobs*Na*(1)),ii)...
                =ObservationsQuad(1:(Nobs*Na*(1)),ii)-...
                    Observations(1:(Nobs*Na*(1))) ; 

            Gy = ObservationsQuad(1:Nobs*Na*(1),ii);

            exponents(ii+(mm-1)*nQs(nn)) = -0.5*Gy'*inv(sigma)*Gy; 
            
            V = sum(QoI(bq,dimtest,delta));
            
            QoIV(ii+(mm-1)*nQs(nn)) = V(2); 
            
            
        end
        
        expmin = min(exponents(1:(nQs(nn)*mm)));
        
        exponents2 = exponents(1:(nQs(nn)*mm))-expmin; 
        
        exponents2 = exp(exponents2); 
        
        QoIV2 = QoIV(1:(nQs(nn)*mm)); 
        
        Integrals{nn}(mm) = QoIV2'*exponents2/sum(exponents2);          
        
    end
    
%     for mm = 1:nLoops(nn)
%         
%         Integrals{nn}(mm) = Integrals{nn}(mm)/(mm*nLoops(nn));
%         
%     end 
%     
end 
