clear all; 
mex mex_funs.cpp IntegrationRules.cpp Adom.cpp kers.cpp prods.cpp geo.cpp in_cells.cpp tchv_cof.cpp greenfunction.cpp CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" -I/home/matlab/boost_1_55_0 -lfftw3 -lm ;


s0 = 50; 
n0 = 512; 

latticeseq_b2('init0');


%points of evaluation of the farfield. 
Nobs = 28; 
%frequencies to try 
kTest = 300; 
%incident angles 
Nangles = 8;
angles = linspace(0,2*pi-2*pi/Nangles,Nangles);


dimtest = 2*s0; 

%%generate the real geometry 
% x= t fixed for now, y2(1) =0; 

Y1r = zeros(s0+1,1);
Y1r(2) =1; 
Y2r = zeros(s0,1); 
Y2r(1:10) = 0.28*(-1).^((1:10)');
Y2r =[0;Y2r]; 
deltar =0.5; 
qr = QoI(Y1r,Y2r,deltar,dimtest); 
% figure;
% plot(qr(:,1),qr(:,2),'r');

%Real farfields 
Delta = Obserbable(Y1r,Y2r,deltar,kTest,angles,Nobs);
Ndelta = Nobs*length(kTest)*Nangles;


%Add Gausian noise (que pasa si no se le añade?) 
sigma = 1;
Delta = Delta+sigma*randn(Ndelta,1);


prevResl = -1; 
result =0; 
iters = 1; 
DnumA =0; 
NnumA =zeros(dimtest,2); 

%quadrature loop. 


xmc = latticeseq_b2(s0,n0)-0.5; 



numerator =zeros(dimtest,2); 
denominator =0;

for ii=1:n0
    Y1 = Y1r;
    Y2 = xmc(:,ii);
    Y2=[0;Y2];
    
    Gy = Obserbable(Y1,Y2,deltar,kTest,angles,Nobs); 
    Gy = Gy-Delta;
    
    
    EPhiy = exp(-0.5*Gy'*inv(sigma)*Gy); 
    
    denominator = denominator+EPhiy; 
    
    numerator = numerator+QoI(Y1,Y2,deltar,dimtest)*EPhiy;    
    
end 

DnumA = DnumA+denominator;

NnumA = NnumA+numerator;

result = NnumA/(DnumA); 

figure
hold on; 
plot(qr(:,1),qr(:,2),'r');
plot(result(:,1),result(:,2)); 


% while( (norm(prevResl-result)/norm(result) > 1e-1)&&(iters < 7))  
%     
%     xmc = latticeseq_b2(s0,n0)-0.5; 
%     
%     prevResul = result;
%     
%     numerator =zeros(dimtest,2); 
%     
%     denominator =0;
%     
%     for ii=1:n0
%     
%         Y1 = Y1r;
%         Y2 = xmc(:,ii);
%         Y2=[0;Y2];
% 
%         Gy = Obserbable(Y1,Y2,deltar,kTest,angles,Nobs); 
%         Gy = Gy-Delta;
% 
% 
%         EPhiy = exp(-0.5*Gy'*inv(sigma)*Gy); 
% 
%         denominator = denominator+EPhiy; 
% 
%         numerator = numerator+QoI(Y1,Y2,deltar,dimtest)*EPhiy;
%     
%     end 
%     
%     DnumA = DnumA+denominator;
% 
%     NnumA = NnumA+numerator;
% 
%     result = NnumA/(DnumA); 
% 
%     iters =iters+1;     
%     
% 
% end 
% 
% figure
% hold on; 
% plot(qr(:,1),qr(:,2),'r');
% plot(result(:,1),result(:,2)); 




