function [freqs,xmin,xmax,ymin,ymax,anchs,amps,phases,xs,ys]...
    = GenerateGeo(Ncurves,freqMin,freqMax,ampMin,ampMax,lmin,lmax)

    if(ampMax < 1e-6) 
        ampMax = 1;
    end
    
    CurvesPerDirX = ceil(sqrt(Ncurves)); 
    
    CurvesPerDirY = ceil(sqrt(Ncurves)); 
    
    divsX =1; 
    
    divsY = 1; 
    
    lx = divsX*CurvesPerDirX*(lmax+1e-6); 
    
    xmin = -lx/2; 
    
    xmax=lx/2;
    
    ly = divsY*CurvesPerDirY*(2*ampMax+1e-6); 
    
    ymin = -ly/2;    
    
    ymax =ly/2;
    
    xs= zeros(Ncurves,1); 
    
    ys= zeros(Ncurves,1); 
    
    freqs = zeros(Ncurves,1); 
    
    amps = zeros(Ncurves,1); 
    
    phases = zeros(Ncurves,1); 
    
    anchs = zeros(Ncurves,1); 
    
    Yselected = randperm(divsY*CurvesPerDirY,CurvesPerDirY);
    
    selected = 0; 
    
    for j=1:CurvesPerDirY
        
        Xselected = randperm(divsX*CurvesPerDirX,CurvesPerDirX);
        
        for i=1:CurvesPerDirX
            
            s=selected+1; 
            
            xs(s) = (Xselected(i)-1)*lmax-xmin; 
            
            ys(s) = (Yselected(j)-1)*(2*ampMax)-ymin; 
            
            freqs(s) = (freqMax-freqMin)*rand+freqMin; 
            
            amps(s) = (ampMax-ampMin)*rand+ampMin; 
            
            phases(s) = 2*pi*rand; 
            
            anchs(s) = (lmax-lmin)*rand+lmin; 
            
            selected = selected+1;
            
            if(selected == Ncurves)
             break; 
            end 
            
        end      
        
        if(selected == Ncurves)
           break; 
        end 
       
    end
    
    xmin = min(xs)-lmax-1; 
    
    xmax = max(xs)+ lmax+1; 
    
    ymin = min(ys) -ampMax -1; 
    
    ymax = max(ys) + ampMax +1;
    
    Muestr = ceil(2*freqMax+2);
    
    ts = linspace(-1,1,Muestr); 
    
    funY=@(t,amp,freq,phase,ymin) amp*sin(freq*t+phase)+ymin; 
    
    funX=@(t,anch,xmin) anch*0.5*(t+1)+xmin; 
    
    C2plotX=zeros(Ncurves,Muestr);
    
    C2plotY=zeros(Ncurves,Muestr);
    
   hold on
   for m=1:Ncurves
       
       C2plotY(m,:) = funY(ts,amps(m),freqs(m),phases(m),ys(m)); 
       
       C2plotX(m,:) = funX(ts,anchs(m),xs(m)); 
       
       plot(C2plotX(m,:) , C2plotY(m,:),'LineWidth',2); 
       
   end 
   hold off
   axis off;


fileID = fopen('Curves.cpp','w');

fprintf(fileID,'#include <math.h> \n');
fprintf(fileID,'#include "Curves.h" \n');
fprintf(fileID,'using namespace std; \n');
fprintf(fileID,'\n');
fprintf(fileID,'void curves(double* x, double* y,double t,int curv){ \n');
fprintf(fileID,'\n');    

for m=1:Ncurves
    
    fprintf(fileID,'if(curv== %d) \n',m);   
    fprintf(fileID,'{ \n'); 
    fprintf(fileID,'    *x = %f*0.5*(t+1)+ %f; \n',anchs(m),xs(m)); 
    fprintf(fileID,'\n');  
    fprintf(fileID,'    *y = %f*sin(%f*t+%f)+%f; \n',amps(m),freqs(m),phases(m),ys(m)); 
    fprintf(fileID,'\n');  
    fprintf(fileID,'} \n'); 
    fprintf(fileID,'\n');     
    
end 

fprintf(fileID,'}\n'); 

fprintf(fileID,'\n');

fprintf(fileID,'void tang_curves(double* x, double* y,double t,int curv){ \n');
fprintf(fileID,'\n');    

for m=1:Ncurves
    
    fprintf(fileID,'if(curv== %d) \n',m);   
    fprintf(fileID,'{ \n'); 
    fprintf(fileID,'    *x = %f*0.5; \n',anchs(m)); 
    fprintf(fileID,'\n');  
    fprintf(fileID,'    *y = %f*%f*cos(%f*t+%f); \n',amps(m),freqs(m),freqs(m),phases(m)); 
    fprintf(fileID,'\n');  
    fprintf(fileID,'} \n'); 
    fprintf(fileID,'\n');     
    
end 

fprintf(fileID,'}\n'); 



fclose(fileID);
    
end