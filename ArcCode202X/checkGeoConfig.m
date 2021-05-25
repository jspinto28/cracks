testself=0;
testcross = 0; 
tTest = linspace(-1,1,128); 
%%check crossing conditions (only depends of the dcoefs and number of terms)!
for ii=1:Narcs
    
    NxCos0 = (1:ceil(Nx(ii)/2))-1;
    NxCos = NxCos0(2:end);
    NxSin = 1:floor(Nx(ii)/2); 
    NyCos0 = (1:ceil(Ny(ii)/2))-1;
    NyCos = NyCos0(2:end);
    NySin = 1:floor(Ny(ii)/2); 
        
    db =@(t) (NxCos.*dcoefs(NxCos+1)')*abs(sin(NxCos'*t))+...
        (NxSin.*dcoefs(NxSin)')*abs(cos(NxSin'*t))+...
        (NyCos.*dcoefs(NyCos+1)')*abs(sin(NyCos'*t))+...
        (NySin.*dcoefs(NySin)')*abs(cos(NySin'*t));
    
    bfun =@(t) (dcoefs(NxCos0+1)')*abs(cos(NxCos0'*t))+...
        (dcoefs(NxSin)')*abs(sin(NxSin'*t))+...
        (dcoefs(NyCos0+1)')*abs(cos(NyCos0'*t))+...
        (dcoefs(NySin)')*abs(sin(NySin'*t));  
  
    
    maxdb = 0; 
    maxb =0; 
    r0p= (sqrt(coefx0{ii}(2)^2+coefy0{ii}(2)^2));
    
    for jj=1:128
       
        testself = testself +(r0p-db(tTest(jj))<0);
        
    end    
    
    for jj=(ii+1):Narcs
        
        
        for kk=1:128
            
            for nn=1:128
                
                crosdist = norm([coefx0{ii}(2),coefy0{ii}(2)]*tTest(kk)+...
                [coefx0{ii}(1),coefy0{ii}(1)]-...
                [coefx0{jj}(2),coefy0{jj}(2)]*tTest(nn)-...
                [coefx0{jj}(1),coefy0{jj}(1)]);
            
%                 testcross = testcross + (crosdist-2*sqrt(2)*maxb<0);


                 testcross = testcross + ...
                     (crosdist-(bfun(tTest(kk))+bfun(tTest(nn)))<0);
           
                
            end 
            
        end 
        
       
        
    end 
    
end 

if(testself ==0) 
    'self crossing ok'
else 
    'possible self crossings'
end 


if(testcross ==0) 
    'cross crossing ok'
else 
    'possible cross crossings'
end 
