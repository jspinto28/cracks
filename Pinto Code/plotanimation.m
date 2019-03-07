figure; 

xlim([-1,1]); 

for i=1:length(bgeo)
   
     zg = GetGeo(ts',bgeo{i},delta);
     
     plot( zr(:,1), zr(:,2), zg(:,1), zg(:,2)); 
     
     ylim([-2,2]);
     
     pause(0.01); 
    
    
end
