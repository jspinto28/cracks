function u=loss(bte,bt,kt,anglet)

    Nsol = 100; 
    
    Nobs = 2*Nsol; 
    
    [thetaobs,wobs] = lgwt(Nobs,0,2*pi); 
    
    ObsExp = Obserbable(bte,kt,anglet,thetaobs,Nsol); 
    
    Obs = Obserbable(bt,kt,anglet,thetaobs,Nsol); 
    
    u= sqrt(wobs*(abs(ObsExp-Obs)).^2);
    

end 