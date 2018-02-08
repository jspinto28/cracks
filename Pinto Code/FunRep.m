function u=FunRep(t,tu,dom_op)
    z=0;
    for m=0:length(tu)-1       
        z=z+tu(m+1)*trial(m,t,0,dom_op);
    end 
    
    u=z;
end 