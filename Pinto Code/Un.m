function u= Un(n,t) 

u = sin((n+1)*acos(t))./sqrt(1-t.*t);

ind = (sqrt(1-t.*t) == 0);  

u(ind) = ((-1)^n)*(n+1); 

end 