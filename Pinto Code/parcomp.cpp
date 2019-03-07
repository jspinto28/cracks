void VdomNsPerSegment4(dcomp** A, double k,int N,int Nq, 
        double* xq,double* wq,int segmentTrial,int segmentTest,
        double tol){
    
    
    int nthread = omp_get_max_threads();
    
    if (nthread == 1) 
    {
            int Nc = 4*N; 
    
            VdomNsPerSegment2(A, k,N,Nc, 0,
            segmentTrial, segmentTest);
            
            return; 
        
        
    }
    
    omp_set_num_threads(2);
  
  
    dcomp** Opd =new dcomp*[Nq];       

    for(int j=0; j<Nq;j++)   
    {
         Opd[j]=new dcomp[Nq];

    }  
    
    int n1=N; 
    
    int n2=N; 
    
    in_cellsVNsPerSegment(k,Opd,segmentTrial,
         segmentTest, Nq, xq,0); 
//     
//     
//     int find =0;  
//     
//     
//     dcomp aux(0.0,0.0); 
//     
//     for(int l=0; l<N; l++)
//     {
//         aux = pVNs(0,l,Opd,Nq,xq,wq); 
//         
//         if(abs(aux)< tol)
//         {
//             find++; 
//         }
//         else
//         {
//             find =0; 
//         }
//         
//         if(find==2) 
//         {
//             n1 = l;             
//             
//             break; 
//         }
//     }
//     
//     find =0; 
    
    //mexPrintf("n1: %d \n",n1);

//     
//     for(int m=1; m<N; m++)
//     {
//         aux = pVNs(m,0,Opd,Nq,xq,wq); 
//         
//         if(abs(aux)< tol)
//         {
//             find++; 
//         }
//         else
//         {
//             find =0; 
//         }
//         
//         if(find==2) 
//         {
//             n2 = m; 
//             
//             break; 
//         }
//     }
//     
// 
//     
    

    
    int Lmax = 3; 
    
    int ns [2] ={N,N};    
    
    
    #pragma omp parallel
    {        
        int lev =0; 
    
        int a = 0; 

        int b= N-1; 

        int ti,tc,td; 

        double vi,vc,vd;

        int mid;        
        
        int id = omp_get_thread_num(); 
        
//         #pragma omp critical
//         {
//             mexPrintf("thread %d \n",id); 
//         }
        
        if(id ==0)
        {
        
            while(lev < Lmax)
            {
                mid = (int)((a+b)/2); 

                ti = mid-1; 

                tc = mid; 

                td = mid+1;

                vi =  abs(pVNs(ti,0,Opd,Nq,xq,wq)); 

                vc =  abs(pVNs(tc,0,Opd,Nq,xq,wq)); 

                vd =  abs(pVNs(td,0,Opd,Nq,xq,wq));

                if( ((vd < 0.5*tol)&&(vc < 0.5*tol))||
                        ((vi < 0.5*tol)&&(vc < 0.5*tol))) 
                {

                    b = mid; 

                }
                else
                {            
                    a = mid; 
                }

                lev++; 

            }   
        
            ns[0] =b;  
        }
        else if(id==1)
        {
            
            while(lev < Lmax)
            {
                mid = (int)((a+b)/2); 

                ti = mid-1; 

                tc = mid; 

                td = mid+1;

                vi =  abs(pVNs(0,ti,Opd,Nq,xq,wq)); 

                vc =  abs(pVNs(0,tc,Opd,Nq,xq,wq)); 

                vd =  abs(pVNs(0,td,Opd,Nq,xq,wq));

                if( ((vd < 0.5*tol)&&(vc < 0.5*tol))||
                        ((vi < 0.5*tol)&&(vc < 0.5*tol))) 
                {

                    b = mid; 

                }
                else
                {            
                    a = mid; 
                }

                lev++; 

            }
         
            ns[1] =b; 
        }


         
         
    }
        
        
    int n = max(ns[0],ns[1]); 
    
//     mexPrintf("%d,",n); 
    
    int Nc = 4*n+128; 
    
    VdomNsPerSegment2(A, k,n,Nc, 0,
         segmentTrial, segmentTest);
    
//     #pragma omp parallel for collapse(2)  
//     for(int l=1; l<n; l++)
//         for(int m=1; m<n; m++)
//         {
//             A[l][m] = pVNs(m,l,Opd,Nq,xq,wq); 
//         }
//     
    
        
    for(int j=0; j<Nq;j++){
        delete [] Opd[j];} 
    delete [] Opd; 




    

}