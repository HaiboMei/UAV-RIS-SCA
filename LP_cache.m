%optimize task cache strategy X by solving problem P2 defined in (14) of the letter, using CVX tool; 
function [x] = LP_cache(D_k,C_c,a_k,H_k,K)  
      cvx_begin
        cvx_quiet true;
        variable x(1,K) nonnegative
        expressions data(1,K)
        
        maximize min(data)
        subject to 
           for k=1:K 
              a_k(1,k)*(1-x(1,k))<=H_k(1,k); 
           end
           
           for k=1:K
               x(1,k)<=1;
               data(1,k)=x(1,k)*a_k(1,k)*D_k(1,k);
           end
           sum(data)<=C_c; 
      cvx_end
      cvx_clear;
end