%optimize task offloading strategy A by solving problem P1 defined in (13) of the letter, using CVX tool; 
function [a_opt, H_k] = LP_offload(C_o,varphi,vartheta,K,R_k,D_k,T_k,F_k,f_k_l,f_k_o)
    E_k_l=zeros(1,K);
    E_k_o=zeros(1,K);
    
    alpha=1;
    H_k=zeros(1,K);
    a_opt=zeros(1,K);
    
    
    for k=1:K
        E_k_l(1,k)=varphi*F_k(1,k)*(f_k_l(1,k)^vartheta);
        E_k_o(1,k)=varphi*F_k(1,k)*(f_k_o(1,k)^vartheta);
    end
    
    for k=1:K
        H_k(1,k)=(T_k(1,k)-F_k(1,k)/f_k_l(1,k))/(D_k(1,k)/R_k(1,k)-((F_k(1,k)/f_k_l(1,k))-(F_k(1,k)/f_k_o(1,k))));
        if (H_k(1,k)<0)
            H_k(1,k)=1;
        end
        if (H_k(1,k)>1)
            H_k(1,k)=1;
        end
       
        a=sqrt(alpha*H_k(1,k)*E_k_l(1,k)/E_k_o(1,k));
        if (a>1)
           a_opt(1,k)=1;
        else
           a_opt(1,k)=a;
        end
    end
    
    if (sum(a_opt(1,:).*f_k_o(1,:))>C_o)     
        H_k=zeros(1,K);
        a_opt=zeros(1,K);
    end
end