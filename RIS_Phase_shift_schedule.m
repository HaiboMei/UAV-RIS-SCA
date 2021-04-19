function [RIS_phase_shift,c_kn,r_kn] = RIS_Phase_shift_schedule(q_n,H,t_n,R_exp,W_k,N,K,W_R,Z_R,p_k,Xi,B,N_0,a,b,M)
    p_los=zeros(K,N);
    r_kn=zeros(K,N);
    c_kn = zeros(K,N);
    c_kn_bl = repmat(eye(K,K),1,round(N/K));
    RIS_phase_shift = zeros(N,M);
    %RIS_phase_shift_bl = zeros(N,M);
    for n=1:N
        d_ur = sqrt((H-Z_R)^2 + (W_R(1,1) - q_n(1,n))^2 + (W_R(1,2) - q_n(2,n))^2);
        for k=1:K
            d_ug=sqrt(H^2+(q_n(1,n)-W_k(1,k))^2+(q_n(2,n)-W_k(2,k))^2);
            d_rg = sqrt(Z_R^2 + (W_R(1,1) - W_k(1,k))^2 + (W_R(1,2) - W_k(2,k))^2);
               
            ratio=H/sqrt((q_n(1,n)-W_k(1,k))^2+(q_n(2,n)-W_k(2,k))^2);
            p_los(k,n)=1+a*exp(1)^(a*b-b*atan(ratio)*(180/pi));
            p_los(k,n)=1/p_los(k,n);
      
            %E = p_los(k,n)*p_k(1,k)*Xi/(B*N_0*d_ug);
            E = 0;
            gain = M;
            %gain =0;
            %for m=1:M
             %   phase_shift = RIS_phase_shift_bl(n,m);
              %  angle_dif = rem(abs(0.5*m*((W_R(1,1)-W_k(1,k))/d_rg-(W_R(1,1)-q_n(1,n))/d_ur)*(2*pi)),pi);
               % angle_dif = angle_dif+phase_shift;
                %ad = sin(angle_dif);
                %gain = gain+ad;
            %end
                
            F = (1-p_los(k,n))*p_k(1,k)*Xi*gain/(B*N_0*d_ur*d_rg);
            r_kn(k,n)=B * log2(1 + E + F);
        end
        %solve the LP knapsack problems;
        f=(r_kn(:,n).*t_n(1,n).*-1)';
        A = [(r_kn(:,n).*-1)';ones(K,1)'];
        bg= [R_exp(1,k).*-1 1]';
        intcon=1:K;
        lb_12=zeros(K,1);
        ub_12=ones(K,1);
        [C_kn,fval]=intlinprog(f,intcon,A,bg,[],[],lb_12,ub_12);
        c_kn(:,n)= C_kn(:);
        %Passive phase-shift of RIS
        for k=1:K
            if (c_kn(k,n)==1)
               d_rg = sqrt(Z_R^2 + (W_R(1,1) - W_k(1,k))^2 + (W_R(1,2) - W_k(2,k))^2);
                for m=1:M
                    RIS_phase_shift(n,m)=pi/2-rem(abs(0.5*m*((W_R(1,1)-W_k(1,k))/d_rg-(W_R(1,1)-q_n(1,n))/d_ur)*(2*pi)),pi);
                end
            end
        end
    end    
end