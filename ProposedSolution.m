function [Energy_service_op,Energy_uav_op,Energy_service_op_no_ps,Energy_uav_op_no_ps,Energy_service_bl,Energy_uav_bl,q_n_result_op,q_n_result_bl,a_k_result_op,x_k_result_op,RIS_phase_shift_op,RIS_phase_shift_bl,r_kn_op,r_kn_op_no_ps,r_kn_bl]...
    = ProposedSolution(q_n,W_k,W_R,Z_R,a,b,B,N_0,Xi,C_o,C_c,D_k,F_k,T_k,f_k_l,f_k_o,v_h_max,Delta_h_max,varphi,...
       vartheta,K,N,H,M,p_k,t_n)
    fprintf('Running proposed solution.\n');
    count=1;
    feasible=1;
    R_k=zeros(1,K);
    RIS_phase_shift = zeros(N,M);
    RIS_phase_shift_bl = zeros(N,M);
    c_kn_bl = repmat(eye(K,K),1,round(N/K)+1);
    c_kn=c_kn_bl;
    q_n_r = q_n;
    q_n_r_bl = q_n;        
    while (feasible)&&(count<=1)
        fprintf('\n');
        fprintf('Iteration:%d.\n',count);
        p_los=zeros(K,N);
        r_kn=zeros(K,N);
        count=count+1;
        for n=1:N
            d_ur = sqrt((H-Z_R)^2 + (W_R(1,1) - q_n_r(1,n))^2 + (W_R(1,2) - q_n_r(2,n))^2);
            for k=1:K
                d_ug=sqrt(H^2+(q_n_r(1,n)-W_k(1,k))^2+(q_n_r(2,n)-W_k(2,k))^2);
                d_rg=sqrt(Z_R^2 + (W_R(1,1) - W_k(1,k))^2 + (W_R(1,2) - W_k(2,k))^2);
               
                ratio=H/sqrt((q_n_r(1,n)-W_k(1,k))^2+(q_n_r(2,n)-W_k(2,k))^2);
                p_los(k,n)=1+a*exp(1)^(a*b-b*atan(ratio)*(180/pi));
                p_los(k,n)=1/p_los(k,n);

                g_los = p_los(k,n)*p_k(1,k)*Xi/(B*N_0*d_ug);
                gain =0;
                for m=1:M
                    if (c_kn(k,n)==0)
                       continue;
                    end
                    phase_shift = RIS_phase_shift(n,m);
                    angle_dif = rem(abs(0.5*m*((W_R(1,1)-W_k(1,k))/d_rg-(W_R(1,1)-q_n_r(1,n))/d_ur)*(2*pi)),pi);
                    angle_dif = angle_dif+phase_shift;
                    ad = sin(angle_dif);
                    gain = gain+ad;
                end
                g_nlos = (1 - p_los(k,n)) * p_k(1,k) * gain * Xi / (B * N_0 * d_ur * d_rg);
                r_kn(k,n)=B * log2(1 + g_los + g_nlos);
            end
        end

        for k=1:K
            R_k(1,k)=sum(r_kn(k,:));
        end
        R_k(1,:)=R_k(1,:)./N;

        %solve P1 by CVX tool;
        fprintf('Solve task offloading.\n');
        a_k_bl = ones(1,K);
        [a_k,H_k] = LP_offload(C_o,varphi,vartheta,K,R_k,D_k,T_k,F_k,f_k_l,f_k_o);
        solve p2 by CVX tool;
        fprintf('Solve task cache.\n');
        x_k_bl = zeros(1,K);
        [x_k] = LP_cache(D_k,C_c,a_k,H_k,K);
        
        R_exp=zeros(1,K);   %R_exp=R_k' as defined in (14)
        for k=1:K         
            R_exp(1,k)=((a_k(1,k)-a_k(1,k)*x_k(1,k))*D_k(1,k))/((a_k(1,k)-a_k(1,k)*x_k(1,k))*(F_k(1,k)/f_k_l(1,k))...
               -a_k(1,k)*(F_k(1,k)/f_k_o(1,k))+x_k(1,k)*(F_k(1,k)/f_k_l(1,k))+T_k(1,k)-(F_k(1,k)/f_k_l(1,k)));
        end
        
        R_exp_bl=zeros(1,K);   %R_exp=R_k' as defined in (14)
        for k=1:K         
            R_exp_bl(1,k)=((a_k_bl(1,k)-a_k_bl(1,k)*x_k_bl(1,k))*D_k(1,k))/((a_k_bl(1,k)-a_k_bl(1,k)*x_k_bl(1,k))*(F_k(1,k)/f_k_l(1,k))...
               -a_k_bl(1,k)*(F_k(1,k)/f_k_o(1,k))+x_k_bl(1,k)*(F_k(1,k)/f_k_l(1,k))+T_k(1,k)-(F_k(1,k)/f_k_l(1,k)));
        end 
        
       fprintf('Optimize UAV trajectory.\n');
       tau =1.6;
       epsilon = 0.8;
       [q_n_r] = QCP_horizontal_rotary_wing(q_n_r,H,B,W_k,t_n,v_h_max,N,K,Delta_h_max,R_exp,epsilon,p_k,Xi,N_0,a,b,tau);
       fprintf('Baseline UAV trajectory.\n'); 
       tau =1.2;
       epsilon = 0.2;
       [q_n_r_bl] = QCP_horizontal_rotary_wing(q_n_r_bl,H,B,W_k,t_n,v_h_max,N,K,Delta_h_max,R_exp_bl,epsilon,p_k,Xi,N_0,a,b,tau);
    
       for n=1:N
           if (isnan(q_n_r_bl(1,n)))||(isnan(q_n_r_bl(2,n)))
              feasible=0;
           end
           if ((q_n_r_bl(1,n)<0)||(q_n_r_bl(2,n)<0))
              feasible=0;
           end
       end
       
       for n=1:N
           if (isnan(q_n_r(1,n)))||(isnan(q_n_r(2,n)))
              feasible=0;
           end
           if ((q_n_r(1,n)<0)||(q_n_r(2,n)<0))
              feasible=0;
           end
       end
       
       if (feasible==1)
          [RIS_phase_shift,c_kn]= RIS_Phase_shift_schedule(q_n_r,H,t_n,R_exp,W_k,N,K,W_R,Z_R,p_k,Xi,B,N_0,a,b,M); 
          %[RIS_phase_shift_bl]= RIS_Phase_shift_schedule(q_n_r_bl,H,t_n,R_exp_bl,W_k,N,K,W_R,Z_R,p_k,Xi,B,N_0,a,b,M); 
       else
          fprintf('Infeasible solution.\n');
          fprintf('\n\n\n');
          return; 
       end  
    end
   
    if (feasible)
       fprintf('UAV solution feasible.\n');
       fprintf('\n\n\n');
       q_n_result_op=q_n_r;
       q_n_result_bl=q_n_r_bl;
       a_k_result_op=a_k;
       x_k_result_op=x_k;
       RIS_phase_shift_op = RIS_phase_shift;

       q_n = q_n_result_bl;
       r_kn_bl=zeros(K,N);
       for n=1:N
            d_ur = sqrt((H-Z_R)^2 + (W_R(1,1) - q_n(1,n))^2 + (W_R(1,2) - q_n(2,n))^2);
            for k=1:K
                d_ug=sqrt(H^2+(q_n(1,n)-W_k(1,k))^2+(q_n(2,n)-W_k(2,k))^2);
                d_rg = sqrt(Z_R^2 + (W_R(1,1) - W_k(1,k))^2 + (W_R(1,2) - W_k(2,k))^2);
               
                ratio=H/sqrt((q_n(1,n)-W_k(1,k))^2+(q_n(2,n)-W_k(2,k))^2);
                p_los(k,n)=1+a*exp(1)^(a*b-b*atan(ratio)*(180/pi));
                p_los(k,n)=1/p_los(k,n);

                g_los = p_los(k,n)*p_k(1,k)*Xi/(B*N_0*d_ug);
                gain =0;
                for m=1:M
                    if (c_kn_bl(k,n)==0)
                        continue;
                    end
                    phase_shift = RIS_phase_shift_bl(n,m);
                    angle_dif = rem(abs(0.5*m*((W_R(1,1)-W_k(1,k))/d_rg-(W_R(1,1)-q_n_r(1,n))/d_ur)*(2*pi)),pi);
                    angle_dif = angle_dif+phase_shift;
                    ad = sin(angle_dif);
                    gain = gain+ad;
                end
                g_nlos = (1 - p_los(k,n)) * p_k(1,k) * gain * Xi/ (B * N_0 * d_ur * d_rg);
                r_kn_bl(k,n)=B * log2(1 + g_los + g_nlos);
            end
        end

       for k=1:K
           R_k(1,k)=sum(r_kn_bl(k,:));
       end
       
       R_k(1,:)=R_k(1,:)./N;  %R_k as defined in (12);
       [Energy_service_bl,Energy_uav_bl]=energy_evaluate(varphi,vartheta,K,N,f_k_l,f_k_o,F_k,D_k,a_k_bl,x_k_bl,R_k,p_k,q_n,t_n);
       
       q_n = q_n_result_op;
       r_kn_op_no_ps=zeros(K,N);
       %without passive phase shift
       for n=1:N
           d_ur = sqrt((H-Z_R)^2 + (W_R(1,1) - q_n(1,n))^2 + (W_R(1,2) - q_n(2,n))^2);
           for k=1:K
                d_ug=sqrt(H^2+(q_n(1,n)-W_k(1,k))^2+(q_n(2,n)-W_k(2,k))^2);
                d_rg = sqrt(Z_R^2 + (W_R(1,1) - W_k(1,k))^2 + (W_R(1,2) - W_k(2,k))^2);
               
                ratio=H/sqrt((q_n(1,n)-W_k(1,k))^2+(q_n(2,n)-W_k(2,k))^2);
                p_los(k,n)=1+a*exp(1)^(a*b-b*atan(ratio)*(180/pi));
                p_los(k,n)=1/p_los(k,n);

                g_los = p_los(k,n)*p_k(1,k)*Xi/(B*N_0*d_ug);
                gain =0;
                for m=1:M
                    if (c_kn_bl(k,n)==0)
                        continue;
                    end
                    phase_shift = RIS_phase_shift_bl(n,m);
                    angle_dif = rem(abs(0.5*m*((W_R(1,1)-W_k(1,k))/d_rg-(W_R(1,1)-q_n_r(1,n))/d_ur)*(2*pi)),pi);
                    angle_dif = angle_dif+phase_shift;
                    ad = sin(angle_dif);
                    gain = gain+ad;
                end
                g_nlos = (1 - p_los(k,n)) * p_k(1,k) * gain * Xi / (B * N_0 * d_ur * d_rg);
                r_kn_op_no_ps(k,n)=B * log2(1 + g_los + g_nlos);
            end
       end
       
       for k=1:K
           R_k(1,k)=sum(r_kn_op_no_ps(k,:));
       end
       
       R_k(1,:)=R_k(1,:)./N;  %R_k as defined in (12);
       [Energy_service_op_no_ps,Energy_uav_op_no_ps]=energy_evaluate(varphi,vartheta,K,N,f_k_l,f_k_o,F_k,D_k,a_k,x_k,R_k,p_k,q_n,t_n);
        
       r_kn_op=zeros(K,N);
       for n=1:N
            d_ur = sqrt((H-Z_R)^2 + (W_R(1,1) - q_n(1,n))^2 + (W_R(1,2) - q_n(2,n))^2);
            for k=1:K
                d_ug=sqrt(H^2+(q_n(1,n)-W_k(1,k))^2+(q_n(2,n)-W_k(2,k))^2);
                d_rg = sqrt(Z_R^2 + (W_R(1,1) - W_k(1,k))^2 + (W_R(1,2) - W_k(2,k))^2);
               
                ratio=H/sqrt((q_n(1,n)-W_k(1,k))^2+(q_n(2,n)-W_k(2,k))^2);
                p_los(k,n)=1+a*exp(1)^(a*b-b*atan(ratio)*(180/pi));
                p_los(k,n)=1/p_los(k,n);

                g_los = p_los(k,n)*p_k(1,k)*Xi/(B*N_0*d_ug);
                gain =0;
                for m=1:M
                    if (c_kn(k,n)==0)
                        continue;
                    end
                    phase_shift = RIS_phase_shift_op(n,m);
                    angle_dif = rem(abs(0.5*m*((W_R(1,1)-W_k(1,k))/d_rg-(W_R(1,1)-q_n_r(1,n))/d_ur)*(2*pi)),pi);
                    angle_dif = angle_dif+phase_shift;
                    ad = sin(angle_dif);
                    gain = gain+ad;
                end
                g_nlos = (1 - p_los(k,n)) * p_k(1,k) * gain * Xi / (B * N_0 * d_ur * d_rg);
                r_kn_op(k,n)=B * log2(1 + g_los + g_nlos);
            end
       end
       
       for k=1:K
           R_k(1,k)=sum(r_kn_op(k,:));
       end
       
       R_k(1,:)=R_k(1,:)./N;  %R_k as defined in (12);
       [Energy_service_op,Energy_uav_op]=energy_evaluate(varphi,vartheta,K,N,f_k_l,f_k_o,F_k,D_k,a_k,x_k,R_k,p_k,q_n,t_n);
    else
       fprintf('Solution infeasible.\n');
       fprintf('\n\n\n');
       return;
    end   
end