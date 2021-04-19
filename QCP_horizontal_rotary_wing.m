function [Q_result] = QCP_horizontal_rotary_wing(q_n,H,B,W_k,t_n,v_h_max,N,K,Delta_h_max,R_exp,epsilon,p_k,Xi,N_0,a,b,tau)
    d_o=0.6;    %fuselage equivalent flat plate area;
    rho=1.225;  %air density in kg/m3;
    s=0.05;     %rotor solidity;
    G=0.503;    %Rotor disc area in m2;
    U_tip=120;  %tip seep of the rotor blade(m/s);
    v_o=4.3;   %mean rotor induced velocity in hover;
    
    omega=300;  %blade angular velocity in radians/second;
    R=0.4;      %rotor radius in meter;
    delta=0.012;%profile drage coefficient;
    k=0.1;      %incremental correction factor to induced power;
    W=20;       %aircraft weight in newton;
    P0=(delta/8)*rho*s*G*(omega^3)*(R^3);
    P1=(1+k)*(W^(3/2)/sqrt(2*rho*G));
     
    cvx_begin
        cvx_quiet false;
        variable Q(2,N+1) nonnegative
        variable phi(1,N) nonnegative
          
        expressions d_l_kn(K,N)
        expressions e_n_phi(1,N)
        expressions rate(K,N)

        minimize max(e_n_phi(1,:))
        %maximize min(rate(1,:))
        subject to 
            Q(1,1)==q_n(1,1);
            Q(2,1)==q_n(2,1);
            Q(1,N+1)==Q(1,1);
            Q(2,N+1)==Q(2,1);
            for n=1:N
                (Q(1,n+1)-Q(1,n))^2+(Q(2,n+1)-Q(2,n))^2<=min((v_h_max*t_n(1,n))^2,Delta_h_max^2);    
                
                v_n_h_l=((q_n(1,n+1)-q_n(1,n))^2+(q_n(2,n+1)-q_n(2,n))^2)/(t_n(1,n))^2;
                v_n_h=((Q(1,n+1)-Q(1,n))^2+(Q(2,n+1)-Q(2,n))^2)/(t_n(1,n))^2;
                phi_l=sqrt(sqrt(1+v_n_h_l^2/(4*v_o^4))-v_n_h_l/(2*v_o^2));
                X_bl=(4*phi_l^3+2*phi_l*v_n_h_l/v_o^2)*phi(1,n)-3*phi_l^4-phi_l^2*v_n_h_l/v_o^2;
                X_bl>=1;
                
                e_n_phi(1,n)=t_n(1,n)*(P0*(1+3*(v_n_h/U_tip^2))+(1/2)*d_o*rho*s*G*pow_pos(v_n_h,3/2)+P1*phi(1,n));
            end

            for k=1:K  
                for n=1:N
                    d_l_kn=sqrt((q_n(1,n)-W_k(1,k))^2+(q_n(2,n)-W_k(2,k))^2);
                    ratio=H/d_l_kn;
                    p_los_l=1+a*exp(1)^(a*b-b*atan(ratio)*(180/pi));
                    p_los_l=1-1/p_los_l;
                    p_los_ad = p_los_l*tau;
                    if p_los_ad>=1
                        p_los_ad = p_los_l;
                    end
                    
                    E = p_k(1,k)*Xi^2/(B*N_0)^2;
                    J_kn_l=(-B*E*(1-p_los_ad))/(log(2)*((d_l_kn^2+H^2)^2+(1-p_los_ad)*E*(d_l_kn^2+H^2)));
                    S_kn_l=B*log2(1+(1-p_los_ad)*E/(d_l_kn^2+H^2));
                    rate(k,n)=J_kn_l*((Q(1,n)-W_k(1,k))^2+(Q(2,n)-W_k(2,k))^2-d_l_kn^2)+S_kn_l;
                    
                    D_threshhod = H/tan((a-(1/b)*log(p_los_ad/(a*(1-p_los_ad))))*pi/180);
                    (Q(1,n)-W_k(1,k))^2+(Q(2,n)-W_k(2,k))^2<= (D_threshhod)^2;
                end
                sum(rate(k,:))>=epsilon*N*R_exp(1,k);
            end
     cvx_end
     Q_result = Q;
     cvx_clear;
end