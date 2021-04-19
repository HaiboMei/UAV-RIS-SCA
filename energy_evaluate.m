function [Energy_service,Energy_uav] =energy_evaluate(varphi,vartheta,K,N,f_k_l,f_k_o,F_k,D_k,a_k,x_k,R_k,p_k,q_n,t_n)  
    E_k_l=zeros(1,K);
    E_k_o=zeros(1,K);
    Energy_service=zeros(1,K);
    Energy_uav=zeros(1,N);
    
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

    for k=1:K
        E_k_l(1,k)=varphi*F_k(1,k)*(f_k_l(1,k)^(vartheta-1));
        E_k_o(1,k)=varphi*F_k(1,k)*(f_k_o(1,k)^(vartheta-1));
        E_gt=a_k(1,k)*(D_k(1,k)/R_k(1,k))*p_k(1,k)+(1-a_k(1,k))*E_k_l(1,k);
        E_uav=a_k(1,k)*E_k_o(1,k);
        Energy_service(1,k)=((1-x_k(1,k))*(E_gt+E_uav)+x_k(1,k)*E_uav)/10^7; 
    end
    for n=1:N
        v_h=sqrt((q_n(1,n+1)-q_n(1,n))^2+(q_n(2,n+1)-q_n(2,n))^2)/t_n(1,n);
        Energy_uav(1,n)=t_n(1,n)*P0*(1+3*(v_h)^2/U_tip^2)+t_n(1,n)*(1/2)*d_o*rho*s*G*(v_h)^3+...
                t_n(1,n)*P1*sqrt(sqrt(1+(v_h)^4/(4*v_o^4))-(v_h)^2/(2*v_o^2));
    end
end