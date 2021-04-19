clc
close all
cvx_clear

baseline=0;
%configuration about genral parameters;
N_0= 10^((-169/3)/10);   %Noise power spectrum density is -169dBm/Hz;
Xi = 10^(3/10);          %the path loss at the reference distance D0 = 1m, 3dB;
N=100;                  %Maximum UAV mission time is 5 minutesd; t_n=1s
t_n=ones(1,N);           %Initial path line duration;

%MEC configuration
C_o=10^5*N;               %MEC server compuation capacity;   
C_c=10^9;                %MEC server cache capacity:1GB; 
varphi =10^(-9);          %configuration about computing and realted energy cost
vartheta=3;

%RIS configuration
W_R = [0,0];
M = 100;  %number of phase element
Z_R = 20; %height of the RIS

%UAV configurations
H = 30;                  %UAV flight height
Delta_h_max=10;          %the UAV can maximum move 5m horizontally, in each path line;
v_h_max=10;
distance = 150;
Dr = 0;  %whether to relase the distance plos constraint or not;
%configuring GT distributionsç
[W_k,K]=GT_distributions(distance);

%Initilize UAV trajectory;
c_ini_x=sum(W_k(1,:))/size(W_k,2);
c_ini_y=sum(W_k(2,:))/size(W_k,2);
r_min=max(sqrt((W_k(1,:)-c_ini_x).^2+(W_k(2,:)-c_ini_y).^2));
r_ini=min(Delta_h_max*N/(2*pi),r_min/2);
%initial horizontal coordinate of the UAV trajectory;
q_n=zeros(2,N+1);
ang=2*pi*((1-1)/N);
q_n(1,1)=c_ini_x+r_ini*cos(ang-pi);
q_n(2,1)=c_ini_y+r_ini*sin(ang-pi);
for n=2:(N+1)
    ang=2*pi*((n-1)/N);
    q_n(1,n)=c_ini_x+r_ini*cos(ang-pi);
    q_n(2,n)=c_ini_y+r_ini*sin(ang-pi);
end

figure(1);
hold on
plot(W_k(1,:),W_k(2,:),'*r');
plot(q_n(1,:),q_n(2,:),'-xb');
grid('on');

%Initialize the wireless environement and some other verables;
a=9.61;
%a=0.961;                          %referenced from paper [Efficient 3-D Placement of an Aerial Base Station in Next Generation Cellular Networks, 2016, ICC]
b=0.16;
%b=1.6;                          %and paper [Optimal LAP Altitude for Maximum Coverage, IEEE WIRELESS COMMUNICATIONS LETTERS, VOL. 3, NO. 6, DECEMBER 2014]
eta_los=1;                       %Loss corresponding to the LoS connections defined in (2)
eta_nlos=20;                     %Loss corresponding to the NLoS connections defined in (2);
A=eta_los-eta_nlos;              %A varable defined in (2)
C=20*log10(4*pi*9/3)+eta_nlos;   %C varable defined in (2), where carrier frequncy is 900Mhz=900*10^6, and light speed is c=3*10^8; then one has f/c=9/3;
B=2*10^6;                       %overall Bandwith is 2Gb; 
P=5*10^7;                        %maximum uplink transimission power of one GT is 5mW;
    
a_k=zeros(1,N+1);
x_k=zeros(1,N+1);

%configuration about GTs, given the initilized r_km, and R_k;
D_k=zeros(1,K);          %Data rate requirment of UE task D_k;
F_k=zeros(1,K);          %Computation Requirement of UE task F_k;
T_k=zeros(1,K);          %Task deadline of GTs;

p_k=P.*ones(1,K);        %transimisson power allocated to a GT;
f_k_o=zeros(1,K);        %computation capacity allocated to a GT from MEC server;
f_k_l=zeros(1,K);        %computation capacity allocated to a GT by GT itself;

right_case=0;
while (~right_case)
    right_case=1;
    for k=1:K
        D_k(1,k)=10^7*randi([6,8]);        %Data requirement of each task: 100~500Mb
        F_k(1,k)=10^5*randi([6,8]);         %Comput cap of each task
        T_k(1,k)=randi([10,20]);           %latency request of each Task: 50~100s
    
        f_k_o(1,k)=10^5*randi([2,4]);      
        f_k_l(1,k)=10^5*randi([1,2]);   

        if (F_k(1,k)/f_k_l(1,k)>=T_k(1,k))
            right_case=0;
        end
    end
end

fprintf('Solving rotary wing.\n');
[Energy_service_op,Energy_uav_op,Energy_service_op_no_ps,Energy_uav_op_no_ps,Energy_service_bl,Energy_uav_bl,q_n_result_op,q_n_result_bl,a_k_result_op,...
    x_k_result_op,RIS_phase_shift_op,RIS_phase_shift_bl,r_kn_op,r_kn_op_no_ps,r_kn_bl]...
    = ProposedSolution(q_n,W_k,W_R,Z_R,a,b,B,N_0,Xi,C_o,C_c,D_k,F_k,T_k,f_k_l,f_k_o,v_h_max,Delta_h_max,varphi,...
       vartheta,K,N,H,M,p_k,t_n);  

%plot the initialize and optimized UAV trajectory and GT distributions for anlysis purpose. 
figure(2);
hold on
plot(W_k(1,:),W_k(2,:),'*r');
plot(q_n(1,:),q_n(2,:),'-xr');
plot(q_n_result_op(1,:),q_n_result_op(2,:),'-xg');
plot(q_n_result_bl(1,:),q_n_result_bl(2,:),'-xb');
grid('on');

fprintf('Energy consumption on service (optimiation):%d\n',sum(Energy_service_op(1,:)));
fprintf('Energy consumption on propulsion (optimiation):%d\n',sum(Energy_uav_op(1,:)));
fprintf('Throughput lead by optimizaiton:%d\n',sum(sum(r_kn_op)));

fprintf('Energy consumption on service (optimiation no RIS phase shift):%d\n',sum(Energy_service_op_no_ps(1,:)));
fprintf('Energy consumption on propulsion (optimiation no RIS phase shift):%d\n',sum(Energy_uav_op_no_ps(1,:)));
fprintf('Throughput lead by optimizaiton no RIS phase shift:%d\n',sum(sum(r_kn_op_no_ps)));

fprintf('Energy consumption on service (baseline):%d\n',sum(Energy_service_bl(1,:)));
fprintf('Energy consumption on propulsion (baseline):%d\n',sum(Energy_uav_bl(1,:)));
fprintf('Throughput lead by baseline:%d\n',sum(sum(r_kn_bl)));

STR=sprintf('Results-%d-GTs-scenario-circle.mat',K);

save(fullfile('./save',STR),'Energy_service_op','Energy_uav_op','Energy_service_op_no_ps','Energy_service_bl','Energy_uav_bl','r_kn_op','r_kn_op_no_ps','r_kn_bl','K','N','D_k','F_k','T_k','f_k_o','f_k_l',...
'q_n_result_op','q_n_result_bl','q_n','W_k','a_k_result_op','x_k_result_op');

fprintf('Simulation Finishes!');

