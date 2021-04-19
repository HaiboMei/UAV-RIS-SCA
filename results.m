clc
close all

STR=sprintf('./save/Results-6-GTs-scenario-circle.mat');
PC_file_exists = exist(STR,'file');
if (PC_file_exists)
   output_results_file= fullfile(STR);
   PC_data = load(output_results_file);
   K=PC_data.K;
   N=PC_data.N;
        
   Energy_service_op=PC_data.Energy_service_op;
   Energy_uav_op =PC_data.Energy_uav_op;
   Energy_service_op_no_ps=PC_data.Energy_service_op_no_ps;
   Energy_service_bl = PC_data.Energy_service_bl;
   Energy_uav_bl = PC_data.Energy_uav_bl;
   r_kn_op = PC_data.r_kn_op;
   r_kn_op_no_ps= PC_data.r_kn_op_no_ps;
   r_kn_bl = PC_data.r_kn_bl;
   D_k = PC_data.D_k;
   F_k = PC_data.F_k;
   T_k = PC_data.T_k;
   q_n_result_op = PC_data.q_n_result_op;
   q_n_result_bl = PC_data.q_n_result_bl;
   q_n = PC_data.q_n;
   W_k = PC_data.W_k;
   a_k_result_op = PC_data.a_k_result_op;
   x_k_result_op = PC_data.x_k_result_op;
   
   f_k_o = PC_data.f_k_o;
   f_k_l = PC_data.f_k_l;
end

x_k_result_bl=zeros(1,K);
a_k_result_bl=ones(1,K);
       
%UAV trajectoryu
figure(1)
%plot(0,0,'x','markersize',25,'markerface','r');
hold on;
plot(W_k(1,:),W_k(2,:),'x','markersize',25,'markerface','k');
hold on;
plot(q_n(1,:),q_n(2,:),'->r','markersize',10,'markerface','r','linewidth',1);
hold on;
plot(q_n_result_op(1,:),q_n_result_op(2,:),'-og','markersize',10,'markerface','g','linewidth',1);
hold on;
plot(q_n_result_bl(1,:),q_n_result_bl(2,:),'-<k','markersize',10,'markerface','k','linewidth',1);
hold on;
%legend('RIS Location','GTs Locations','Initial','Optimal','Benchmark');
legend('GTs Locations','Initial Trajectory','Proposed Solution','Benchmark');
grid('on');
set(gca,'FontSize',40)
xlabel('x(m)','fontsize',40);
ylabel('y(m)','fontsize',40);

%plot data rate 
R_k_bl=zeros(1,K);
R_k_op=zeros(1,K);
R_k_op_no_ps=zeros(1,K);
for k=1:K
    R_k_bl(1,k)=sum(r_kn_bl(k,:));
    R_k_op(1,k) =sum(r_kn_op(k,:));
    R_k_op_no_ps(1,k) =sum(r_kn_op_no_ps(k,:));
end

figure(2);
hold on
plot(1:K,R_k_bl(1,:),'r-.s','LineWidth',5,'markersize',15,'markerface','r'); 
plot(1:K,R_k_op(1,:),'g--d','LineWidth',5,'markersize',15,'markerface','g');
plot(1:K,R_k_op_no_ps(1,:),'k:>','LineWidth',5,'markersize',15,'markerface','k'); 

legend('Benchmark Solution','Optimal Solution','No RIS phase-shift');
grid('on');
set(gca,'FontSize',40)
ylabel('Data Rate (bps)','fontsize',40);
xlabel('GTs','fontsize',40);

%data rate CDF
endValue1=ceil(max(r_kn_bl(:,:)));
endValue2=ceil(max(r_kn_op(:,:)));
endValue3=ceil(max(r_kn_op_no_ps(:,:)));
    
[xTime1,yPercentage1]=funcCDF(K*N,0,endValue1,r_kn_bl(:,:)');
[xTime2,yPercentage2]=funcCDF(K*N,0,endValue2,r_kn_op(:,:)');
[xTime3,yPercentage3]=funcCDF(K*N,0,endValue3,r_kn_op_no_ps(:,:)');
    
figure(3);
hold on
plot(xTime1(450:10:K*N),yPercentage1(450:10:K*N),'r-.s','LineWidth',5,'markersize',15,'markerface','r'); 
plot(xTime2(450:10:K*N),yPercentage2(450:10:K*N),'g--d','LineWidth',5,'markersize',15,'markerface','g');
plot(xTime3(450:10:K*N),yPercentage3(450:10:K*N),'k:>','LineWidth',5,'markersize',15,'markerface','k'); 
 
legend('Benchmark Solution','Proposede Solution','RIS no passive phase-shift');
grid('on');
set(gca,'FontSize',40)
xlabel('Data Rate (bps)','fontsize',40);
ylabel('CDF','fontsize',40);

%energy consumption on service
figure(4);
hold on
plot(1:K,Energy_service_bl(1,:),'r-.s','LineWidth',5,'markersize',15,'markerface','r'); 
plot(1:K,Energy_service_op(1,:),'g--d','LineWidth',5,'markersize',15,'markerface','g');
plot(1:K,Energy_service_op_no_ps(1,:),'k:>','LineWidth',5,'markersize',15,'markerface','k'); 

legend('Benchmark Solution','Proposed Solution','RIS No passive phase-shift');
grid('on');
set(gca,'FontSize',40)
ylabel('Energy Consumption on Services (J)','fontsize',40);
xlabel('GTs','fontsize',40);

%energy consumption on UAV propulsion
endValue1=ceil(max(Energy_uav_bl(1,:)));
endValue2=ceil(max(Energy_uav_op(1,:)));
    
[xTime1,yPercentage1]=funcCDF(N,0,endValue1,Energy_uav_bl(1,:)');
[xTime2,yPercentage2]=funcCDF(N,0,endValue2,Energy_uav_op(1,:)');

figure(5);
hold on
plot(xTime1(70:3:N),yPercentage1(70:3:N),'r-.s','LineWidth',5,'markersize',15,'markerface','r'); 
plot(xTime2(70:3:N),yPercentage2(70:3:N),'g--d','LineWidth',5,'markersize',15,'markerface','g');

legend('Benchmark Solution','Proposed Solution');
grid('on');
set(gca,'FontSize',40)
xlabel('Energy Consumption on UAV Propulsion(J)','fontsize',40);
ylabel('CDF','fontsize',40);

%plot latency
R_k_bl=zeros(1,K);
R_k_op=zeros(1,K);
R_k_op_no_ps=zeros(1,K);

for k=1:K
    R_k_bl(1,k)=sum(r_kn_bl(k,:))/N;
    R_k_op(1,k)=sum(r_kn_op(k,:))/N;
    R_k_op_no_ps(1,k)=sum(r_kn_op_no_ps(k,:))/N;
end

result1=zeros(1,K);
result2=zeros(1,K);
result3=zeros(1,K);

%result3=T_k;
for k=1:K
    result1(1,k)=(1-x_k_result_bl(1,k))*(a_k_result_bl(1,k)*(F_k(1,k)/f_k_o(1,k)+D_k(1,k)/R_k_bl(1,k))+(1-a_k_result_bl(1,k))*(F_k(1,k)/f_k_l(1,k)))+x_k_result_bl(1,k)*a_k_result_bl(1,k)*(F_k(1,k)/f_k_o(1,k));
    result2(1,k)=(1-x_k_result_op(1,k))*(a_k_result_op(1,k)*(F_k(1,k)/f_k_o(1,k)+D_k(1,k)/R_k_op(1,k))+(1-a_k_result_op(1,k))*(F_k(1,k)/f_k_l(1,k)))+x_k_result_op(1,k)*a_k_result_op(1,k)*(F_k(1,k)/f_k_o(1,k));
    result3(1,k)=(1-x_k_result_op(1,k))*(a_k_result_op(1,k)*(F_k(1,k)/f_k_o(1,k)+D_k(1,k)/R_k_op_no_ps(1,k))+(1-a_k_result_op(1,k))*(F_k(1,k)/f_k_l(1,k)))+x_k_result_op(1,k)*a_k_result_op(1,k)*(F_k(1,k)/f_k_o(1,k));
end
figure(6);
hold on
plot(1:K,result1(1,:),'r-.s','LineWidth',5,'markersize',15,'markerface','r'); 
plot(1:K,result2(1,:),'g--d','LineWidth',5,'markersize',15,'markerface','g');
plot(1:K,result3(1,:),'k:>','LineWidth',5,'markersize',15,'markerface','k'); 

legend('Benchmark Solution','Proposed Solution','RIS no passive phase-shift');
grid('on');
set(gca,'FontSize',40)
ylabel('Task Latency(s)','fontsize',40);
xlabel('GTs','fontsize',40);
return;

