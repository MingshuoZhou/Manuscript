clear 
clc
fluid1='oxygen';
T0=60;%K
T1=400;%K
delta_T=5;
delta_P=2000;%kpa;
times=(T1-T0)/delta_T+1;

for j=1:3
    P(1,j)=1000+(j-1)*delta_P;
    for i=1:times
        P(i,j)=P(1,j)
        T(i,j)=T0+(i-1)*delta_T;
        rho_liq(i,j)=refpropm('D','T',T(i,j),'P',P(i,j),fluid1);
    end
end
xlswrite('density_predictions.csv',rho_liq,1,'A1')
%%
%alpha PR & SRK for oxygen
omega = 0.0222; %Acentric factor
Tc = refpropm('T','C',0,'',0,fluid1); %K
Pc=refpropm('P','C',0,'',0,fluid1)*1000; %Pa
R = 8.31446; %J/K*mol
Zc = refpropm('Z','C',0,'',0,fluid1)*1.168;%0.327
M = 31.9988/1000; %kg/mol
P=P*1E3;%Pa
b1 = 0.07780*(R*Tc/Pc);
a1 = 0.45724*(((R^2)*(Tc^2))/Pc);
b2 = 0.08664*(R*Tc/Pc);
a2 = 0.42747*(((R^2)*(Tc^2))/Pc);

omega_c = 1-3*Zc;
omega_b = 0.08446114781;
omega_a = 3*Zc^2+3*(1-2*Zc)*omega_b+omega_b^2+1-3*Zc;
a3 = omega_a*(((R^2)*(Tc^2))/Pc);
b3 = omega_b*(R*Tc/Pc);
c3 = omega_c*(R*Tc/Pc);
S1=0.37464+1.54226*omega-0.26992*omega^2;
S2=0.48508+1.55171*omega-0.176*omega^2;
%%
for j=1:3
    for i = 1:1:times
    Pr = P(i,j)/Pc; 
    Tr = T(i,j)/Tc;
    Tr_list(i,j)=Tr;
    Pr_list(i,j)=Pr;
    Vm_liq(i,j) = M/rho_liq(i,j);
    %    Vm_vap(i) = M/rho_vap(i);
    alpha_PT_liq(i,j) = ((R*T(i,j))/(Vm_liq(i,j)-b3)-P(i,j))*((Vm_liq(i,j)^2+b3*Vm_liq(i,j)+c3*Vm_liq(i,j)-c3*b3)/a3);
    %    alpha_PT_vap(i) = ((R*T)/(Vm_vap(i)-b3)-P(i))*((Vm_vap(i)^2+b3*Vm_vap(i)+c3*Vm_vap(i)-c3*b3)/a3);
    alpha_PR_liq(i,j) = ((R*T(i,j))/(Vm_liq(i,j)-b1)-P(i,j))*((Vm_liq(i,j)^2+2*b1*Vm_liq(i,j)-b1^2)/a1);
    alpha_SRK_liq(i) = ((R*T(i))/(Vm_liq(i)-b2)-P(i))*((Vm_liq(i)^2+b2*Vm_liq(i))/a2);
    %    alpha_PR_vap(i) = ((R*T)/(Vm_vap(i)-b1)-P(i))*((Vm_vap(i)^2+2*b1*Vm_vap(i)-b1^2)/a1);
    %    alpha_SRK_vap(i) = ((R*T)/(Vm_vap(i)-b2)-P(i))*((Vm_vap(i)^2+b2*Vm_vap(i))/a2);
    alpha_soave_PR(i,j) =(1+S1*(1-sqrt(Tr)))^2;
    alpha_soave_SRK(i,j) =(1+S2*(1-sqrt(Tr)))^2;
    alpha_Gasem(i,j) = exp((1.943+0.926*Tr)*(1-(Tr)^0.1441));
    %    alpha_twu(i) = (Tr)^(0.081043)*exp(0.915696*(1-(Tr)^(2.61622)));
    end
end
%%
%读写工作
for j=1:9
    if mod(j,3)==1
        column=(j+2)/3;
        data_write(:,j)=Tr_list(:,column);
    elseif mod(j,3)==2
        column=(j+1)/3;
        data_write(:,j)=Pr_list(:,column);
    elseif mod(j,3)==0
        column=j/3;
        data_write(:,j)=alpha_PR_liq(:,column);
    end
end
data_1=data_write(:,1:3);
data_2=data_write(:,4:6);
data_3=data_write(:,7:9);
%%        
% 列名称
title={'Tr','Pr','alpha'}

%生成表格，按列生成
% VariableNames 参数用于设置列头
result_table=table(data_1(:,1),data_1(:,2),data_1(:,3),'VariableNames',title)

% 保存数据
writetable(result_table, 'alpha.csv');
writetable(result_table, 'raw alpha for exp_data.csv');
% xlswrite('alpha.csv',header,1,'A1')
% xlswrite('alpha.csv',data_1,1,'A2')
% xlswrite('raw alpha for exp_data.csv',header,1,'A1')
% xlswrite('raw alpha for exp_data.csv',data_1,1,'A2')
%%
 plot(Tr_list,alpha_PR_liq,'h','MarkerSize',8,'MarkerFaceColor','r','LineWidth',1.2)
% plot(Tr_list,alpha_PT_liq,'h','MarkerSize',8,'MarkerFaceColor','g','LineWidth',1.2)
% plot(Tr_list,alpha_SRK_vap,'h','MarkerSize',8,'MarkerFaceColor','m','LineWidth',1.2)
%  plot(Tr_list,alpha_PR_vap,'^','MarkerSize',8,'MarkerFaceColor','r','LineWidth',1.2)
% plot(Tr_list,alpha_PT_vap,'^','MarkerSize',8,'MarkerFaceColor','g','LineWidth',1.2)
% plot(Tr_list,alpha_SRK_vap,'^','MarkerSize',8,'MarkerFaceColor','m','LineWidth',1.2)
 plot(Tr_list,alpha_soave_PR,'r','LineWidth',2)
% plot(Tr_list,alpha_soave_SRK,'g','LineWidth',2)
% plot(Tr_list,alpha_Gasem,'b','LineWidth',2)

hold off
%legend('PR NIST liq','PT NIST liq','SRK NIST liq','PR NIST vap','PT NIST vap','SRK NIST vapor','Soave PR','Soave SRK','Gasem')
xlabel('Pr')
ylabel('\alpha')
title('EoS Alpha functions value vs Pr for oxygen')
% axis([0.6,1,0.8,2.4])