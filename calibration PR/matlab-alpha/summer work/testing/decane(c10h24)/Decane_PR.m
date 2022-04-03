clear 
clc
load decane_rho_different_pressure_data.mat
%PR Decane
P = [0.1E6,0.3E6,0.5E6,1E6,3E6,5E6,10E6,20E6]; %pa
omega = 0.4884; %Acentric factor
Tc = 617.7; %K
Pc = 2103000; %Pa
R = 8.31446 %J/K*mol
b = 0.07780*(R*Tc/Pc);
a = 0.45724*(((R^2)*(Tc^2))/Pc);
T = 300; %k
M = 142.2817/1000; %kg/mol
S1=0.37464+1.54226*omega-0.26992*omega^2;
for i = 1:1:36
   Tr=T/Tc;
   Tr_list(i)=Tr;
   for j = 1:1:8
       Vm(i,j) = M/rho(i,j);
   end
   alpha_0(i) = (Tr)^(-0.171813)*exp(0.125283*(1-(Tr)^(1.77634)));
   alpha_1(i) = (Tr)^(-0.607352)*exp(0.511614*(1-(Tr)^(2.20517)));
   alpha_PR_01(i) = ((R*T)/(Vm(i,1)-b)-P(1))*((Vm(i,1)^2+2*b*Vm(i,1)-b^2)/a);
   alpha_PR_03(i) = ((R*T)/(Vm(i,2)-b)-P(2))*((Vm(i,2)^2+2*b*Vm(i,2)-b^2)/a);
   alpha_PR_05(i) = ((R*T)/(Vm(i,3)-b)-P(3))*((Vm(i,3)^2+2*b*Vm(i,3)-b^2)/a);
   alpha_PR_1(i) = ((R*T)/(Vm(i,4)-b)-P(4))*((Vm(i,4)^2+2*b*Vm(i,4)-b^2)/a);
   alpha_PR_3(i) = ((R*T)/(Vm(i,5)-b)-P(5))*((Vm(i,5)^2+2*b*Vm(i,5)-b^2)/a);
   alpha_PR_5(i) = ((R*T)/(Vm(i,6)-b)-P(6))*((Vm(i,6)^2+2*b*Vm(i,6)-b^2)/a);
   alpha_PR_10(i) = ((R*T)/(Vm(i,7)-b)-P(7))*((Vm(i,7)^2+2*b*Vm(i,7)-b^2)/a);
   alpha_PR_20(i) = ((R*T)/(Vm(i,8)-b)-P(8))*((Vm(i,8)^2+2*b*Vm(i,8)-b^2)/a);
   alpha_soave_PR(i) =(1+S1*(1-sqrt(Tr)))^2;
   alpha_Gasem(i) = exp((1.943+0.926*Tr)*(1-(Tr)^0.1441));
   alpha_twu(i) = alpha_0(i)+omega*(alpha_1(i)-alpha_0(i));
   Temp(i)=T;
   T = T+10;
end
hold on
plot(Tr_list,alpha_PR_01,'h','MarkerSize',8,'MarkerFaceColor','r','LineWidth',1.2)
plot(Tr_list,alpha_PR_03,'^','MarkerSize',8,'MarkerFaceColor','g','LineWidth',1.2)
plot(Tr_list,alpha_PR_05,'p','MarkerSize',8,'MarkerFaceColor','m','LineWidth',1.2)
plot(Tr_list,alpha_PR_1,'o','MarkerSize',8,'MarkerFaceColor','c','LineWidth',1.2)
plot(Tr_list,alpha_PR_3,'>','MarkerSize',8,'MarkerFaceColor','g','LineWidth',1.2)
plot(Tr_list,alpha_PR_5,'<','MarkerSize',8,'MarkerFaceColor','b','LineWidth',1.2)
plot(Tr_list,alpha_PR_10,'d','MarkerSize',8,'MarkerFaceColor','k','LineWidth',1.2)
plot(Tr_list,alpha_PR_20,'^','MarkerSize',8,'MarkerFaceColor','y','LineWidth',1.2)
plot(Tr_list,alpha_soave_PR,'r','LineWidth',2)
plot(Tr_list,alpha_Gasem,'b','LineWidth',2)
plot(Tr_list,alpha_twu,'g','LineWidth',2)
hold off
legend('PR NIST P0.1','PR NIST P0.3','PR NIST P0.5','PR NIST P1','PR NIST P3','PR NIST P5','PR NIST P10','PR NIST P20','Soave PR','Gasem','Twu')
xlabel('Tr')
ylabel('Alpha')
title('PR EoS Alpha functions value comparison for decane')