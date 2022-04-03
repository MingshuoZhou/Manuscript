clear 
clc
load dodecane_rho_different_pressure_data_300_400.mat
%PR Decane
P = [0.1E6,0.3E6,0.5E6,1E6,3E6,5E6,10E6,20E6]; %pa
omega = 0.5742; %Acentric factor
Tc = 658.1; %K
Pc = 1817000; %Pa
R = 8.31446 %J/K*mol
b = 0.07780*(R*Tc/Pc);
a = 0.45724*(((R^2)*(Tc^2))/Pc);
T = 300; %k
M = 170.33/1000; %kg/mol
S1=0.37464+1.54226*omega-0.26992*omega^2;
L = 0.411574
m = 0.802000
N = 2.25086
for i = 1:1:101
   Tr=T/Tc;
   for j = 1:1:8
       Vm(i,j) = M/rho(i,j);
   end
%    alpha_0(i) = (Tr)^(-0.171813)*exp(0.125283*(1-(Tr)^(1.77634)));
%    alpha_1(i) = (Tr)^(-0.607352)*exp(0.511614*(1-(Tr)^(2.20517)));
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
   alpha_twu(i) = (Tr)^(N*(m-1))*exp(L*(1-(Tr)^(N*m)));
   Temp(i)=T;
   T = T+1;
end
hold on
plot(Temp,alpha_PR_01,'h','MarkerSize',8,'MarkerFaceColor','r','LineWidth',1.2)
plot(Temp,alpha_PR_03,'^','MarkerSize',8,'MarkerFaceColor','g','LineWidth',1.2)
plot(Temp,alpha_PR_05,'p','MarkerSize',8,'MarkerFaceColor','m','LineWidth',1.2)
plot(Temp,alpha_PR_1,'o','MarkerSize',8,'MarkerFaceColor','c','LineWidth',1.2)
plot(Temp,alpha_PR_3,'>','MarkerSize',8,'MarkerFaceColor','g','LineWidth',1.2)
plot(Temp,alpha_PR_5,'<','MarkerSize',8,'MarkerFaceColor','b','LineWidth',1.2)
plot(Temp,alpha_PR_10,'d','MarkerSize',8,'MarkerFaceColor','k','LineWidth',1.2)
plot(Temp,alpha_PR_20,'^','MarkerSize',8,'MarkerFaceColor','y','LineWidth',1.2)
plot(Temp,alpha_soave_PR,'r','LineWidth',2)
plot(Temp,alpha_Gasem,'b','LineWidth',2)
plot(Temp,alpha_twu,'g','LineWidth',2)
hold off
legend('PR NIST P0.1','PR NIST P0.3','PR NIST P0.5','PR NIST P1','PR NIST P3','PR NIST P5','PR NIST P10','PR NIST P20','Soave PR','Gasem','Twu')
xlabel('Temperature(k)')
ylabel('Alpha')
axis([315 360 -5 5])
title('PR EoS Alpha functions value comparison for dodecane under 300-400K')