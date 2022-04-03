clear 
clc
load nitrogen_rho_different_pressure_data1.mat
%PR Nitrogen
P = [0.01E6,0.05E6,0.1E6,0.2E6,0.5E6]; %pa
omega = 0.038; %Acentric factor
Tc = 126.2; %K
Pc = 3.395E6; %Pa
R = 8.31446 %J/K*mol
b = 0.07780*(R*Tc/Pc);
a = 0.45724*(((R^2)*(Tc^2))/Pc);
T = 70; %k
M = 28/1000; %kg/mol
S1=0.37464+1.54226*omega-0.26992*omega^2;
for i = 1:1:34
   Tr=T/Tc;
   for j = 1:1:5
       Vm(i,j) = M/rho(i,j);
   end
   alpha_0(i) = (Tr)^(-0.171813)*exp(0.125283*(1-(Tr)^(1.77634)));
   alpha_1(i) = (Tr)^(-0.607352)*exp(0.511614*(1-(Tr)^(2.20517)));
   alpha_PR_001(i) = ((R*T)/(Vm(i,1)-b)-P(1))*((Vm(i,1)^2+2*b*Vm(i,1)-b^2)/a);
   alpha_PR_005(i) = ((R*T)/(Vm(i,2)-b)-P(2))*((Vm(i,2)^2+2*b*Vm(i,2)-b^2)/a);
   alpha_PR_01(i) = ((R*T)/(Vm(i,3)-b)-P(3))*((Vm(i,3)^2+2*b*Vm(i,3)-b^2)/a);
   alpha_PR_02(i) = ((R*T)/(Vm(i,4)-b)-P(4))*((Vm(i,4)^2+2*b*Vm(i,4)-b^2)/a);
   alpha_PR_05(i) = ((R*T)/(Vm(i,5)-b)-P(5))*((Vm(i,5)^2+2*b*Vm(i,5)-b^2)/a);
   alpha_soave_PR(i) =(1+S1*(1-sqrt(Tr)))^2;
   alpha_Gasem(i) = exp((1.943+0.926*Tr)*(1-(Tr)^0.1441));
   alpha_twu(i) = alpha_0(i)+omega*(alpha_1(i)-alpha_0(i));
   Temp(i)=T;
   T = T+10;
end
hold on
plot(Temp,alpha_PR_001,'h','MarkerSize',10,'MarkerFaceColor','r')
plot(Temp,alpha_PR_005,'h')
plot(Temp,alpha_PR_01,'h')
plot(Temp,alpha_PR_02,'h')
plot(Temp,alpha_PR_05,'p')
% plot(Temp,alpha_PR_3,'o')
% plot(Temp,alpha_PR_5,'x')
% plot(Temp,alpha_PR_10,'*')
% plot(Temp,alpha_PR_15,'d')
% plot(Temp,alpha_PR_25,'^')
plot(Temp,alpha_soave_PR)
plot(Temp,alpha_Gasem)
plot(Temp,alpha_twu)
hold off
legend('PR NIST P00.1','PR NIST P00.5','PR NIST P0.1','PR NIST P0.2','PR NIST P0.5','Soave PR','Gasem','Twu')
xlabel('Temperature(k)')
ylabel('Alpha')
title('PR EoS Alpha functions value comparison for Nitrogen')