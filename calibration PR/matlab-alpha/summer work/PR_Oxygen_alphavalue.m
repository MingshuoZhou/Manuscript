clear all
clc
load oxygen_rho_data.mat
%PR
P = 3E6; %pa
omega = 0.022; %Acentric factor
Tc = 154.58; %K
Pc = 5082125.6; %Pa
R = 8.31446 %J/K*mol
b1 = 0.07780*(R*Tc/Pc);
a1 = 0.45724*(((R^2)*(Tc^2))/Pc);
b2 = 0.08664*(R*Tc/Pc);%SRK
a2 = 0.42747*(((R^2)*(Tc^2))/Pc);
T = 70; %k
M = 31.9988/1000; %kg/mol
S1=0.37464+1.54226*omega-0.26992*omega^2;
S2=0.48508+1.55171*omega-0.15613*omega^2;
for i = 1:1:166
   Tr=T/Tc;
   Vm(i) = M/rho(i); 
   alpha_0(i) = (Tr)^(-0.171813)*exp(0.125283*(1-(Tr)^(1.77634)));
   alpha_1(i) = (Tr)^(-0.607352)*exp(0.511614*(1-(Tr)^(2.20517)));
   alpha_PR(i) = ((R*T)/(Vm(i)-b1)-P)*((Vm(i)^2+2*b1*Vm(i)-b1^2)/a1);
   alpha_SRK(i) = ((R*T)/(Vm(i)-b2)-P)*((Vm(i)^2+b2*Vm(i))/a2);
   alpha_soave_PR(i) =(1+S1*(1-sqrt(Tr)))^2;
   alpha_soave_SRK(i) =(1+S2*(1-sqrt(Tr)))^2;
   alpha_Gasem(i) = exp((1.943+0.926*Tr)*(1-(Tr)^0.1441));
   alpha_twu(i) = alpha_0(i)+omega*(alpha_1(i)-alpha_0(i));
   Temp(i)=T;
   T = T+2;
end
hold on
plot(Temp,alpha_PR)
plot(Temp,alpha_SRK)
plot(Temp,alpha_soave_PR)
plot(Temp,alpha_soave_SRK)
plot(Temp,alpha_Gasem)
plot(Temp,alpha_twu)
hold off
legend('alpha PR','alpha SRK','alpha Soave PR','alpha Soave SRK','alpha Gasem','alpha Twu')
xlabel('Temperature(k)')
ylabel('Alpha')
title('alpha functions value comparison for oxygen')


