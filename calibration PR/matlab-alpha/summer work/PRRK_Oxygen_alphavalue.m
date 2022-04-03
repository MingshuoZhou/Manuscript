% clc
% clear all
%PR
P = 3E6; %pa
omega = 0.022; %Acentric factor
Zc = 0.288;
Tc = 154.58; %K
Pc = 5082125.6; %Pa
R = 8.31446 %J/K*mol
d1 = 0.428363;
d2 = 18.496215;
d3 = 0.338426;
d4 = 0.66;
d5 = 789.723105;
d6 = 2.512392;
A1 = -2.4407;
A0 = 0.0017;
B1 = 7.4513;
B0 = 1.9681;
C1 = 12.5040;
C0 = -2.7238;
delta1 = d1+d2*(d3-1.168*Zc)^d4+d5*(d3-1.168*Zc)^d6;
delta2 = (1-delta1)/(1+delta1);
d = (1+delta1^2)/(1+delta1);
y = 1 +(2*(1+delta1))^(1/3)+(4/(1+delta1))^(1/3);
a = ((3*y^2+3*y*d+d^2+d-1)/((3*y+d-1)^2))*(R^2*Tc^2/Pc);
b = (1/(3*y+d-1))*(R*Tc/Pc);
T = 70; %k
M = 31.9988/1000; %kg/mol
S=0.37464+1.54226*omega-0.26992*omega^2;
for i = 1:1:166
   alpha(i) = (rho(i)*R*T/(M-b*rho(i))-P)*((M+delta1*b*rho(i))*(M+delta2*b*rho(i))/(a*rho(i)^2))
   alpha2(i) =(1+S*(1-sqrt(T/Tc)))^2;
   Temp(i)=T;
   T = T+2;
end
hold on
plot(Temp,alpha)
plot(Temp,alpha1)
legend('alpha','alpha1')