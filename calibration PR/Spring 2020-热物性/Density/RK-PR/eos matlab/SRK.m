%SRK
clc
clear all
omega = 0.0222; %Acentric factor
Tc = 154.581; %K
Pc = 5.0430E6; %Pa
S=0.48508+1.55171*omega-0.15613*omega^2;
T = 56;
R = 8.31446 %J/Kmol
b = 0.08664*(R*Tc/Pc)
a = 0.42747*(((R^2)*(Tc^2))/Pc)
M = 31.998;

for i = 1:1:73;
alpha =(1+S*(1-sqrt(T/Tc)))^2;
eqn = @(rho) ((rho*R*T)/(M-b*rho))-((a*alpha*rho^2)/(M^2+M*b*rho))-3E6
J=fzero(eqn,1000)
time{i}=T;
RHO{i}=J;
T = T+2; 
end