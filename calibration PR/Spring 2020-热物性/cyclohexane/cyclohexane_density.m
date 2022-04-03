%SRK cyclohexane
clc
clear all
P = 50E6; %pa
omega = 0.212; %Acentric factor
Tc = 553.8; %K
Pc = 4070000; %Pa
R = 8.31446; %J/Kmol
b = 0.08664*(R*Tc/Pc);
a = 0.42747*(((R^2)*(Tc^2))/Pc);
S= 0.48508+1.55171*omega-0.15613*omega^2;
T = 300; %k
M = 84.1595;

for i = 1:1:101
alpha =(1+S*(1-sqrt(T/Tc)))^2;
A = (a*alpha*P)/(R^2*T^2);
B = (b*P)/(R*T);
z = roots([1 -1 (A-B-B^2) (-A*B)]);
% decide the state whether it is vap or liq
if imag(z(1,1))==0
    z_vap = z(1,1)
    rho_vap = P*M/(z_vap*R*T*1000);
    temp_vap{i}=T;
    ComR_vap{i}=z(1,1);
    RHO_vap{i}=rho_vap;
end
if imag(z(3,1))==0
    z_liq = z(3,1) 
    rho_liq = P*M/(z_liq*R*T*1000);
    temp_liq{i}=T;
    ComR_liq{i}=z(3,1);
    RHO_liq{i}=rho_liq;
end
T = T+5; 
end

