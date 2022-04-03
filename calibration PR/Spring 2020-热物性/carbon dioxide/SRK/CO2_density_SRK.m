%SRK decane
clc
clear all
P = 0.1E6; %pa
omega = 0.224; %Acentric factor
Tc = 304.1; %K
Pc = 7380000; %Pa
R = 8.31446; %J/Kmol
b = 0.08664*(R*Tc/Pc);
a = 0.42747*(((R^2)*(Tc^2))/Pc);
S= 0.48508+1.55171*omega-0.15613*omega^2;
T = 220; %k
M = 44.1;

for i = 1:1:59
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
    RHO_vap(i,:)=rho_vap;
end
if imag(z(3,1))==0
    z_liq = z(3,1) 
    rho_liq = P*M/(z_liq*R*T*1000);
    temp_liq{i}=T;
    ComR_liq{i}=z(3,1);
    RHO_liq(i,:)=rho_liq;
end
T = T+10; 
end

