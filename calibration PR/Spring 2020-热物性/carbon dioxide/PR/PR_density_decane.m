clc
clear all
%PR
P = 0.1E6; %pa
% omega = 0.4884; %Acentric factor
% Tc = 617.7; %K
% Pc = 2103000; %Pa
% R = 8.31446 %J/Kmol
% b = 0.07780*(R*Tc/Pc)
% a = 0.45724*(((R^2)*(Tc^2))/Pc)
% S=0.37464+1.54226*omega-0.26992*omega^2;
% T = 300; %k
% M = 142.2817;
omega = 0.22394; %Acentric factor
Tc = 304.1; %K
Pc = 7380000; %Pa
R = 8.31446; %J/(mol*K)
b = 0.07780*(R*Tc/Pc); % m^3/mol
a = 0.45724*(((R^2)*(Tc^2))/Pc); % J*m^3/mol
M = 44.10/1000; %kg/mol
S=0.37464+1.54226*omega-0.26992*omega^2;
T = 220; %k

for i = 1:1:101
alpha =(1+S*(1-sqrt(T/Tc)))^2;
ALPHA(i) = alpha;
A = (a*alpha*P)/(R^2*T^2);
B = (b*P)/(R*T);
z = roots([1 -(1-B) (A-2*B-3*B^2) -(A*B-B^2-B^3)])
% decide the state whether it is vap or liq
if imag(z(1,1))==0
    z_vap = z(1,1);
    rho_vap = P*M/(z_vap*R*T);
    temp_vap{i}=T;
    ComR_vap{i}=z(1,1);
    RHO_vap{i}=rho_vap;
end
if imag(z(3,1))==0
    z_liq = z(3,1)
    rho_liq = P*M/(z_liq*R*T);
    temp_liq{i}=T;
    ComR_liq{i}=z(3,1);
    RHO_liq{i}=rho_liq;
end
T = T+5; 
end