clc
clear all
%PR for kerosene
%kerosene properties calculation
Tc_d = 617.7;
rho_d = 233;
Pc_d = 2.103E6;
omega_d = 0.488;
Tc_c = 553.64;
rho_c = 273;
Pc_c = 4.075E6;
omega_c = 0.20926;
Tc_t = 591.75;
rho_t = 292;
Pc_t = 4.1263E6;
omega_t = 0.266;
Tc = 0.78*Tc_d+0.098*Tc_c+0.122*Tc_t;
rho = 0.78*rho_d+0.098*rho_c+0.122*rho_t;
omega = 0.78*omega_d+0.098*omega_c+0.122*omega_t;
Pc = 0.78*Pc_d+0.098*Pc_c+0.122*Pc_t;

% PR
P = 3E6; %pa
R = 8.31446; %J/Kmol
b = 0.07780*(R*Tc/Pc);
a = 0.45724*(((R^2)*(Tc^2))/Pc);
S=0.37464+1.54226*omega-0.26992*omega^2;
T = 300; %k

for i = 1:1:101
alpha =(1+S*(1-sqrt(T/Tc)))^2;
A = (a*alpha*P)/(R^2*T^2);
B = (b*P)/(R*T);
z = roots([1 -(1-B) (A-2*B-3*B^2) -(A*B-B^2-B^3)])
% decide the state whether it is vap or liq
if imag(z(1,1))==0
    z_vap = z(1,1)
    rho_vap = P*32/(z_vap*R*T*1000);
    temp_vap{i}=T;
    ComR_vap{i}=z(1,1);
    RHO_vap{i}=rho_vap;
end
if imag(z(3,1))==0
    z_liq = z(3,1)
    rho_liq = P*32/(z_liq*R*T*1000);
    temp_liq{i}=T;
    ComR_liq{i}=z(3,1);
    RHO_liq{i}=rho_liq;
end
T = T+5; 
end
