clear 
clc
load oxygen_3D_data.mat
%PR 3D Decane
omega = 0.022; %Acentric factor
Tc = 154.58; %K
Pc = 5082125.6; %Pa
R = 8.31446 %J/K*mol
Zc = 0.327;
b1 = 0.07780*(R*Tc/Pc);
a1 = 0.45724*(((R^2)*(Tc^2))/Pc);
b2 = 0.08664*(R*Tc/Pc);
a2 = 0.42747*(((R^2)*(Tc^2))/Pc);
omega_c = 1-3*Zc;
omega_b = 0.07430872445;
omega_a = 3*Zc^2+3*(1-2*Zc)*omega_b+omega_b^2+1-3*Zc;
a3 = omega_a*(((R^2)*(Tc^2))/Pc);
b3 = omega_b*(R*Tc/Pc);
c3 = omega_c*(R*Tc/Pc);
M = 31.9988/1000; %kg/mol
S1=0.37464+1.54226*omega-0.26992*omega^2;
% S2=0.48508+1.55171*omega-0.176*omega^2;
Pr_list = P/Pc
for j = 1:1:50
   T = 56; %k
   for i = 1:1:50
       Vm(i,j) = M / rho(i,j);
       alpha_PT(i,j) = ((R*T)/(Vm(i,j)-b3)-P(j))*((Vm(i,j)^2+b3*Vm(i,j)+c3*Vm(i,j)-c3*b3)/a3);
       alpha_PR(i,j) = ((R*T)/(Vm(i,j)-b1)-P(j))*((Vm(i,j)^2+2*b1*Vm(i,j)-b1^2)/a1);
       alpha_SRK(i,j) = ((R*T)/(Vm(i,j)-b2)-P(j))*((Vm(i,j)^2+b2*Vm(i,j))/a2);
       Tr = T/Tc;
       Tr_list(i) = Tr;
       alpha_soave_PR(i,j) =(1+S1*(1-sqrt(Tr)))^2;
       T = T + 2;
   end
end
subplot(4,1,1);
surfc(Tr_list,Pr_list,alpha_PT);
xlabel('Tr')
ylabel('Pr')
zlabel('alpha')
colorbar
title('PT EoS')


subplot(4,1,2);
surfc(Tr_list,Pr_list,alpha_PR);
xlabel('Tr')
ylabel('Pr')
zlabel('alpha')
title('PR EoS')
colorbar


subplot(4,1,3);
surfc(Tr_list,Pr_list,alpha_SRK);
xlabel('Tr')
ylabel('Pr')
zlabel('alpha')
title('SRK EoS')
colorbar


subplot(4,1,4);
surfc(Tr_list,Pr_list,alpha_soave_PR);
xlabel('Tr')
ylabel('Pr')
zlabel('alpha')
title('Soave alpha function')
colorbar