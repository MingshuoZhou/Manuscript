clear 
clc
load decane_3D_data.mat
%PR 3D Decane
omega = 0.484; %Acentric factor
Tc = 617.7; %K
Pc = 2103000; %Pa
R = 8.31446 %J/K*mol
Zc = 0.297;
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
M = 142.2817/1000; %kg/mol
S1=0.37464+1.54226*omega-0.26992*omega^2;
% S2=0.48508+1.55171*omega-0.176*omega^2;
Pr_list = P/Pc
for j = 1:1:73
   T = 250; %k
   for i = 1:1:73
       Vm(i,j) = M / rho(i,j);
       alpha_PT(i,j) = ((R*T)/(Vm(i,j)-b3)-P(j))*((Vm(i,j)^2+b3*Vm(i,j)+c3*Vm(i,j)-c3*b3)/a3);
       alpha_PR(i,j) = ((R*T)/(Vm(i,j)-b1)-P(j))*((Vm(i,j)^2+2*b1*Vm(i,j)-b1^2)/a1);
       alpha_SRK(i,j) = ((R*T)/(Vm(i,j)-b2)-P(j))*((Vm(i,j)^2+b2*Vm(i,j))/a2);
       Tr = T/Tc;
       Tr_list(i) = Tr;
       alpha_soave_PR(i,j) =(1+S1*(1-sqrt(Tr)))^2;
       T = T + 5;
   end
end
subplot(4,1,1);
surfc(Tr_list,Pr_list,alpha_PT);
xlabel('Tr')
ylabel('Pr')
zlabel('alpha')
colorbar
title('PT EoS')
axis([0.4 1 0 0.9 0 10])

subplot(4,1,2);
surfc(Tr_list,Pr_list,alpha_PR);
xlabel('Tr')
ylabel('Pr')
zlabel('alpha')
title('PR EoS')
colorbar
axis([0.4 1 0 0.9 0 12])

subplot(4,1,3);
surfc(Tr_list,Pr_list,alpha_SRK);
xlabel('Tr')
ylabel('Pr')
zlabel('alpha')
title('SRK EoS')
colorbar
axis([0.4 1 0 0.9 0 12])

subplot(4,1,4);
surfc(Tr_list,Pr_list,alpha_soave_PR);
xlabel('Tr')
ylabel('Pr')
zlabel('alpha')
title('Soave alpha function')
colorbar

%    Vm_liq(i) = M/rho_liq(i);
%    Vm_vap(i) = M/rho_vap(i);
%    alpha_PT_liq(i) = ((R*T)/(Vm_liq(i)-b3)-P(i))*((Vm_liq(i)^2+b3*Vm_liq(i)+c3*Vm_liq(i)-c3*b3)/a3);
%    alpha_PT_vap(i) = ((R*T)/(Vm_vap(i)-b3)-P(i))*((Vm_vap(i)^2+b3*Vm_vap(i)+c3*Vm_vap(i)-c3*b3)/a3);
%   
% %    alpha_SRK_liq(i) = ((R*T)/(Vm_liq(i)-b2)-P(i))*((Vm_liq(i)^2+b2*Vm_liq(i))/a2);
%    alpha_PR_vap(i) = ((R*T)/(Vm_vap(i)-b1)-P(i))*((Vm_vap(i)^2+2*b1*Vm_vap(i)-b1^2)/a1);
% %    alpha_SRK_vap(i) = ((R*T)/(Vm_vap(i)-b2)-P(i))*((Vm_vap(i)^2+b2*Vm_vap(i))/a2);
%    alpha_soave_PR(i) =(1+S1*(1-sqrt(Tr)))^2;
% %    alpha_soave_SRK(i) =(1+S2*(1-sqrt(Tr)))^2;
%    alpha_Gasem(i) = exp((1.997+0.934*Tr)*(1-(Tr)^0.3645));
% %    alpha_twu(i) = (Tr)^(0.081043)*exp(0.915696*(1-(Tr)^(2.61622)));
%    Temp(i)=T;
%    T = T+5;
% hold on
% plot(sqrt(Pr),sqrt(alpha_PR_liq),'h','MarkerSize',8,'MarkerFaceColor','r','LineWidth',1.2)
% plot(sqrt(Pr),sqrt(alpha_PT_liq),'^','MarkerSize',8,'MarkerFaceColor','g','LineWidth',1.2)
% plot(sqrt(Pr),sqrt(alpha_PR_vap),'h','MarkerSize',8,'MarkerFaceColor','k','LineWidth',1.2)
% plot(sqrt(Pr),sqrt(alpha_PT_vap),'^','MarkerSize',8,'MarkerFaceColor','b','LineWidth',1.2)
% plot(sqrt(Pr),sqrt(alpha_soave_PR),'r','LineWidth',2)
% % plot(sqrt(Tr_list),sqrt(alpha_soave_SRK),'g','LineWidth',2)
% plot(sqrt(Pr),sqrt(alpha_Gasem),'b','LineWidth',2)
% plot(Tr_list,Pr)
% hold off
% legend('PR NIST liq','PT NIST liq','PR NIST vap','PT NIST vap','Soave PR','Gasem')
% xlabel('Tr^{0.5}')
% ylabel('Alpha^{0.5}')
% title('EoS Alpha functions value comparison for decane')
