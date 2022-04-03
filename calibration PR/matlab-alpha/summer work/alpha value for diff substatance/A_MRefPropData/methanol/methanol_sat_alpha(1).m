clear 
clc
load methanol_sat_data.mat
%alpha PR & SRK for heptafluoropropane
P = P*10^6;
omega = 0.565; %Acentric factor
Tc = 512.6; %K
Pc = 8103500; %Pa
R = 8.31446 %J/K*mol
Zc = 0.222;
b1 = 0.07780*(R*Tc/Pc);
a1 = 0.45724*(((R^2)*(Tc^2))/Pc);
b2 = 0.08664*(R*Tc/Pc);
a2 = 0.42747*(((R^2)*(Tc^2))/Pc);
omega_c = 1-3*Zc;
omega_b = 0.08446114781;
omega_a = 3*Zc^2+3*(1-2*Zc)*omega_b+omega_b^2+1-3*Zc;
a3 = omega_a*(((R^2)*(Tc^2))/Pc);
b3 = omega_b*(R*Tc/Pc);
c3 = omega_c*(R*Tc/Pc);
T = 252.0; %k
M = 32.04/1000; %kg/mol
S1=0.37464+1.54226*omega-0.26992*omega^2;
S2=0.48508+1.55171*omega-0.176*omega^2;

for i = 1:1:128
   Pr = P/Pc; 
   Tr = T/Tc;
   Tr_list(i)=Tr;
   Vm_liq(i) = M/rho_liq(i);
   Vm_vap(i) = M/rho_vap(i);
   alpha_PT_liq(i) = ((R*T)/(Vm_liq(i)-b3)-P(i))*((Vm_liq(i)^2+b3*Vm_liq(i)+c3*Vm_liq(i)-c3*b3)/a3);
   alpha_PT_vap(i) = ((R*T)/(Vm_vap(i)-b3)-P(i))*((Vm_vap(i)^2+b3*Vm_vap(i)+c3*Vm_vap(i)-c3*b3)/a3);
   alpha_PR_liq(i) = ((R*T)/(Vm_liq(i)-b1)-P(i))*((Vm_liq(i)^2+2*b1*Vm_liq(i)-b1^2)/a1)
   alpha_SRK_liq(i) = ((R*T)/(Vm_liq(i)-b2)-P(i))*((Vm_liq(i)^2+b2*Vm_liq(i))/a2);
   alpha_PR_vap(i) = ((R*T)/(Vm_vap(i)-b1)-P(i))*((Vm_vap(i)^2+2*b1*Vm_vap(i)-b1^2)/a1);
   alpha_SRK_vap(i) = ((R*T)/(Vm_vap(i)-b2)-P(i))*((Vm_vap(i)^2+b2*Vm_vap(i))/a2);
   alpha_soave_PR(i) =(1+S1*(1-sqrt(Tr)))^2;
   alpha_soave_SRK(i) =(1+S2*(1-sqrt(Tr)))^2;
   alpha_Gasem(i) = exp((1.943+0.926*Tr)*(1-(Tr)^0.1441));
%    alpha_twu(i) = (Tr)^(0.081043)*exp(0.915696*(1-(Tr)^(2.61622)));
   Temp(i)=T;
   T = T+2;
end
hold on

 plot(Tr_list,alpha_PR_liq,'h','MarkerSize',8,'MarkerFaceColor','r','LineWidth',1.2)
% plot(Tr_list,alpha_PT_liq,'h','MarkerSize',8,'MarkerFaceColor','g','LineWidth',1.2)
% plot(Tr_list,alpha_SRK_vap,'h','MarkerSize',8,'MarkerFaceColor','m','LineWidth',1.2)
 plot(Tr_list,alpha_PR_vap,'^','MarkerSize',8,'MarkerFaceColor','r','LineWidth',1.2)
% plot(Tr_list,alpha_PT_vap,'^','MarkerSize',8,'MarkerFaceColor','g','LineWidth',1.2)
% plot(Tr_list,alpha_SRK_vap,'^','MarkerSize',8,'MarkerFaceColor','m','LineWidth',1.2)
 plot(Tr_list,alpha_soave_PR,'r','LineWidth',2)
% plot(Tr_list,alpha_soave_SRK,'g','LineWidth',2)
% plot(Tr_list,alpha_Gasem,'b','LineWidth',2)

hold off
legend('PR NIST liq','PT NIST liq','SRK NIST liq','PR NIST vap','PT NIST vap','SRK NIST vapor','Soave PR','Soave SRK','Gasem')
xlabel('Tr')
ylabel('\alpha')
title('EoS Alpha functions value vs Pr for methanol')
% axis([0.6,1,0.8,2.4])