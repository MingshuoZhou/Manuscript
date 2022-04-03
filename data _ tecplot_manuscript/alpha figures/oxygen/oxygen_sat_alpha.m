clear 
clc
load oxygen_density.mat
%alpha PR & SRK for oxygen
omega = 0.022; %Acentric factor
Tc = 154.58; %K
Pc = 5043125.6; %Pa
R = 8.31446 %J/K*mol
Zc = 0.327;
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
M = 31.9988/1000; %kg/mol
S1=0.37464+1.54226*omega-0.26992*omega^2;
S2=0.48508+1.55171*omega-0.176*omega^2;
P = [0.1*10^6, 1*10^6, 2*10^6, 3*10^6, 4*10^6, 5.0*10^6, 7*10^6, 8*10^6, 10*10^6]
for j = 1:1:9
    T = 60; %k
    for i = 1:1:341
        Pr = P/Pc;
        Tr = T/Tc;
        Tr_list(i)=Tr;
        Vm_liq(i,j) = M/rho_liq(i,j);
        alpha_PR_liq(i,j) = ((R*T)/(Vm_liq(i,j)-b1)-P(j))*((Vm_liq(i,j)^2+2*b1*Vm_liq(i,j)-b1^2)/a1)
        alpha_soave_PR(i,j) =(1+S1*(1-sqrt(Tr)))^2;
        alpha_soave_SRK(i,j) =(1+S2*(1-sqrt(Tr)))^2;
        alpha_Gasem(i,j) = exp((1.943+0.926*Tr)*(1-(Tr)^0.1441));
        Temp(i,:)=T;
        T = T+1;
    end
end
hold on

 plot(Tr_list,alpha_PR_liq,'h','MarkerSize',8,'MarkerFaceColor','r','LineWidth',1.2)
 plot(Tr_list,alpha_soave_PR,'r','LineWidth',2)


hold off
legend('PR NIST liq','PT NIST liq','SRK NIST liq','PR NIST vap','PT NIST vap','SRK NIST vapor','Soave PR','Soave SRK','Gasem')
xlabel('Pr')
ylabel('\alpha')
title('EoS Alpha functions value vs Pr for oxygen')
 axis([0.6,1,0.8,2.4])
