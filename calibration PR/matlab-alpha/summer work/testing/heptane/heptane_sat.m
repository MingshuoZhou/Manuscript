clear 
clc
load heptane_data.mat
omega = 0.349; %Acentric factor
Tc = 540.61; %K
Pc = 2736000; %Pa
R = 8.31446 %J/K*mol
b1 = 0.07780*(R*Tc/Pc);
a1 = 0.45724*(((R^2)*(Tc^2))/Pc);
b2 = 0.08664*(R*Tc/Pc);
a2 = 0.42747*(((R^2)*(Tc^2))/Pc);
T = 185; %k
M = 100.21/1000; %kg/mol
S1=0.37464+1.54226*omega-0.26992*omega^2;
S2=0.48508+1.55171*omega-0.176*omega^2;
for i = 1:1:71
   Tr=T/Tc;
   Tr_list(i)=Tr;
   Vm(i) = M/rho(i);
   alpha_PR(i) = ((R*T)/(Vm(i)-b1)-P(i))*((Vm(i)^2+2*b1*Vm(i)-b1^2)/a1);
   alpha_SRK(i) = ((R*T)/(Vm(i)-b2)-P(i))*((Vm(i)^2+b2*Vm(i))/a2);
   alpha_soave_PR(i) =(1+S1*(1-sqrt(Tr)))^2;
   alpha_soave_SRK(i) =(1+S2*(1-sqrt(Tr)))^2;
   alpha_Gasem(i) = exp((1.927+1.122*Tr)*(1-(Tr)^0.1324));
%    alpha_twu(i) = (Tr)^(0.081043)*exp(0.915696*(1-(Tr)^(2.61622)));
   Temp(i)=T;
   T = T+5;
end
hold on
plot(Tr_list,alpha_PR,'h','MarkerSize',8,'MarkerFaceColor','r','LineWidth',1.2)
plot(Tr_list,alpha_SRK,'^','MarkerSize',8,'MarkerFaceColor','g','LineWidth',1.2)
plot(Tr_list,alpha_soave_PR,'r','LineWidth',2)
plot(Tr_list,alpha_soave_SRK,'g','LineWidth',2)
plot(Tr_list,alpha_Gasem,'b','LineWidth',2)
% plot(Tr_list,alpha_twu,'g','LineWidth',2)
hold off
legend('PR NIST','SRK NIST','Soave PR','Soave SRK','Gasem')
xlabel('Tr^{0.5}')
ylabel('Alpha^{0.5}')
title('EoS Alpha functions value comparison for heptane')
axis([0.471,0.999 -5 5])