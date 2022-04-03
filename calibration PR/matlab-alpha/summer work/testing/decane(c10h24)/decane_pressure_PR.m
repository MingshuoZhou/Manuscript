%decane pressure prediction for PR
load('decane presure data.mat')
omega = 0.484; %Acentric factor
Tc = 617.7; %K
Pc = 2103000; %Pa
R = 8.31446 %J/K*mol
b1 = 0.07780*(R*Tc/Pc);
a1 = 0.45724*(((R^2)*(Tc^2))/Pc);
T = 250; %k
M = 142.2817/1000; %kg/mol
S1=0.37464+1.54226*omega-0.26992*omega^2;
for i = 1:1:72
   Tr=T/Tc;
   Tr_list(i)=Tr;
   Vm_vap(i) = M/rho_vap(i);
   Vm_liq(i) = M/rho_liq(i);
   alpha_soave_PR(i) =(1+S1*(1-sqrt(Tr)))^2;
   P_alpha_liq(i) = (R*T/(Vm_liq(i)-b1))-(a1*alpha_liq(i))/(Vm_liq(i)^2+2*b1*Vm_liq(i)-b1^2);
   P_alpha_vap(i) = (R*T/(Vm_vap(i)-b1))-(a1*alpha_vap(i))/(Vm_vap(i)^2+2*b1*Vm_vap(i)-b1^2);
   P_soave_vap(i) = (R*T/(Vm_vap(i)-b1))-(a1*alpha_soave_PR(i))/(Vm_vap(i)^2+2*b1*Vm_vap(i)-b1^2);
   P_soave_liq(i) = (R*T/(Vm_liq(i)-b1))-(a1*alpha_soave_PR(i))/(Vm_liq(i)^2+2*b1*Vm_liq(i)-b1^2);
   Temp(i)=T;
   T = T+5;
end
figure(1)
hold on
plot(Temp,P_alpha_liq,'h','MarkerSize',8,'MarkerFaceColor','r','LineWidth',1.2)
plot(Temp,P_alpha_vap,'^','MarkerSize',8,'MarkerFaceColor','g','LineWidth',1.2)
plot(Temp,P_soave_liq,'r','LineWidth',2)
plot(Temp,P_soave_vap,'g','LineWidth',2)
legend('NIST P liq','NIST P vap','P soave liq','P soave vap')
xlabel('Temperture(K)')
ylabel('Pressure(Pa)')
title('Decane presure prediction')
% axis([445,620,0,2E6])
hold off