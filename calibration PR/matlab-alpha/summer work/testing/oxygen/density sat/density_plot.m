load density_oxygen_sat_data.mat;
Tc = 154.581;
Tr = T/Tc;
hold on 
plot(Tr,NIST_denisty_liq,'r','LineWidth',2);
plot(Tr,NIST_denisty_vap,'r','LineWidth',2);
plot(Tr,PR_denisty_liq,'b','LineWidth',2);
plot(Tr,PR_denisty_vap,'b','LineWidth',2);
plot(Tr,SRK_denisty_liq,'g','LineWidth',2);
plot(Tr,SRK_denisty_vap,'g','LineWidth',2);
xlabel('Tr')
ylabel('Density(kg/m^3)')
title('Oxygen density for different EoS')
legend('NIST liquid density','NIST vapor density','PR liquid density','PR vapor density','SRK liquid density','SRK vapor density')