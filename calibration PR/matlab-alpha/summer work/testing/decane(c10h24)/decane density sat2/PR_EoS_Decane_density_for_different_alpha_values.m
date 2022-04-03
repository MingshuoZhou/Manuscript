load density_decane_sat_data1.mat;
Tc = 617.7;
Tr = T/Tc;
Tr1 = T1/Tc;
hold on 
plot(Tr,NIST_denisty_liq,'r','LineWidth',2);
plot(Tr,NIST_denisty_vap,'r','LineWidth',2);
plot(Tr,Soave_denisty_liq,'o','MarkerSize',8,'MarkerFaceColor','c','MarkerEdgeColor','k','LineWidth',1.2);
plot(Tr,Soave_denisty_vap,'o','MarkerSize',8,'MarkerFaceColor','c','MarkerEdgeColor','k','LineWidth',1.2);
plot(Tr,liq_alpha_denisty_liq,'s','MarkerSize',8,'MarkerFaceColor','g','MarkerEdgeColor','k','LineWidth',1.2);
plot(Tr1,liq_alpha_denisty_vap,'s','MarkerSize',8,'MarkerFaceColor','g','MarkerEdgeColor','k','LineWidth',1.2);
plot(Tr,vap_alpha_denisty_liq,'p','MarkerSize',8,'MarkerFaceColor','y','MarkerEdgeColor','k','LineWidth',1.2);
plot(Tr,vap_alpha_denisty_vap,'p','MarkerSize',8,'MarkerFaceColor','y','MarkerEdgeColor','k','LineWidth',1.2);
xlabel('Tr')
ylabel('Density(kg/m^3)')
title('SRK EoS Decane density for different alpha values')
legend('NIST liquid density','NIST vapor density','Soave liquid density','Soave vapor density','liq alpha liquid density','liq alpha vapor density','vap alpha liquid density','vap alpha vapor density')
