load Water_plot.mat
hold on
plot(Tr,alpha1,'b--o')
plot(Tr,alpha_soave,'r-*')
hold off
legend('\alpha','soave \alpha function')
xlabel('T_r')
ylabel('\alpha')
title('Water alpha vs soave alpha function for PR EoS')