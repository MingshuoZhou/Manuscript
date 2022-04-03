load Oxygen_plot.mat
hold on
plot(Tr,Alpha,'b--o')
plot(Tr,Alpha_soave,'r-*')
hold off
legend('\alpha','soave \alpha function')
xlabel('T_r')
ylabel('\alpha')
title('Oxygen alpha vs soave alpha function for PR EoS')