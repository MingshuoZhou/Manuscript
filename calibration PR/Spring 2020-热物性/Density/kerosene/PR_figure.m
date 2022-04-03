hold on
plot(BWRTemp3,BWRDensity3,'ro')
plot(BWRTemp3,BWRDensity6,'m+')
plot(BWRTemp3,BWRDensity12,'g*')
plot(BWRTemp3,BWRDensity25,'bsquare')
plot(BWRTemp3,BWRDensity50,'kx')
plot(PRTemp6,PRDensity3,'r')
plot(PRTemp6,PRDensity6,'m')
plot(PRTemp6,PRDensity12,'g')
plot(PRTemp6,PRDensity25,'b')
plot(PRTemp6,PRDensity50,'k')
legend('NIST webbook (3Mpa)','NIST webbook (6Mpa)','NIST webbook (12Mpa)','NIST webbook (25Mpa)','NIST webbook (50Mpa)','PR(3Mpa)','PR(6Mpa)','PR(12Mpa)','PR(25Mpa)','PR(50Mpa)')
xlabel('Temperature(K)');
ylabel('Density(kg/m^3)');
axis([300 800 0 800])
title('(b) Peng-Robinson EoS')
save('PR_data.mat')
