hold on
plot(BWRTemp3,BWRDensity3,'ro')
plot(BWRTemp6,BWRDensity6,'m+')
plot(BWRTemp12,BWRDensity12,'g*')
plot(BWRTemp25,BWRDensity25,'bsquare')
plot(BWRTemp50,BWRDensity50,'kx')
plot(PRTemp3,PRDensity3,'r')
plot(PRTemp3,PRDensity6,'m')
plot(PRTemp3,PRDensity12,'g')
plot(PRTemp3,PRDensity25,'b')
plot(PRTemp3,PRDensity50,'k')
legend('NIST webbook (3Mpa)','NIST webbook (6Mpa)','NIST webbook (12Mpa)','NIST webbook (25Mpa)','NIST webbook (50Mpa)','PR(3Mpa)','PR(6Mpa)','PR(12Mpa)','PR(25Mpa)','PR(50Mpa)')
xlabel('Temperature(K)');
ylabel('Density(kg/m^3)');
axis([60 400 0 1600])
title('(b) Peng-Robinson EoS')
save('PR_data.mat')
