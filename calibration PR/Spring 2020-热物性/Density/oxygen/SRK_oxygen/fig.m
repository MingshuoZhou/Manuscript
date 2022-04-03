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
legend('NIST webbook (3Mpa)','NIST webbook (6Mpa)','NIST webbook (12Mpa)','NIST webbook (25Mpa)','NIST webbook (50Mpa)','SRK(3Mpa)','SRK(6Mpa)','SRK(12Mpa)','SRK(25Mpa)','SRK(50Mpa)')
xlabel('Temperature(K)');
ylabel('Density(kg/m^3)');
axis([60 400 0 1600])
title('(a) Soave-Redlich-kwong EoS')
save('SRK_data.mat')
