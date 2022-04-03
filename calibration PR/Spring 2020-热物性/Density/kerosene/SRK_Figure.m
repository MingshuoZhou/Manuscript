hold on
plot(BWRTemp3,BWRDensity3,'ro')
plot(BWRTemp3,BWRDensity6,'m+')
plot(BWRTemp3,BWRDensity12,'g*')
plot(BWRTemp3,BWRDensity25,'bsquare')
plot(BWRTemp3,BWRDensity50,'kx')
plot(SRKTemp3,SRKDensity3,'r')
plot(SRKTemp3,SRKDensity6,'m')
plot(SRKTemp3,SRKDensity12,'g')
plot(SRKTemp3,SRKDensity25,'b')
plot(SRKTemp3,SRKDensity50,'k')
legend('NIST webbook (3Mpa)','NIST webbook (6Mpa)','NIST webbook (12Mpa)','NIST webbook (25Mpa)','NIST webbook (50Mpa)','SRK(3Mpa)','SRK(6Mpa)','SRK(12Mpa)','SRK(25Mpa)','SRK(50Mpa)')
xlabel('Temperature(K)');
ylabel('Density(kg/m^3)');
axis([300 800 0 800])
title('(a) Soave-Redlich-Kwong EoS')
save('SRK_data.mat')
