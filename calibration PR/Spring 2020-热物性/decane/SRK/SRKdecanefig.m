hold on
plot(NISTTemp,NISTDensity3,'ro')
plot(NISTTemp,NISTDensity6,'m+')
plot(NISTTemp,NISTDensity12,'g*')
plot(NISTTemp,NISTDensity25,'bsquare')
plot(NISTTemp,NISTDensity50,'kx')
plot(SRKTemp,SRKDensity3,'r')
plot(SRKTemp,SRKDensity6,'m')
plot(SRKTemp,SRKDensity12,'g')
plot(SRKTemp,SRKDensity25,'b')
plot(SRKTemp,SRKDensity50,'k')
legend('NIST webbook (3Mpa)','NIST webbook (6Mpa)','NIST webbook (12Mpa)','NIST webbook (25Mpa)','NIST webbook (50Mpa)','SRK(3Mpa)','SRK(6Mpa)','SRK(12Mpa)','SRK(25Mpa)','SRK(50Mpa)')
xlabel('Temperature(K)');
ylabel('Density(kg/m^3)');
axis([300 800 0 800])
title('SRK EoS Density for Decane')

