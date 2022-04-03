load PT_data.mat
hold on
plot(NISTTemp,NISTDensity3,'ro')
plot(NISTTemp,NISTDensity6,'m+')
plot(NISTTemp,NISTDensity12,'g*')
plot(NISTTemp,NISTDensity25,'bsquare')
plot(NISTTemp,NISTDensity50,'kx')
plot(PRRKTemp,PRRKDensity3,'r')
plot(PRRKTemp,PRRKDensity6,'m')
plot(PRRKTemp,PRRKDensity12,'g')
plot(PRRKTemp,PRRKDensity25,'b')
plot(PRRKTemp,PRRKDensity50,'k')
legend('NIST webbook (3Mpa)','NIST webbook (6Mpa)','NIST webbook (12Mpa)','NIST webbook (25Mpa)','NIST webbook (50Mpa)','PT(3Mpa)','PT(6Mpa)','PT(12Mpa)','PT(25Mpa)','PT(50Mpa)')
xlabel('Temperature(K)');
ylabel('Density(kg/m^3)');
axis([300 800 0 800])
title('PT EoS Density for Decane')

