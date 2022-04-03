load PR_data.mat
hold on
plot(NISTTemp,NISTDensity3,'ro')
plot(NISTTemp,NISTDensity6,'m+')
plot(NISTTemp,NISTDensity12,'g*')
plot(NISTTemp,NISTDensity25,'bsquare')
plot(NISTTemp,NISTDensity50,'kx')
plot(PRTemp,PRDensity3,'r')
plot(PRTemp,PRDensity6,'m')
plot(PRTemp,PRDensity12,'g')
plot(PRTemp,PRDensity25,'b')
plot(PRTemp,PRDensity50,'k')
legend('NIST webbook (3Mpa)','NIST webbook (6Mpa)','NIST webbook (12Mpa)','NIST webbook (25Mpa)','NIST webbook (50Mpa)','PR(3Mpa)','PR(6Mpa)','PR(12Mpa)','PR(25Mpa)','PR(50Mpa)')
xlabel('Temperature(K)');
ylabel('Density(kg/m^3)');
axis([300 800 0 800])
title('PR EoS Density for Decane')

