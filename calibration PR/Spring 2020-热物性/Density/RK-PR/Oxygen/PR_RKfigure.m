hold on
plot(BWRTemp,BWRDensity3,'ro')
plot(BWRTemp,BWRDensity6,'m+')
plot(BWRTemp,BWRDensity12,'g*')
plot(BWRTemp,BWRDensity25,'bsquare')
plot(BWRTemp,BWRDensity50,'kx')
plot(PRRKTemp,PRRKDensity3,'r')
plot(PRRKTemp,PRRKDensity6,'m')
plot(PRRKTemp,PRRKDensity12,'g')
plot(PRRKTemp,PRRKDensity25,'b')
plot(PRRKTemp,PRRKDensity50,'k')
legend('NIST webbook (3Mpa)','NIST webbook (6Mpa)','NIST webbook (12Mpa)','NIST webbook (25Mpa)','NIST webbook (50Mpa)','PR(3Mpa)','PR(6Mpa)','PR(12Mpa)','PR(25Mpa)','PR(50Mpa)')
xlabel('Temperature(K)');
ylabel('Density(kg/m^3)');
axis([60 400 0 1600])
title('(C) Generalized RK-PR EoS')
save('PR_RK_data.mat')
