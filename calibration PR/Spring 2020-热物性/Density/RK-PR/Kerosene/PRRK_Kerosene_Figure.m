hold on
plot(BWRTemp3,BWRDensity3,'ro')
plot(BWRTemp3,BWRDensity6,'m+')
plot(BWRTemp3,BWRDensity12,'g*')
plot(BWRTemp3,BWRDensity25,'bsquare')
plot(BWRTemp3,BWRDensity50,'kx')
plot(PRRKTemp,PRRKDensity3,'r')
plot(PRRKTemp,PRRKDensity6,'m')
plot(PRRKTemp,PRRKDensity12,'g')
plot(PRRKTemp,PRRKDensity25,'b')
plot(PRRKTemp,PRRKDensity50,'k')
legend('NIST webbook (3Mpa)','NIST webbook (6Mpa)','NIST webbook (12Mpa)','NIST webbook (25Mpa)','NIST webbook (50Mpa)','RK-PR(3Mpa)','RK-PR(6Mpa)','RK-PR(12Mpa)','RK-PR(25Mpa)','RK-PR(50Mpa)')
xlabel('Temperature(K)');
ylabel('Density(kg/m^3)');
axis([300 800 0 800])
title('(c) Generalized RK-PR EoS')
save('PRRK_data.mat')
