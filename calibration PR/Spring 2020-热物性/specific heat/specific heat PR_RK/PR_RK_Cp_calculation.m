%RK-PR cp calculation
clear all
clc
load PRRK_Cpdata.mat
omega = 0.0222; %Acentric factor
Zc = 0.288;
Tc = 154.581; %K
Pc = 5.0430E6; %Pa
R = 8.31446; %J/(mol*K)
d1 = 0.428363;
d2 = 18.496215;
d3 = 0.338426;
d4 = 0.66;
d5 = 789.723105;
d6 = 2.512392;
A1 = -2.4407;
A0 = 0.0017;
B1 = 7.4513;
B0 = 1.9681;
C1 = 12.5040;
C0 = -2.7238;
delta1 = d1+d2*(d3-1.168*Zc)^d4+d5*(d3-1.168*Zc)^d6;
delta2 = (1-delta1)/(1+delta1);
d = (1+delta1^2)/(1+delta1);
y = 1 +(2*(1+delta1))^(1/3)+(4/(1+delta1))^(1/3);
k = (1.168*Zc*A1+A0)*omega^2+(1.1168*Zc*B1+B0)*omega+(1.168*Zc*C1+C0);
a = ((3*y^2+3*y*d+d^2+d-1)/((3*y+d-1)^2))*(R^2*Tc^2/Pc);
b = (1/(3*y+d-1))*(R*Tc/Pc);
M = 32/1000; %kg/mol
Cp_3_list = zeros(69,1);
Cp_6_list = zeros(69,1);
Cp_12_list = zeros(69,1);
Cp_25_list = zeros(69,1);
Cp_50_list = zeros(69,1);
Cp_0_list = zeros(69,1);
T = 60; %k
Sos_3_list = zeros(171,1);
Sos_6_list = zeros(171,1);
Sos_12_list = zeros(171,1);
Sos_25_list = zeros(171,1);
Sos_50_list = zeros(171,1);
h_3_list = zeros(171,1);
h_6_list = zeros(171,1);
h_12_list = zeros(171,1);
h_25_list = zeros(171,1);
h_50_list = zeros(171,1);
for i = 1:1:171
    
    alpha = (3/(2+T/Tc))^k;% no unit
    
    Cp_0 = (25.48+1.52E-2*T-0.7155E-5*T^2+1.312E-9*T^3)/M; %J/(kg*k)
    
    Cp_0_list(i,1) = Cp_0/1000;%kJ/(kg*k)
    
    Cv_0 = Cp_0 - R/M;  %J/(kg*k)
    
    e_0 = Cv_0*T;%J
    
    daaldT = a*(-3^k*k/(Tc*(2+T/Tc)^(k+1)));%1 deriative aalpha to T;
    
    d2aaldT = a*(3^k*k*(k+1)/(Tc^2*(2+T/Tc)^(k+2)));%2 derivative aalpha to temp
    %% Cp_3 calculation
    
    dpdT_3 = getdpdT(PRRKDensity3(i,1),R,M,b,daaldT,delta1,delta2);
    
    dpdrho_3 = getdpdrho(PRRKDensity3(i,1),R,M,b,T,a,alpha,delta1,delta2);
    
    Cv_3 = getCv(PRRKDensity3(i,1),Cv_0,d2aaldT,M,b,T,delta1,delta2);
    
    Cp_3 = getCp(PRRKDensity3(i,1),Cv_3,T,dpdT_3,dpdrho_3);
    
    Cp_3_list(i,1) = Cp_3/1000;%kJ/(kg*k)
    
    Sos_3 = getSos(Cp_3,Cv_3,dpdrho_3);%m/s
    
    Sos_3_list(i,1) = Sos_3;
          
    e_3 = gete(e_0,T,b,M,daaldT,a,alpha,PRRKDensity3(i,1));
    
    h_3=geth(e_3,3E6,PRRKDensity3(i,1));
    
    h_3_list(i,1)=h_3;
    
    
    %% Cp_6 calculation
    
    dpdT_6 = getdpdT(PRRKDensity6(i,1),R,M,b,daaldT,delta1,delta2);
    
    dpdrho_6 = getdpdrho(PRRKDensity6(i,1),R,M,b,T,a,alpha,delta1,delta2);
    
    Cv_6 = getCv(PRRKDensity6(i,1),Cv_0,d2aaldT,M,b,T,delta1,delta2);
    
    Cp_6 = getCp(PRRKDensity6(i,1),Cv_6,T,dpdT_6,dpdrho_6);
    
    Cp_6_list(i,1) = Cp_6/1000;%kJ/(kg*k)
    
    Sos_6 = getSos(Cp_6,Cv_6,dpdrho_6)%m/s
    
    Sos_6_list(i,1) = Sos_6;
    
    e_6 = gete(e_0,T,b,M,daaldT,a,alpha,PRRKDensity6(i,1));
    
    h_6=geth(e_6,6E6,PRRKDensity6(i,1));
    
    h_6_list(i,1)=h_6;
    
    %% Cp_12 calculation
    
    dpdT_12 = getdpdT(PRRKDensity12(i,1),R,M,b,daaldT,delta1,delta2);
    
    dpdrho_12 = getdpdrho(PRRKDensity12(i,1),R,M,b,T,a,alpha,delta1,delta2);
    
    Cv_12 = getCv(PRRKDensity12(i,1),Cv_0,d2aaldT,M,b,T,delta1,delta2);
    
    Cp_12 = getCp(PRRKDensity12(i,1),Cv_12,T,dpdT_12,dpdrho_12);
    
    Cp_12_list(i,1) = Cp_12/1000;%kJ/(kg*k)
    
    Sos_12 = getSos(Cp_12,Cv_12,dpdrho_12)%m/s
        
    Sos_12_list(i,1) = Sos_12;
    
    e_12 = gete(e_0,T,b,M,daaldT,a,alpha,PRRKDensity12(i,1));
    
    h_12=geth(e_12,12E6,PRRKDensity12(i,1));
     
    h_12_list(i,1)=h_12;
        
    %% Cp_25 calculation
    
    dpdT_25 = getdpdT(PRRKDensity25(i,1),R,M,b,daaldT,delta1,delta2);
    
    dpdrho_25 = getdpdrho(PRRKDensity25(i,1),R,M,b,T,a,alpha,delta1,delta2);
    
    Cv_25 = getCv(PRRKDensity25(i,1),Cv_0,d2aaldT,M,b,T,delta1,delta2);
    
    Cp_25 = getCp(PRRKDensity25(i,1),Cv_25,T,dpdT_25,dpdrho_25);
    
    Cp_25_list(i,1) = Cp_25/1000;%kJ/(kg*k)
    
    Sos_25 = getSos(Cp_25,Cv_25,dpdrho_25)%m/s
    
    Sos_25_list(i,1) = Sos_25;
    
    e_25 = gete(e_0,T,b,M,daaldT,a,alpha,PRRKDensity25(i,1));
    
    h_25=geth(e_25,25E6,PRRKDensity25(i,1));
    
    h_25_list(i,1)=h_25;

    %% Cp_50 calculation
    
    dpdT_50 = getdpdT(PRRKDensity50(i,1),R,M,b,daaldT,delta1,delta2);
    
    dpdrho_50 = getdpdrho(PRRKDensity50(i,1),R,M,b,T,a,alpha,delta1,delta2);
    
    Cv_50 = getCv(PRRKDensity50(i,1),Cv_0,d2aaldT,M,b,T,delta1,delta2);
    
    Cp_50 = getCp(PRRKDensity50(i,1),Cv_50,T,dpdT_50,dpdrho_50);
    
    Cp_50_list(i,1) = Cp_50/1000;%kJ/(kg*k)
    
    Sos_50 = getSos(Cp_50,Cv_50,dpdrho_50)%m/s
    
    Sos_50_list(i,1) = Sos_50;
    
    e_50 = gete(e_0,T,b,M,daaldT,a,alpha,PRRKDensity50(i,1));
    
    h_50=geth(e_50,50E6,PRRKDensity50(i,1));
     
    h_50_list(i,1)=h_50;
    
    T = T+2;
    
   
end

%%Cp vs Temp figure
figure(1);
hold on
plot(Temp,Cp_0_list,'--k');
plot(Temp,Cp_3_list,'-m');
plot(Temp,Cp_6_list,'-g');
plot(Temp,Cp_12_list,'-b');
plot(Temp,Cp_25_list,'-k');
plot(Temp,Cp_50_list,'-r');
plot(Temp_3,NIST_Cp_3,'om');
plot(Temp_6,NIST_Cp_6,'*g');
plot(Temp_6,NIST_Cp_12,'xb');
plot(Temp_6,NIST_Cp_25,'sk');
plot(Temp_6,NIST_Cp_50,'+r');
legend('Ideal gas eos','RK-PR @ 3Mpa','RK-PR @ 6Mpa','RK-PR @ 12Mpa','RK-PR @ 25Mpa','RK-PR @ 50Mpa','NIST webbook (3MPa)','NIST webbook (6MPa)','NIST webbook (12MPa)','NIST webbook (25MPa)','NIST webbook (50MPa)' )
xlabel('Temperature(K)')
ylabel('Constant pressure specific heat(kJ/kg/K)')
title('RK-PR EoS Constant pressure specific heat')
axis([50 400 0 6]);
hold off

%%speed of sound vs Temp figure
figure(2);
hold on
plot(Temp,Sos_3_list,'-m');
plot(Temp,Sos_6_list,'-g');
plot(Temp,Sos_12_list,'-b');
plot(Temp,Sos_25_list,'-k');
plot(Temp,Sos_50_list,'-r');
plot(NIST_temp_3,NIST_Sos_3,'om');
plot(NIST_temp_6,NIST_Sos_6,'*g');
plot(NIST_temp_6,NIST_Sos_12,'xb');
plot(NIST_temp_6,NIST_Sos_25,'sk');
plot(NIST_temp_6,NIST_Sos_50,'+r');
legend('RK-PR @ 3Mpa','RK-PR @ 6Mpa','RK-PR @ 12Mpa','RK-PR @ 25Mpa','RK-PR @ 50Mpa','NIST webbook (3MPa)','NIST webbook (6MPa)','NIST webbook (12MPa)','NIST webbook (25MPa)','NIST webbook (50MPa)' )
xlabel('Temperature(K)')
ylabel('Speed of sound (m/s)')
title('RK-PR Speed of Sound')
axis([50 400 0 1800])
hold off

%%h vs Temp figure
figure(3);
hold on
plot(Temp,h_3_list,'-m');
plot(Temp,h_6_list,'-g');
plot(Temp,h_12_list,'-b');
plot(Temp,h_25_list,'-k');
plot(Temp,h_50_list,'-r');
plot(NIST_temp_3,NIST_h_3,'om');
plot(NIST_temp_6,NIST_h_6,'*g');
plot(NIST_temp_6,NIST_h_12,'xb');
plot(NIST_temp_6,NIST_h_25,'sk');
plot(NIST_temp_6,NIST_h_50,'+r');
legend('PRRK @ 3Mpa','PRRK @ 6Mpa','PRRK @ 12Mpa','PRRK @ 25Mpa','PRRK @ 50Mpa','NIST webbook (3MPa)','NIST webbook (6MPa)','NIST webbook (12MPa)','NIST webbook (25MPa)','NIST webbook (50MPa)' )
xlabel('Temperature(K)')
ylabel('Specific enthalpy (kJ/kg)')
title('RK-PR Specific enthalpy')
axis([60 400 -200 400]);
hold off

function dpdT = getdpdT(rho,R,M,b,daaldT,delta1,delta2)
dpdT = rho*R/(M-b*rho)-daaldT*rho^2/((M+delta1*b*rho)*(M+delta2*b*rho)); %partical derivative p to Temp;
end

function Cv = getCv(rho,Cv_0,d2aaldT,M,b,T,delta1,delta2)
Cv = Cv_0 + T/((delta1-delta2)*b*M)*d2aaldT*log((M+delta1*b*rho)/(M+delta2*b*rho));%J/(kg*k)
end

function dpdrho = getdpdrho(rho,R,M,b,T,a,alpha,delta1,delta2)
dpdrho = -(M*(a*alpha*rho*(b*rho-M)^2*(b*rho*(delta1+delta2)+2*M)- R*T*(b*rho*delta1+M)^2*(b*rho*delta2+M)^2))/((b*rho-M)^2*(b*rho*delta1+M)^2*(b*rho*delta2+M)^2);
end

function Cp = getCp(rho,Cv,T,Q,O)
Cp = Cv + T/rho^2*Q^2/O;%J/(kg*k)
end

function Sos = getSos(Cp,Cv,dpdrho)
Sos = sqrt(Cp/Cv*dpdrho)%m/s
end

function e=gete(e_0,T,b,M,daaldT,a,alpha,rho)
e = e_0+1/(b*M)*(T*daaldT-a*alpha)*log((M+b*rho)/(M)); %J
end

function h=geth(e,P,rho)
h = (e + P/rho)/1000;
end