%RK-PR cp calculation
clear all
clc
load RKPR_Cpdata.mat
omega = 0.4884; %Acentric factor
Zc = 0.247;
Tc = 617.7; %K
Pc = 2103000; %Pa
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
C1 = 12.540;
C0 = -2.7238;
delta1 = d1+d2*(d3-1.168*Zc)^d4+d5*(d3-1.168*Zc)^d6;
delta2 = (1-delta1)/(1+delta1);
d = (1+delta1^2)/(1+delta1);
y = 1 +(2*(1+delta1))^(1/3)+(4/(1+delta1))^(1/3);
k = (1.168*Zc*A1+A0)*omega^2+(1.1168*Zc*B1+B0)*omega+(1.168*Zc*C1+C0);
a = ((3*y^2+3*y*d+d^2+d-1)/((3*y+d-1)^2))*(R^2*Tc^2/Pc);
b = (1/(3*y+d-1))*(R*Tc/Pc);
M = 142/1000; %kg/mol
Cp_3_list = zeros(71,1);
Cp_6_list = zeros(71,1);
Cp_12_list = zeros(71,1);
Cp_25_list = zeros(71,1);
Cp_10_list = zeros(71,1);
Cp_0_list = zeros(71,1);
T = 300; %k
Sos_3_list = zeros(71,1);
Sos_6_list = zeros(71,1);
Sos_12_list = zeros(71,1);
Sos_25_list = zeros(71,1);
Sos_10_list = zeros(71,1);
h_3_list = zeros(71,1);
h_6_list = zeros(71,1);
h_12_list = zeros(71,1);
h_25_list = zeros(71,1);
h_5_list = zeros(71,1);
Tr = 447;
for i = 1:1:71
    
    alpha = (3/(2+T/Tc))^k;% no unit66
    
    Cp_0 = 2.407E2+5.099*T-(6.29E-4)*T^2-(1.07E-6)*T^3; %J/(kg*k)
    
    Cp_0_list(i,1) = Cp_0/1000;%kJ/(kg*k)
    
    Cv_0 = Cp_0 - R/M;  %J/(kg*k)
    
    e_0 = 2.407E2*T+5.099/2*T^2-(6.29E-4)/3*T^3-(1.07E-6)/4*T^4-R/M*T-(2.407E2*Tr+5.099/2*Tr^2-(6.29E-4)/3*Tr^3-(1.07E-6)/4*Tr^4-R/M*Tr);
    
    daaldT = a*(-3^k*k/(Tc*(2+T/Tc)^(k+1)));%1 deriative aalpha to T;
    
    d2aaldT = a*(3^k*k*(k+1)/(Tc^2*(2+T/Tc)^(k+2)));%2 derivative aalpha to temp
    %% Cp_3 calculation
    
    dpdT_3 = getdpdT(RKPRDensity3(i,1),R,M,b,daaldT,delta1,delta2);
    
    dpdrho_3 = getdpdrho(RKPRDensity3(i,1),R,M,b,T,a,alpha,delta1,delta2);
    
    Sos_0 = sqrt((Cp_0/Cv_0)*R*T/M);
    
    Sos_0_list(i) = Sos_0;
    
    Cv_3 = getCv(RKPRDensity3(i,1),Cv_0,d2aaldT,M,b,T,delta1,delta2);
    
    Cp_3 = getCp(RKPRDensity3(i,1),Cv_3,T,dpdT_3,dpdrho_3);
    
    Cp_3_list(i,1) = Cp_3/1000;%kJ/(kg*k)
    
    Sos_3 = getSos(Cp_3,Cv_3,dpdrho_3);%m/s
    
    Sos_3_list(i,1) = Sos_3;
          
    e_3 = gete(e_0,T,b,M,daaldT,a,alpha,RKPRDensity3(i,1));
    
    h_3=geth(e_3,3E6,RKPRDensity3(i,1));
    
    h_3_list(i,1)=h_3;
    
    
       %% Cp_01 calculation
    
    dpdT_01 = getdpdT(RKPRDensity01(i,1),R,M,b,daaldT,delta1,delta2);
    
    dpdrho_01 = getdpdrho(RKPRDensity01(i,1),R,M,b,T,a,alpha,delta1,delta2);
    
    Cv_01 = getCv(RKPRDensity01(i,1),Cv_0,d2aaldT,M,b,T,delta1,delta2);
    
    Cp_01 = getCp(RKPRDensity01(i,1),Cv_01,T,dpdT_01,dpdrho_01);
    
    Cp_01_list(i,1) = Cp_01/1000;%kJ/(kg*k)
    
    Sos_01 = getSos(Cp_01,Cv_01,dpdrho_01);%m/s
    
    Sos_01_list(i,1) = Sos_01;
    
    e_01 = gete(e_0,T,b,M,daaldT,a,alpha,RKPRDensity01(i,1));
    
    h_01=geth(e_01,0.1E6,RKPRDensity01(i,1));
    
    h_01_list(i,1)=h_01;
    
    %% Cp_1 calculation
    
    dpdT_1 = getdpdT(RKPRDensity1(i,1),R,M,b,daaldT,delta1,delta2);
    
    dpdrho_1 = getdpdrho(RKPRDensity1(i,1),R,M,b,T,a,alpha,delta1,delta2);
    
    Cv_1 = getCv(RKPRDensity1(i,1),Cv_0,d2aaldT,M,b,T,delta1,delta2);
    
    Cp_1 = getCp(RKPRDensity1(i,1),Cv_1,T,dpdT_1,dpdrho_1);
    
    Cp_1_list(i,1) = Cp_1/1000;%kJ/(kg*k)
    
    Sos_1 = getSos(Cp_1,Cv_1,dpdrho_1)%m/s
        
    Sos_1_list(i,1) = Sos_1;
    
    e_1 = gete(e_0,T,b,M,daaldT,a,alpha,RKPRDensity1(i,1));
    
    h_1=geth(e_1,1E6,RKPRDensity1(i,1));
     
    h_1_list(i,1)=h_1;
        
    %% Cp_2 calculation
    
    dpdT_2 = getdpdT(RKPRDensity2(i,1),R,M,b,daaldT,delta1,delta2);
    
    dpdrho_2 = getdpdrho(RKPRDensity2(i,1),R,M,b,T,a,alpha,delta1,delta2);
    
    Cv_2 = getCv(RKPRDensity2(i,1),Cv_0,d2aaldT,M,b,T,delta1,delta2);
    
    Cp_2 = getCp(RKPRDensity2(i,1),Cv_2,T,dpdT_2,dpdrho_2);
    
    Cp_2_list(i,1) = Cp_2/1000;%kJ/(kg*k)
    
    Sos_2 = getSos(Cp_2,Cv_2,dpdrho_2)%m/s
    
    Sos_2_list(i,1) = Sos_2;
    
    e_2 = gete(e_0,T,b,M,daaldT,a,alpha,RKPRDensity2(i,1));
    
    h_2=geth(e_2,2E6,RKPRDensity2(i,1));
    
    h_2_list(i,1)=h_2;

    %% Cp_5 calculation
    
    dpdT_5 = getdpdT(RKPRDensity5(i,1),R,M,b,daaldT,delta1,delta2);
    
    dpdrho_5 = getdpdrho(RKPRDensity5(i,1),R,M,b,T,a,alpha,delta1,delta2);
    
    Cv_5 = getCv(RKPRDensity5(i,1),Cv_0,d2aaldT,M,b,T,delta1,delta2);
    
    Cp_5 = getCp(RKPRDensity5(i,1),Cv_5,T,dpdT_5,dpdrho_5);
    
    Cp_5_list(i,1) = Cp_5/1000;%kJ/(kg*k)
    
    Sos_5 = getSos(Cp_5,Cv_5,dpdrho_5)%m/s
    
    Sos_5_list(i,1) = Sos_5;
    
    e_5 = gete(e_0,T,b,M,daaldT,a,alpha,RKPRDensity5(i,1));
    
    h_5=geth(e_5,5E6,RKPRDensity5(i,1));
     
    h_5_list(i,1)=h_5;
    
      %% Cp_10 calculation
    
    dpdT_10 = getdpdT(RKPRDensity10(i,1),R,M,b,daaldT,delta1,delta2);
    
    dpdrho_10 = getdpdrho(RKPRDensity10(i,1),R,M,b,T,a,alpha,delta1,delta2);
    
    Cv_10 = getCv(RKPRDensity10(i,1),Cv_0,d2aaldT,M,b,T,delta1,delta2);
    
    Cp_10 = getCp(RKPRDensity10(i,1),Cv_10,T,dpdT_10,dpdrho_10);
    
    Cp_10_list(i,1) = Cp_10/1000;%kJ/(kg*k)
    
    Sos_10 = getSos(Cp_10,Cv_10,dpdrho_10)%m/s
    
    Sos_10_list(i,1) = Sos_10;
    
    e_10 = gete(e_0,T,b,M,daaldT,a,alpha,RKPRDensity10(i,1));
    
    h_10=geth(e_10,5E6,RKPRDensity5(i,1));
     
    h_10_list(i,1)=h_10;
    
    T = T+10;
   
end

%%Cp vs Temp figure
figure(1);
hold on
plot(Temp,Cp_01_list,'--k');
plot(Temp,Cp_3_list,'-m');
plot(Temp,Cp_5_list,'-g');
plot(Temp,Cp_1_list,'-b');
plot(Temp,Cp_2_list,'-k');
plot(Temp,Cp_10_list,'-r');
% plot(Temp,NIST_Cp_3,'om');
% plot(Temp,NIST_Cp_6,'*g');
% plot(Temp,NIST_Cp_1,'xb');
% plot(Temp,NIST_Cp_2,'sk');
% plot(Temp,NIST_Cp_5,'+r');
legend('Ideal gas eos','PT @ 3Mpa','PT @ 6Mpa','PT @ 1Mpa','PT @ 2Mpa','PT @ 5Mpa','NIST webbook (3MPa)','NIST webbook (6MPa)','NIST webbook (1MPa)','NIST webbook (2MPa)','NIST webbook (5MPa)' )
xlabel('Temperature(K)')
ylabel('Constant pressure specific heat(kJ/kg/K)')
title('PT EoS Constant Pressure Specific Heat for Decane')
axis([300 800 1.5 6]);
hold off

%%speed of sound vs Temp figure
figure(2);
hold on
plot(Temp,Sos_01_list,'--k');
plot(Temp,Sos_3_list,'-m');
plot(Temp,Sos_5_list,'-g');
plot(Temp,Sos_1_list,'-b');
plot(Temp,Sos_2_list,'-k');
plot(Temp,Sos_10_list,'-r');
% plot(Temp,NIST_Sos_3,'om');
% plot(Temp,NIST_Sos_6,'*g');
% plot(Temp,NIST_Sos_1,'xb');
% plot(Temp,NIST_Sos_2,'sk');
% plot(Temp,NIST_Sos_5,'+r');
legend('Ideal gas eos','PT @ 3Mpa','PT @ 6Mpa','PT @ 1Mpa','PT @ 2Mpa','PT @ 5Mpa','NIST webbook (3MPa)','NIST webbook (6MPa)','NIST webbook (1MPa)','NIST webbook (2MPa)','NIST webbook (5MPa)' )
xlabel('Temperature(K)')
ylabel('Speed of sound (m/s)')
title('PT EoS Speed of Sound for Decane')
axis([300 800 0 2000])
hold off

%%h vs Temp figure
% figure(3);
% hold on
% plot(Temp,h_3_list,'-m');
% plot(Temp,h_6_list,'-g');
% plot(Temp,h_1_list,'-b');
% plot(Temp,h_2_list,'-k');
% plot(Temp,h_5_list,'-r');
% plot(NIST_temp_3,NIST_h_3,'om');
% plot(NIST_temp_6,NIST_h_6,'*g');
% plot(NIST_temp_6,NIST_h_1,'xb');
% plot(NIST_temp_6,NIST_h_2,'sk');
% plot(NIST_temp_6,NIST_h_5,'+r');
% legend('RKPR @ 3Mpa','RKPR @ 6Mpa','RKPR @ 1Mpa','RKPR @ 2Mpa','RKPR @ 5Mpa','NIST webbook (3MPa)','NIST webbook (6MPa)','NIST webbook (1MPa)','NIST webbook (2MPa)','NIST webbook (5MPa)' )
% xlabel('Temperature(K)')
% ylabel('Specific enthalpy (kJ/kg)')
% title('RK-PR Specific enthalpy')
% axis([60 400 -200 400]);
% hold off

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
h = (e + P/rho)/1000+300;
end