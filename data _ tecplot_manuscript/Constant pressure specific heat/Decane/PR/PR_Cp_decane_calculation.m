%PR cp calculation
clear all
clc
load PR_Cpdata.mat
omega = 0.4884; %Acentric factor
Tc = 617.7; %K
Pc = 2103000; %Pa
R = 8.31446; %J/(mol*K)
b = 0.07780*(R*Tc/Pc); % m^3/mol
a = 0.45724*(((R^2)*(Tc^2))/Pc); % J*m^3/mol
S = 0.37464+1.54226*omega-0.26992*omega^2;% no unit
M = 142/1000; %kg/mol
Cp_01_list = zeros(71,1);
Cp_1_list = zeros(71,1);
Cp_3_list = zeros(71,1);
Cp_5_list = zeros(71,1);
Cp_2_list = zeros(71,1);
Cp_10_list = zeros(71,1);
T = 300; %k
Sos_01_list = zeros(71,1);
Sos_1_list = zeros(71,1);
Sos_3_list = zeros(71,1);
Sos_5_list = zeros(71,1);
Sos_2_list = zeros(71,1);
Sos_10_list = zeros(71,1);
h_01_list = zeros(71,1);
h_1_list = zeros(71,1);
h_3_list = zeros(71,1);
h_5_list = zeros(71,1);
h_2_list = zeros(71,1);
h_10_list = zeros(71,1);
Tr = 447;
for i = 1:1:71
    
    
    alpha =(1+S*(1-sqrt(T/Tc)))^2;% no unit
    
    Cp_0 = 2.407E2+5.099*T-(6.29E-4)*T^2-(1.07E-6)*T^3; %J/(kg*k)
    
    Cp_0_list(i,1) = Cp_0/1000;%kJ/(kg*k)
    
    Cv_0 = Cp_0 - R/M;  %J/(kg*k)
    
    e_0 = 2.407E2*T+5.099/2*T^2-(6.29E-4)/3*T^3-(1.07E-6)/4*T^4-R/M*T-(2.407E2*Tr+5.099/2*Tr^2-(6.29E-4)/3*Tr^3-(1.07E-6)/4*Tr^4-R/M*Tr);
    
    daaldT = a*(-S/sqrt(T*Tc)*(1+S*(1-sqrt(T/Tc))));%1 deriative aalpha to T;
    
    d2aaldT = a/2*(S^2/(T*Tc)+S/sqrt(T^3*Tc)*(1+S*(1-sqrt(T/Tc))));%2 derivative aalpha to temp
    %% Cp_01 calculation
    
    dpdT_01 = getdpdT(PRDensity01(i,1),R,M,b,daaldT);
    
    dpdrho_01 = getdpdrho(PRDensity01(i,1),R,M,b,T,a,alpha);
    
    Cv_01 = getCv(PRDensity01(i,1),Cv_0,d2aaldT,M,b,T);
    
    Cp_01 = getCp(PRDensity01(i,1),Cv_01,T,dpdT_01,dpdrho_01);
    
    Cp_01_list(i,1) = Cp_01/1000;%kJ/(kg*k)
    
    Sos_01 = getSos(Cp_01,Cv_01,dpdrho_01);%m/s
    
    Sos_01_list(i,1) = Sos_01;
    
    e_01 = gete(e_0,T,b,M,daaldT,a,alpha,PRDensity01(i,1));
    
    h_01=geth(e_01,0.1E6,PRDensity01(i,1));
    
    h_01_list(i,1)=h_01;
    
      %% Cp_1 calculation
    
    dpdT_1 = getdpdT(PRDensity1(i,1),R,M,b,daaldT);
    
    dpdrho_1 = getdpdrho(PRDensity1(i,1),R,M,b,T,a,alpha);
    
    Cv_1 = getCv(PRDensity1(i,1),Cv_0,d2aaldT,M,b,T);
    
    Cp_1 = getCp(PRDensity1(i,1),Cv_1,T,dpdT_1,dpdrho_1);
    
    Cp_1_list(i,1) = Cp_1/1000;%kJ/(kg*k)
    
    Sos_1 = getSos(Cp_1,Cv_1,dpdrho_1);%m/s
    
    Sos_1_list(i,1) = Sos_1;
    
    e_1 = gete(e_0,T,b,M,daaldT,a,alpha,PRDensity1(i,1));
    
    h_1=geth(e_1,1E6,PRDensity1(i,1));
    
    h_1_list(i,1)=h_1;
    %% Cp_3 calculation
    
    dpdT_3 = getdpdT(PRDensity3(i,1),R,M,b,daaldT);
    
    dpdrho_3 = getdpdrho(PRDensity3(i,1),R,M,b,T,a,alpha);
    
    Sos_0 = sqrt((Cp_0/Cv_0)*R*T/M);
    
    Sos_0_list(i) = Sos_0;
    
    Cv_3 = getCv(PRDensity3(i,1),Cv_0,d2aaldT,M,b,T);
    
    Cp_3 = getCp(PRDensity3(i,1),Cv_3,T,dpdT_3,dpdrho_3);
    
    Cp_3_list(i,1) = Cp_3/1000;%kJ/(kg*k)
    
    Sos_3 = getSos(Cp_3,Cv_3,dpdrho_3);%m/s
    
    Sos_3_list(i,1) = Sos_3;
    
    e_3 = gete(e_0,T,b,M,daaldT,a,alpha,PRDensity3(i,1));
    
    h_3=geth(e_3,3E6,PRDensity3(i,1));
    
    h_3_list(i,1)=h_3;
        
    %% Cp_5 calculation
    
    dpdT_5 = getdpdT(PRDensity5(i,1),R,M,b,daaldT);
    
    dpdrho_5 = getdpdrho(PRDensity5(i,1),R,M,b,T,a,alpha);
    
    Cv_5 = getCv(PRDensity5(i,1),Cv_0,d2aaldT,M,b,T);
    
    Cp_5 = getCp(PRDensity5(i,1),Cv_5,T,dpdT_5,dpdrho_5);
    
    Cp_5_list(i,1) = Cp_5/1000;%kJ/(kg*k)
    
    Sos_5 = getSos(Cp_5,Cv_5,dpdrho_5)%m/s
    
    Sos_5_list(i,1) = Sos_5;
    
    e_5 = gete(e_0,T,b,M,daaldT,a,alpha,PRDensity5(i,1));
    
    h_5=geth(e_5,5E6,PRDensity5(i,1));
    
    h_5_list(i,1)=h_5;

    %% Cp_2 calculation
    
    dpdT_2 = getdpdT(PRDensity2(i,1),R,M,b,daaldT);
    
    dpdrho_2 = getdpdrho(PRDensity2(i,1),R,M,b,T,a,alpha);
    
    Cv_2 = getCv(PRDensity2(i,1),Cv_0,d2aaldT,M,b,T);
    
    Cp_2 = getCp(PRDensity2(i,1),Cv_2,T,dpdT_2,dpdrho_2);
    
    Cp_2_list(i,1) = Cp_2/1000;%kJ/(kg*k)
    
    Sos_2 = getSos(Cp_2,Cv_2,dpdrho_2)%m/s
    
    Sos_2_list(i,1) = Sos_2;
    
    e_2 = gete(e_0,T,b,M,daaldT,a,alpha,PRDensity2(i,1));
    
    h_2=geth(e_2,2E6,PRDensity2(i,1));
     
    h_2_list(i,1)=h_2;
    

        %% Cp_10 calculation
    
    dpdT_10 = getdpdT(PRDensity10(i,1),R,M,b,daaldT);
    
    dpdrho_10 = getdpdrho(PRDensity10(i,1),R,M,b,T,a,alpha);
    
    Cv_10 = getCv(PRDensity10(i,1),Cv_0,d2aaldT,M,b,T);
    
    Cp_10 = getCp(PRDensity10(i,1),Cv_10,T,dpdT_10,dpdrho_10);
    
    Cp_10_list(i,1) = Cp_10/1000;%kJ/(kg*k)
    
    Sos_10 = getSos(Cp_10,Cv_10,dpdrho_10)%m/s
    
    Sos_10_list(i,1) = Sos_10;
    
    e_10 = gete(e_0,T,b,M,daaldT,a,alpha,PRDensity10(i,1));
    
    h_10=geth(e_10,10E6,PRDensity10(i,1));
     
    h_10_list(i,1)=h_10;
    
    T = T+10;
end

%%Cp vs Temp figure
figure(1);
hold on
plot(Temp,Cp_10_list,'--k');
plot(Temp,Cp_01_list,'-m');
plot(Temp,Cp_1_list,'-g');
plot(Temp,Cp_3_list,'-b');
plot(Temp,Cp_5_list,'-k');
plot(Temp,Cp_2_list,'-r');
plot(Temp,Cp_10_list,'-y');
% plot(Temp_3,NIST_Cp_3,'om');
% plot(Temp_6,NIST_Cp_6,'*g');
% plot(Temp_6,NIST_Cp_12,'xb');
% plot(Temp_6,NIST_Cp_25,'sk');
% plot(Temp_6,NIST_Cp_50,'+r');
legend('Ideal gas eos','PR @ 3Mpa','PR @ 6Mpa','PR @ 12Mpa','PR @ 25Mpa','PR @ 50Mpa','NIST webbook (3MPa)','NIST webbook (6MPa)','NIST webbook (12MPa)','NIST webbook (25MPa)','NIST webbook (50MPa)' )
xlabel('Temperature(K)')
ylabel('Constant pressure specific heat(kJ/kg/K)')
title('PR EoS Constant pressure specific heat')
% axis([50 400 0 6]);
hold off

%%speed of sound vs Temp figure
figure(2);
hold on
plot(Temp,Sos_01_list,'-m');
plot(Temp,Sos_1_list,'-g');
plot(Temp,Sos_3_list,'-b');
plot(Temp,Sos_5_list,'-k');
plot(Temp,Sos_2_list,'-r');
plot(Temp,Sos_10_list,'-y');
% plot(NIST_temp_3,NIST_Sos_3,'om');
% plot(NIST_temp_6,NIST_Sos_6,'*g');
% plot(NIST_temp_6,NIST_Sos_12,'xb');
% plot(NIST_temp_6,NIST_Sos_25,'sk');
% plot(NIST_temp_6,NIST_Sos_50,'+r');
legend('PR @ 3Mpa','PR @ 6Mpa','PR @ 12Mpa','PR @ 25Mpa','PR @ 50Mpa','NIST webbook (3MPa)','NIST webbook (6MPa)','NIST webbook (12MPa)','NIST webbook (25MPa)','NIST webbook (50MPa)' )
xlabel('Temperature(K)')
ylabel('Speed of sound (m/s)')
title('PR Speed of Sound')
hold off

%%h vs Temp figure
figure(3);
hold on
plot(Temp,h_01_list,'-m');
plot(Temp,h_1_list,'-g');
plot(Temp,h_3_list,'-b');
plot(Temp,h_5_list,'-k');
plot(Temp,h_2_list,'-r');
plot(Temp,h_10_list,'-y');
% plot(NIST_temp_3,NIST_h_3,'om');
% plot(NIST_temp_6,NIST_h_6,'*g');
% plot(NIST_temp_6,NIST_h_12,'xb');
% plot(NIST_temp_6,NIST_h_25,'sk');
% plot(NIST_temp_6,NIST_h_50,'+r');
legend('PR @ 3Mpa','PR @ 6Mpa','PR @ 12Mpa','PR @ 25Mpa','PR @ 50Mpa','NIST webbook (3MPa)','NIST webbook (6MPa)','NIST webbook (12MPa)','NIST webbook (25MPa)','NIST webbook (50MPa)' )
xlabel('Temperature(K)')
ylabel('Specific enthalpy (kJ/kg)')
title('PR Specific enthalpy')
% axis([60 400 -200 400]);
hold off

function dpdT = getdpdT(rho,R,M,b,daaldT)
dpdT = rho*R/(M-b*rho)-daaldT*rho^2/(M^2+2*b*M*rho-b^2*rho^2); %partical derivative p to Temp;
end

function Cv = getCv(rho,Cv_0,d2aaldT,M,b,T)
Cv = Cv_0 + T/(sqrt(8)*b*M)*d2aaldT*log((M+(1+sqrt(2)*b*rho))/(M+(1-sqrt(2)*b*rho)));%J/(kg*k)
end

function dpdrho = getdpdrho(rho,R,M,b,T,a,alpha)
% dpdrho = M*R*T/(M-b*rho)^2-a*alpha*rho*(2*M+b*rho)/(M*(M+b*rho)^2);%partical derivative p to Rho;
dpdrho = M*(R*T/(b*rho-M)^2-(2*a*rho*alpha*(b*rho+M))/(b^2*rho^2-2*b*M*rho-M^2)^2);
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
h = (e + P/rho)/1000+290;
end