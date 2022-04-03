%PR cp calculation
clear all
clc
load PR_data.mat
omega = 0.4884; %Acentric factor
Tc = 617.7; %K
Pc = 2103000; %Pa
R = 8.31446; %J/(mol*K)
b = 0.07780*(R*Tc/Pc); % m^3/mol
a = 0.45724*(((R^2)*(Tc^2))/Pc); % J*m^3/mol
S = 0.37464+1.54226*omega-0.26992*omega^2;% no unit
M = 142/1000; %kg/mol
Cp_3_list = zeros(101,1);
Cp_6_list = zeros(101,1);
Cp_12_list = zeros(101,1);
Cp_25_list = zeros(101,1);
Cp_50_list = zeros(101,1);
Cp_0_list = zeros(101,1);
T = 300; %k
Sos_3_list = zeros(101,1);
Sos_6_list = zeros(101,1);
Sos_12_list = zeros(101,1);
Sos_25_list = zeros(101,1);
Sos_50_list = zeros(101,1);
h_3_list = zeros(101,1);
h_6_list = zeros(101,1);
h_12_list = zeros(101,1);
h_25_list = zeros(101,1);
h_50_list = zeros(101,1);
Tr = 447;
for i = 1:1:101
    
    alpha =(1+S*(1-sqrt(T/Tc)))^2;% no unit
    
    Cp_0 = 2.407E2+5.099*T-(6.29E-4)*T^2-(1.07E-6)*T^3; %J/(kg*k) 这个表哪里查完整的
    
    Cp_0_list(i,1) = Cp_0/1000;%kJ/(kg*k)
    
    Cv_0 = Cp_0 - R/M;  %J/(kg*k)
    
    e_0 = 2.407E2*T+5.099/2*T^2-(6.29E-4)/3*T^3-(1.07E-6)/4*T^4-R/M*T-(2.407E2*Tr+5.099/2*Tr^2-(6.29E-4)/3*Tr^3-(1.07E-6)/4*Tr^4-R/M*Tr);%Tr？
    
    daaldT = a*(-S/sqrt(T*Tc)*(1+S*(1-sqrt(T/Tc))));%1 deriative aalpha to T;
    
    d2aaldT = a/2*(S^2/(T*Tc)+S/sqrt(T^3*Tc)*(1+S*(1-sqrt(T/Tc))));%2 derivative aalpha to temp
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
    
    %% Cp_6 calculation
    
    dpdT_6 = getdpdT(PRDensity6(i,1),R,M,b,daaldT);
    
    dpdrho_6 = getdpdrho(PRDensity6(i,1),R,M,b,T,a,alpha);
    
    Cv_6 = getCv(PRDensity6(i,1),Cv_0,d2aaldT,M,b,T);
    
    Cp_6 = getCp(PRDensity6(i,1),Cv_6,T,dpdT_6,dpdrho_6);
    
    Cp_6_list(i,1) = Cp_6/1000;%kJ/(kg*k)
    
    Sos_6 = getSos(Cp_6,Cv_6,dpdrho_6)%m/s
    
    Sos_6_list(i,1) = Sos_6;
    
    e_6 = gete(e_0,T,b,M,daaldT,a,alpha,PRDensity6(i,1));
    
    h_6=geth(e_6,6E6,PRDensity6(i,1));
    
    h_6_list(i,1)=h_6;
    %% Cp_12 calculation
    
    dpdT_12 = getdpdT(PRDensity12(i,1),R,M,b,daaldT);
    
    dpdrho_12 = getdpdrho(PRDensity12(i,1),R,M,b,T,a,alpha);
    
    Cv_12 = getCv(PRDensity12(i,1),Cv_0,d2aaldT,M,b,T);
    
    Cp_12 = getCp(PRDensity12(i,1),Cv_12,T,dpdT_12,dpdrho_12);
    
    Cp_12_list(i,1) = Cp_12/1000;%kJ/(kg*k)
    
    Sos_12 = getSos(Cp_12,Cv_12,dpdrho_12)%m/s
        
    Sos_12_list(i,1) = Sos_12;
    
    e_12 = gete(e_0,T,b,M,daaldT,a,alpha,PRDensity12(i,1));
    
    h_12=geth(e_12,12E6,PRDensity12(i,1));
     
    h_12_list(i,1)=h_12;
        
    %% Cp_25 calculation
    
    dpdT_25 = getdpdT(PRDensity25(i,1),R,M,b,daaldT);
    
    dpdrho_25 = getdpdrho(PRDensity25(i,1),R,M,b,T,a,alpha);
    
    Cv_25 = getCv(PRDensity25(i,1),Cv_0,d2aaldT,M,b,T);
    
    Cp_25 = getCp(PRDensity25(i,1),Cv_25,T,dpdT_25,dpdrho_25);
    
    Cp_25_list(i,1) = Cp_25/1000;%kJ/(kg*k)
    
    Sos_25 = getSos(Cp_25,Cv_25,dpdrho_25)%m/s
    
    Sos_25_list(i,1) = Sos_25;
    
    e_25 = gete(e_0,T,b,M,daaldT,a,alpha,PRDensity25(i,1));
    
    h_25=geth(e_25,25E6,PRDensity25(i,1));
    
    h_25_list(i,1)=h_25;

    %% Cp_50 calculation
    
    dpdT_50 = getdpdT(PRDensity50(i,1),R,M,b,daaldT);
    
    dpdrho_50 = getdpdrho(PRDensity50(i,1),R,M,b,T,a,alpha);
    
    Cv_50 = getCv(PRDensity50(i,1),Cv_0,d2aaldT,M,b,T);
    
    Cp_50 = getCp(PRDensity50(i,1),Cv_50,T,dpdT_50,dpdrho_50);
    
    Cp_50_list(i,1) = Cp_50/1000;%kJ/(kg*k)
    
    Sos_50 = getSos(Cp_50,Cv_50,dpdrho_50)%m/s
    
    Sos_50_list(i,1) = Sos_50;
    
    e_50 = gete(e_0,T,b,M,daaldT,a,alpha,PRDensity50(i,1));
    
    h_50=geth(e_50,50E6,PRDensity50(i,1));
     
    h_50_list(i,1)=h_50;
    
    T = T+5;
end

%%Cp vs Temp figure
figure(1);
hold on
plot(PRTemp,Cp_0_list,'--k');
plot(PRTemp,Cp_3_list,'-m');
plot(PRTemp,Cp_6_list,'-g');
plot(PRTemp,Cp_12_list,'-b');
plot(PRTemp,Cp_25_list,'-k');
plot(PRTemp,Cp_50_list,'-r');
plot(NISTTemp,NIST_Cp_3,'om');
plot(NISTTemp,NIST_Cp_6,'*g');
plot(NISTTemp,NIST_Cp_12,'xb');
plot(NISTTemp,NIST_Cp_25,'sk');
plot(NISTTemp,NIST_Cp_50,'+r');
legend('Ideal gas eos','PR @ 3Mpa','PR @ 6Mpa','PR @ 12Mpa','PR @ 25Mpa','PR @ 50Mpa','NIST webbook (3MPa)','NIST webbook (6MPa)','NIST webbook (12MPa)','NIST webbook (25MPa)','NIST webbook (50MPa)' )
xlabel('Temperature(K)')
ylabel('Constant pressure specific heat(kJ/kg/K)')
title('PR EoS Constant Pressure Specific Heat for Decane')
axis([300 800 1.5 6]);
hold off

%%speed of sound vs Temp figure
figure(2);
hold on
plot(PRTemp,Sos_0_list,'--k');
plot(PRTemp,Sos_3_list,'-m');
plot(PRTemp,Sos_6_list,'-g');
plot(PRTemp,Sos_12_list,'-b');
plot(PRTemp,Sos_25_list,'-k');
plot(PRTemp,Sos_50_list,'-r');
plot(NISTTemp,NIST_Sos_3,'om');
plot(NISTTemp,NIST_Sos_6,'*g');
plot(NISTTemp,NIST_Sos_12,'xb');
plot(NISTTemp,NIST_Sos_25,'sk');
plot(NISTTemp,NIST_Sos_50,'+r');
legend('Ideal gas eos','PR @ 3Mpa','PR @ 6Mpa','PR @ 12Mpa','PR @ 25Mpa','PR @ 50Mpa','NIST webbook (3MPa)','NIST webbook (6MPa)','NIST webbook (12MPa)','NIST webbook (25MPa)','NIST webbook (50MPa)' )
xlabel('Temperature(K)')
ylabel('Speed of sound (m/s)')
title('PR Speed of Sound for Decane')
hold off

%h vs Temp figure
figure(3);
hold on
plot(PRTemp,h_3_list,'-m');
plot(PRTemp,h_6_list,'-g');
plot(PRTemp,h_12_list,'-b');
plot(PRTemp,h_25_list,'-k');
plot(PRTemp,h_50_list,'-r');
plot(NISTTemp,NIST_h_3,'om');
plot(NISTTemp,NIST_h_6,'*g');
plot(NISTTemp,NIST_h_12,'xb');
plot(NISTTemp,NIST_h_25,'sk');
plot(NISTTemp,NIST_h_50,'+r');
legend('PR @ 3Mpa','PR @ 6Mpa','PR @ 12Mpa','PR @ 25Mpa','PR @ 50Mpa','NIST webbook (3MPa)','NIST webbook (6MPa)','NIST webbook (12MPa)','NIST webbook (25MPa)','NIST webbook (50MPa)' )
xlabel('Temperature(K)')
ylabel('Specific enthalpy (kJ/kg)')
title('PR Specific enthalpy')
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
h = (e + P/rho)/1000;
end