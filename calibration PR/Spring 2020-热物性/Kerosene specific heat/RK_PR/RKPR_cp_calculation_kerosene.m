%PR cp calculation for kerosene
clear all
clc
load RKPR_Kerosene_Cpdata.mat
%% for decane
Tc_d = 617.7; %K
rho_d = 233; %kg/m^3
Pc_d = 2.103E6; %Pa
omega_d = 0.488;  
MM(1) = 0.142;
Tr(1) = 447;
%% for cyclohexane
Tc_c = 553.64;
rho_c = 273;
Pc_c = 4.075E6;
omega_c = 0.20926;
MM(2) = 0.084;
Tr(2) = 356.15;
%% for toluene
Tc_t = 591.75;
rho_t = 292;
Pc_t = 4.1263E6;
omega_t = 0.266;
MM(3) = 0.092;
Tr(3) = 383.75;
%% general properties
Zc = [0.2564 0.2729 0.2638];
R = 8.31446;%J/Kmol
M = 130.216/1000;  %g/mol              moleculor mass 129.72 kg/mol
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
%% delta1
delta1s(1)=getdelta1(Zc(1),d1,d2,d3,d4,d5,d6);
delta1s(2)=getdelta1(Zc(2),d1,d2,d3,d4,d5,d6);
delta1s(3)=getdelta1(Zc(3),d1,d2,d3,d4,d5,d6);
%% delta2
delta2s(1)=getdelta2(delta1s(1));
delta2s(2)=getdelta2(delta1s(2));
delta2s(3)=getdelta2(delta1s(3));
%% d
d(1)=getd(delta1s(1));
d(2)=getd(delta1s(2));
d(3)=getd(delta1s(3));
%% y
y(1)=gety(delta1s(1));
y(2)=gety(delta1s(2));
y(3)=gety(delta1s(3));
%% a
a(1)=geta(Tc_d,Pc_d,R,d(1),y(1));
a(2)=geta(Tc_c,Pc_c,R,d(2),y(2));
a(3)=geta(Tc_t,Pc_t,R,d(3),y(3));
%% b
b(1) = getb(Tc_d,Pc_d,R,d(1),y(1));
b(2) = getb(Tc_c,Pc_c,R,d(2),y(2));
b(3) = getb(Tc_t,Pc_t,R,d(3),y(3));
X(1)=0.78;
X(2)=0.098;
X(3)=0.122;
bmix = X(1)*b(1)+X(2)*b(2)+X(3)*b(3);
%% delta1&2
delta1 = X(1)*delta1s(1)+X(2)*delta1s(2)+X(3)*delta1s(3);
delta2 = X(1)*delta2s(1)+X(2)*delta2s(2)+X(3)*delta2s(3);
%% k
k(1) =getk(A0,A1,B0,B1,C0,C1,Zc(1),omega_d);
k(2) =getk(A0,A1,B0,B1,C0,C1,Zc(2),omega_c);
k(3) =getk(A0,A1,B0,B1,C0,C1,Zc(3),omega_t);


T = 300; %k


%Tc = 683K 
%Pc = 2.3MPa
Cp_3_list = zeros(101,1);
Cp_6_list = zeros(101,1);
Cp_12_list = zeros(101,1);
Cp_25_list = zeros(101,1);
Cp_50_list = zeros(101,1);


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

  for i = 1:1:101
     daaldT = 0;
     d2aaldT2 = 0;  
     aalpha = 0;
     Sum1(1)=0;
     Sum1(2)=0;
     Sum1(3)=0;
      

    
     Alpha(1) = (3/(2+T/Tc_d))^k(1);
     Alpha(2) = (3/(2+T/Tc_c))^k(2);
     Alpha(3) = (3/(2+T/Tc_t))^k(3);

    
    for y =1:1:3
        for z=1:1:3
            aalpha = aalpha + X(y)*X(z)*sqrt(a(y)*a(z)*Alpha(y)*Alpha(z));
        end
    end
    
    dalidT_d = -3^k(1)*k(1)/(Tc_d*(2+T/Tc_d)^(k(1)+1));%1 deriative alpha to T for decane; 
    dalidT_c = -3^k(2)*k(2)/(Tc_c*(2+T/Tc_c)^(k(2)+1));%1 deriative alpha to T for cyclohexane  
    dalidT_t = -3^k(3)*k(3)/(Tc_t*(2+T/Tc_t)^(k(3)+1));%1 deriative alpha to T for toluene   
    dalidT(1) = dalidT_d;
    dalidT(2) = dalidT_c;
    dalidT(3) = dalidT_t;
    
     for y =1:1:3
         for z=1:1:3
             daaldT = daaldT + X(y)*X(z)*sqrt(a(y)*a(z))*0.5*(sqrt(Alpha(y)/Alpha(z))*dalidT(z)+sqrt(Alpha(z)/Alpha(y))*dalidT(y));
         end
     end

     
        d2alidT_d = a(1)*(3^k(1)*k(1)*(k(1)+1)/(Tc_d^2*(2+T/Tc_d)^(k(1)+2)));
        d2alidT_c = a(2)*(3^k(2)*k(2)*(k(2)+1)/(Tc_c^2*(2+T/Tc_c)^(k(2)+2)));
        d2alidT_t = a(3)*(3^k(3)*k(3)*(k(3)+1)/(Tc_t^2*(2+T/Tc_t)^(k(3)+2)));
        d2alidT(1) = d2alidT_d;
        d2alidT(2) = d2alidT_c;
        d2alidT(3) = d2alidT_t;
       
     for y =1:1:3
         for z=1:1:3
             d2aaldT2 = d2aaldT2 + X(y)*X(z)*sqrt(a(y)*a(z))*(1/(2*sqrt(Alpha(y)*Alpha(z)))*dalidT(y)*dalidT(z)-1/4*sqrt(Alpha(y)/Alpha(z)^3)*dalidT(z)^2-1/4*sqrt(Alpha(z)/Alpha(y)^3)*dalidT(y)^2+1/2*sqrt(Alpha(y)/Alpha(z))*d2alidT(z)+1/2*sqrt(Alpha(z)/Alpha(y))*d2alidT(y));
         end
     end

     
     Cp_0(1) = 2.407E2+5.099*T-(6.29E-4)*T^2-(1.07E-6)*T^3;%decane
     Cp_0(2) = -143.687 + 2.2338877*T+(1.5757957E-2)*T^2-(3.322767E-5)*T^3+(2.992193E-8)*T^4-(1.3118055E-11)*T^5+(2.278605E-15)*T^6;%cyclehexane
     Cp_0(3) = -310.29+5.640685*T-2.81410224E-3*T^2+2.4913887E-7*T^3;%toluene
     
     Cp0mix = X(1)*Cp_0(1)+X(2)*Cp_0(2)+X(3)*Cp_0(3);
     Cp0(i) = Cp0mix/1000;    
     Cv_0 = Cp0*1000 - R/M;  %J/(kg*k)
     
     e_0(1) = 2.407E2*T+5.099/2*T^2-(6.29E-4)/3*T^3-(1.07E-6)/4*T^4-R/MM(1)*T-(2.407E2*Tr(1)+5.099/2*Tr(1)^2-(6.29E-4)/3*Tr(1)^3-(1.07E-6)/4*Tr(1)^4-R/MM(1)*Tr(1));
     e_0(2) = -143.687*T + 2.2338877/2*T^2+(1.5757957E-2)/3*T^3-(3.322767E-5)/4*T^4+(2.992193E-8)/5*T^5-(1.3118055E-11)/6*T^6+(2.278605E-15)/7*T^7-R/MM(2)*T-(-143.687*Tr(2) + 2.2338877/2*Tr(2)^2+(1.5757957E-2)/3*Tr(2)^3-(3.322767E-5)/4*Tr(2)^4+(2.992193E-8)/5*Tr(2)^5-(1.3118055E-11)/6*Tr(2)^6+(2.278605E-15)/7*Tr(2)^7-R/MM(2)*Tr(2));
     e_0(3) = -310.29*T+5.640685/2*T^2-(2.81410224E-3)/3*T^3+(2.4913887E-7)/4*T^4-R/MM(3)*T-(-310.29*Tr(3)+5.640685/2*Tr(3)^2-(2.81410224E-3)/3*Tr(3)^3+(2.4913887E-7)/4*Tr(3)^4-R/MM(3)*Tr(3));

     e_0mix = X(1)*e_0(1)+X(2)*e_0(2)+X(3)*e_0(3);
%      for y = 1:1:3;
%          for z = 1:1:3;
%              if z ~= y
%                  Sum1(y) = Sum1(y) +(T*sqrt(a(y)*a(z))*0.5*(sqrt(Alpha(y)/Alpha(z))*dalidT(z)+sqrt(Alpha(z)/Alpha(y))*dalidT(y))-sqrt(a(y)*a(z)*Alpha(y)*Alpha(z)));
%              end
%          end
%      end
%        
%      for y = 1:1:3;
%          e_tilde(y) = e_0(y) + 2 * Sum1(y)/((delta1-delta2)*bmix*MM(y))*log((M+delta1*bmix*rho)/(M+delta2*bmix*rho)+b(y)*(T*daaldT-aalpha)/((delta1-delta2)*bmix*MM(y))*((delta1-delta2)*rho*M/((M+delta1*b*rho)*(M+delta2*b*rho))-1/bmix*log((M+delta1*bmix*rho)/(M+delta2*bmix*rho);
%      end
% %     Cp_0 = (25.48+1.52E-2*T-0.7155E-5*T^2+1.312E-9*T^3)/M; %J/(kg*k)
% %     
% %     Cp_0_list(i,1) = Cp_0/1000;%kJ/(kg*k)
%      

%      
%     e_0 = Cv_0*T;
%     
% 
%     
% %     d2aaldT = a/2*(S^2/(T*Tc)+S/sqrt(T^3*Tc)*(1+S*(1-sqrt(T/Tc))));%2 derivative aalpha to temp
% %        
     %% Cp_3 calculation
     
     dpdT_3 = getdpdT(PRRKDensity3(i,1),R,M,bmix,daaldT);
     
     dpdrho_3 = getdpdrho(PRRKDensity3(i,1),R,M,bmix,T,aalpha);
    
     Cv_3 = getCv(PRRKDensity3(i,1),Cv_0(i,1),d2aaldT2,M,bmix,T);                                                     
    
     Cp_3 = getCp(PRRKDensity3(i,1),Cv_3,T,dpdT_3,dpdrho_3);
    
     Cp_3_list(i,1) = Cp_3/1000;%kJ/(kg*k)
    
     Sos_3 = getSos(Cp_3,Cv_3,dpdrho_3);%m/s
    
     Sos_3_list(i,1) = Sos_3;
    
    e_3 = gete(e_0mix,T,bmix,M,daaldT,aalpha,PRRKDensity3(i,1),delta1,delta2);
    
    h_3=geth(e_3,3E6,PRRKDensity3(i,1));
    
    h_3_list(i,1)=h_3;
     %% Cp_6 calculation
    
    dpdT_6 = getdpdT(PRRKDensity6(i,1),R,M,bmix,daaldT);
    
    dpdrho_6 = getdpdrho(PRRKDensity6(i,1),R,M,bmix,T,aalpha);
    
    Cv_6 = getCv(PRRKDensity6(i,1),Cv_0(i,1),d2aaldT2,M,bmix,T);
    
    Cp_6 = getCp(PRRKDensity6(i,1),Cv_6,T,dpdT_6,dpdrho_6);
    
    Cp_6_list(i,1) = Cp_6/1000;%kJ/(kg*k)
    
    Sos_6 = getSos(Cp_6,Cv_6,dpdrho_6);%m/s
    
    Sos_6_list(i,1) = Sos_6;
    
    e_6 = gete(e_0mix,T,bmix,M,daaldT,aalpha,PRRKDensity6(i,1),delta1,delta2);
    
    h_6=geth(e_6,6E6,PRRKDensity6(i,1));
    
    h_6_list(i,1)=h_6;
     %% Cp_12 calculation
    
    dpdT_12 = getdpdT(PRRKDensity12(i,1),R,M,bmix,daaldT);
    
    dpdrho_12 = getdpdrho(PRRKDensity12(i,1),R,M,bmix,T,aalpha);
    
    Cv_12 = getCv(PRRKDensity12(i,1),Cv_0(i,1),d2aaldT2,M,bmix,T);
    
    Cp_12 = getCp(PRRKDensity12(i,1),Cv_12,T,dpdT_12,dpdrho_12);
    
    Cp_12_list(i,1) = Cp_12/1000;%kJ/(kg*k)
    
    Sos_12 = getSos(Cp_12,Cv_12,dpdrho_12);%m/s
        
    Sos_12_list(i,1) = Sos_12;

    e_12 = gete(e_0mix,T,bmix,M,daaldT,aalpha,PRRKDensity12(i,1),delta1,delta2);
    
    h_12=geth(e_12,12E6,PRRKDensity12(i,1));
    
    h_12_list(i,1)=h_12;
     %% Cp_25 calculation
    
    dpdT_25 = getdpdT(PRRKDensity25(i,1),R,M,bmix,daaldT);
    
    dpdrho_25 = getdpdrho(PRRKDensity25(i,1),R,M,bmix,T,aalpha);
    
    Cv_25 = getCv(PRRKDensity25(i,1),Cv_0(i,1),d2aaldT2,M,bmix,T);
    
    Cp_25 = getCp(PRRKDensity25(i,1),Cv_25,T,dpdT_25,dpdrho_25);
    
    Cp_25_list(i,1) = Cp_25/1000;%kJ/(kg*k)
    
    Sos_25 = getSos(Cp_25,Cv_25,dpdrho_25);%m/s
    
    Sos_25_list(i,1) = Sos_25;

    e_25 = gete(e_0mix,T,bmix,M,daaldT,aalpha,PRRKDensity25(i,1),delta1,delta2);
    
    h_25=geth(e_25,25E6,PRRKDensity25(i,1));
    
    h_25_list(i,1)=h_25;
     %% Cp_50 calculation
    
    dpdT_50 = getdpdT(PRRKDensity50(i,1),R,M,bmix,daaldT);
    
    dpdrho_50 = getdpdrho(PRRKDensity50(i,1),R,M,bmix,T,aalpha);
    
    Cv_50 = getCv(PRRKDensity50(i,1),Cv_0(i,1),d2aaldT2,M,bmix,T);
    
    Cp_50 = getCp(PRRKDensity50(i,1),Cv_50,T,dpdT_50,dpdrho_50);
    
    Cp_50_list(i,1) = Cp_50/1000;%kJ/(kg*k)
    
    Sos_50 = getSos(Cp_50,Cv_50,dpdrho_50);%m/s
    
    Sos_50_list(i,1) = Sos_50;
   
    e_50 = gete(e_0mix,T,bmix,M,daaldT,aalpha,PRRKDensity50(i,1),delta1,delta2);
    
    h_50=geth(e_50,50E6,PRRKDensity50(i,1));
    
    h_50_list(i,1)=h_50;
%     
    T = T+5;
 end

%%Cp vs Temp figure
figure(1);
hold on
plot(PRRKTemp,Cp0,'--k');
plot(PRRKTemp,Cp_3_list,'-m');
plot(PRRKTemp,Cp_6_list,'-g');
plot(PRRKTemp,Cp_12_list,'-b');
plot(PRRKTemp,Cp_25_list,'-k');
plot(PRRKTemp,Cp_50_list,'-r');
plot(NISTTemp,NIST_Cp_3,'om');
plot(NISTTemp,NIST_Cp_6,'*g');
plot(NISTTemp,NIST_Cp_12,'xb');
plot(NISTTemp,NIST_Cp_25,'sk');
plot(NISTTemp,NIST_Cp_50,'+r');
legend('Ideal gas eos','RKPR @ 3Mpa','RKPR @ 6Mpa','RKPR @ 12Mpa','RKPR @ 25Mpa','RKPR @ 50Mpa','NIST webbook (3MPa)','NIST webbook (6MPa)','NIST webbook (12MPa)','NIST webbook (25MPa)','NIST webbook (50MPa)' )
xlabel('Temperature(K)')
ylabel('Constant pressure specific heat(kJ/kg/K)')
title('RKPR EoS Constant pressure specific heat')
axis([300 800 1.5 6]);
hold off

%%speed of sound vs Temp figure
figure(2);
hold on
plot(PRRKTemp,Sos_3_list,'-m');
plot(PRRKTemp,Sos_6_list,'-g');
plot(PRRKTemp,Sos_12_list,'-b');
plot(PRRKTemp,Sos_25_list,'-k');
plot(PRRKTemp,Sos_50_list,'-r');
plot(NISTTemp,NIST_Sos_3,'om');
plot(NISTTemp,NIST_Sos_6,'*g');
plot(NISTTemp,NIST_Sos_12,'xb');
plot(NISTTemp,NIST_Sos_25,'sk');
plot(NISTTemp,NIST_Sos_50,'+r');
legend('PRRK @ 3Mpa','PRRK @ 6Mpa','PRRK @ 12Mpa','PRRK @ 25Mpa','PRRK @ 50Mpa','NIST webbook (3MPa)','NIST webbook (6MPa)','NIST webbook (12MPa)','NIST webbook (25MPa)','NIST webbook (50MPa)' )
xlabel('Temperature(K)')
ylabel('Speed of sound (m/s)')
title('PRRK Speed of Sound')
axis([300 800 0 2000]);
hold off

%%h vs Temp figure
figure(3);
hold on
plot(PRRKTemp,h_3_list,'-m');
plot(PRRKTemp,h_6_list,'-g');
plot(PRRKTemp,h_12_list,'-b');
plot(PRRKTemp,h_25_list,'-k');
plot(PRRKTemp,h_50_list,'-r');
plot(NISTTemp,NIST_h_3,'om');
plot(NISTTemp,NIST_h_6,'*g');
plot(NISTTemp,NIST_h_12,'xb');
plot(NISTTemp,NIST_h_25,'sk');
plot(NISTTemp,NIST_h_50,'+r');
legend('PRRK @ 3Mpa','PRRK @ 6Mpa','PRRK @ 12Mpa','PRRK @ 25Mpa','PRRK @ 50Mpa','NIST webbook (3MPa)','NIST webbook (6MPa)','NIST webbook (12MPa)','NIST webbook (25MPa)','NIST webbook (50MPa)' )
xlabel('Temperature(K)')
ylabel('Specific enthalpy (kJ/kg)')
title('PRRK Specific Enthalpy')
% axis([60 400 -200 400]);
hold off

function delta1s = getdelta1(Zc,d1,d2,d3,d4,d5,d6)
delta1s = d1+d2*(d3-1.168*Zc)^d4+d5*(d3-1.168*Zc)^d6;
end
function delta2s = getdelta2(delta1)
delta2s = (1-delta1)/(1+delta1);
end
function d = getd(delta1)
d = (1+delta1^2)/(1+delta1);
end
function y = gety(delta1)
y = 1 +(2*(1+delta1))^(1/3)+(4/(1+delta1))^(1/3);
end
function a = geta(Tc,Pc,R,d,y)
a = ((3*y^2+3*y*d+d^2+d-1)/((3*y+d-1)^2))*(R^2*Tc^2/Pc);
end
function b = getb(Tc,Pc,R,d,y)
b = (1/(3*y+d-1))*(R*Tc/Pc);
end
function k =getk(A0,A1,B0,B1,C0,C1,Zc,omega)
k = (1.168*Zc*A1+A0)*omega^2+(1.1168*Zc*B1+B0)*omega+(1.168*Zc*C1+C0);
end
function alpha = getalpha(k,T,Tc)
alpha = (3/(2+T/Tc))^k;
end

function dpdT = getdpdT(rho,R,M,b,daaldT)
dpdT = rho*R/(M-b*rho)-daaldT*rho^2/(M^2+2*b*M*rho-b^2*rho^2); %partical derivative p to Temp;
end

function Cv = getCv(rho,Cv_0,d2aaldT,M,b,T)
Cv = Cv_0 + T/(sqrt(8)*b*M)*d2aaldT*log((M+(1+sqrt(2)*b*rho))/(M+(1-sqrt(2)*b*rho)));%J/(kg*k)
end

function dpdrho = getdpdrho(rho,R,M,b,T,aalpha)
dpdrho = M*(R*T/(b*rho-M)^2-(2*rho*aalpha*(b*rho+M))/(b^2*rho^2-2*b*M*rho-M^2)^2);%partical derivative p to Rho;
end

function Cp = getCp(rho,Cv,T,Q,O)
Cp = Cv + T/rho^2*Q^2/O;%J/(kg*k)
end

function Sos = getSos(Cp,Cv,dpdrho)
Sos = sqrt(Cp/Cv*dpdrho);%m/s
end

function e=gete(e_0,T,b,M,daaldT,aalpha,rho,delta1,delta2)
e = e_0+1/((delta1-delta2)*b*M)*(T*daaldT-aalpha)*log((M+delta1*b*rho)/(M+delta2*b*rho)); %J
end

function h=geth(e,P,rho)
h = (e + P/rho)/1000+330;
end