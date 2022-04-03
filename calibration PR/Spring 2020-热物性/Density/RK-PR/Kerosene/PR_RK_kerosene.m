% RK-PR EoS for kerosene
%test 1
clear all
clc
%% for decane
Tc_d = 617.7;
rho_d = 233;
Pc_d = 2.103E6;
omega_d = 0.488;
%% for cyclohexane
Tc_c = 553.64;
rho_c = 273;
Pc_c = 4.075E6;
omega_c = 0.20926;
%% for toluene
Tc_t = 591.75;
rho_t = 292;
Pc_t = 4.1263E6;
omega_t = 0.266;
%% general properties
Zc = [0.2564 0.2729 0.2638];
R = 8.31446;%J/Kmol
P = 50E6;
M = 129.72;%moleculor mass 130.216
T = 300; %k
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
b_d = getb(Tc_d,Pc_d,R,d(1),y(1));
b_c = getb(Tc_c,Pc_c,R,d(2),y(2));
b_t = getb(Tc_t,Pc_t,R,d(3),y(3));
X(1)=0.78;
X(2)=0.098;
X(3)=0.122;
b = X(1)*b_d+X(2)*b_c+X(3)*b_t;
%% delta1&2
delta1 = X(1)*delta1s(1)+X(2)*delta1s(2)+X(3)*delta1s(3);
delta2 = X(1)*delta2s(1)+X(2)*delta2s(2)+X(3)*delta2s(3);
%% k
k(1) =getk(A0,A1,B0,B1,C0,C1,Zc(1),omega_d);
k(2) =getk(A0,A1,B0,B1,C0,C1,Zc(2),omega_c);
k(3) =getk(A0,A1,B0,B1,C0,C1,Zc(3),omega_t);
%% main function
for i = 1:1:101
    
     alpha(1) = getalpha(k(1),T,Tc_d);
     alpha(2) = getalpha(k(2),T,Tc_c);
     alpha(3) = getalpha(k(3),T,Tc_t);
     
     aalpha = 0;
     
     for j =1:1:3
         for z=1:1:3
             aalpha = aalpha + X(j)*X(z)*sqrt(a(j)*a(z)*alpha(j)*alpha(z));
         end
     end

    V = roots([(-P) (R*T-P*delta1*b-P*delta2*b+b*P) (b*R*T*delta1+b*R*T*delta2-aalpha-delta1*delta2*b*b*P+delta1*b*b*P+delta2*b*b*P) (delta1*delta2*b*b*R*T+aalpha*b+delta1*delta2*b*b*b*P)]); 
    if (imag(V(1,1))==0 & V(1,1)>0)
        v1{i} = V(1,1);
        rho1{i} = M/(1000*V(1,1));
    end
    if (imag(V(2,1))==0 & V(2,1)>0)
        v2{i} = V(2,1);
        rho2{i} = M/(1000*V(2,1));
    end
    if (imag(V(3,1))==0 & V(3,1)>0) 
        v3{i} = V(3,1);
        rho3{i} = M/(1000*V(3,1));
    end
T = T+5; 
end

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

