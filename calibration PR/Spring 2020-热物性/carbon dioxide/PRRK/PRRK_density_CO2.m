%test 1
%RK-PR eos for CO2
clear all
clc
Zc = 0.274;
omega = 0.224;
R = 8.31446;%J/Kmol
Tc = 304.1; %K
Pc = 7.38E6; %Pa
P = 10E6;
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
M = 44.1;%moleculor mass
T = 220; %k
for i = 1:1:59
    alpha = (3/(2+T/Tc))^k
    V = roots([(-P) (R*T-P*delta1*b-P*delta2*b+b*P) (b*R*T*delta1+b*R*T*delta2-a*alpha-delta1*delta2*b*b*P+delta1*b*b*P+delta2*b*b*P) (delta1*delta2*b*b*R*T+a*alpha*b+delta1*delta2*b*b*b*P)])
    V_write(:,i)=V;
    if imag(V(1,1))==0
        v1{i} = V(1,1);
        rho1(i,:) = M/(1000*V(1,1));
    end
    if imag(V(2,1))==0
        v2{i} = V(2,1);
        rho2{i} = M/(1000*V(2,1));
    end
    if imag(V(3,1))==0
        v3{i} = V(3,1);
        rho3(i,:) = M/(1000*V(3,1));
    end
    T=T+10;
end
