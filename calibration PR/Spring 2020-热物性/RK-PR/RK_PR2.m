%PR_RK 2
clc
clear all
P = 3E6; %pa
omega = 0.0222; %Acentric factor
Tc = 154.581; %K
Pc = 5.0430E6; %Pa
R = 8.31446 %J/Kmol
Zc = 0.289;
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
M = 32;%moleculor mass for oxygen
delta1 = d1+d2*(d3-1.168*Zc)^d4+d5*(d3-1.168*Zc)^d6;
delta2 = (1-delta1)/(1+delta1);
u = delta1+delta2;
w = delta1*delta2;
d = (1+delta1^2)/(1+delta1);
y = 1+(2*(1+delta1))^(1/3)+(4/(1+delta1))^(1/3);
a = ((3*y^2+3*y*d+d^2+d-1)/(3*y+d-1)^2)*(R^2*Tc^2/Pc);
b = 1/(3*y+d-1)*(R^2*Tc^2/Pc);
k = (1.168*Zc*A1+A0)*omega^2+(1.168*Zc*B1+B0)*omega+(1.168*Zc*C1+C0);
T = 60; %k
for i = 1:1:170
    alpha = (3/(2+T/Tc))^k
    rho=fsolve(@(rho) P-(rho*R*T)/(M-b*rho)-(a*alpha*rho^2)/((M+delta1*b*rho)*(M+delta2*b*rho)),500)
    T=T+2;
end
%     rho = roots([(-P*delta1*delta2*b^3-a*alpha*b-delta1*delta2*b^2*R*T) 
%         (delta1*delta2*b^2*M*P-b*M^2*P-M*delta1*b^2*P-P*M*delta2*b^2-M*delta1*b*R*T-M*delta2*b*R*T+M*a*alpha) 
%         (M^2*delta1*b*P+M^2*delta2*b*P-R*T*M^2)
%         (M^3*P)]);
%     if imag(rho(1,1))==0
%         rho1{i} =(rho(1,1));
%     end
%     if imag(rho(2,1))==0
%         rho2{i} = (rho(2,1));
%     end
%     if imag(rho(3,1))==0
%         rho3{i} = (rho(3,1));
%     end



