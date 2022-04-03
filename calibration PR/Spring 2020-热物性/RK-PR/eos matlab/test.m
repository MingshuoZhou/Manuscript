%SRK
clc
clear all
P = (150+14.74)/14.74; %pa
omega = 0.227; %Acentric factor
Tc = 314.2; %K
Pc = 70.9; %Pa
R = 0.08205 %J/Kmol
b = 0.08664*(R*Tc/Pc)
a = 0.42747*(((R^2)*(Tc^2))/Pc)
S = 0.48508+1.55171*omega-0.15613*omega^2;
T = 298; %k

for i = 1:1:1
alpha =(1+S*(1-sqrt(T/Tc)))^2
A = (a*alpha*P)/(R^2*T^2)
B = (b*P)/(R*T)
z = roots([1 -1 (A-B-B^2) (-A*B)])
v =(z(1,1)*R*T)/P
time{i}=T;
CR{i}=z;
RHO{i}=v;
T = T+2; 
end