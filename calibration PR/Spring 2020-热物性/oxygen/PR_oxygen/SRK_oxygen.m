clc
clear all
%SRK
omega = 0.0222; %Acentric factor
Tc = 154.581; %K
Pc = 5043000; %Pa
R = 8.31446 %J/Kmol
b = 0.08664*(R*Tc/Pc);
a = 0.42747*(((R^2)*(Tc^2))/Pc);
S=0.48508+1.55171*omega-0.15613*omega^2;
T = 55; %k
M = 32/1000;
% delta_P = [0.01,0.01,0.1,0.1,1,1,1,1,1,1,1,1,1,10,10,10,10,10,10,10]
P = 50;
for i = 1:1:20
    
    delta_P(i) = P/100000
    div = 10;
    
    while div > 0.0001
        P = P + delta_P(i)
        alpha =(1+S*(1-sqrt(T/Tc)))^2;
        A = (a*alpha*P)/(R^2*T^2);
        B = (b*P)/(R*T);
        z = roots([1 -1 (A-B-B^2) -(A*B)]);
        % decide the state whether it is vap or liq
        if imag(z(1,1))==0
            z_vap(i) = z(1,1);
            rho_vap(i) = P*M/(z_vap(i)*R*T);
            temp_vap(i)=T;
            
        end
        if imag(z(3,1))==0
            z_liq(i) = z(3,1);
            rho_liq(i) = P*M/(z_liq(i)*R*T);
            temp_liq(i)=T;
            
        end
        f_L = P*exp(z_liq(i)-1-log(z_liq(i)-B)-A/B*log((z_liq(i)+B)/z_liq(i)));
        f_V = P*exp(z_vap(i)-1-log(z_vap(i)-B)-A/B*log((z_vap(i)+B)/z_vap(i)));
        div = f_L-f_V;
        % m = m +1;
    end
    P_sat(i) = P*10^(-6)
    T = T+5;
end
