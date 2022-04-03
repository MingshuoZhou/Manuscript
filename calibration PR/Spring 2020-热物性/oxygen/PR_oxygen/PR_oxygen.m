clc
clear all
%PR
omega = 0.0222; %Acentric factor
Tc = 154.581; %K
Pc = 5.0430E6; %Pa
R = 8.31446 %J/Kmol
b = 0.07780*(R*Tc/Pc)
a = 0.45724*(((R^2)*(Tc^2))/Pc)
S = 0.37464+1.54226*omega-0.26992*omega^2;
T = 55; %k
M = 32/1000;

P = 100;
for i = 1:1:20
    
    delta_P(i) = P/10000
    div = 10;
    
    while div > 0.001
        P = P + delta_P(i)
        alpha =(1+S*(1-sqrt(T/Tc)))^2;
        A = (a*alpha*P)/(R^2*T^2);
        B = (b*P)/(R*T);
        z = roots([1 -(1-B) (A-2*B-3*B^2) -(A*B-B^2-B^3)]);
        % decide the state whether it is vap or liq
        if imag(z(1,1))==0
            z_vap(i) = z(1,1);
            rho_vap = P*M/(z_vap(i)*R*T);
            temp_vap(i)=T;
            
        end
        if imag(z(3,1))==0
            z_liq(i) = z(3,1);
            rho_liq(i) = P*M/(z_liq(i)*R*T);
            temp_liq(i)=T;
            
        end
        f_L = P*exp(z_liq(i)-1-log(z_liq(i)-B)-A/(2.8284*B)*log((z_liq(i)+2.4142*B)/(z_liq(i)-0.4142*B)));
        f_V = P*exp(z_vap(i)-1-log(z_vap(i)-B)-A/(2.8284*B)*log((z_vap(i)+2.414*B)/(z_vap(i)-0.4142*B)));
        div = f_L-f_V;
        % m = m +1;
    end
    P_sat(i) = P*10^(-6)
    T = T+5;
end
