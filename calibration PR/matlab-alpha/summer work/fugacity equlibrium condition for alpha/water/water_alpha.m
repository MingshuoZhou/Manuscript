clc
clear all
%PR
load Water_P_sat_data.mat
omega = 0.344; %Acentric factor
Tc = 647.096; %K
Pc = 22.064E6; %Pa
R = 8.31446 %J/Kmol
b = 0.07780*(R*Tc/Pc)
a = 0.45724*(((R^2)*(Tc^2))/Pc)
S = 0.37464+1.54226*omega-0.26992*omega^2;
T = 275; %k
M = 18/1000;
P = P*10^6;
alpha = 1.8;
 for i = 1:1:19
    
    div1 = 10000000;
    div2 = 1000000;
    
    
    while div1 > div2
        div1 = div2;
        alpha = alpha - 0.00000001;
        A = (a*alpha*P(i))/(R^2*T^2);
        B = (b*P(i))/(R*T);
        z = roots([1 -(1-B) (A-2*B-3*B^2) -(A*B-B^2-B^3)]);
        % decide the state whether it is vap or liq
        if imag(z(1,1))==0
            z_vap(i) = z(1,1);
            rho_vap(i) = P(i)*M/(z_vap(i)*R*T);
            temp_vap(i)=T;
            
        end
        if imag(z(3,1))==0
            z_liq(i) = z(3,1);
            rho_liq(i) = P(i)*M/(z_liq(i)*R*T);
            temp_liq(i)=T;
            
        end
        f_L = P(i)*exp(z_liq(i)-1-log(z_liq(i)-B)-A/(2.8284*B)*log((z_liq(i)+2.4142*B)/(z_liq(i)-0.4142*B)));
        f_V = P(i)*exp(z_vap(i)-1-log(z_vap(i)-B)-A/(2.8284*B)*log((z_vap(i)+2.4142*B)/(z_vap(i)-0.4142*B)));
        div2 = abs(f_L-f_V);
        
        % m = m +1;
    end
    alpha_value(i) = alpha;
    div_value(i) = div2;
    alpha_soave(i) = (1+S*(1-sqrt(T/Tc)))^2;
    T = T+20;
end