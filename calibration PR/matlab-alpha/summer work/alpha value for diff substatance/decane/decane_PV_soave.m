clear all
clc
% % load oxygen_sat_data.mat
omega = 0.489; %Acentric factor
Tc = 617.7; %K
Pc = 2.12E6; %Pa
R = 8.31446; %J/(mol*K)
b = 0.07780*(R*Tc/Pc); % m^3/mol
a = 0.45724*(((R^2)*(Tc^2))/Pc); % J*m^3/mol
M = 142.29/1000; %kg/mol
S = 0.37464+1.54226*omega-0.26992*omega^2;% no unit
T = [300,440,580,720,860,1000];
delta_P = 2E5;
for j = 1:1:6
    P = 0.5E6;
    for i = 1:1:51
        alpha = sqrt(1+S*(1-(T(j)/Tc)^0.5))
        A = (a*alpha*P)/(R^2*T(j)^2);
        B = (b*P)/(R*T(j));
        z = roots([1 -(1-B) (A-2*B-3*B^2) -(A*B-B^2-B^3)]);
        
        %decide the state whether it is vap or liq
        if imag(z(1,1)) == 0
            z_vap(i,j) = z(1,1);
            rho_vap(i,j) = P*M/(z_vap(i,j)*R*T(j));
            v_vap(i,j) = 1/rho_vap(i,j);
            pres_vap(i,j)=P;
        end
        if imag(z(2,1)) == 0
            z_mid(i,j) = z(2,1);
            rho_mid(i,j) = P*M/(z_mid(i,j)*R*T(j));
            v_mid(i,j) = 1/rho_mid(i,j);
            pres_mid(i,j)=P;
            
        end
        if imag(z(3,1))==0
            z_liq(i,j) = z(3,1);
            rho_liq(i,j) = P*M/(z_liq(i,j)*R*T(j));
            v_liq(i,j) = 1/rho_liq(i,j);
            pres_liq(i,j)=P;
        end
        P_save(i)= P;
        P = P + delta_P;
    end
end

for j = 1:1:6
    for i = 1:1:51
        if pres_vap(i,j)== 0
           pres_vap(i,j) = NaN;
           v_vap(i,j) = NaN;
           rho_vap(i,j)=NaN;
        end
    end
end
for j = 1:1:4
    for i = 1:1:51
        if pres_liq(i,j)== 0
           pres_liq(i,j) = NaN;
           v_liq(i,j) = NaN;
           rho_liq(i,j) =NaN;
        end
    end
end
for j = 1:1:3
    for i = 1:1:6
        if pres_mid(i,j)== 0
           pres_mid(i,j) = NaN;
           v_mid(i,j) = NaN;
           rho_mid(i,j) =NaN;
        end
    end
end

color = ["m","c","r","g","b","k"];

figure(1)
hold on 
plot([NaN;v_liq(:,1);NaN;v_mid(:,1);NaN;v_vap(:,1);NaN],[NaN;pres_liq(:,1);NaN;pres_mid(:,1);NaN;pres_vap(:,1);NaN],color(1),'linewidth',1.5);
plot([NaN;v_liq(:,2);NaN;v_mid(:,2);NaN;v_vap(:,2);NaN],[NaN;pres_liq(:,2);NaN;pres_mid(:,2);NaN;pres_vap(:,2);NaN],color(2),'linewidth',1.5);
plot([NaN;v_liq(:,3);NaN;v_mid(:,3);NaN;v_vap(:,3);NaN],[NaN;pres_liq(:,3);NaN;pres_mid(:,3);NaN;pres_vap(:,3);NaN],color(3),'linewidth',1.5);
plot([NaN;v_liq(:,4);NaN;v_mid(:,4);NaN;v_vap(:,4);NaN],[NaN;pres_liq(:,4);NaN;pres_mid(:,4);NaN;pres_vap(:,4);NaN],color(4),'linewidth',1.5);
plot([NaN;v_liq(:,5);NaN;v_mid(:,5);NaN;v_vap(:,5);NaN],[NaN;pres_liq(:,5);NaN;pres_mid(:,5);NaN;pres_vap(:,5);NaN],color(5),'linewidth',1.5);
plot([NaN;v_liq(:,6);NaN;v_mid(:,6);NaN;v_vap(:,6);NaN],[NaN;pres_liq(:,6);NaN;pres_mid(:,6);NaN;pres_vap(:,6);NaN],color(6),'linewidth',1.5);
% plot([NaN;v_liq(:,7);NaN;v_mid(:,7);NaN;v_vap(:,7);NaN],[NaN;pres_liq(:,7);NaN;pres_mid(:,7);NaN;pres_vap(:,7);NaN],'linewidth',1.5);
% plot(v_liq_sat,P_sat,"rp");
% plot(v_vapor_sat,P_sat,"rp");
legend("60k","75k","90k","105k","120k","135k","150K","saturation line")
title("Constant Temperture Molar Volume VS Pressure for oxygen using Soave PR EoS")
xlabel("Molar Volume (m^3/mol)")
ylabel("Pressure(Pa)")
axis([0 4*10^(-3) 0.5*10^6 5*10^6])
hold off

figure(2)
hold on 
plot([NaN;rho_liq(:,1);NaN;rho_mid(:,1);NaN;rho_vap(:,1);NaN],[NaN;pres_liq(:,1);NaN;pres_mid(:,1);NaN;pres_vap(:,1);NaN],color(1),'linewidth',1.5);
plot([NaN;rho_liq(:,2);NaN;rho_mid(:,2);NaN;rho_vap(:,2);NaN],[NaN;pres_liq(:,2);NaN;pres_mid(:,2);NaN;pres_vap(:,2);NaN],color(2),'linewidth',1.5);
plot([NaN;rho_liq(:,3);NaN;rho_mid(:,3);NaN;rho_vap(:,3);NaN],[NaN;pres_liq(:,3);NaN;pres_mid(:,3);NaN;pres_vap(:,3);NaN],color(3),'linewidth',1.5);
plot([NaN;rho_liq(:,4);NaN;rho_mid(:,4);NaN;rho_vap(:,4);NaN],[NaN;pres_liq(:,4);NaN;pres_mid(:,4);NaN;pres_vap(:,4);NaN],color(4),'linewidth',1.5);
plot([NaN;rho_liq(:,5);NaN;rho_mid(:,5);NaN;rho_vap(:,5);NaN],[NaN;pres_liq(:,5);NaN;pres_mid(:,5);NaN;pres_vap(:,5);NaN],color(5),'linewidth',1.5);
plot([NaN;rho_liq(:,6);NaN;rho_mid(:,6);NaN;rho_vap(:,6);NaN],[NaN;pres_liq(:,6);NaN;pres_mid(:,6);NaN;pres_vap(:,6);NaN],color(6),'linewidth',1.5);
% plot([NaN;rho_liq(:,7);NaN;rho_mid(:,7);NaN;rho_vap(:,7);NaN],[NaN;pres_liq(:,7);NaN;pres_mid(:,7);NaN;pres_vap(:,7);NaN],'linewidth',1.5);
% plot(rho_liq_sat,P_sat,"rp");
% plot(rho_vapor_sat,P_sat,"rp");
legend("60k","75k","90k","105k","120k","135k","150K","saturation line")
title("Constant Temperture Density VS Pressure for oxygen using Soave PR EoS")
xlabel("Density (kg/m^3)")
ylabel("Pressure(Pa)")
hold off