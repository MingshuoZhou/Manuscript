clear
clc
load('oxygen_alpha_R.mat')

%%
% load oxygen_sat_data.mat
omega = 0.0222; %Acentric factor
Tc = 154.581; %K
Pc = 5.0430E6; %Pa
R = 8.31446; %J/(mol*K)
b = 0.07780*(R*Tc/Pc); % m^3/mol
a = 0.45724*(((R^2)*(Tc^2))/Pc); % J*m^3/mol
M = 32/1000; %kg/mol
P = 5E6;
delta_T = 5;
for j = 1:1:1
    T = 60;
    for i = 1:1:69
        A = (a*alpha(i,j)*P(j))/(R^2*T^2);
        B = (b*P(j))/(R*T);
        z = roots([1 -(1-B) (A-2*B-3*B^2) -(A*B-B^2-B^3)]);
        
        %decide the state whether it is vap or liq
        if imag(z(3,1))==0
            z_liq(i,j) = z(3,1);
            rho_liq(i,j) = P(j)*M/(z_liq(i,j)*R*T);
            v_liq(i,j) = 1/rho_liq(i,j);
            temp_liq(i,j)=T;
        end
        if imag(z(2,1)) == 0
            z_mid(i,j) = z(2,1);
            rho_mid(i,j) = P(j)*M/(z_mid(i,j)*R*T);
            v_mid(i,j) = 1/rho_mid(i,j);
            temp_mid(i,j)=T;
        end
        if imag(z(1,1)) == 0
            z_vap(i,j) = z(1,1);
            rho_vap(i,j) = P(j)*M/(z_vap(i,j)*R*T);
            v_vap(i,j) = 1/rho_vap(i,j);
            temp_vap(i,j)=T;
            temp_liq(i,j)=0;
            tem_mid(i,j)=0;
        end
        T_save(i)= T;
        T = T + delta_T;
    end
end

for j = 1:1:1
    for i = 1:1:69
        if temp_vap(i,j)== 0
           temp_vap(i,j) = NaN;
           v_vap(i,j) = NaN;
           rho_vap(i,j)=NaN;
        end
    end
end
for j = 1:1:1
    for i = 1:1:33
        if temp_liq(i,j)== 0
           temp_liq(i,j) = NaN;
           v_liq(i,j) = NaN;
           rho_liq(i,j) =NaN;
        end
    end
end
for j = 1:1:1
    for i = 1:1:27
        if temp_mid(i,j)== 0
           temp_mid(i,j) = NaN;
           v_mid(i,j) = NaN;
           rho_mid(i,j) =NaN;
        end
    end
end

color = ["m","c","r","g","b","k"];

figure(1)
hold on 
plot([NaN;v_liq(:,1);NaN;v_mid(:,1);NaN;v_vap(:,1);NaN],[NaN;temp_liq(:,1);NaN;temp_mid(:,1);NaN;temp_vap(:,1);NaN],color(1),'linewidth',1.5);
plot([NaN;v_liq(:,2);NaN;v_mid(:,2);NaN;v_vap(:,2);NaN],[NaN;temp_liq(:,2);NaN;temp_mid(:,2);NaN;temp_vap(:,2);NaN],color(2),'linewidth',1.5);
plot([NaN;v_liq(:,3);NaN;v_mid(:,3);NaN;v_vap(:,3);NaN],[NaN;temp_liq(:,3);NaN;temp_mid(:,3);NaN;temp_vap(:,3);NaN],color(3),'linewidth',1.5);
plot([NaN;v_liq(:,4);NaN],[NaN;temp_liq(:,4);NaN],color(4),'linewidth',1.5);
plot([NaN;v_liq(:,5);NaN],[NaN;temp_liq(:,5);NaN],color(5),'linewidth',1.5);
plot([v_liq(:,6)],[temp_liq(:,6)],color(6),'linewidth',1.5);
% plot(v_liq_sat,T_sat,"rp");
% plot(v_vapor_sat,T_sat,"rp");
legend("1Mpa","3Mpa","4Mpa","5Mpa","6Mpa","10Mpa","Saturation line")
title("Constant pressure Molar Volume VS Temperture for oxygen using new model PR EoS")
xlabel("Molar Volume (m^3/mol)")
ylabel("Temperature (K)")
% axis([0 0.8*10^-3 40 160])
hold off

figure(2)
hold on 
plot([NaN;rho_liq(:,1);NaN;rho_mid(:,1);NaN;rho_vap(:,1);NaN],[NaN;temp_liq(:,1);NaN;temp_mid(:,1);NaN;temp_vap(:,1);NaN],color(1),'linewidth',1.5);
plot([NaN;rho_liq(:,2);NaN;rho_mid(:,2);NaN;rho_vap(:,2);NaN],[NaN;temp_liq(:,2);NaN;temp_mid(:,2);NaN;temp_vap(:,2);NaN],color(2),'linewidth',1.5);
plot([NaN;rho_liq(:,3);NaN;rho_mid(:,3);NaN;rho_vap(:,3);NaN],[NaN;temp_liq(:,3);NaN;temp_mid(:,3);NaN;temp_vap(:,3);NaN],color(3),'linewidth',1.5);
plot([NaN;rho_liq(:,4);NaN],[NaN;temp_liq(:,4);NaN],color(4),'linewidth',1.5);
plot([NaN;rho_liq(:,5);NaN],[NaN;temp_liq(:,5);NaN],color(5),'linewidth',1.5);
plot([rho_liq(:,6)],[temp_liq(:,6)],color(6),'linewidth',1.5);
% plot(rho_liq_sat,T_sat,"rp");
% plot(rho_vapor_sat,T_sat,"rp");
legend("1Mpa","3Mpa","4Mpa","5Mpa","6Mpa","10Mpa","Saturation line")
title("Constant pressure Density VS Temperture for oxygen using new model PR EoS")
xlabel("Density (kg/m^3)")
ylabel("Temperature (K)")
hold off