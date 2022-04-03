clear
clc
fluid1='CO2'
T0=250;%K
T1=780;%K
delta_T=5;
times=(T1-T0)/delta_T+1;
alpha=xlsread('predictions.csv',1,'C2:C214');
P=xlsread('pt.csv',1,'A1:A213');
T=xlsread('pt.csv',1,'B1:B213');
times=length(P);
%%
omega = 0.22394; %Acentric factor
Tc = refpropm('T','C',0,'',0,fluid1); %K
Pc = refpropm('P','C',0,'',0,fluid1)*1000; %Pa
R = 8.31446; %J/(mol*K)
b = 0.07780*(R*Tc/Pc); % m^3/mol
a = 0.45724*(((R^2)*(Tc^2))/Pc); % J*m^3/mol
M = 44.0095/1000; %kg/mol
for j = 1:1:1
    for i = 1:1:times
        A = (a*alpha(i,j)*P(i,j))/(R^2*T(i,j)^2);
        B = (b*P(i,j))/(R*T(i,j));
        z = roots([1 -(1-B) (A-2*B-3*B^2) -(A*B-B^2-B^3)]);
        zwrite(:,i)=z;
        
        %decide the state whether it is vap or liq
        if imag(z(3,1))==0
            z_liq(i,j) = z(3,1);
            rho_liq(i,j) = P(i,j)*M/(z_liq(i,j)*R*T(i,j));
            v_liq(i,j) = 1/rho_liq(i,j);
            temp_liq(i,j)=T(i,j);
        end
        if imag(z(2,1)) == 0
            z_mid(i,j) = z(1,1);
            rho_mid(i,j) = P(i,j)*M/(z_mid(i,j)*R*T(i,j));
            v_mid(i,j) = 1/rho_mid(i,j);
            temp_mid(i,j)=T(i,j);
        end
        if imag(z(1,1)) == 0
            z_vap(i,j) = z(1,1);
            rho_vap(i,j) = P(i,j)*M/(z_vap(i,j)*R*T(i,j));
            v_vap(i,j) = 1/rho_vap(i,j);
            temp_vap(i,j)=T(i,j);
            temp_liq(i,j)=0;
            temp_mid(i,j)=0;
        end
    end
end
%%
%存储
k=1;
for j = 1:1:1
    for i = 1:1:times
        if temp_liq(i,j)== 0
           temp_liq(i,j) = NaN;
           v_liq(i,j) = NaN;
           rho_liq(i,j) =NaN;
        else
            density_prediction(i,1)=rho_liq(i,j);
        end
    end
end
for j = 1:1:1
    for i = 1:times
        if temp_mid(i,j)== 0
           temp_mid(i,j) = NaN;
           v_mid(i,j) = NaN;
           rho_mid(i,j) =NaN;
        else
           density_prediction(i,1)=rho_mid(i,j);
        end
    end
end
for j = 1:1:1
    for i = 1:1:times
        if temp_vap(i,j)== 0
           temp_vap(i,j) = NaN;
           v_vap(i,j) = NaN;
           rho_vap(i,j)=NaN;
        else
           density_prediction(i,1)=rho_vap(i,j);
        end
    end
end
xlswrite('density_predictions.csv',density_prediction,1,'B1')
%%
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