clear
clc
load('decane8.mat')

[numRows,numCols] = size(x)
j = linspace(low_lim_x,up_lim_x,71);

for m = 1:1:71
% h_theta xi calculation 
for i = 1:1:numRows
h_theta_1 = @(xi) exp(-(x(i)-xi).^2/gamma);
h_theta1 = integral(h_theta_1,low_lim_x,up_lim_x);
h_theta_2 = @(xi) xi.*exp(-(x(i)-xi).^2/gamma);
h_theta2 = integral(h_theta_2,low_lim_x,up_lim_x);
h_theta(1,i) = h_theta1;
h_theta(2,i) = h_theta2;
end

%h_theta x calculation which j is the input I am interested
h_theta_1j = @(xi) exp(-(j(m)-xi).^2/gamma);
h_theta1j = integral(h_theta_1j,low_lim_x,up_lim_x);
h_theta_2j = @(xi) xi.*exp(-(j(m)-xi).^2/gamma);
h_theta2j = integral(h_theta_2j,low_lim_x,up_lim_x);
h_thetaj(1) = h_theta1j;
h_thetaj(2) = h_theta2j;


%H_theta
 H_theta_1 = @(x,y) exp(-(x-y).^2./gamma);
 H_theta1 = integral2(H_theta_1,low_lim_x,up_lim_x,low_lim_x,up_lim_x);
 H_theta_2 = @(x,y) y.*exp(-(x-y).^2./gamma);
 H_theta2 = integral2(H_theta_2,low_lim_x,up_lim_x,low_lim_x,up_lim_x);
 H_theta_3 = @(x,y) x.*exp(-(x-y).^2./gamma);
 H_theta3 = integral2(H_theta_3,low_lim_x,up_lim_x,low_lim_x,up_lim_x) ;
 H_theta_4 = @(x,y) x.*y.*exp(-(x-y).^2./gamma);
 H_theta4 = integral2(H_theta_4,low_lim_x,up_lim_x,low_lim_x,up_lim_x);
 H_theta = [H_theta1,H_theta2
            H_theta3,H_theta4];
        
%Kn matrix calculation
 for i = 1:1:numRows
     for h = 1:1:numRows
         Kn(i,h) =exp(-(x(i)-x(h))^2/gamma) -[h_theta(1,i),h_theta(2,i)]*inv(H_theta)*[h_theta(1,h);h_theta(2,h)];
     end
 end
 
 for i = 1:1:numRows
     k(i) =exp(-(j(m)-x(i))^2/gamma) - [h_thetaj(1),h_thetaj(2)]*inv(H_theta)*[h_theta(1,i);h_theta(2,i)];
 end
 
 
 k_xx = exp(-(j(m)-j(m))^2/gamma) - [h_thetaj(1),h_thetaj(2)]*inv(H_theta)*[h_thetaj(1);h_thetaj(2)];
 

 v = r*ones(1,numRows);%model$g
 V = diag(v);
 f_theta = theta2*x+theta1;
 mu(m) = k*inv(Kn+V)*(Y-f_theta)

 sd(m) = nu * k_xx + r - nu * k *inv(Kn+V)*k';
end
model = theta2*j +theta1
alpha = model+mu
%figure 1 model vs experimental data
figure(3)
hold on
plot(x,Y,'bo')
plot(j,model,'r--')
plot(j,alpha,'k-')
hold off
title('Bias corrected alpha function model for n-decane')
ylabel('alpha value')
xlabel('T_r')
%%
% % load oxygen_sat_data.mat
omega = 0.489; %Acentric factor
Tc = 617.7; %K
Pc = 2.12E6; %Pa
R = 8.31446; %J/(mol*K)
b = 0.07780*(R*Tc/Pc); % m^3/mol
a = 0.45724*(((R^2)*(Tc^2))/Pc); % J*m^3/mol
M = 142.29/1000; %kg/mol
P = [0.1E6,1E6,2E6,3E6,5E6,10E6];
delta_T = 10;
for j = 1:1:6
    T = 350;
    for i = 1:1:66
        A = (a*alpha(i)*P(j))/(R^2*T^2);
        B = (b*P(j))/(R*T);
        z = roots([1 -(1-B) (A-2*B-3*B^2) -(A*B-B^2-B^3)]);
        
        %decide the state whether it is vap or liq
        if imag(z(1,1)) == 0
            z_vap(i,j) = z(1,1);
            rho_vap(i,j) = P(j)*M/(z_vap(i,j)*R*T);
            v_vap(i,j) = 1/rho_vap(i,j);
            temp_vap(i,j)=T;
        end
        if imag(z(2,1)) == 0
            z_mid(i,j) = z(2,1);
            rho_mid(i,j) = P(j)*M/(z_mid(i,j)*R*T);
            v_mid(i,j) = 1/rho_mid(i,j);
            temp_mid(i,j)=T;
            
        end
        if imag(z(3,1))==0
            z_liq(i,j) = z(3,1);
            rho_liq(i,j) = P(j)*M/(z_liq(i,j)*R*T);
            v_liq(i,j) = 1/rho_liq(i,j);
            temp_liq(i,j)=T;
        end
        T_save(i)= T;
        T = T + delta_T;
    end
end

for j = 1:1:6
    for i = 1:1:71
        if temp_vap(i,j)== 0
           temp_vap(i,j) = NaN;
           v_vap(i,j) = NaN;
           rho_vap(i,j)=NaN;
        end
    end
end
for j = 1:1:6
    for i = 1:1:39
        if temp_liq(i,j)== 0
           temp_liq(i,j) = NaN;
           v_liq(i,j) = NaN;
           rho_liq(i,j) =NaN;
        end
    end
end
for j = 1:1:2
    for i = 1:1:31
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