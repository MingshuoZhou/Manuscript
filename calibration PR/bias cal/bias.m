clear
clc
fluid1='CO2';
load('carbon dioxide.mat')

[numRows,numCols] = size(x)
n=30;
j = linspace(low_lim_x,up_lim_x,n);

% h_theta xi calculation 
for i = 1:1:numRows
    h_theta_1 = @(xi) exp(-(x(i)-xi).^2/gamma);
    h_theta1 = integral(h_theta_1,low_lim_x,up_lim_x)
    h_theta_2 = @(xi) xi.*exp(-(x(i)-xi).^2/gamma);
    h_theta2 = integral(h_theta_2,low_lim_x,up_lim_x)
    h_theta(1,i) = h_theta1;
    h_theta(2,i) = h_theta2;
end
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
for m = 1:1:n
    %h_theta x calculation which j is the input I am interested
    h_theta_1j = @(xi) exp(-(j(m)-xi).^2/gamma);
    h_theta1j = integral(h_theta_1j,low_lim_x,up_lim_x);
    h_theta_2j = @(xi) xi.*exp(-(j(m)-xi).^2/gamma);
    h_theta2j = integral(h_theta_2j,low_lim_x,up_lim_x);
    h_thetaj(1) = h_theta1j;
    h_thetaj(2) = h_theta2j;
 
    for i = 1:1:numRows
        k(i) =exp(-(j(m)-x(i))^2/gamma) - [h_thetaj(1),h_thetaj(2)]*inv(H_theta)*[h_theta(1,i);h_theta(2,i)];
    end
 
    k_xx = exp(-(j(m)-j(m))^2/gamma) - [h_thetaj(1),h_thetaj(2)]*inv(H_theta)*[h_thetaj(1);h_thetaj(2)];
 

    v = r*ones(1,numRows);%model$g
    V = diag(v);
    f_theta = theta2*x+theta1;
    mu(m) = k*inv(Kn+V)*(Y-f_theta);
    sd(m) = nu * k_xx + r - nu * k *inv(Kn+V)*k';
end
model = theta2*j +theta1
alpha = model+mu
%figure 1 model vs experimental data
figure(1)
hold on
plot(x,Y,'bo')
plot(j,model,'r--')
plot(j,alpha,'k-')
hold off
title('Bias corrected alpha function model for CO_2')
ylabel('alpha value')
xlabel('T_r')
 
% %% PR EoS Saturation Pressure Calcaulation
% % omega = 0.224; %Acentric factor
% % Tc = 304.1; %K
% % Pc = 7380000; %Pa
% % R = 8.31446 %J/Kmol
% % M = 44.01/1000;
 delta_T = (up_lim_x - low_lim_x)*Tc/29;
% 
b = 0.07780*(R*Tc/Pc);
a = 0.45724*(((R^2)*(Tc^2))/Pc);
% T = low_lim_x * Tc; %k
% 
% 
% P = 5; %initial P
% for i = 1:1:28
%     
%     
%     div = 6;
%     
%     while div > 0.01 
%         if div > 2
%             delta_P(i) = P/1000000;
%         else
%             delta_P(i) = P/10000000000;
%         end
%         
%         P = P + delta_P(i);
%         A = (a*alpha(i)*P)/(R^2*T^2);
%         B = (b*P)/(R*T);
%         z = roots([1 -(1-B) (A-2*B-3*B^2) -(A*B-B^2-B^3)]);
%         % decide the state whether it is vap or liq
%                 if imag(z(1,1))==0
%                     z_vap(i) = z(1,1);
%         %             rho_vap(i) = P*M/(z_vap(i)*R*T);
%         %             temp_vap(i)=T;
%         %
%                 end
%                 if imag(z(3,1))==0
%                     z_liq(i) = z(3,1);
%         %             rho_liq(i) = P*M/(z_liq(i)*R*T);
%         %             temp_liq(i)=T;
%                 end
%         
%         f_L = P*exp(z_liq(i)-1-log(z_liq(i)-B)-A/(2.8284*B)*log((z_liq(i)+2.4142*B)/(z_liq(i)-0.4142*B)));
%         f_V = P*exp(z_vap(i)-1-log(z_vap(i)-B)-A/(2.8284*B)*log((z_vap(i)+2.4142*B)/(z_vap(i)-0.4142*B)));
%         div = abs(f_L-f_V)
%         
%     end
%     P_sat(i) = P*10^(-6)
%     diviation(i) = div
%     Tr_psat(i) = T;   
%     T = T + delta_T;
% 
% end
% figure(2)
% hold on
% plot(Tr_psat,P_sat,'r--')
% plot(j,P_exp/10^6,'b-')
% hold off
% title('Saturation pressure prediction')
% ylabel('Saturation Pressure(Pa)')
% xlabel('T_r')
%% density caluclation
T = low_lim_x * Tc; %k
        
for i = 1:1:30
    
    
    A = (a*alpha(i)*P_exp(i))/(R^2*T^2);
    B = (b*P_exp(i))/(R*T);
    z = roots([1 -(1-B) (A-2*B-3*B^2) -(A*B-B^2-B^3)]);
    zwrite(:,i)=z;
    
    rho_trueliq(i)=refpropm('D','T',T,'P',P_exp(i)/1000,fluid1);
        
    %decide the state whether it is vap or liq
    if imag(z(1,1)) == 0
        z_vap(i) = z(1,1);
        rho_vap(i) = P_exp(i)*M/(z_vap(i)*R*T);
        temp_vap(i)=T;
    else
        z_vap(i) =0;
        rho_vap(i) = 0;
        temp_vap(i)=0;
    end
    if imag(z(3,1))==0
        z_liq(i) = z(3,1);
        rho_liq(i) = P_exp(i)*M/(z_liq(i)*R*T);
        temp_liq(i)=T;
    else
        z_liq(i) = 0;
        rho_liq(i) = 0;
        temp_liq(i)=0;
    end
    T_save(i)= T;
    T = T + delta_T;
end
    figure(3)
    hold on
    plot(j,rho_liq,'r--')
    plot(j,rho_vap,'b-')
    hold off
    title('Density prediction')
    ylabel('Density (kg/m^3)')
    xlabel('T_r')
