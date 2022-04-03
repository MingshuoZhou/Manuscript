clear all
clc
%kerosene properties calculation SRK
%%for decane
Tc_d = 617.7;
rho_d = 233;
Pc_d = 2.103E6;
omega_d = 0.488;
%% for cyclohexane
Tc_c = 553.64;
rho_c = 273;
Pc_c = 4.075E6;
omega_c = 0.20926;
%% for toluene
Tc_t = 591.75;
rho_t = 292;
Pc_t = 4.1263E6;
omega_t = 0.266;
%% general property
P = 50E6; %pa
%Tc = 683K 
%Pc = 2.3MPa
R = 8.31446 %J/Kmol
a(1) = geta(Tc_d,Pc_d,R);
b_d = getb(Tc_d,Pc_d,R);
a(2) = geta(Tc_c,Pc_c,R);
b_c = getb(Tc_c,Pc_c,R);
a(3) = geta(Tc_t,Pc_t,R);
b_t = getb(Tc_t,Pc_t,R);
X(1)=0.78
X(2)=0.098;
X(3)=0.122;
b = X(1)*b_d+X(2)*b_c+X(3)*b_t;
%b = 0.07780*(R*Tc/Pc)
%a = 0.45724*(((R^2)*(Tc^2))/Pc)
S_d=gets(omega_d);
S_c=gets(omega_c);
S_t=gets(omega_t);
T = 300;

for i = 1:1:101
    Alpha(1)= getalpha(S_d,T,Tc_d)
    Alpha(2) = getalpha(S_c,T,Tc_c)
    Alpha(3)= getalpha(S_d,T,Tc_t)
    aalpha = 0;

    for j =1:1:3
        for k=1:1:3
            aalpha = aalpha + X(j)*X(k)*sqrt(a(j)*a(k)*Alpha(j)*Alpha(k));
        end
    end
    A = (aalpha*P)/(R^2*T^2);
    B = (b*P)/(R*T);
    z = roots([1 -(1-B) (A-2*B-3*B^2) -(A*B-B^2-B^3)]);
    % decide the state whether it is vap or liq
    if (z(1,1)>0)
        if imag(z(1,1))==0
            z_vap = z(1,1);
            rho_vap = P*130.216/(z_vap*R*T*1000);
            temp_vap{i}=T;
            ComR_vap{i}=z(1,1);
            RHO_vap{i}=rho_vap;
        end
    end
   if (z(2,1)>0)
        if imag(z(2,1))==0
            z_u = z(2,1);
            rho_u = P*130.216/(z_u*R*T*1000);
            temp_u{i}=T;
            ComR_u{i}=z(2,1);
            RHO_u{i}=rho_u;
        end
   end
    if (z(3,1)>0)
        if imag(z(3,1))==0
        z_liq = z(3,1)
        rho_liq = P*130.216/(z_liq*R*T*1000);
        temp_liq{i}=T;
        ComR_liq{i}=z(3,1);
        RHO_liq{i}=rho_liq;
        end
    end
T = T+5; 
end

function a = geta(Tc,Pc,R)
a = 0.42747*(((R^2)*(Tc^2))/Pc);
end
function b = getb(Tc,Pc,R)
b = 0.08664*(R*Tc/Pc);
end
function S =gets(omega)
S=0.48508+1.55171*omega-0.15613*omega^2;
end
function alpha = getalpha(S,T,Tc)
alpha =(1+S*(1-sqrt(T/Tc)))^2;
end
 