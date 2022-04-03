from tarfile import PAX_FIELDS
from utils import *

# settings for C12
# fluid = "C12"
# name = "c12h26"
# M = 170.33 / 1000
# Pcrit = CP.PropsSI(fluid, 'pcrit')
# Tcrit = CP.PropsSI(fluid, 'Tcrit')
# a, b = PR(Tcrit, Pcrit)
# T_step = 5
# D_step = 50
# T_lo, T_hi = 400, 800
# # P_arr = np.array([0.8,1.0,1.5,2.0]) * Pcrit
# # # P_arr = np.array([0.5, 1, 1.817, 3, 5, 7.5, 10, 15, 20, 25, 30, 35, 40, 45, 50]) * 1e6
# P_arr = np.linspace(0.8, 4, 10) * Pcrit

# settings for CH4
fluid = "CH4"
name = "ch4"
M = 16.043 / 1000
Pcrit = CP.PropsSI(fluid, 'pcrit')
Tcrit = CP.PropsSI(fluid, 'Tcrit')
a, b = PR(Tcrit, Pcrit)
T_step = 5
D_step = 10
T_lo, T_hi = 100, 600
# P_arr = np.linspace(2, 20, 10) * 1e6
P_arr = np.array([0.8,0.9,1.0]) * Pcrit

# # settings for O2
# fluid = "oxygen"
# name = "o2"
# M = 31.999 / 1000
# Pcrit = CP.PropsSI(fluid, 'pcrit')
# Tcrit = CP.PropsSI(fluid, 'Tcrit')
# a, b = PR(Tcrit, Pcrit)
# T_step = 5
# D_step = 40
# T_lo, T_hi = 60, 400
# P_arr = np.array([0.8,0.9,1.0])*Pcrit
# P_arr = np.array([0.01984,0.1,0.1984,0.3968,0.5952,0.79365,0.992,1.1,1.3889,1.5873,1.78,1.9841])*Pcrit
# P_arr = np.array([0.1,1,3,5,7,10]) * 1e6
# print(P_arr/Pcrit)

# settings for C7H8
# fluid = "toluene"
# name = "C6H5CH3"
# M = 92.138 / 1000
# Pcrit = CP.PropsSI(fluid, 'pcrit')
# Tcrit = CP.PropsSI(fluid, 'Tcrit')
# a, b = PR(Tcrit, Pcrit)
# T_step = 5
# D_step = 50
# T_lo, T_hi = 400, 1000
# P_arr = np.array([0.1,1,2,3,4,5,7,9,10]) * 1e6
# P_arr = np.linspace(0.1, 10, 20) * 1e6

# # settings for decane
# fluid = "decane"
# name = "nc10h22"
# M = 142.28 / 1000
# Pcrit = CP.PropsSI(fluid, 'pcrit')
# Tcrit = CP.PropsSI(fluid, 'Tcrit')
# a, b = PR(Tcrit, Pcrit)
# T_step = 5
# D_step = 50
# T_lo, T_hi = 350, 1000
# P_arr = np.array([0.5,0.71,0.95])*Pcrit
# P_arr = np.array([0.04717,0.4717,0.9434,1.4151,1.8868,2.3595,3.3019,3.7736,4.7170,])*Pcrit
# P_arr = np.array([3,4,5,6,7,8,9,10]) * 1e6
# P_arr = np.linspace(0.1,1.5,2) * 1e6

# settings for N2
# fluid = "nitrogen"
# name = "n2"
# M = 28.013 / 1000
# Pcrit = CP.PropsSI(fluid, 'pcrit')
# Tcrit = CP.PropsSI(fluid, 'Tcrit')
# a, b = PR(Tcrit, Pcrit)
# T_step = 5
# D_step = 20
# T_lo, T_hi = 100, 500
# P_arr = np.array([0.5,1.0,1.5,2.0]) * Pcrit
# P_arr = np.linspace(0.5,4,10)*Pcrit

# # settings for c6h12
# fluid = "Cyclohexane"
# name = "c6h12"
# M = 84.161 / 1000
# Pcrit = CP.PropsSI(fluid, 'pcrit')
# Tcrit = CP.PropsSI(fluid, 'Tcrit')
# a, b = PR(Tcrit, Pcrit)
# T_step = 5
# D_step = 50
# T_lo, T_hi = 300,1000
# P_arr = np.array([0.1,1,2,3,5,7,8,10]) * 1e6

# settings for h2o
# fluid = "H2O"
# name = "h2o"
# M = 18.015 / 1000
# Pcrit = CP.PropsSI(fluid, 'pcrit')
# Tcrit = CP.PropsSI(fluid, 'Tcrit')
# a, b = PR(Tcrit, Pcrit)
# T_step = 5
# D_step = 20
# T_lo, T_hi = 600,1000
# P_arr = np.array([0.1,1,2,3,5,10]) * 1e6

# # settings for co2
# fluid = "CO2"
# name = "co2"
# M = 44.01 / 1000
# Pcrit = CP.PropsSI(fluid, 'pcrit')
# Tcrit = CP.PropsSI(fluid, 'Tcrit')
# a, b = PR(Tcrit, Pcrit)
# T_step = 5
# D_step = 20
# T_lo, T_hi = 220,800
# P_arr = np.array([0.82, 0.89, 0.95]) * Pcrit
# P_arr = np.array([0.01355,0.1355,0.271,0.54201,0.67751,0.94851,1.08401,1.35501]) * Pcrit

# load omega
AS = CP.AbstractState("HEOS", fluid)
omega = AS.acentric_factor()

print("Critical Properties: \r\nTc: %f \r\nPc: %f"%(Tcrit, Pcrit))
print("YAML \r\na: %.4e \r\nb: %.4e"%PR(Tcrit, Pcrit, R_u=R*1e6))
print("acentric-factor:", omega)
TPD_arr = []
for P in P_arr:
    TPD_arr += get_TPD_under_P(fluid, P, T_lo, T_hi, T_step, D_step)

TPD_arr = np.array(TPD_arr)

T = TPD_arr[:,0]
P = TPD_arr[:,1]
D = TPD_arr[:,2]
V = M / D
Alpha = (R*T/(V-b) - P) / a * (V*(V+b) + b*(V-b))
PR_Alpha  = PR_alpha(T, P, Tcrit, Pcrit, omega)

# ===============
# # show results
plt.figure()
plt.plot(T, D, 's', label="NIST Desity", alpha=0.5)

# plt.figure()
# plt.plot(T, V, 'p', label="NIST Desity", alpha=0.5)

plt.figure()
plt.plot(T, Alpha, 's', label="NIST alpha", alpha=0.5)
plt.plot(T, PR_Alpha, 'o', label="PR alpha", alpha=0.5)
plt.xlabel("T")
plt.ylabel("Alpha")
plt.legend()


# ===============
# # save results
Data = np.zeros(TPD_arr.shape)
Data[:,0] = T / Tcrit
Data[:,1] = P / Pcrit
Data[:,2] = Alpha
np.savetxt("mech_model3/Alpha/%s.csv"%name, Data, delimiter=', ')

# only save Temperature
# Data = np.zeros((len(TPD_arr), 2))
# Data[:,0] = T / Tcrit
# Data[:,1] = Alpha
# np.savetxt("mech_model2/Alpha/%s.csv"%name, Data, delimiter=', ')

plt.show()
