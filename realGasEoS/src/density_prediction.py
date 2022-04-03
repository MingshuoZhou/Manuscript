from pyparsing import col
from utils import *
from copy import deepcopy
import pandas as pd

mechs = ['mech_model2/nDodecane_temp.yaml','mech_model2/nDodecane_AlphaGP.yaml']#, 'mech_model2/nDodecane_Reitz.yaml', 'mech_model2/nDodecane_Reitz.yaml'
names = ['nDodecane_PR','nDodecane_PR_ALphaGP']#, 'nDodecane_RK', 'nDodecane_IG'
lines = ['o', 'd', 'v', '.']
colors = ['r', 'b', 'g', 'c']
columns = ['T', 'P', 'NIST','PR', 'New' , 'PR_rate', 'New_rate']

# ================================================
# settings

#settings for C12
fluid = "C12"
X = {'c12h26':1}
T_step = 5
D_step = 2000
T_lo, T_hi = 400, 800
Tc = CP.PropsSI(fluid, 'Tcrit')
Pc = CP.PropsSI(fluid, 'pcrit')
P_arr = np.array([0.6,1.0,2,3]) * Pc

# settings for o2
# fluid = "oxygen"
# X = {'o2':1}
# T_step = 5
# D_step = 2000
# T_lo, T_hi = 60, 400
# P_arr = np.array([0.1,1,3,5,7,10]) * 1e6

#settings for C7H8
# fluid = "toluene"
# X = {'C6H5CH3':1}
# T_step = 5
# D_step = 50
# T_lo, T_hi = 400, 1000
# P_arr = np.array([0.1,1,3,5,7,10]) * 1e6

#settings for decane
# fluid = "decane"
# X = {'nc10h22':1}
# T_step = 5
# D_step = 1000
# T_lo, T_hi = 360, 1000
# P_arr = np.array([0.1,1,2,3,5,10]) * 1e6

# settings for N2
# fluid = "nitrogen"
# X = {'n2':1}
# T_step = 5
# D_step = 500
# T_lo, T_hi = 100, 400
# Tc = CP.PropsSI(fluid, 'Tcrit')
# Pc = CP.PropsSI(fluid, 'pcrit')
# P_arr = np.array([0.6,1,2,3]) * Pc

# # settings for c6h12
# fluid = "Cyclohexane"
# X = {'c6h12':1}
# T_step = 5
# D_step = 50
# T_lo, T_hi = 400,1000
# P_arr = np.array([0.1,1,3,5,10]) * 1e6

# # settings for ch4
# fluid = "CH4"
# X = {'ch4':1}
# T_step = 5
# D_step = 10000
# T_lo, T_hi = 100,600
# Tc = CP.PropsSI(fluid, 'Tcrit')
# Pc = CP.PropsSI(fluid, 'pcrit')
# P_arr = np.array([0.6,1.0,2.0,3.0]) * Pc

# settings for h2o
# fluid = "H2O"
# X = {'h2o':1}
# T_step = 5
# D_step = 20
# T_lo, T_hi = 600,1000
# P_arr = np.array([0.1,1,2,3,5,10]) * 1e6

# # settings for co2
# fluid = "CO2"
# X = {'co2':1}
# T_step = 5
# D_step = 2000
# T_lo, T_hi = 220,800
# P_arr = np.array([0.1,1,5,7,8,10]) * 1e6

# ================================================
# get adaptive TP list and NIST data
TPD_arr = []
column = 3
for P in P_arr:
    TPD_arr += get_TPD_under_P(fluid, P, T_lo, T_hi, T_step, D_step)
TPD_arr = np.array(TPD_arr)

TPD_data = np.zeros([TPD_arr.shape[0],len(columns)])
TPD_data[:,0] = TPD_arr[:,0]
TPD_data[:,1] = TPD_arr[:,1]
TPD_data[:,2] = TPD_arr[:,2]

# ================================================
# Density
fig = plt.figure()
plt.plot(TPD_arr[:,0], TPD_arr[:,2], 'ks', label="NIST", fillstyle='none')

TPD_rate = []
for k,name in enumerate(names):
    gas = ct.Solution(mechs[k], name)
    TPD_calc = deepcopy(TPD_arr)
    
    t0 = time.time()
    for i,(T,P,_) in enumerate(TPD_calc):
        gas.TPX = T,P, X
        TPD_calc[i,2] = gas.density
    print("Density cost of %-20s = %.5f s"%(name, time.time()-t0))
    TPD_rate = np.abs((TPD_calc[:,2]-TPD_arr[:,2])/TPD_arr[:,2]*100)
    plt.plot(TPD_calc[:,0], TPD_calc[:,2], colors[k]+lines[k], label=name, alpha=0.8, fillstyle='none')
    TPD_data[:,column] = TPD_calc[:,2]
    TPD_data[:,column+2] = TPD_rate
    column = column + 1

dt = pd.DataFrame(TPD_data, columns=columns)
dt.to_csv("results/%s_model3.csv"%fluid, index=0)

plt.xlabel("Temperature [K]")
plt.ylabel("Density [kg/m^3]")
plt.legend()
plt.savefig("figs/model3/PRAlphaGP_%s_Density.png"%fluid)


plt.show()
