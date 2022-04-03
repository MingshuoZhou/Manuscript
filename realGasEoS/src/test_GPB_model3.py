from utils import *
from scipy import stats
from copy import deepcopy
import pandas as pd
np.random.seed(0)

columns = ['Tr','Pr','bias corrected model', 'Nist', 'bias uncorrected model', ' upper', 'lower', 'bias','PR_alpha']

## settings for C12
Mole = 170.33 / 1000
fluid = "C12"
name = 'c12h26'

# # settings for O2
# Mole = 31.999 / 1000
# fluid = "oxygen"
# name = "o2"

# # settings for c10
# Mole = 44.01 / 1000
# fluid = "CO2"
# name = 'co2'

# # settings for c10
# Mole = 142.28 / 1000
# fluid = "decane"
# name = 'nc10h22'

# # settings for n2
# fluid = "nitrogen"
# name = 'n2'

# # settings for CH4
# Mole = 16.043 / 1000
# fluid = "CH4"
# name = 'ch4'
# ================================================
# load data
data = np.loadtxt("mech_model3/Alpha/%s.csv"%name, delimiter=',')
dim = data.shape[1]-1
X = data[:,:dim]
y = data[:,dim]
N = len(y)

para = np.loadtxt('mech_model3/Alpha/%s_para.csv'%name, delimiter=',')
𝛾 = para[0,:dim] # kernel size
σ = para[0,dim]  # kernel multiplier
θ = para[1,:]    # basis function's parameters

AS = CP.AbstractState("HEOS", fluid)
Tc = CP.PropsSI(fluid, 'Tcrit')
Pc = CP.PropsSI(fluid, 'pcrit')
a, b = PR(Tc, Pc)
omega = AS.acentric_factor()

print(𝛾, σ)
# ================================================
# define covariance
def k0(X1, X2, 𝛾=𝛾, σ=σ):
    cov = np.zeros((len(X1), len(X2)))
    for i in range(len(X1)):
        for j in range(len(X2)):
            cov[i,j] = σ**2 * np.exp(-np.sum((X1[i] - X2[j])**2/ 2 / 𝛾**2))
    return cov

# define basis function
def f(x, θ):
    return θ[0] + θ[1]*x[:,0] + θ[2]*x[:,1]

def dfdθ(x, θ):
    if len(x.shape)<=1:
        return np.array([1, x[0], x[1]]) #, x[1]
    else:
        N = len(x)
        return np.vstack([np.ones(N), x[:,0], x[:,1]]).T #, x[:,1]

K0  = k0(X, X)
df = dfdθ(X,θ)
Hm = [np.outer(df[i], df[j]) * K0[i,j] for i in range(len(X)) for j in range(len(X))]
H = np.mean(Hm, axis=0)

def h(X1, X, θ):
    K0x = k0(X1, X)
    df = dfdθ(X1,θ)
    return df.T * np.mean(K0x, axis=1)

def k(X1, X2, X, H, θ):
    return k0(X1,X2) + h(X1,X,θ).T @ np.linalg.inv(H) @ h(X2,X,θ)


# ================================================
# random queries
tool = np.arange(400,800,5)
M = len(tool)
# M = 100

# pr = np.linspace(0.047,4.7,10)
pr = np.array([0.6,1,2,3])
# pr = np.array([0.1,1,5,7,8,10])*1e6/Pc
Xnew = np.empty([M*len(pr),dim])
alpha_data = np.zeros([M*len(pr),9])
for i in range(len(pr)): 
    # Xnew[i*M:M*(i+1),0] = np.linspace(np.min(X[:,0]), np.max(X[:,0]),M)
    Xnew[i*M:M*(i+1),0] = tool/Tc
    Xnew[i*M:M*(i+1),1] = np.ones(M)*pr[i] 
alpha_data[:,0] = Xnew[:,0]
alpha_data[:,1] = Xnew[:,1]
for i in range(M*len(pr)):
    density = CP.PropsSI("D", "T", alpha_data[i,0]*Tc, "P", alpha_data[i,1]*Pc, fluid)
    V = Mole/density
    alpha_data[i,3] = (R*alpha_data[i,0]*Tc/(V-b) - alpha_data[i,1]*Pc) / a * (V*(V+b) + b*(V-b))
    alpha_data[i,8] = PR_alpha(alpha_data[i,0]*Tc, alpha_data[i,1]*Pc, Tc, Pc, omega)
    # print(density,alpha_data[i,3])

# x1lim = [np.min(X[:,0]), np.max(X[:,0])]
# x2lim = [np.min(X[:,1]), np.max(X[:,1])]
# Xnew = np.random.rand(M, dim)
# Xnew[:,0] = Xnew[:,0] * (x1lim[1] - x1lim[0]) + x1lim[0]
# Xnew[:,1] = Xnew[:,1] * (x2lim[1] - x2lim[0]) + x2lim[0]
# Xnew = Xnew[np.argsort(Xnew[:,1]),:]
# Xnew = Xnew[np.argsort(Xnew[:,0]),:]
#print(Xnewfront)

# ================================================
# GP predict
K = k(X, X, X, H, θ)
Ki = np.linalg.inv(K + 1e-4*np.diag(np.ones(len(X))))
Kxx = k(Xnew, Xnew, X, H, θ)
Kx = k(X, Xnew, X, H, θ)

mu = f(Xnew,θ) + Kx.T @ Ki @ (y-f(X,θ))
cov = Kxx - Kx.T @ Ki @ Kx
sigma = np.diag(cov)
sd = np.sqrt(np.max(sigma))
print(sd)
upper = mu + stats.norm.isf(0.05)*sd
lower = mu + stats.norm.ppf(0.05)*sd

Ypred = f(X,θ) + K.T @ Ki @ (y-f(X,θ))

alpha_data[:,2] = mu
alpha_data[:,4] = f(Xnew,θ)
alpha_data[:,5] = upper
alpha_data[:,6] = lower
alpha_data[:,7] = Kx.T @ Ki @ (y-f(X,θ))

dt = pd.DataFrame(alpha_data,columns=columns)
dt.to_csv("results/density/%s_model3.csv"%fluid, index=0)
# ================================================
# 2d plot
plt.plot(y, y, 'k--')
plt.plot(y, Ypred, 'gp', fillstyle='none')
plt.xlabel("Groundtruth")
plt.ylabel("Prediction")

plt.figure()
plt.plot(Xnew[:,0], mu, 'gp')
plt.plot(Xnew[:,0], lower, 'k--')
plt.plot(Xnew[:,0], upper, 'b--')
plt.xlabel("Groundtruth")
plt.ylabel("Prediction")
# # 3d plot
# fig, ax = plt.subplots(subplot_kw={"projection":"3d"})
# ax.scatter(X[:,0], X[:,1], y, 'r')
# ax.scatter(Xnew[:,0], Xnew[:,1], mu, 'g')
# ax.set_xlabel("Tr")
# ax.set_ylabel("Pr")
# ax.set_zlabel("Alpha")
# plt.legend(["Groundtruth ", "Prediction"])

plt.show()
