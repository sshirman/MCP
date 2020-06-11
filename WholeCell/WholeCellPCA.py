import numpy as np
from scipy.integrate import odeint
import sys
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
from colour import Color
from sklearn.decomposition import PCA as sklearnPCA
from scipy.constants import *


print('Loading Paramters\n')
kcat0 = 3.e2 #Reaction rate for S->I for a single reaction site (reaction/s)
N0 = 1.5e3 #Number of reaction sites for S->I reaction
K0 = 0.5e3 #Half max concentration for S->I reaction (microM)
kcat1 = 55. #Reaction rate for I->P for a single reaction site (reaction/s)
N1 = 2.5e3 #Number of reaction sites for I->P reaction
K1 = 15.e3 #Half max concentration for I->P reaction(microM)
perm0c = 1.e-5 #Permeability of S through compartment (cm/s)
perm1c = 1.e-5 #Permeability of I through compartment (cm/s)
perm2c = 1.e-5 #Permeability of P through compartment (cm/s)
perm0b = 1.e-3 #Permeability of S through cell (cm/s)
perm1b = 1.e-3 #Permeability of I through cell (cm/s)
perm2b = 1.e-3 #Permeability of P through cell (cm/s)
Rc = 1.e-5 #Radius of compartment (cm)
Rb = 5.e-5 #Effective Radius of cell (cm)
Diff = 1.e-5 #Diffusion coefficient (cm^2/s)
SInit = 70. #Initial concentration of external substrate (mM)

V0T = kcat0*N0*1.e9/(Avogadro*(4.*pi/3.)*Rc**3) #Max reaction rate S->I (microM/s)
V1T = kcat1*N1*1.e9/(Avogadro*(4.*pi/3.)*Rc**3) #Max reaction rate I->P (microM/s)


ngrid = 10 #number of desired grid points for the integration inside the cell (but outside the MCP)
Vratio = 100. #This is the multiplier between the external and cell volumes in the solution (Vratio = Vexternal/Vcell)


print('Defining Derivatives\n')
def SDeriv(x,t):
    Mb = (Rb/Rc)**3/3.
    Vext = 4*pi*Mb*Vratio
    Mc = 1./3.
    DeltaM = np.divide((Mb-Mc),(ngrid))
    M = np.linspace(Mc-DeltaM/2.,Mb+DeltaM/2.,ngrid+2)
    D = Diff/(Rc**2)
    k0c = perm0c/Rc
    k1c = perm1c/Rc
    k2c = perm2c/Rc
    k0b = perm0b/Rc
    k1b = perm1b/Rc
    k2b = perm2b/Rc
    rc = 1.
    rb = Rb/Rc
    M[0] = Mc
    M[-1] = Mb
    v0 = V0/K0
    v1 = V1/K1
    kappa01 = K0/K1
    assert len(x) == 3*(2+(ngrid))
    d = np.zeros((len(x)))
    d[0] = -(x[0]*v0)/(1. + 1000.*x[0]) + 3.*k0c*(x[3] - x[0])/rc
    d[1] = (x[0]*v0*kappa01)/(1. + 1000.*x[0]) -(x[1]*v1)/(1. + 1000.*x[1]) + 3.*k1c*(x[4] - x[1])/rc
    d[2] = (x[1]*v1)/(2.*(1. + 1000.*x[1])) + 3.*k2c*(x[5] - x[2])/rc
    
    d[-3] = 12.*pi*Mb*k0b*(x[-6] - x[-3])/(rb*Vext)
    d[-2] = 12.*pi*Mb*k1b*(x[-5] - x[-2])/(rb*Vext)
    d[-1] = 12.*pi*Mb*k2b*(x[-4] - x[-1])/(rb*Vext)
    
    d[6:-6:3] = ((3**(4./3.))*D/(DeltaM)**2)*((M[2:-2])**(4./3.)*(x[9:-3:3]-x[6:-6:3])- (M[1:-3])**(4./3.)*(x[6:-6:3]-x[3:-9:3]))
    d[7:-5:3] = ((3**(4./3.))*D/(DeltaM)**2)*((M[2:-2])**(4./3.)*(x[10:-2:3]-x[7:-5:3])- (M[1:-3])**(4./3.)*(x[7:-5:3]-x[4:-8:3]))
    d[8:-4:3] = ((3**(4./3.))*D/(DeltaM)**2)*((M[2:-2])**(4./3.)*(x[11:-1:3]-x[8:-4:3])- (M[1:-3])**(4./3.)*(x[8:-4:3]-x[5:-7:3]))
    
    '''for k in range(2,(ngrid+1)):
        i = 3*k
        d[i] =   ((3**(4./3.))*D/(DeltaM)**2)*((M[k])**(4./3.)*(x[i+3]-x[i])- (M[k-1])**(4./3.)*(x[i]-x[i-3]))
        d[i+1] = ((3**(4./3.))*D/(DeltaM)**2)*((M[k])**(4./3.)*(x[i+4]-x[i+1])- (M[k-1])**(4./3.)*(x[i+1]-x[i-2]))
        d[i+2] = ((3**(4./3.))*D/(DeltaM)**2)*((M[k])**(4./3.)*(x[i+5]-x[i+2])- (M[k-1])**(4./3.)*(x[i+2]-x[i-1]))''' 
    
    d[3] = ((3**(4./3.))*D/(DeltaM)**2)*((M[1])**(4./3.)*(x[6]-x[3])) - 3.*k0c*(x[3] - x[0])*Mc/(rc*DeltaM)
    d[4] = ((3**(4./3.))*D/(DeltaM)**2)*((M[1])**(4./3.)*(x[7]-x[4])) - 3.*k1c*(x[4] - x[1])*Mc/(rc*DeltaM)
    d[5] = ((3**(4./3.))*D/(DeltaM)**2)*((M[1])**(4./3.)*(x[8]-x[5])) - 3.*k2c*(x[5] - x[2])*Mc/(rc*DeltaM)
    
    d[-6] = -3.*k0b*(x[-6] - x[-3])*Mb/(rb*DeltaM) - ((3**(4./3.))*D/(DeltaM)**2)*((M[ngrid])**(4./3.)*(x[-6]-x[-9]))
    d[-5] = -3.*k1b*(x[-5] - x[-2])*Mb/(rb*DeltaM) - ((3**(4./3.))*D/(DeltaM)**2)*((M[ngrid])**(4./3.)*(x[-5]-x[-8]))
    d[-4] = -3.*k2b*(x[-4] - x[-1])*Mb/(rb*DeltaM) - ((3**(4./3.))*D/(DeltaM)**2)*((M[ngrid])**(4./3.)*(x[-4]-x[-7]))
    return d


print('Defining Initial Conditions\n')
y0 = np.zeros((ngrid*3+6))
y0[-3] = SInit/K0 #y0[-3] gives the initial state of the external substrate. The /K0 turns the value into a dimensionless quantity
Kscale = np.array([K0,K1,K1]) #Creates a vector of metabolite scaling values
Equilibrium = Vratio*SInit/(Vratio+1.) #This takes into account that there is no substrate inside the cell at the start and finds the equilibrium concentration when it has diffused into the cell


print('Loading Sweep Simulation Data\n')
EqTimes = np.loadtxt('Times.txt')
Pars = np.loadtxt('Parameters.txt')
Substrate = np.loadtxt('Substrate.txt')
Intermediate = np.loadtxt('Intermediate.txt')
Product = np.loadtxt('Product.txt')


print('Defining Phase Space\n')
S = Substrate/SInit
P = 2.*Product/(Equilibrium)
I = Intermediate/(Equilibrium)

FinalState = np.array([0.,0.,1.])
InitialState = np.array([1.,0.,0.])
Trajectory = FinalState-InitialState
NormTraj = Trajectory/np.linalg.norm(Trajectory)

States = np.transpose(np.array([S - InitialState[0],I - InitialState[1],P-InitialState[2]]),(1,2,0))/np.linalg.norm(Trajectory)

print('Converting Simulation to Phase Space variables\n')
pSweep = np.dot(States,NormTraj)
ParVec = np.reshape(pSweep,(len(pSweep),500,1))*np.reshape(NormTraj,(1,-1))
perpendicular = States - ParVec

print('Running PCA\n')
Components = np.empty((100,3,3))
Values = np.empty((100,3))
ValueRatio = np.empty((100,3))
count = 0
for i in np.arange(0,1,0.01):
    testset= perpendicular[np.logical_and(pSweep>i, pSweep<i+0.01)]
    testset_centered = testset - np.mean(testset,axis = 0)
    pca = sklearnPCA()
    pca.fit(testset_centered)
    Components[count] = pca.components_
    Values[count] = pca.explained_variance_
    ValueRatio[count] = pca.explained_variance_ratio_
    count +=1

print('Plotting Eigenvalues\n')
plt.plot(Values[1:95,0],'.')
plt.plot(Values[1:95,1],'.')
plt.plot(Values[1:95,2],'.')
plt.show()

print('Saving PCA Results\n')
np.save('Eigenvectors.npy', Components)
np.save('Eigenvalues.npy', Values)
np.save('Eigenratios.npy', ValueRatio)