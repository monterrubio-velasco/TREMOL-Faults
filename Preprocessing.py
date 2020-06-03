import numpy as np
from InitialLoadConfiguration import *
from FBMaftershocsMAIN import *
from pytictoc import TicToc
t = TicToc() #create instance of class
#def Preprocessi(B)ng:

t.tic() #Start timer

print('TREMOL active Faults and aftershocks generation')
print('Marisol Monterrubio-Velasco')
print('marisol.monterrubio@bsc.es')

VectorP = [0.0]  # [0.0, 0.08, 0.16, 0.24, 0.32, 0.38]
VectorPiFrac = [0.9] #, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0]
VectorPiBG = [0.65]


VectorN = [180] #, 60]
rho =  30
FaultConfigurations = ('Northridge') #, Landers 'Northridge', 'HectorMine', 'NewZealand')
NumStatistics = 1  # Number of Statistical  realizations each (VectorN, VectorP, VectorPiFrac,)
ResultsTot = np.zeros((NumStatistics*len(VectorPiFrac), 103))

  # Creation of all the input parameters that goes to MPI's
def cartesian(*arrays):
      mesh = np.meshgrid(*arrays)  # standard numpy meshgrid
      dim = len(mesh)  # number of dimensions
      elements = mesh[0].size  # number of elements, any index will do
      flat = np.concatenate(mesh).ravel()  # flatten the whole meshgrid
      reshape = np.reshape(flat, (dim, elements)).T  # reshape and transpose
      return reshape



A = np.arange(NumStatistics)
b = cartesian(VectorN, VectorP, A)
a = cartesian(VectorN, VectorP, VectorPiFrac, A) # The matriz that gives the number of nodes to used. Each element goes in a MPI process
b = np.asarray(b)
print('b',b)
contador = 0
fileString = "Northridge-"+str(VectorN[0])+"Test-FBMfaults"
for i in enumerate(b):
  print(i[1])
  print(contador)
  Nbox=int(i[1][0])
  print(Nbox)
  P = i[1][1]
  print(P)
  VecPosi = InitialLoadConfiguration(Nbox,P) # Here we create the inital load matrix, of a N size and P order particular value.
  MatFaults = np.loadtxt('Matriz'+FaultConfigurations+'Faults'+str(int(i[1][0]))+'.dat') # esta matriz ha de entrar al MPI
  np.savetxt(fileString+'VecPosi.dat',VecPosi)
  
  print('Sali')
  smin =  3*pow(i[1][0],2)/4
  print(smin)
  smin = int(smin)
  print('smin(Int)',smin)
  for j in enumerate(VectorPiFrac):
   
    print(j)
    
    Results = FBMaftershocsMAIN(MatFaults,i[1][0],i[1][1],j[1],VectorPiBG[0],rho,smin,VecPosi,contador)  
    ResultsTot[contador,:] = Results
    contador += 1
    

np.savetxt('CI_'+fileString+'.csv', ResultsTot, delimiter=',')
t.toc('End Program No Optimizados') #Time elapsed since t.tic()
