# Con este script se podrán determinar las características más destacables de las series de réplicas b-valor, exponente de Husrt, Dimensión fractal.
# Se leen los datos partiendo de P y se coloca en el eje abcisas a fhiFracture
#include("Hurst.jl")

import numpy as np
import scipy
from scipy import optimize
from scipy.optimize import minimize


def ParametrosResultantes(XX,file):   
 
#  XX[l,1] = VecDatos[ij,8]  #x
#  XX[l,2] = VecDatos[ij,9]  #y
#  XX[l,3] = VecDatos[ij,2]  #tiempo 
#  XX[l,5] = VecDatos[ij,1]  #number position in original database
#  XX[l,4] = NUM+1         #numero evenctos avalanchas
    print(len(XX[:,0]))

    VectorDistancia = np.zeros(len(XX[:]))
    VectorTiempoInter = np.zeros(len(XX[:]))
    VectorTiempo = np.zeros(len(XX[:]))
    VecAreaAval = np.zeros(len(XX[:]))
    MagniResults = np.zeros(3)
    Results = 0.0
    contador = -1
    NUM = 0

# Se crea los vectores de distancia y tiempo entre sismos consecutivos

    for jk,_ in enumerate(XX[:,0]):
        contador += 1
        VectorDistancia[jk-1] = np.sqrt(np.power((XX[jk,0]-XX[jk-1,0]),2) + np.power((XX[jk,1]-XX[jk-1,1]),2))
        VectorTiempoInter[jk-1] = (XX[jk,2]) - (XX[jk-1,2])
        VectorTiempo[jk-1] = XX[jk,2]
        VecAreaAval[jk-1] = XX[jk,3]


    paramDistan = escalat(VectorDistancia[0:contador]) # paramDistan = [a0,b0,db0,rho02,res0]
    paramInterevent= escalat(VectorTiempoInter[0:contador])
    paramMagni= escalat(VecAreaAval[0:contador])

    if len(VecAreaAval[0:contador]) != 0:
        Results = np.max(VecAreaAval[0:contador])
    else:
        Results = 0

    VecResultados = np.zeros(16)
    VecResultados[0:5] = paramDistan
    VecResultados[5:10] = paramInterevent
    VecResultados[10:15] = paramMagni
    VecResultados[15] = Results

    return VecResultados


def escalat(data):

    nanys = len(data)  # 261
    nlim = 30000  # 18000
    nmaxr = 200  # 500

    ntope = nanys
    VecPas = np.arange(2, 11, 1, dtype=np.int)

    if nanys < 100:
        VecPas = np.append(VecPas,  np.arange(10, nanys, 10, dtype=np.int))
    elif nanys > 100:
        VecPas = np.append(VecPas, np.arange(10, 100, 10,dtype=np.int))
        VecPas = np.append(VecPas, np.arange(100, nanys, 100, dtype=np.int))
        
        
    VecSerie = np.zeros((2, len(VecPas)))

    if nlim > len(data):
      nlim = len(data)


    for ik,_ in enumerate(VecPas):
      it = VecPas[ik]
      if it > nanys:
        ikk = ik - 1

      xrs = HurstFor(nmaxr, it, nlim, data)
      VecSerie[0, ik] = np.log10(it)
      VecSerie[1, ik] = np.log10(xrs)

    ikk = len(VecPas)

    a0, b0, db0, rho02, res0 = regres(VecSerie, ikk)

    VecResults = np.zeros(5)
    VecResults[0] = a0
    VecResults[1] = b0
    VecResults[2] = db0
    VecResults[3] = rho02
    VecResults[4] = res0
    
    return VecResults

def HurstFor(nmaxr, it, nlim, data):

  x = data

  xd = np.zeros(30000)
  xrs = 0
  nmaxrf = 0
  Vecnmaxr = np.arange(1,nmaxr,1)


  for kk,_ in enumerate(Vecnmaxr):

    it1 = it * kk

    it2 = it * (kk+1)

    if it2 > nlim or it1 > nlim:
        break

    nmaxrf += 1

    total = 0.0
    total = np.sum(x[it1:it2+1])
    xm = total/it

    for ii in range(it1,it2+1):
      if it1 > len(x) or  it2 > len(x) or ii > len(x):
        break

      total = 0
      for iii in range(it1,ii):
        total = total + (x[iii] - xm)

    xd[ii] = total
    dmax = -1.0e+12
    dmin = 1.0e+12

    for ii in range(it1,it2+1):
      if xd[ii] >= dmax:
        dmax = xd[ii]
      if xd[ii] < dmin:
        dmin = xd[ii]

    xr = dmax - dmin
    total = 0
    xs = 0

    for ii in range(it1,it2):
      if it1 > len(x) or  it2 > len(x) or ii > len(x):
        break
      total = total + np.power((x[ii]-xm), 2)
      xs = np.sqrt(total/it)

    #if xr == 0 or xs == 0:
    #  nmaxrf = nmaxrf - 1
    #  break
    xrs = xrs + xr/xs

  xrs = xrs/nmaxrf
  return xrs

def regres(serie, npunts):

  xm = 0.
  ym = 0.
  xym = 0.
  sx = 0.
  sy = 0.

  for j in range(0,npunts):
    xm += serie[0, j]
    ym += serie[1, j]

  xm = xm / npunts
  ym = ym / npunts

  for j in range(1,npunts):
    xym += (serie[0, j] - xm) * (serie[1, j] - ym)
    sx += np.power((serie[0, j] - xm), 2)
    sy += np.power((serie[1, j] - ym), 2)

  xym = xym / npunts
  sx = sx / npunts
  sy = sy / npunts
  b0 = xym / sx
  a0 = ym - b0 * xm
  rho0 = b0 * np.sqrt(sx) / np.sqrt(sy)
  rho02 = np.power(rho0, 2)
  res0 = 0.

  for j in range(0,npunts):
    res0 += np.power((serie[1, j] - (a0 + b0 * serie[0, j])) , 2)
  

  res0 = res0 / (npunts - 2)
  db0 = res0 / (npunts * sx * sx)

  return a0, b0, db0, rho02, res0











