# coding=utf-8
import numpy as np
import scipy
from scipy.optimize import curve_fit

def OmoriFit(time_as, numEventos): 

 
  def Omori(x,p,c,k):
   
    B = pow(c,(1-p))
    C = pow((c+x),(1-p)) 
    A = B - C
    param = A/(p-1)
   # cumnr_model = p[2]*param 
   #i = np.arange(0,len(x))
    return k*param    #  k*(pow((fTend+c),(1-p))-pow((fTstart+c),(1-p)))/(1-p); 
   
  p0 = [1.1, 0.1, 20]
  popt1, pcov1 = curve_fit(Omori, time_as, numEventos, p0, bounds = ((0.2, 0.001, 1), (2.9, 30,5000)), max_nfev= 1000)
 # popt1 = minimize(Omori, p0, method='SLSQP', bounds = ((0.2, 2.9), (0.0001,30), (1,4000)))
#  plt.plot(time_as, numEventos, 'ob-', label='data')
#  plt.plot(time_as, Omori(time_as, *popt1), 'g--',label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt1))
#  plt.plot(time_as, result.init_fit, 'k--')
#  plt.plot(time_as, result.best_fit, 'r-')
#  plt.legend()
#  plt.show()
  #pause()
    
  perr1 = np.sqrt(np.diag(pcov1))
  pv = popt1[0]
  cv = popt1[1]
  kv = popt1[2]
  minRes_p = perr1[0]
  #minRes_c = 
  #minRes_k = 
  return pv,cv,kv,minRes_p

