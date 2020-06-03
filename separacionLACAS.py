# coding=utf-8
import numpy as np
def separacionLACAS(nsismo,tiempoacu,tiempo,AVANOR,kmax):
  #print(vectk1)
  #nsismo = vectk1[:,0];
  #tiempoacu = vectk1[:,1];
  #tiempo = vectk1[:,2]; 
  #AVANOR = vectk1[:,3];
  nserie = len(nsismo);
    #!! las dos primeras replicas son leadings
  ncas = -1;
  num = 0;
  nlead = 0;
  nflag = 0;
  k = -1;
  npor = 0;
  countCas = 0
  tiempox = np.zeros(nserie)
  nsismolead = np.zeros(nserie)
  tiempoacul = np.zeros(nserie)
  vectLead = np.zeros([nserie,6]) 
  vecCas = np.zeros([nserie,5])
  nsismoca = np.zeros(nserie)
  tiempocas = np.zeros(nserie);
  tiempoacuca = np.zeros(nserie);
  vectorProbab = np.zeros(nserie)
  vecNumCas = np.zeros([nserie,2])
  vectorNumCasca = np.zeros(nserie)

  i = 0;
  deltatempl = 0.0;
  tiempox[nlead] = 0.0;	 # ESTE VECTOR GUARDA LA DIFERENCIA ENTRE DOS LEADING AFTERSHOCKS SUCESIVOS
  nsismolead[nlead] = nlead;	 # ESTE VECTOR GUARDA EL NUMERO SECUENCIAL DE LEADINGS
  tiempoacul[nlead] = tiempoacu[i]  # ESTE VECTOR GUARDA EL TIEMPO ACUMULADO EN EL QUE OCURRE CADA LEADING
  vectLead[nlead,0] = nsismo[i];
  vectLead[nlead,1] = nsismolead[nlead];
  vectLead[nlead,2] = tiempox[nlead];
  vectLead[nlead,3] = deltatempl;
  vectLead[nlead,4] = tiempoacul[nlead];
  vectLead[nlead,5] = AVANOR[i];
  nlead = nlead+1;
  i = 1;
  tiempol = tiempoacu[i];
  deltatempl = tiempol - tiempoacu[i-1];
  tiempox[nlead] = deltatempl;
  nsismolead[nlead] = nlead;
  tiempoacul[nlead] = tiempoacu[i];
  vectLead[nlead,0] = nsismo[i];
  vectLead[nlead,1] = nsismolead[nlead];
  vectLead[nlead,2] = tiempox[nlead];
  vectLead[nlead,3] = deltatempl;
  vectLead[nlead,4] = tiempoacul[nlead];
  vectLead[nlead,5] = AVANOR[i];

  for ii in np.arange(2,nserie,dtype=np.int):
    if tiempoacu[ii]-tiempoacu[ii-1] > (tiempoacu[ii-1]-tiempoacu[ii-2]) and (tiempoacu[ii]-tiempol) > deltatempl:    
      if nflag == 1:
        nflag = 0;
        num = 0;
        ncas = ncas+1;
        nsismoca[ncas] = num;
        tiempoacuca[ncas] = tiempoacul[nlead];
      
      nlead = nlead+1;
      deltatempl = tiempoacu[ii] - tiempol;
      tiempol = tiempoacu[ii];
      tiempox[nlead] = deltatempl+tiempox[nlead-1];
      nsismolead[nlead] = nlead;
      tiempoacul[nlead] = tiempoacu[ii];      
      vectLead[nlead,0] = nsismo[ii];
      vectLead[nlead,1] = nsismolead[nlead];
      vectLead[nlead,2] = tiempox[nlead];
      vectLead[nlead,3] = deltatempl;
      vectLead[nlead,4] = tiempoacul[nlead];
      vectLead[nlead,5] = AVANOR[ii];          
      k=k+1;
      vecCas[k,0] = nlead;
      vecCas[k,1] = 0;
      vecCas[k,2] = tiempo[ii];
      vecCas[k,3] = tiempoacul[nlead];
      vecCas[k,4] = AVANOR[ii];
      
    else:    
      nflag = 1;
      num = num+1;
      ncas = ncas+1;
      nsismoca[ncas] = num;
      tiempocas[ncas] = tiempo[ii];
      tiempoacuca[ncas] = tiempoacu[ii];
    
      k=k+1;
      vecCas[k,0] = num;
      vecCas[k,1] = nsismoca[ncas];
      vecCas[k,2] = tiempocas[ncas];
      vecCas[k,3] = tiempoacuca[ncas];
      vecCas[k,4] = AVANOR[ii];      
     
  return vectLead[1:nlead,:]  
