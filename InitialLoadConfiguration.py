# coding=utf-8
import numpy as np
# This function createas the initial load configuration that depends on N y and P
def InitialLoadConfiguration(Nbox,P):
  np.random.seed(1)    
    #The initial load matrix will be ordered in respect to the P value    
# 1. First we ordered for P=1
 ## Nbox = int(Vector[0][0])
 ## P = Vector[0][1]
  OriginalVec = np.random.rand(pow(Nbox,2)) # matrix of random numbers P = 0 
  numCentral= int((Nbox/3)*(Nbox/3))   # number of boxes which is the matrix dividing 
  VecCentr = OriginalVec[0:numCentral]
  VecExtern = OriginalVec[numCentral:len(OriginalVec)]
  VecNEW = sorted(VecExtern, reverse = True)
  sizeCuad = int(len(OriginalVec)/9)
  Part = int(len(VecExtern)/2)
  VecOrdCouPosi = VecNEW[0:Part]
  VecOrdCouNega = VecNEW[Part:len(VecNEW)]
  Vec1Posi = VecCentr 
  sizebox = int(Nbox/3)
  
  Matriz3 = distInicialOrdenadaV0Comp(Vec1Posi,sizebox)  
 # np.savetxt('CentralMatriz.dat',Matriz3) 
  
  Vec2Posi = VecOrdCouPosi[0:len(VecOrdCouPosi):4]
  Vec3Posi = VecOrdCouPosi[1:len(VecOrdCouPosi):4]
  Vec4Posi = VecOrdCouPosi[2:len(VecOrdCouPosi):4]
  Vec5Posi = VecOrdCouPosi[3:len(VecOrdCouPosi):4]
  Vec6Posi = VecOrdCouNega[0:len(VecOrdCouNega):4]
  Vec7Posi = VecOrdCouNega[1:len(VecOrdCouNega):4]
  Vec8Posi = VecOrdCouNega[2:len(VecOrdCouNega):4]
  Vec9Posi = VecOrdCouNega[3:len(VecOrdCouNega):4]
  
  sizeNeg = int(pow(len(Vec6Posi),0.5))
  sizePos = int(pow(len(Vec2Posi),0.5))

  #Vec6Posi=sorted(Vec6Posi)
  #Vec7Posi=sorted(Vec7Posi)
  #Vec8Posi=sorted(Vec8Posi)
  #Vec9Posi=sorted(Vec9Posi)

  Matriz6=matcoulUniforV5(Vec6Posi,sizeNeg,6)
  Matriz7=matcoulUniforV5(Vec7Posi,sizeNeg,7)
  Matriz8=matcoulUniforV5(Vec8Posi,sizeNeg,8)
  Matriz9=matcoulUniforV5(Vec9Posi,sizeNeg,9)
 
  Matriz1=matcoulUniforV5(Vec2Posi,sizePos,2)
  Matriz2=matcoulUniforV5(Vec3Posi,sizePos,3)
  Matriz4=matcoulUniforV5(Vec4Posi,sizePos,4)
  Matriz5=matcoulUniforV5(Vec5Posi,sizePos,5)

  L1= np.concatenate((Matriz6[0:len(Matriz6),0:len(Matriz6)], Matriz1[0:len(Matriz1),0:len(Matriz1)], Matriz9[0:len(Matriz9),0:len(Matriz9)]),axis = 1)
  L2= np.concatenate((Matriz4[0:len(Matriz4),0:len(Matriz4)], Matriz3[0:len(Matriz3),0:len(Matriz3)], Matriz5[0:len(Matriz5),0:len(Matriz5)]),axis = 1)
  L3= np.concatenate((Matriz8[0:len(Matriz8),0:len(Matriz8)], Matriz2[0:len(Matriz2),0:len(Matriz2)], Matriz7[0:len(Matriz7),0:len(Matriz7)]),axis = 1)
  MAT = np.concatenate((L1, L2, L3),axis= 0)   # MAT contains the initial loads of each cell following a 100% order

 #PyPlot.figure(79)
	 #PyPlot.surf(C,rstride=1, cstride=1, cmap="seismic")
	 #PyPlot.title("CoulombUniformNbox90 P=100%")
	#  PyPlot.savefig("CoulombUniformNbox90 P=100%.eps",dpi=800)
  VecPosi = disorderAlgorithm(P,MAT,Nbox)
  
 # file12 = ('InitialLoadCon-P'+str(P)+'-'+str(Nbox)+'.dat')
  #  file12 = ("VecAsig100orderNbox"*string(Nbox)*".dat")
 # np.savetxt(file12,VecPosi)
  return VecPosi

  
# 2. We will disorder the matrix MAT according to The P value

def disorderAlgorithm(P,YY,Nboxx):
  
  contador = 0
  contadori = -1
  contadorj = -1
  alfa = 0
  Y = np.copy(YY)
 
  VecPosi = np.zeros((Nboxx,Nboxx))
  VecMAP = np.zeros((pow(Nboxx,2),3))
  VecResha = Y.reshape(pow(len(Y[:,1]),2))

 # VecSort=sort(VecResha,rev=true)
  POrd = int(len(VecResha)*P)
  PDes = int(len(VecResha*(1.0-P)))  
 
  VecNflagOrd = np.zeros(len(VecResha))

  if PDes==0:
    VecNflag = np.zeros(1)
  else:
    VecNflag = np.zeros(len(VecResha))
    
  contikAlfa = -1
  contikk = -1

  for i in enumerate(VecResha):
  
    alfa = np.random.rand()
    
    if 0 < alfa <= P:      
      contikAlfa += 1
      VecNflagOrd[contikAlfa] = i[1]
      
    elif P < alfa <= 1:
      contikk += 1
      VecNflag[contikk] = i[1]

  
  if contikk == 0:
    VecDesor = VecNflag[len(VecNflag)]
  else:
    VecDesor = np.random.permutation(VecNflag[0:contikk])

  
  #Aqui se acomodan los elementos con Probabilidad P
  for i in enumerate(Y[:,1]):
    contadori += 1
    for j in enumerate(Y[:,1]):
      contadorj += 1
      for k in VecNflagOrd[0:contikAlfa]:
        if Y[contadori,contadorj] == k:
          VecPosi[contadori,contadorj] = Y[contadori,contadorj]
          break    
    contadorj = -1
  
  if P == 1.0:
    VecPosi= VecPosi
  else:
    contadori = -1
    contadorj = -1
    for i in enumerate(Y[:,1]):
      contadori += 1
      for j in enumerate(Y[:,1]):
        contadorj += 1
        if len(VecDesor) == 0:
          break
        if VecPosi[contadori,contadorj] == 0:
          VecPosi[contadori,contadorj] = VecDesor[0]
          VecDesor = np.delete(VecDesor,0)
        else:
          continue
      contadorj = -1

  return VecPosi

 
def matcoulUniforV5(Vec, Nbox, nflag):
  
  A = int(Nbox/2)-1
  B = A
  X = np.copy(Vec)
  Matriz = np.zeros((Nbox,Nbox))
  numVec = np.zeros(pow(Nbox,2))
  #Matriz[A,B] = X[0]
  #X = np.delete(X,0)
  numPas = 0
  k = 0
  ranginf = 0
  if nflag == 6 or nflag == 7 or nflag == 8 or nflag == 9:
    Xrev = np.copy(X)
  else:
    Xrev = np.sort(X)

  for i in np.arange(A+1):
    Matriz[i,i:Nbox-i] = Xrev[0:len(Matriz[i,i:Nbox-i])]
    Xrev = np.delete(Xrev,np.arange(len(Xrev[0:len(Matriz[i,i:Nbox-i])])))
    Matriz[Nbox-i-1,i:Nbox-i] = Xrev[0:len(Matriz[Nbox-i-1,i:Nbox-i])]
    Xrev = np.delete(Xrev,np.arange(len(Xrev[0:len(Matriz[Nbox-i-1,i:Nbox-i])])))
    Matriz[i+1:Nbox-i-1,i] = Xrev[0:len(Matriz[i+1:Nbox-i-1,i])]
    Xrev = np.delete(Xrev,np.arange(len(Xrev[0:len(Matriz[i+1:Nbox-i-1,i])])))
    Matriz[i+1:Nbox-i-1,Nbox-1-i] = Xrev[0:len(Matriz[i+1:Nbox-i-1,Nbox-1-i])]
    Xrev = np.delete(Xrev,np.arange(len(Xrev[0:len(Matriz[i+1:Nbox-i-1,Nbox-1-i])])))

  return Matriz

def distInicialOrdenadaV0Comp(X,Nboxx):

  VecPosi = np.zeros((Nboxx,Nboxx))
  VecMAP = np.zeros((pow(Nboxx,2),3))
  Vec = np.zeros(pow(Nboxx,2))

  for cont in enumerate(Vec):
    Vec[cont[0]]=X[cont[0]] 

  fracc=1.0     ##porcentaje de Ordenación con respecto a la malla central. Dejarla en este valor 
  
  VecNEW = np.sort(Vec)
  k = -1
  contador = -1
  VecMax = VecNEW[len(VecNEW)-1] # el elemento mas grande
  VecNEW = np.delete(VecNEW,len(VecNEW)-1)
  numAle = 1
  p = np.random.rand()
   
  VecMax1 = VecNEW[len(VecNEW)-1]  # en la grafica ORIENTACION NS    
  VecNEW = np.delete(VecNEW,len(VecNEW)-1) # borra el elemento que se acaba de asignar

  for j in np.arange(int(Nboxx/2)):
    for i in np.arange(int(Nboxx)):
      if len(VecNEW) == 0:
          break
      
      if j == int(Nboxx/2)-1 and i == int(Nboxx/2)-1:
        
        VecPosi[i,j] = VecMax
        VecPosi[i,(Nboxx-1)-j]=VecMax1  
        
      else:
        VecPosi[i,1+(j-1)] = VecNEW[0]  # en la grafica ORIENTACION NS
        VecNEW = np.delete(VecNEW,0)

        VecPosi[i,(Nboxx-1)-j]=VecNEW[0]
        VecNEW = np.delete(VecNEW,0)
                 

 #PyPlot.figure(718)
 #PyPlot.surf(VecPosi,rstride=1, cstride=1, cmap="seismic")
  return VecPosi




  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  




    
    
    






