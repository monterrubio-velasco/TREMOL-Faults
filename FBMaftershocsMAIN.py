# coding=utf-8
import numpy as np
from numpy.random import rand
from FBMV12FracturasOPT import FBMV12FracturasOPT
from CorreDimensionAvalanchesV11 import CorreDimensionAvalanchesV11
from calcuMagniInKM import calcuMagniInKM
from ParametrosResultantes import ParametrosResultantes
# # include("CreacionFractReal.jl")
from OmoriFit import OmoriFit
from separacionLACAS import separacionLACAS
#import SpatialCoordenadas

def FBMaftershocsMAIN(MatFracturasCompleta,Nboxi,P,fhiFractura,fhi,rho,smin,VecPosi,contExperimentos):

    Nbox = len(VecPosi[:,1])
    n_celdas = np.power(Nbox,2)
    MatFracturas = MatFracturasCompleta #[0:Nboxi-3,0:Nboxi-3]
    #     MatFracturas = MatFracturasCompleta[0:Nboxi-2,0:Nboxi-2]
    za=str(fhi)
    xa=str(rho)
    ResulEstad = np.zeros(103)
    fileString= 'NR-'+str(Nboxi)+"-"+str(P)+'-'+str(fhiFractura)+'-'+str(fhi)
    print('Test-ENTRE')
    VecDatos,kfin =  FBMV12FracturasOPT(MatFracturas,VecPosi,Nbox,smin,fhi,fhiFractura,fileString,rho,contExperimentos)
    
    np.savetxt('VecDatos'+fileString+'.dat',VecDatos)

    #        !!!!!!!!!!!!!!!Análisis Post-proceso !!!!!!!!!!!!!!!!
    #        se genera una serie que contenga la serie filtrada para las avalanchas de un solo elemento independiente del tamaño de la red

    
    XX = np.zeros((len(VecDatos)-1,5))
    print(len(VecDatos))
    print(len(VecDatos[:,1]))
    l = -1
    NUM = 0

    for ij,_ in enumerate(VecDatos):
        if VecDatos[ij,3] == 1 and VecDatos[ij-1,3] == 0:
            NUM=0
            XX[l,0]=VecDatos[ij,7]  #x
            XX[l,1]=VecDatos[ij,8]  #y
            XX[l,2]=VecDatos[ij,1]  #tiempo
            XX[l,4]=VecDatos[ij,0]  #number position in original database
            l += 1
        elif VecDatos[ij,3] == 1 and VecDatos[ij-1,3] == 1:
            NUM=NUM+1
        elif VecDatos[ij,3] == 0 and VecDatos[ij-1,3] == 1:
            XX[l,3] = NUM+1

    print('XX',XX[0,:],'l',l,len(VecDatos))
    # se genera una serie que contenga la serie filtrada para un valor minimo de magnitud definido
    # Landers earthquake (50 km x 50 km) = area 2500 km²
    # Darfield earthquake test0 (200 km x 200 km) = area 40000 km²

    numCeldasRef = 177^2
    AreaCeldaRef = 4300/numCeldasRef
    numCeldas = Nbox^2
    AreaCelda = 4300/numCeldas

    VecArea = AreaCelda * XX[1:l,4]
    Mag1cell = round(4/3*np.log10(AreaCelda) + 3.07,1)
    Mag2cell = round(4/3*np.log10(AreaCelda*2) + 3.07,1)
    XXmin = np.zeros((l,5))
    ll = -1
    for i,_ in enumerate(XXmin):
        if XX[i,3] >= 2:
            ll += 1
            XXmin[ll,:] = XX[i,:]
    XXRef = np.zeros((l,5))
    llRef = -1
    ValRef = int(AreaCeldaRef/AreaCelda+1)

    for i ,_ in enumerate(XXmin):
     if XX[i,3] >= ValRef:
         llRef += 1
         XXRef[llRef,:] = XX[i,:]

     Mag1cellRef = round(4/3*np.log10(AreaCeldaRef) + 3.07,1)
     XXRef2 = np.zeros((l,5))
     llRef2 = -1
     for i, _ in enumerate(XXmin):
        if XX[i,3] >= ValRef*2:
            llRef2 += 1
            XXRef2[llRef2,:] = XX[i,:]

    Mag2cellRef = round(4/3*np.log10(AreaCeldaRef*2) + 3.07,1)
    #
    # #     # función para graficar el espacio
    # #       SpatialCoordenadas(VecDatos,MatFracturas,Nbox,AreaCelda,fileString)
    #
         # función que grafica las dimensiones fractales

    fract0 = MatFracturas[MatFracturas==1]
    NumCeldasFrac = len(fract0)
    fract1 = MatFracturas[MatFracturas==0]
    NumCeldasNonFrac = len(fract1)

    ResulEstad[0] = P
    ResulEstad[1] = fhiFractura
    ResulEstad[2] = Nbox
    ResulEstad[3] = NumCeldasFrac
    ResulEstad[4] = NumCeldasNonFrac

    D0magRef2, D1magRef2, D2magRef2 = CorreDimensionAvalanchesV11(XXRef2[0:llRef2,:])
    
    ResulEstad[5] = D0magRef2[0]
    ResulEstad[6] = D0magRef2[1]
    ResulEstad[7] = D0magRef2[2]
    ResulEstad[8] = D0magRef2[3]
    ResulEstad[9] = D1magRef2[0]
    ResulEstad[10] = D1magRef2[1]
    ResulEstad[11] = D1magRef2[2]
    ResulEstad[12] = D1magRef2[3]
    ResulEstad[13] = D2magRef2[0]
    ResulEstad[14] = D2magRef2[1] 
    ResulEstad[15] = D2magRef2[2]
    ResulEstad[16] = D2magRef2[3]
    
    meanmHB1Ref2, bHB1Ref2, sigHB1Ref2, av2HB1Ref2,paramGRHB11Ref2,paramGRHB12Ref2,paramGRHB13Ref2,paramGRHB14Ref2,paramGRHB15Ref2,aVLSSilvaRef2, qVLSSilvaRef2, minVLSSilvaRef2, aVLSTelescaRef2, qVLSTelescaRef2, minVLSTelescaRef2,numSerieMagRef2,cuenrepHB1Ref2 =calcuMagniInKM(XXRef2[1:llRef2-1,3],AreaCelda)
    
    
    ResulEstad[17] = meanmHB1Ref2  # FH
    ResulEstad[18] = bHB1Ref2
    ResulEstad[19] = sigHB1Ref2
    ResulEstad[20] = av2HB1Ref2
    ResulEstad[21] = cuenrepHB1Ref2
    ResulEstad[22] = 0
    ResulEstad[23] = numSerieMagRef2
    ResulEstad[24] = paramGRHB11Ref2   #FP
    ResulEstad[25] = paramGRHB12Ref2
    ResulEstad[26] = paramGRHB13Ref2
    ResulEstad[27] = paramGRHB14Ref2
    ResulEstad[28] = paramGRHB15Ref2               
    ResulEstad[29] = aVLSSilvaRef2     #GL               
    ResulEstad[30] = qVLSSilvaRef2
    ResulEstad[31] = minVLSSilvaRef2
    ResulEstad[32] = aVLSTelescaRef2  #GO
    ResulEstad[33] = qVLSTelescaRef2
    ResulEstad[34] = minVLSTelescaRef2
    ResulEstad[35] = 0
    
    # función en la que se calcula el parámetro de Hurst y la distribucion estadística del vector de distancias

    Vec_Hursr_DistProb_Ref2 = ParametrosResultantes(XXRef2[0:llRef2,:],fileString)
            #         paramDistan = [a0,b0,db0,rho02,res0]
            #  VecResultados = zeros(26)
            #     Vec_Hursr_DistProb [1:5] = paramDistan
            #    Vec_Hursr_DistProb [6:10] = paramInterevent
            #     Vec_Hursr_DistProb [11:15] = paramMagni
            #     Vec_Hursr_DistProb [16] = Results
            #     Vec_Hursr_DistProb [17] = ncont
            #     Vec_Hursr_DistProb [18] = rl1
            #     Vec_Hursr_DistProb [19] = rl2
            #     Vec_Hursr_DistProb [20] = rl3
            #     Vec_Hursr_DistProb [21] = rl4
            #     Vec_Hursr_DistProb [22] = t3
            #     Vec_Hursr_DistProb [23] = t4
            #     Vec_Hursr_DistProb [24] = q1_CDF_fit
            #     Vec_Hursr_DistProb [25] = q2_CDF_fit
            #     Vec_Hursr_DistProb [26] = rhocoefcorr_CDF_fit
   
    ResulEstad[36] = Vec_Hursr_DistProb_Ref2[0]
    ResulEstad[37] = Vec_Hursr_DistProb_Ref2[1] # JS
    ResulEstad[38] = Vec_Hursr_DistProb_Ref2[2]
    ResulEstad[39] = Vec_Hursr_DistProb_Ref2[3]
    ResulEstad[40] = Vec_Hursr_DistProb_Ref2[4]
    ResulEstad[41] = Vec_Hursr_DistProb_Ref2[5]
    ResulEstad[42] = Vec_Hursr_DistProb_Ref2[6] #JX
    ResulEstad[43] = Vec_Hursr_DistProb_Ref2[7]
    ResulEstad[44] = Vec_Hursr_DistProb_Ref2[8]
    ResulEstad[45] = Vec_Hursr_DistProb_Ref2[9]
    ResulEstad[46] = Vec_Hursr_DistProb_Ref2[10]
    ResulEstad[47] = Vec_Hursr_DistProb_Ref2[11]
    ResulEstad[48] = Vec_Hursr_DistProb_Ref2[12]
    ResulEstad[49] = Vec_Hursr_DistProb_Ref2[13]
    ResulEstad[50] = Vec_Hursr_DistProb_Ref2[14]
    ResulEstad[51] = Vec_Hursr_DistProb_Ref2[15] #KG
    #ResulEstad[52] = Vec_Hursr_DistProb_Ref2[16]
    #ResulEstad[53] = Vec_Hursr_DistProb_Ref2[17]
    #ResulEstad[54] = Vec_Hursr_DistProb_Ref2[18]
    #ResulEstad[55] = Vec_Hursr_DistProb_Ref2[19]
    #ResulEstad[56] = Vec_Hursr_DistProb_Ref2[20]
    #ResulEstad[57] = Vec_Hursr_DistProb_Ref2[21]
    #ResulEstad[58] = Vec_Hursr_DistProb_Ref2[22]
    #ResulEstad[59] = Vec_Hursr_DistProb_Ref2[23]
    #ResulEstad[60] = Vec_Hursr_DistProb_Ref2[24] #KP
    #ResulEstad[61] = Vec_Hursr_DistProb_Ref2[25]
   
    pv = np.zeros(6)
    cv = np.zeros(6)
    kv = np.zeros(6)
    minerror = np.zeros(6)
    numOmoriData = np.zeros(6)
    
    interTime = np.zeros(llRef)
    interTime[0] = 0.0
    interTime[1:llRef] = XXRef[1:llRef,2]-XXRef[0:llRef-1,2]
    VecLA = separacionLACAS(np.arange(0, llRef,dtype = np.int), XXRef[0:llRef,2], interTime, np.ones(llRef),llRef)
      # np.savetxt("LApython.dat",VecLA )
        
        #numFit = llRef2 - (llRef2*(0.05*np.power(2,(0)))).astype(np.int)
        #numOmoriData[0] = (llRef2*(0.05*np.power(2,(0)))).astype(np.int)     
        #if llRef2 - numFit > 10:
            #pv[0],cv[0],kv[0],minerror[0] = OmoriFit(XXRef2[numFit:llRef2,2],np.arange(llRef2-numFit))
            #ResulEstad[97] = len(np.arange(llRef2-numFit)) 
        ##    print( pv[0],cv[0],kv[0],minerror[0])
        
        #numFit = llRef2 - (llRef2*(0.05*np.power(2,(1)))).astype(np.int)
        #numOmoriData[1] = (llRef2*(0.05*np.power(2,(1)))).astype(np.int)    
        #if llRef2 - numFit > 10:
            #pv[1],cv[1],kv[1],minerror[1] = OmoriFit(XXRef2[numFit:llRef2,2],np.arange(llRef2-numFit))
            #ResulEstad[98] = len(np.arange(llRef2-numFit))
        ##    print( pv[1],cv[1],kv[1],minerror[1])

        #numFit = llRef2 - (llRef2*(0.05*np.power(2,(2)))).astype(np.int)
        #numOmoriData[2] = (llRef2*(0.05*np.power(2,(2)))).astype(np.int)    
        #if llRef2 - numFit > 10:
            #pv[2],cv[2],kv[2],minerror[2] = OmoriFit(XXRef2[numFit:llRef2,2],np.arange(llRef2-numFit))
            #ResulEstad[99] = len(np.arange(llRef2-numFit))
        ##    print( pv[2],cv[2],kv[2],minerror[2])

        #numFit = llRef2 - (llRef2*(0.05*np.power(2,(3)))).astype(np.int)
        #numOmoriData[3] = (llRef2*(0.05*np.power(2,(3)))).astype(np.int)    
        #if llRef2 - numFit > 10:
            #pv[3],cv[3],kv[3],minerror[3] = OmoriFit(XXRef2[numFit:llRef2,2],np.arange(llRef2-numFit))
            #ResulEstad[100] = len(np.arange(llRef2-numFit))
        ##    print(pv[3],cv[3],kv[3],minerror[3])

        #numFit = llRef2 - (llRef2*(0.05*np.power(2,(4)))).astype(np.int)
        #numOmoriData[4] = (llRef2*(0.05*np.power(2,(4)))).astype(np.int)    
        #if llRef2 - numFit > 10:
            #pv[4],cv[4],kv[4],minerror[4] = OmoriFit(XXRef2[numFit:llRef2,2],np.arange(llRef2-numFit)) 
            #ResulEstad[101] = len(np.arange(llRef2-numFit))
        ##    print(pv[4],cv[4],kv[4],minerror[4])

        
        #YY = [np.arange(0, llRef2,dtype = np.int), XXRef2[0:llRef2,2], interTime, np.ones(llRef2)]
        
    numOmoriData[5] = len(VecLA[:,4])
    if  numOmoriData[5] > 10:
        print('numOmoriData[5]')
        print(numOmoriData[5])
        pv[5],cv[5],kv[5],minerror[5] = OmoriFit(VecLA[:,4],np.arange(len(VecLA[:,4])))  
        print(pv[5],cv[5],kv[5],minerror[5])
    else:
        pv[5] = 0
        cv[5] = 0
        kv[5] = 0
        minerror[5] = 0
        
    ResulEstad[62] = pv[0] #BK
    ResulEstad[63] = pv[1] 
    ResulEstad[64] = pv[2]
    ResulEstad[65] = pv[3]
    ResulEstad[66] = pv[4]
    ResulEstad[67] = pv[5]
    ResulEstad[68] = cv[0] #BQ
    ResulEstad[69] = cv[1]
    ResulEstad[70] = cv[2]
    ResulEstad[71] = cv[3]
    ResulEstad[72] = cv[4]
    ResulEstad[73] = cv[5]
    ResulEstad[74] = kv[0]  # BW
    ResulEstad[75] = kv[1]
    ResulEstad[76] = kv[2]
    ResulEstad[77] = kv[3]
    ResulEstad[78] = kv[4]
    ResulEstad[79] = kv[5]
    ResulEstad[80] = minerror[0]
    ResulEstad[81] = minerror[1]
    ResulEstad[82] = minerror[2]
    ResulEstad[83] = minerror[3]
    ResulEstad[84] = minerror[4]
    ResulEstad[85] = minerror[5]
    ResulEstad[86] = numSerieMagRef2
    ResulEstad[87] = VecDatos[-1,4]
    ResulEstad[88] = VecDatos[-1,5] #LT
    ResulEstad[89] = AreaCelda
    ResulEstad[90] = AreaCeldaRef
    ResulEstad[91] = ValRef
    ResulEstad[92] = smin #LX
    ResulEstad[93] = l
    ResulEstad[94] = ll
    ResulEstad[95] = llRef
    ResulEstad[96] = llRef2
    ResulEstad[102] = len(VecLA[:,4])  
               
   
    return ResulEstad[:]
              
