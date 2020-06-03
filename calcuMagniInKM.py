# coding=utf-8
import numpy as np
import scipy
import sympy as sp
from scipy import optimize
from scipy.misc import derivative
from scipy.optimize import minimize
from scipy.optimize import differential_evolution

def calcuMagniInKM(NumEventos, AreaCelda):
    print('NUmEventos',NumEventos)
    print('AreaCelda', AreaCelda)
# Darfield earthquake test0 (200 km x 200 km) = area 40000 km²
# Landers earthquake (50 km x 50 km) = area 2500 km²
    VecArea = AreaCelda * NumEventos
    print('VecArea',VecArea)
    VecMagLE = np.log10(VecArea) + 3.99
    VecMagHB = np.log10(VecArea) + 3.98
    VecMagWY = np.log10(VecArea) + 4.1
    VecMagHB1 = 4 / 3 * np.log10(VecArea) + 3.07
    VecMagSTRA = 4.441 + (0.846 * np.log10(VecArea))
    VecMagVIL = 3.39 + (1.33 * np.log10(VecArea))
    VecMagWECO = 4.07 + (0.98 * np.log10(VecArea))

    nboot = 50
    npor = 1

##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1Cálculo de b mediante Minimos cuadrados!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    paramGRHB1, cuenrepHB1, vecWriteMagniHB1 = gutenberRich(VecMagHB1, 1)
    # paramGRVIL, cuenrepVIL, vecWriteMagniVIL = gutenberRich(VecMagVIL, 2)
    # paramGRWYS, cuenrepWYS, vecWriteMagniWYS = gutenberRich(VecMagWY, 3)
    # paramGRHB, cuenrepHB, vecWriteMagniHB = gutenberRich(VecMagHB, 4)
    # paramGRLE, cuenrepLE, vecWriteMagniLE = gutenberRich(VecMagLE, 5)
    # paramGRSTRA, cuenrepSTRA, vecWriteMagniSTRA = gutenberRich(VecMagSTRA, 6)
    # paramGRWECO, cuenrepWECO, vecWriteMagniWECO = gutenberRich(VecMagWECO, 7)

##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Cálculo de b MLE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    meanmHB1, bHB1, sigHB1, av2HB1 = bmemag(VecMagHB1)
    # meanmHB, bHB, sigHB, av2HB = bmemag(VecMagHB)
    # meanmVIL, bVIL, sigVIL, av2VIL = bmemag(VecMagVIL)
    # meanmWYS, bWYS, sigWYS, av2WYS = bmemag(VecMagWY)
    # meanmLE, bLE, sigLE, av2LE = bmemag(VecMagLE)
    # meanmWECO, bWECO, sigWECO, av2WECO = bmemag(VecMagWECO)
    # meanmSTRA, bSTRA, sigSTRA, av2STRA = bmemag(VecMagSTRA)

    # println(meanmHB1, "_", bHB1, "_", sigHB1, "_", av2HB1, "", paramGRHB1[4], "", paramGRHB1[5])
    # println(meanmHB, "_", bHB, "_", sigHB, "_", av2HB, "", paramGRHB[4], "", paramGRHB[5])
    # println(meanmVIL, "_", bVIL, "_", sigVIL, "_", av2VIL, "", paramGRVIL[4], "", paramGRVIL[5])
    # println(meanmWYS, "_", bWYS, "_", sigWYS, "_", av2WYS, "", paramGRWYS[4], "", paramGRWYS[5])
    # println(meanmSTRA, "", bSTRA, "", sigSTRA, "", av2STRA, "", paramGRSTRA[4], "", paramGRSTRA[5])
    # println(meanmLE, "", bLE, "", sigLE, "", av2LE, "", paramGRLE[4], "", paramGRLE[5])
    # println(meanmWECO, "", bWECO, "", sigWECO, "", av2WECO, "", paramGRWECO[4], "", paramGRWECO[5])

#     figure(88,figsize = (10,5))
#     xlabel(L"$Magnitude$",fontsize = 14,weight="demi")
#     ylabel(L"$Cumulative \hspace{.3} number \hspace{.3} of \hspace{.3} aftershocks$",fontsize=14,weight="demi")
#     PyPlot.minorticks_off()
#     hold("True")
#     xticks(fontsize=13,weight="demi")
#     yticks(fontsize=13,weight="demi")
#     semilogy(vecWriteMagniHB1[:,1],(vecWriteMagniHB1[:,2]),"ok",lw=2,label=L"$Simulated \hspace{.3} Series$")
#     Yfit =((-bHB1*vecWriteMagniHB1[:,1])+ av2HB1)
#     PyPlot.plot(vecWriteMagniHB1[:,1],10.^Yfit,"--r",label = L"$Gutenberg-Richter \hspace{.3} fit$")
#     xlim(0,8)
#     annotate("a =" * string(round(av2HB1,2)), fontsize=12,weight="bold",[1,3])
#     annotate("b =" * string(round(bHB1,2)), fontsize=12,weight="bold",[1,2])
#     annotate("rho =" * string(round(1-sigHB1,2)), fontsize=12,weight="bold",[1,1])
#     #annotate(L"$Landers \hspace{.5} M_{min}=1.6$", fontsize=14,weight="bold",[6,10])
#     legend(fontsize="medium", bbox_to_anchor=(0.8,0.3))
#     savefig("GR-Simulated.eps",dpi = 500)


 #   NumberValueB = cuenrepHB1 / len(VecMagHB1)
    nflag = 1

## Non-extensive formulation
    np.savetxt('vecWriteMagniHB1.dat',vecWriteMagniHB1)
    if len(vecWriteMagniHB1) != 0:
       aVLSSilva, qVLSSilva, minVLSSilva, aVLSTelesca, qVLSTelesca, minVLSTelesca = valorq1(vecWriteMagniHB1, nflag)
    else:
       aVLSSilva = 0
       qVLSSilva = 0
       minVLSSilva = 0
       aVLSTelesca = 0
       qVLSTelesca = 0
       minVLSTelesca = 0

    return meanmHB1, bHB1, sigHB1, av2HB1, paramGRHB1[0], paramGRHB1[1], paramGRHB1[2], paramGRHB1[3], paramGRHB1[4],aVLSSilva, qVLSSilva, minVLSSilva, aVLSTelesca, qVLSTelesca, minVLSTelesca, len(vecWriteMagniHB1), cuenrepHB1

def bmemag(b):
    # function calculates the mean magnitute, the b value based
    # on the mean and the standart deviation
    # Stefan Wiemer 03/95
    newcat = np.copy(b)
    maxmag = np.max(newcat)
    mima = np.min(newcat)

    if mima > 0:
        mima = 0

    # calculate the mean magnitude, b(mean) and std
    n = len(newcat)
    meanm1 = np.mean(newcat)
    b1 = (1 / (meanm1 - np.min(newcat - 0.05))) * np.log10(np.exp(1))
    sig1 = (sum(np.power((newcat - meanm1),2))) / (n * (n - 1))
    sig1 = np.sqrt(sig1)
    sig1 = 2.30 * sig1 * np.power(b1,2)  # standard deviation
    # disp ([' b-value segment 1 = ' num2str(b1) ]);
    # disp ([' standard dev b_val_1 = ' num2str(sig1) ]);
    av2 = np.log10(len(newcat)) + b1 * np.min(newcat)

    return meanm1, b1, sig1, av2

def valorq1(VecGR, nflag):
     nflag = 0
     # VecGR = readdlm("actopan.dat")
     # VecGR = readdlm("VecGutenberg-Rich.dat")
     X = VecGR[:, 0]  # magnitude de maxima a mínima
     Y = VecGR[:, 1]  # numero acumulado de eventos de mínimo a máximo
     N = np.max(Y)  # normaliza la curva
     maxMag = np.max(X)
     minMag = np.min(X)
     interval = 0.1
     interMag = int((maxMag - minMag) / interval)
     Y1 = np.log10(Y / (N - 1))
    # Y1 = (Y / (N - 1))
     Y1t = Y1[1:len(Y1)]
     Xt = X[1:len(Y1)]
     VecBorrar = np.zeros(len(Xt))
     contBorrar = -1
     print(Xt)
     for i,_ in  enumerate(Xt):
         if Xt[i] == float('-inf') or Xt[i] == float('inf') or Xt[i] == 0 or Y1t[i] == 0:
             print("##################", i, "-INF", Xt[i], "____", len(Xt))
             contBorrar += 1
             VecBorrar[contBorrar] = i
     Xt = np.delete(Xt,VecBorrar[0:contBorrar],0)
     Y1t = np.delete(Y1t,VecBorrar[0:contBorrar],0)
     if len(Y1t) > 3:
         print('HOLA####len(Y1t)')
         print(Y1t)
         print(Xt)
        # fun = lambda p:  (p[0] - 1)**2 + (p[1] - 2.5)**2
        # bnds = ((0.25, 0.75), (0, 2.0))  (2, 0)
         p = np.array([0,0])
         def fun(p):
             A = np.power(10, np.multiply(2, Xt))
             B = np.power(p[1],(2/3))
             C = 1-((1-p[0])/(2-p[0]))
             D = (2-p[0])/(1-p[0])
             E = Y1t - (D*(np.log10(C*(A/B))))
             return  np.sqrt(sum(np.power(E,2)))
         bnds = ((1.01, 1.99), (1e+1, 1e+10))
         pguess = (1.6, 1e+5)
       #  jacobia = derivative(fun, pguess)
         res = minimize(fun, pguess, bounds = bnds)
         print(res.x)
      #   res = minimize(NonExten_Silva, , bounds= (1.0, 2.0))  # ,optimizer=GradientDescent)
         print(res)
        # resultados
         pres1 = res.x[0]
         pres2 = res.x[1]
         rmse_Silva = 0
         
         pT = np.array([0,0])
         def funT(pT):
             A = np.power(10,Xt)
             B = np.power(pT[1],(2/3))
             C = 1-((1-pT[0])/(2-pT[0]))
             D = (2-pT[0])/(1-pT[0])
             E = Y1t - (D*(np.log10(C*(A/B))))
             return  np.sqrt(sum(np.power(E,2)))
         bnds = ((1.01, 1.99), (1e+1, 1e+10))
         pguess = (1.6, 1e+5)
       #  jacobia = derivative(fun, pguess)
         resT = minimize(funT, pguess, bounds = bnds)
         print(resT.x)
      #   res = minimize(NonExten_Silva, , bounds= (1.0, 2.0))  # ,optimizer=GradientDescent)
         print(resT)
         pres1_1 = resT.x[0]
         pres1_2 = resT.x[1]
         rmse_Telesca = 0     
         
     elif len(Y1t) <= 3:
         pres1 = 0
         pres2 = 0
         rmse_Silva = 0
         pres1_1 = 0
         pres1_2 = 0
         rmse_Telesca = 0
         
     return pres1, pres2, rmse_Silva, pres1_1, pres1_2, rmse_Telesca

def gutenberRich(VecMagHB1,num):

    VecMagni = np.copy(VecMagHB1)
    print(VecMagni)
    contador = 0
    paramGRAreaWY = np.zeros(5)
    magMax = np.max(VecMagni)
    magMin = np.min(VecMagni)
    magRang = magMax-magMin
    print("magMax=",magMax)
    print("magMin=",magMin)
    print("magRang=",magRang)
    numMag = int(magRang/0.10)
    limMin = np.zeros(numMag)
    NumSismos = np.zeros(numMag)
    vecWrite = np.zeros((numMag,2))

    for i,_ in enumerate(limMin):
        limMin[i]=magMin+((i*0.10)-0.10)
        for j,_ in enumerate(VecMagni):
            if VecMagni[j] >= limMin[i]:
                contador = contador+1
        NumSismos[i] = contador
        contador = 0

    cuenrep = 0
    for i,_ in enumerate(limMin[0:numMag-1]):
        if NumSismos[i+1] == NumSismos[i]:
            cuenrep += 1

    if magRang > 0:

        vecWrite[0:numMag,0] = limMin[0:numMag]
        vecWrite[0:numMag,1] = NumSismos[0:numMag]
        p1,p2 = np.polyfit(limMin[0:numMag],np.log10(NumSismos[0:numMag]),1)
        rho = np.corrcoef(limMin[0:numMag],np.log10(NumSismos[0:numMag]))

        paramGRAreaWY[0] = magMax
        paramGRAreaWY[1] = magMin
        paramGRAreaWY[2] = p1
        paramGRAreaWY[3] = p2
        paramGRAreaWY[4] = rho[0,0]

    elif magRang == 0:

        paramGRAreaWY[0] = 0
        paramGRAreaWY[1] = 0
        paramGRAreaWY[2] = 0
        paramGRAreaWY[3] = 0
        paramGRAreaWY[4] = 0


    return  paramGRAreaWY, cuenrep, vecWrite[0:numMag,:]
