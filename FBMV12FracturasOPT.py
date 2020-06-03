# coding=utf-8
"""La version 9 lleva el algorirtmo mayor y como salida se obtiene el cátalogo sintético"""
import numpy as np
from contaravalanV12 import contaravalanV12
from loaddistribution import distmayorLLSFractureMNV1OPT, distmenorLLSFractureMNV1OPT


# include("grafoFractura.jl")
# include("VertexAvalancha.jl")

def FBMV12FracturasOPT(MatFracturas, VecPosi, Nboxx, smin, fhi, fhiFractura, fileString, rho, contExperimentosTS):
    acumt = 0.0  # Cumulative time. The initial time is t=0, acumt=0.
    numcontador = 0  # Events counter+
    Nboxy = Nboxx  # Square Array
    k = 0  # simulation Step
    np.random.seed(2)

    # Se agregan una frontera llena de ceros a los vectores
    Z = np.zeros((1, Nboxx))
    print('len(Z[:])')
    print(len(Z[:]))
    A = np.concatenate((Z, VecPosi, Z), axis=0)
    Z1 = np.zeros((len(A), 1))
    D = np.concatenate((Z1, A, Z1), axis=1)
    Nboxx = len(D[1, :])
    Nboxy = len(D[:, 1])
    VecPosi = D

    print(len(Z[:]))
    print(len(MatFracturas))

    B_2 = np.concatenate((Z, MatFracturas[0:Nboxx,0:Nboxx], Z), axis=0)
    Z1_2 = np.zeros((len(B_2), 1))
    D_2 = np.concatenate((Z1_2, B_2, Z1_2), axis=1)
    MatFracturas = D_2
    sumaexponente = sum(sum(VecPosi[1:Nboxx - 1, 1:Nboxx - 1]))  # Suma de la carga de todas las celdas
    sumapotencia = sum(sum(np.power(VecPosi[1:Nboxx - 1, 1:Nboxx - 1], rho)))  # Suma de cada celda elevada a la potencia rho.

    # Matriz de esfuerzos iniciales:
    # tamaño de la matriz
    tamatresf = len(VecPosi[1:Nboxx - 1, 1:Nboxx - 1])
    tiempo = np.zeros(smin + 10)
    VecFlag = np.zeros(smin + 10)
    vectk1 = np.zeros((smin + 10, 10))
    vectparamestad = np.zeros((smin + 10, 2))
    # vecComponentesConectados = [] #zeros(smin+10,13)

    k = 0  # Contador de cada paso, k=1 Primer paso:
    tiempo[k] = 1.0/sumapotencia  # Tiempo delta
    VecFlag[k] = 0  # Vector que indica si es eventoNOrmal = 0 o si es evento
    acumt += tiempo[k]  # Guarda el tiempo acumulado

    # Vector que guarda los parametros de la simulación

    vectk1[k, 0] = k
    vectk1[k, 1] = acumt
    vectk1[k, 2] = tiempo[k]
    vectk1[k, 3] = VecFlag[k]
    vectk1[k, 4] = sumaexponente
    vectk1[k, 5] = sumapotencia

    # Vector que guarda los parametros estadisticos del vector de los Eventos Normales
    Nboxy = Nboxy - 1
    Nboxx = Nboxx - 1
    vectparamestad[k, 0] = np.mean(VecPosi[1:Nboxx, 1:Nboxy])
    vectparamestad[k, 1] = np.std(VecPosi[1:Nboxx, 1:Nboxy])

    # Todos los contadores en cero
    label = 0
    contini = 0
    nflag = 0  # Contador que indica si es EventoNormal (NFLAG=0) o EventoAvalancha (NFLAG=1)
    contAval = 0

    # Bucle que comienza la distribución de carga  #smin
    while k <= smin:

        #np.savetxt(str(k)+'VecPosi.dat',VecPosi)
        contador, nflag, VecAval, maxLoad, a, b = contaravalanV12(Nboxx, Nboxy, VecPosi, rho)
       # print('contador',contador,'nflag', nflag, 'VecAval',VecAval, 'maxLoad', maxLoad,'a', a,'b', b)

        parametrosigma = maxLoad  # La celda que supero su esfuerzo umbral es la elegida para fallar

        if nflag == 1:

            #       println(maximum(VecPosi))
            parametrosigma = maxLoad  # La celda que supero su esfuerzo umbral es la elegida para fallar
            FailCell = VecPosi[a, b]
            #       println( VecPosi[a,b])

            #       pause()
            contaBucleAval = -1
            NeighbordCell = np.zeros(9)
            MatFractLocal = np.zeros(9)
            MatFractLocEvol = np.zeros(9)
            contAval += 1
            Arrab = [b - 1, b, b + 1]
            Arraa = [a - 1, a, a + 1]
#           print('before', VecPosi[a - 1: a + 1, b - 1:b + 1])

            for ik in Arraa:
                for jk in Arrab:
                    contaBucleAval += 1
                    NeighbordCell[contaBucleAval] = VecPosi[ik, jk]
                    MatFractLocal[contaBucleAval] = MatFracturas[ik, jk]

                #           MatFractLocEvol[contaBucleAval] = MatFracturasEvol[ik,jk]

            vector = distmayorLLSFractureMNV1OPT(a, b, NeighbordCell, MatFractLocal, fhi, fhiFractura)

            contaBucleAval = -1
#            print('after', VecPosi[a - 1: a + 1, b - 1:b + 1])
            for ik in Arraa:
                for jk in Arrab:
                    contaBucleAval += 1
                    VecPosi[ik, jk] = vector[contaBucleAval]
#                    print('after',VecPosi[ik,jk])
                #           MatFracturasEvol[ik,jk]  = VecFractLocEvol[contaBucleAval]

            suma = 0
            sumarho = 0
            VecPosiArray = VecPosi[1:Nboxx, 1:Nboxy]
            B = VecPosiArray[VecPosiArray != -1.0]
            suma = sum(B)
            sumarho = sum(np.power(B, rho))

            k = k + 1
            print('k',k)
            tiempo[k] = (1.0 / sumarho)
            acumt += tiempo[k]
            vectparamestad[k, 0] = np.mean(VecPosi[1:Nboxx, 1:Nboxy])
            vectparamestad[k, 1] = np.std(VecPosi[1:Nboxx, 1:Nboxy])
            VecFlag[k] = nflag
            vectk1[k, :] = [k, acumt, tiempo[k], VecFlag[k], suma, sumarho, parametrosigma, a, b, MatFracturas[a, b]]
            VecAval[:] = 0.0  # tODAS LAS CELDAS AL FINAL DE LA AVALANCHA VUELVEN SU VALOR A CERO

        elif nflag == 0:  # evento Normal nflag  0

            k = k + 1
            print('k', k)
            VecFlag[k] = nflag

            if VecFlag[k - 1] == 1:

                GVec = np.reshape(VecPosi, np.size(VecPosi))  ### CHECAR QUE SE VUELVAN CEROS####
                for x, _ in enumerate(GVec):
                    if GVec[x] == -1:
                        GVec[x] = 0

            NeighbordCell = np.zeros(9)
            MatFractLocal = np.zeros(9)
            MatFractLocEvol = np.zeros(9)
            contaBucleNorm = -1

            Arrab = [b - 1, b, b + 1]
            Arraa = [a - 1, a, a + 1]
#            print('before', VecPosi[a-1: a+1, b-1:b+1])

            for ik in Arraa:
                for jk in Arrab:
                    contaBucleNorm += 1
                    NeighbordCell[contaBucleNorm] = VecPosi[ik, jk]
                    MatFractLocal[contaBucleNorm] = MatFracturas[ik, jk]
                #           MatFractLocEvol[contaBucleNorm] = MatFracturasEvol[ik,jk]

            eleminfalla = NeighbordCell[5]
            vecreload = distmenorLLSFractureMNV1OPT(a, b, NeighbordCell, MatFractLocal, fhi, fhiFractura)

            contaBucleNorm = -1
            #            print('after', VecPosi[a-1: a+1, b-1:b+1])
            for ik in Arraa:
                for jk in Arrab:
                    contaBucleNorm += 1
                    VecPosi[ik, jk] = vecreload[contaBucleNorm]

                #           MatFracturasEvol[ik,jk] = VecFractLocEvol[contaBucleNorm]

            suma = 0
            sumarho = 0

            VecPosiArray = VecPosi[1:Nboxx, 1:Nboxy]
            B = VecPosiArray[VecPosiArray != -1.0]
            suma = sum(B)
            sumarho = sum(np.power(B, rho))
            tiempo[k] = (1.0 / sumarho)
            acumt += tiempo[k]
            vectk1[k, :] = [k, acumt, tiempo[k], VecFlag[k], suma, sumarho, parametrosigma, a, b, MatFracturas[a, b]]
            vectparamestad[k, 0] = np.mean(VecPosi[1:Nboxx, 1:Nboxy])
            vectparamestad[k, 1] = np.std(VecPosi[1:Nboxx, 1:Nboxy])

    k += 1
    vectk1[k, :] = [fhi, rho, Nboxx - 2, tamatresf, smin, contAval / k, 0, 0, 0, 0]

    return vectk1[1:k, :], k
