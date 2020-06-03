# coding=utf-8
import numpy as np

def contaravalanV12(nx,ny,VectorP,rho):

    contador = 0
    Nflag = 0
    VecAval = np.zeros(nx * ny)
    iout = 0
    jout = 0
    maxtempo = 0.0
    maxtempo1 = 0.0

    # print('nx',nx,'ny',ny,'length(VectorP[2,:])',len(VectorP[2,:]), 'length(VectorP[:,2])',len(VectorP[:,2]))

    VecPNotBord = VectorP[1:nx, 1:ny]
    VecAval = VecPNotBord[VecPNotBord >= 1]
    contador = len(VecAval)

    if contador == 0:
        Nflag = 0
    else:
        Nflag = 1

    if Nflag == 0:
        rN = 0
        lN = 0
        suma = 0
        sumarho = 0
        tempsmig = 0
        vecNorm = np.zeros(len(VecPNotBord))
        contNorm = 0
        MenosNega = VecPNotBord[VecPNotBord >= 0]
        suma = sum(MenosNega)
        sumarho = sum(np.power(MenosNega, rho))
        contNorm = len(MenosNega)
        vecNorm = np.copy(MenosNega)

        tempsmig = 1.0 / sumarho
        probkFuera = np.zeros(nx * ny)
        limiteinferior = 0
        probab = np.random.rand()  # numero aleatorio entre 0 y 1 para saber que elemento fallará se utiliza un número aleatorio entre 0 y 1,se recorre la lista
        # de elementos (no importa el orden) y se van sumando uno a uno. Cuando con la suma de un elemento x se sobrepasa el valor del
        # número aleatorio que se había tomado, el elemento x será el indicado para fallar

        contador2 = 0
        ioutN = 0
        joutN = 0
        contFueraX = 0
        eleminfalla = 0

        for i,_ in enumerate(vecNorm):
            limiteinferior = contador2
            probkFuera[i] = np.power(vecNorm[i], rho) * tempsmig  # se calcula la probabilidad pk de cada elemento k-ésim
            contador2 += probkFuera[i]  # se van sumando acumulada de la probabilidad de cada una de las celdas.
            if limiteinferior < probab < contador2:  # se compara si el numero aleatorio elegido se encuentra dentro del intervalo
                eleminfalla = vecNorm[i]
                C = np.where(VectorP == eleminfalla)
                ioutN = C[0][0]
                joutN = C[1][0]
                maxtempo1 = VectorP[ioutN, joutN]
                break

        contFuera = 0
        iout = ioutN
        jout = joutN

        return contador, Nflag, VecAval[1:contador], maxtempo1, iout, jout

    elif Nflag == 1:

        B =  np.where(VectorP==np.max(VectorP))
              #     println(size(VectorP),"size(VectorP)",indmax(VectorP),"indmax(VectorP)")
        iout = B[0][0]
        jout = B[1][0]
        maxtempo = VectorP[iout, jout]

        return contador, Nflag, VecAval[1:contador], maxtempo, iout, jout
