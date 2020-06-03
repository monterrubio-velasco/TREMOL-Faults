# coding=utf-8
#This function computes the load distribution in an Avalanche event.
import numpy as np

def distmenorLLSFractureMNV1OPT(a,b,vector1,MatFracturas,fhi,Fraction):

    if MatFracturas[4] == 0:  # la celda esta fuera de la fractura
        vecEsfTot = vector1[4] * fhi
        numProb = np.random.rand()
        # Se reparte en diagonal (NW-NE-SW-SE) un porcentaje pequeño de Esfuerzo
        vecEsfDiag = vecEsfTot * 0.02  # 2% de esfuerzo a las vecinas NE-NW-SE-SW
        factDiag = vecEsfDiag / 4  # las cuatro direcciones
        # Se reparte en N-S-E-W un porcentaje mayor de Esfuerzo
        vecEsf = vecEsfTot * 0.98  # 98% de esfuerzo a las vecinas N-S-E-W
        factPerpe = vecEsf / 4
        vector1[4] = 0.0
        vector1[7] += factPerpe
        vector1[1] += factPerpe
        vector1[3] += factPerpe
        vector1[5] += factPerpe

    # vectores orientados SE-SW-NE-NW

        vector1[8] += factDiag
        vector1[2] += factDiag
        vector1[6] += factDiag
        vector1[0] += factDiag


    elif MatFracturas[4] == 1:

        contVecinosFracturaPerpend = 0
        contVecinosFracturaDiago = 0
        contVecinosFueraPerpend = 0
        contVecinosFueraDiago = 0

        LoadTotal = vector1[4]
        LoadTotalPercentage = vector1[4] * Fraction
        factPerpe = LoadTotalPercentage * 0.98 / 4
        factDiag = LoadTotalPercentage * 0.02 / 4

        vector1[4] = 0.0
        vector1[7] += factPerpe
        vector1[1] += factPerpe
        vector1[3] += factPerpe
        vector1[5] += factPerpe

# vectores orientados SE-SW-NE-NW

        vector1[8] += factDiag
        vector1[2] += factDiag
        vector1[6] += factDiag
        vector1[0] += factDiag

    return vector1  # ,MatFractLocEvol


def distmayorLLSFractureMNV1OPT(a, b, vector, MatFractLocal, fhi, Fraction):

    suma = 0.0
    sumarho = 0.0
    nflag1 = 0

    # srand(3)  #NOTA:Aqui no se puede colocar una semilla porque sino cada que entra al programa se repiten los mismos valores
    # probab=X[k,6]

    if MatFractLocal[4] == 0:  # la celda esta fuera de la fractura
        vecEsfTot = vector[4] * fhi
        numProb = np.random.rand()
        # Se reparte a las celdas diagonales (NW-NE-SW-SE) un porcentaje pequeño de Esfuerzo
        vecEsfDiag = vecEsfTot * 0.02  # 2% de esfuerzo a las vecinas NE-NW-SE-SW
        factDiag = vecEsfDiag / 4  # las cuatro direcciones

        # Se reparte en N-S-E-W un porcentaje mayor de Esfuerzo
        vecEsf = vecEsfTot * 0.98  # 98 % de esfuerzo a las vecinas N-S-E-W
        factPerpe = vecEsf / 4

    elif MatFractLocal[4] == 1:
        numProb = np.random.rand()
        vecEsfTot = vector[4] * Fraction
        vecEsfDiag = vecEsfTot * 0.02  # 2% de esfuerzo a las vecinas NE-NW-SE-SW
        factDiag = vecEsfDiag / 4  # las cuatro direcciones
        # Se reparte en N-S-E-W un porcentaje mayor de Esfuerzo
        vecEsf = vecEsfTot * 0.98  # 98 % de esfuerzo a las vecinas N-S-E-W
        factPerpe = vecEsf / 4

# Se buscan vecinos prohibidos de entre los ocho que están alrededor de la celda que falla

    if vector[0] == -1.0 or vector[1] == -1.0 or vector[2] == -1.0 or vector[3] == -1.0 or vector[5] == -1.0 or vector[6] == -1.0 or vector[7] == -1.0 or vector[8] == -1:
        vectclon = vecinosLLS(vector, vecEsfTot, a, b, numProb)
        nflag1 = 1

    else:

        # Se reparte la carga

        vector[4] = -1.0
        # vectores orientados N-S-E-W
        vector[7] += factPerpe
        vector[1] += factPerpe
        vector[3] += factPerpe
        vector[5] += factPerpe

        # vectores orientados SE-SW-NE-NW
        vector[8] += factDiag
        vector[2] += factDiag
        vector[6] += factDiag
        vector[0] += factDiag

    # Checar si alguno de los vecinos a los que se les comparte carga tienen etiqueta de fractura = 1, si es así entonces el elemento vector[a,b] también se le cambiará a 1 su valor de MatFractLocEvol[a,b]

    #         for i = 1:9
    #             if MatFractLocEvol[i] == 1 && i !=5
    #                 MatFractLocEvol[i] = 1
    #                 break
    #             end
    #         end

    if nflag1 == 1:
        vector = vectclon

    return vector


def vecinosLLS(vector, vecEsfTot, a, b, numProb):


    vecEsfDiag = vecEsfTot * 0.02  # 2% de esfuerzo a las vecinas NE-NW-SE-SW
    factDiag = vecEsfDiag / 4.0  # las cuatro direcciones

    # Se reparte en N-S-E-W un porcentaje mayor de Esfuerzo
    vecEsf = vecEsfTot * 0.98  # 98% de esfuerzo a las vecinas N-S-E-W
    factPerpe = vecEsf / 4.0

    i = a
    j = b

# Bucle que mira si los vecinos son prohibidos. Si lo son el esfuerzo
# se distribuye entre los vecinos no prohibidos, segun el número de vecinos a repartir

# Opción 1: TODOS LOS VECINOS N-S-E-W SON PROHIBIDOS
    if vector[1] == -1.0 and vector[3] == -1.0 and vector[5] == -1.0 and vector[7] == -1.0:
        vector[4] = -1.0

# Opcion 2: TRES VECINOS SON PROHIBIDOS

    if vector[3] == -1.0 and vector[5] == -1.0 and vector[7] == -1.0 and vector[1] != -1.0:
        vector[1] = vector[1] + vecEsf
        vector[4] = -1.0
    
    elif vector[5] == -1.0 and vector[7] == -1.0 and vector[1] == -1.0 and vector[3] != -1.0:
        vector[3] = vector[3] + vecEsf
        vector[4] = -1.0
     
    elif vector[7] == -1.0 and vector[1] == -1.0 and vector[3] == -1.0 and vector[5] != -1.0:
         vector[5] = vector[5] + vecEsf
         vector[4] = -1.0
    
    elif vector[1] == -1.0 and vector[3] == -1.0 and vector[5] == -1.0 and vector[7] != -1.0:
        vector[7] = vector[7] + vecEsf
        vector[4] = -1.0
    

# Opcion 3: DOS VECINOS SON PROHIBIDOS


    if vector[7] == -1.0 and vector[1] == -1.0 and vector[3] != -1.0 and vector[5] != -1.0:
        vector[3] = vector[3] + (vecEsf / 2)
        vector[5] = vector[5] + (vecEsf / 2)
        vector[4] = -1.0

    elif vector[1] == -1.0 and vector[3] == -1.0 and vector[7] != -1.0 and vector[5] != -1.0:
        vector[7] = vector[7] + (vecEsf / 2)
        vector[5] = vector[5] + (vecEsf / 2)
        vector[4] = -1.0

    elif vector[1] == -1.0 and vector[5] == -1.0 and vector[7] != -1.0 and vector[3] != -1.0:
        vector[6] = vector[6] + (vecEsf / 2)
        vector[2] = vector[2] + (vecEsf / 2)
        vector[4] = -1.0

    elif vector[3] == -1.0 and vector[5] == -1.0 and vector[7] != -1.0 and vector[1] != -1.0:
        vector[7] = vector[7] + (vecEsf / 2)
        vector[1] = vector[1] + (vecEsf / 2)
        vector[4] = -1.0

    elif vector[7] == -1.0 and vector[3] == -1.0 and vector[1] != -1.0 and vector[5] != -1.0:
        vector[1] = vector[1] + (vecEsf / 2)
        vector[5] = vector[5] + (vecEsf / 2)
        vector[4] = -1.0

    elif vector[7] == -1.0 and vector[5] == -1.0 and vector[1] != -1.0 and vector[3] != -1.0:
        vector[1] = vector[1] + (vecEsf / 2)
        vector[3] = vector[3] + (vecEsf / 2)
        vector[4] = -1.0


# Opcion 4: UN VECINO ES PROHIBIDO

    if vector[7] == -1.0 and vector[1] != -1.0 and vector[3] != -1.0 and vector[5] != -1.0:
        vector[1] = vector[1] + (vecEsf / 3)
        vector[3] = vector[3] + (vecEsf / 3)
        vector[5] = vector[5] + (vecEsf / 3)
        vector[4] = -1.0

    elif vector[1] == -1.0 and vector[7] != -1.0 and vector[3] != -1.0 and vector[5] != -1.0:
        vector[7] = vector[7] + (vecEsf / 3)
        vector[3] = vector[3] + (vecEsf / 3)
        vector[5] = vector[5] + (vecEsf / 3)
        vector[4] = -1.0

    elif vector[5] == -1.0 and vector[1] != -1.0 and vector[7] != -1.0 and vector[3] != -1.0:
        vector[1] = vector[1] + (vecEsf / 3)
        vector[7] = vector[7] + (vecEsf / 3)
        vector[3] = vector[3] + (vecEsf / 3)
        vector[4] = -1.0

    elif vector[3] == -1.0 and vector[1] != -1.0 and vector[7] != -1.0 and vector[5] != -1.0:
        vector[1] = vector[1] + (vecEsf / 3)
        vector[7] = vector[7] + (vecEsf / 3)
        vector[5] = vector[5] + (vecEsf / 3)
        vector[4] = -1.0


# Los vecinos diagonales
    D = np.zeros(4)
    D[0] = vector[8]
    D[1] = vector[6]
    D[2] = vector[0]
    D[3] = vector[2]

    vectorXconf, B = veciProhiDiag(a, b, vector, vecEsfDiag, D)

    vector[8] = B[0]
    vector[6] = B[1]
    vector[0] = B[2]
    vector[2] = B[3]

    vector[4] = vectorXconf
    vectclon = np.copy(vector)

    return vectclon



def veciProhiDiag(a,b,vector,vecEsfDiag,A):


    contador = 0
    for i,_ in enumerate(A):
        if A[i] == -1:
            contador = contador + 1

    if contador != len(A):
        num = len(A) - contador
        newsigma = vecEsfDiag / num

        for i,_ in enumerate(A):
            if A[i] != -1:
                A[i] = A[i] + newsigma

        vector[4] = -1.0
    else:
        vector[4] = -1.0

    vectorXconf = vector[4]

    return vectorXconf, A
