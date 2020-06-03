import numpy as np
import math
# calcula las dimensiones fractales D0, D1 y D2

def CorreDimensionAvalanchesV11(XX):
    print(XX[:, 0])
    print(XX[:, 1])
    xmin = np.min(XX[:, 0])
    xmax = np.max(XX[:, 0])
    ymin = np.min(XX[:, 1])
    ymax = np.max(XX[:, 1])
    rmax = 300  # maximum size of window  300
    rmin = 15  # minimum size of window  15
    nr = 50  # Number of rs to fit
    n = len(XX[:, 1])
    nm = n - 1
    x = XX[:, 0]
    y = XX[:, 1]
    ds = np.zeros((n, n))  # matrix of squared interevent distances
    conti = -1
    contj = -1
    for i,_ in enumerate(x[0:len(x)-1]):
        conti = conti + 1
        xi = x[i]
        yi = y[i]
        for j,_  in enumerate(x[i+1:len(x)]):
            contj = contj + 1
            ds[i, j+1] = np.power((xi - x[j+1]),2) + np.power((yi - y[j+1]),2)
            ds[j+1, i] = ds[i, j+1]
    lgrmin = np.log10(rmin)
    lgrmax = np.log10(rmax)
    dlgr = (lgrmax - lgrmin) / (nr - 1)
    lgr = np.zeros(nr)
    for i,_ in enumerate(lgr):
        lgr[i] = lgrmax - (((i+1) - 1) * dlgr)  # equispaced log r array

    rs = np.power(10,(2*lgr))  # r^2  array
    c = np.zeros((3, nr))  # Correlation integrals 0, 1, & 2 array
    i = 0
    j = 0
    k = 0
    ni = np.zeros(n)
    nii = np.zeros(n)
    for k,_ in enumerate(lgr):  # ==========  Loop for each length ============================= #
        value = rs[k]
        ni = np.zeros(n)
        nii = np.zeros(n)
        nii2 = np.zeros(n)
        nii3 = np.zeros(n)
        for i,_ in enumerate(ni):
            cont = 0
            for j,_ in enumerate(ni):
                if ds[i, j] < value:
                    cont += 1
            ni[i] = (float(cont) - 1) / nm
        nii = ni[ni != 0]
        contix = -1
        for ij,_ in enumerate(nii):
            contix += 1
            nii2[contix] = ni[int(nii[ij])]
        contnotnan = -1
        for jj,_ in enumerate((nii2[0:contix])):
            if nii2[jj] != float('nan'):
                contnotnan += 1
                nii3[contnotnan] = nii2[jj]
            else:
                println("NaN number")
                println(nii2[jj])
        if contnotnan == -1:
            nr = k - 1
            lgr = lgr[0:nr]
            c = c[:, 0:nr]
            break
        c[0, k] = sum(np.power(nii3[0:contnotnan],-1)) / n
        c[1, k] = sum(np.log10(nii3[0:contnotnan])) / n  # C1
        c[2, k] = sum(nii3[0:contnotnan]) / n  # C2
    c[0, :] = -np.log10(c[0, :])
    c[2, :] = np.log10(c[2, :])

# ======================================================================  Function dimfit #

    def dimfit(lgr, lgc):
        nr = len(lgr)
        ix1 = 0
        ix2 = nr-1
        lgcm = abs(lgc[ix1] + lgc[ix2]) / 2
        print('lgcm',lgcm)
        p1, p2 = np.polyfit(lgr[ix1:ix2], lgc[ix1:ix2], 1)
        print('p1',p1,'p2',p2)
        rho = np.corrcoef(lgr[ix1:ix2], lgc[ix1:ix2])
        if lgcm != 0:
            f = 100 * rho[0,0] / lgcm
        else:
            f = 0
        ix1a = ix1
        ix2a = ix2
        return p1, p2, rho, f, len(lgr), len(lgc)


# ====================================================================== END Function dimfit #

    D0 = np.zeros(6)
    D1 = np.zeros(6)
    D2 = np.zeros(6)
    D = np.zeros((3, 6))

    for i,_ in enumerate(c):  # ----------  Presents each C for D fitting --------------------#
        if len(c[i, :]) != 0:
            p1, p2, rho, f, lenlgr, lenlgc = dimfit(lgr, c[i, :])
            D[i, 0] = p1
            D[i, 1] = p2
            D[i, 2] = rho[0,0]
            D[i, 3] = f
            D[i, 4] = lenlgr
            D[i, 5] = lenlgc
        else:
            D[i, 0] = 0
            D[i, 1] = 0
            D[i, 2] = 0
            D[i, 3] = 0
            D[i, 4] = 0
            D[i, 5] = 0

    D0[:] = D[0, :]
    D1[:] = D[1, :]
    D2[:] = D[2, :]
    return D0, D1, D2

