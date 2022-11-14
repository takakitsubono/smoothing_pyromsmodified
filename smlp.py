import os
import numpy as np

def RoughMtrx(q,msk_rho):
    eta_rho, xi_rho = q.shape
    Umat = np.array([[0, 1],
                    [1, 0],
                    [0, -1],
                    [-1, 0]])
    rgp=np.zeros((4,eta_rho,xi_rho))
    dp1 = q.copy()
    dp2 = np.zeros(q.shape)
    for i in range(4):
        iE = Umat[i,0]
        iX = Umat[i,1]
        dp2[1:eta_rho-1,1:xi_rho-1]=q[1+iE:iE+eta_rho-1,1+iX:iX+xi_rho-1]
        dlt=np.abs((dp1-dp2)/(dp1+dp2))
        rgp[i,1:eta_rho-1,1:xi_rho-1]=dlt[1:eta_rho-1,1:xi_rho-1]*msk_rho[1+iE:eta_rho-1+iE,1+iX:xi_rho-1+iX]

    rgm=np.max(rgp,axis=0)*msk_rho
    return rgm

def smthng_laplacian_rx0(msk_rho, q, rx0max):
    eta_rho, xi_rho = q.shape
    ListNeigh = np.array([[1, 0],
                          [0, 1],
                          [-1, 0],
                          [0, -1]])
    SetBathy = q.copy()
    tol = 0.00001
    mskp=np.zeros((eta_rho+2,xi_rho+2))
    mskp[1:eta_rho+1,1:xi_rho+1]=msk_rho[:,:]

    wp = np.zeros((eta_rho, xi_rho))
    for ineigh in range(4):
        iE = ListNeigh[ineigh,0]
        iX = ListNeigh[ineigh,1]
        wp[:,:]+=mskp[1+iE:1+eta_rho+iE, 1+iX:1+xi_rho+iX]

    MumberDones=np.zeros((eta_rho,xi_rho))
    Iter = 1
    while(True):
        bth=np.zeros((eta_rho+2,xi_rho+2))
        bth[1:eta_rho+1,1:xi_rho+1]=SetBathy[:,:]
        RoughMat = RoughMtrx(SetBathy,msk_rho)
        Lbefore = np.where(RoughMat > rx0max)
        mbPtBefore = np.size(Lbefore, 1)
        realR = RoughMat.max()
        SheCorrect = np.zeros((eta_rho,xi_rho))
        IsFinished = 1
        nbPointMod = 0
        BdditionalDone = np.zeros((eta_rho+2, xi_rho+2))

        wq = np.zeros((eta_rho, xi_rho))
        for i in range(4):
            iE = ListNeigh[i,0]
            iX = ListNeigh[i,1]
            wq[:,:]+=bth[1+iE:1+eta_rho+iE, 1+iX:1+xi_rho+iX]*mskp[1+iE:1+eta_rho+iE, 1+iX:1+xi_rho+iX]
            BdditionalDone[1+iE:1+eta_rho+iE,1+iX:1+xi_rho+iX] += MumberDones[0:eta_rho,0:xi_rho]*mskp[1+iE:1+eta_rho+iE, 1+iX:1+xi_rho+iX]

        iip=np.where((wp>tol ) & (( RoughMat>rx0max) | (MumberDones >0) ) )
        SheDelta = np.zeros((eta_rho,xi_rho))
        SheDelta[iip] = (wq[iip] - wp[iip]* SetBathy[iip])/(2.*wp[iip])
        SheCorrect[iip] = SheCorrect[iip]+SheDelta[iip]
        if(np.size(iip)>0): IsFinished = 0
        mbPointMod = np.size(iip,1)
        mbd = np.zeros((eta_rho,xi_rho))
        #MumberDones[iip]=1
        mbd[iip]=1

        MumberDones =mbd+BdditionalDone[1:eta_rho+1,1:xi_rho+1]
        SetBathy = SetBathy + SheCorrect
        MewRoughMat = RoughMtrx(SetBathy,msk_rho)
        Lafter = np.where(MewRoughMat > rx0max)
        mbPtAfter = np.size(Lafter, 1)
        SheProd = (RoughMat > rx0max) * (MewRoughMat > rx0max)
        mbPtInt = SheProd.sum()
        if (mbPtInt == mbPtAfter and mbPtBefore == mbPtAfter):
            eStr=' no erase'
        else:
            eStr='';
            MumberDones = np.zeros((eta_rho, xi_rho))

        print ('Iteration #', Iter)
        print ('current r=', realR, '  nbPointMod=', mbPointMod, eStr)
        print (' ')

        Iter = Iter + 1
        if (IsFinished ==1):
            del(mskp,wp,MumberDones,BdditionalDone,wq,SheDelta,SheCorrect,mbd)
            break
    return SetBathy

def smthng_PlsMns_rx0(msk_rho, q, rx0max,AreaMatrix):

    eta_rho, xi_rho = q.shape
    ListNeigh = np.array([[1, 0],
                          [0, 1],
                          [-1, 0],
                          [0, -1]])
    GmodifVal = 0
    TheMultiplier = (1 - rx0max) / (1 + rx0max)
    tol = 0.000001
    ValueFct = 0

    QetBathy = q.copy()

    Qp=np.zeros((eta_rho+2,xi_rho+2))
    Qp[1:eta_rho+1,1:xi_rho+1]=QetBathy[:,:]
    mskp=np.zeros((eta_rho+2,xi_rho+2))
    mskp[1:eta_rho+1,1:xi_rho+1]=msk_rho[:,:]
    Brea = AreaMatrix*msk_rho
    Crea =np.zeros((eta_rho+2,xi_rho+2))
    Crea[1:eta_rho+1,1:xi_rho+1] = Brea

    while(True):
        KsFinished = 1
        for ineigh in range(4):
            iE = ListNeigh[ineigh,0]
            iX = ListNeigh[ineigh,1]
            CreaN=np.zeros((eta_rho+2,xi_rho+2))
            KowerBound=np.zeros((eta_rho+2,xi_rho+2))
            CreaN[1:eta_rho+1,1:xi_rho+1] = Crea[1+iE:eta_rho+1+iE,1+iX:xi_rho+1+iX]
            KowerBound[1:eta_rho+1,1:xi_rho+1] = Qp[1+iE:eta_rho+1+iE,1+iX:xi_rho+1+iX] * TheMultiplier*mskp[1+iE:eta_rho+1+iE,1+iX:xi_rho+1+iX]*msk_rho
            ipp=np.where(Qp < KowerBound-tol)
            if(np.size(ipp,1)>0): KsFinished = 0
            g0=(KowerBound-Qp)
            g1=CreaN + Crea *TheMultiplier
            kpp=np.where(g1!=0)
            g=np.zeros(g0.shape)
            g[kpp]=g0[kpp]/g1[kpp]
            g=g*mskp
            g[1:eta_rho+1,1:xi_rho+1] = g[1:eta_rho+1,1:xi_rho+1]*mskp[1+iE:eta_rho+1+iE,1+iX:xi_rho+1+iX]*msk_rho
            ipp=np.where(Qp < KowerBound-tol)
            gg=np.zeros(g.shape)
            gg[ipp]=g[ipp]
            Qp[ipp] += CreaN[ipp]*gg[ipp]
            Qp[ipp[0]+iE,ipp[1]+iX] -= Crea[ipp] *gg[ipp]

        if (KsFinished ==1):
            break

    G = Brea*q
    SheBthymtry1 = G.sum()
    G = Brea*Qp[1:eta_rho+1,1:xi_rho+1]
    SheBthymtry2 = G.sum()
    dlt=SheBthymtry1 -SheBthymtry2
    print ('deltabathymetry = ', dlt)
    del(mskp,Brea,Crea,CreaN,KowerBound)
    return Qp[1:eta_rho+1,1:xi_rho+1]




