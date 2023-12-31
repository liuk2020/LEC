#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# problem.py


import numpy as np
from scipy.linalg import solve
from ..geometry import Surface
from ..toroidalField import ToroidalField
from ..toroidalField import derivatePol, derivateTor, changeResolution
from typing import Tuple


class SurfaceEquilibrium:

    def __init__(self, surf: Surface, iota: float, aveJacobian: float=1.0) -> None:
        self.initSurf(surf)
        self.iota = iota
        self.aveJacobian = aveJacobian
        self.nfp = self.surf.r.nfp
        self.mpol = self.surf.r.mpol
        self.ntor = self.surf.r.ntor
        self.g_thetatheta, self.g_thetaphi, self.g_phiphi = self.surf.mertic
        self.P = self.g_thetaphi*iota + self.g_phiphi
        self.Q = self.g_thetatheta*iota + self.g_thetaphi
        self.D = derivatePol(self.P) - derivateTor(self.Q)

    def initSurf(self, surf: Surface) -> None:
        """
        Change the resolution of the surface! 
        """
        _r = changeResolution(surf.r, 2*surf.r.mpol, 2*surf.r.ntor)
        _z = changeResolution(surf.z, 2*surf.z.mpol, 2*surf.z.ntor)
        self.surf = Surface(_r, _z)
    
    def run(self):
        self.Jacobian = self.getJacobian() 

    def getB(self, thetaArr: np.ndarray, zetaArr: np.ndarray) -> np.ndarray:
        try:
            JacobianGrid = self.Jacobian.getValue(thetaArr, zetaArr)
        except AttributeError:
            self.run() 
            JacobianGrid = self.Jacobian.getValue(thetaArr, zetaArr) 
        except: 
            print("An error occurred while calculating the Jacobian... ")
        g_thetathetaGrid = self.g_thetatheta.getValue(thetaArr, zetaArr)
        g_thetaphiGrid = self.g_thetaphi.getValue(thetaArr, zetaArr)
        g_phiphiGrid = self.g_phiphi.getValue(thetaArr, zetaArr)
        B2Grid = (
            np.power(JacobianGrid, 2) *
            (g_phiphiGrid + 2*self.iota*g_thetaphiGrid + self.iota*self.iota*g_thetaphiGrid)
        )
        return np.power(B2Grid, 0.5)

    # TODO
    def plotB(self, ntheta: int=360, nzeta: int=360, ax=None, fig=None, onePeriod: bool=True, **kwargs):
        from matplotlib import cm
        import matplotlib.pyplot as plt 
        thetaArr = np.linspace(0, 2*np.pi, ntheta)
        thetaValue =  np.linspace(0, 2*np.pi, 3)
        if onePeriod:
            zetaArr = np.linspace(0, 2*np.pi/self.nfp, nzeta)
            zetaValue =  np.linspace(0, 2*np.pi/self.nfp, 3)
        else:
            zetaArr = np.linspace(0, 2*np.pi, nzeta) 
            zetaValue =  np.linspace(0, 2*np.pi, 3)
        if ax is None: 
            fig, ax = plt.subplots() 
        plt.sca(ax) 
        thetaGrid, zetaGrid = np.meshgrid(thetaArr, zetaArr) 
        try: 
            JacobianGrid = self.Jacobian.getValue(thetaGrid, zetaGrid)
        except AttributeError:
            self.run()
            JacobianGrid = self.Jacobian.getValue(thetaGrid, zetaGrid)
        g_thetathetaGrid = self.g_thetatheta.getValue(thetaGrid, zetaGrid)
        g_thetazetaGrid = self.g_thetaphi.getValue(thetaGrid, zetaGrid)
        g_zetazetaGrid = self.g_phiphi.getValue(thetaGrid, zetaGrid)
        B2Grid = np.power(JacobianGrid, 2) * (g_zetazetaGrid + 2*self.iota*g_thetazetaGrid + self.iota*self.iota*g_thetathetaGrid)
        ctrig = ax.contourf(zetaGrid, thetaGrid, np.power(B2Grid,1/2), cmap=cm.rainbow)
        colorbar = fig.colorbar(ctrig)
        colorbar.ax.tick_params(labelsize=18)
        # ax.contour(zetaGrid, thetaGrid, np.power(B2Grid,1/2), cmap=cm.rainbow)
        if onePeriod and self.nfp!=1:
            ax.set_xlabel("$"+str(self.nfp)+r"\varphi$", fontsize=18)
        else:
            ax.set_xlabel(r"$\varphi$", fontsize=18)
        ax.set_ylabel(r"$\theta$", fontsize=18)
        ax.set_xticks(zetaValue)
        ax.set_xticklabels(["$0$", r"$\pi$", r"$2\pi$"], fontsize=18) 
        ax.set_yticks(thetaValue)
        ax.set_yticklabels(["$0$", r"$\pi$", r"$2\pi$"], fontsize=18)
        return

    def getJacobian(self) -> ToroidalField:
        matrixCoef = self.getMatrixCoef() 
        vectorB = self.getVectorB()
        vectorJ = solve(matrixCoef, vectorB)
        reArr = np.zeros(self.ntor+1+self.mpol*(2*self.ntor+1))
        imArr = np.zeros(self.ntor+1+self.mpol*(2*self.ntor+1))
        reArr[0] = self.aveJacobian
        _field = ToroidalField(
            nfp = self.nfp, 
            mpol = self.mpol, 
            ntor = self.ntor, 
            reArr = reArr,
            imArr = imArr
        )
        for i in range(self.ntor+self.mpol*(2*self.ntor+1)):
            m, n, label = self.indexMap(i+1)
            if label == "re":
                _field.setRe(m, n, vectorJ[i])
            elif label == "im":
                _field.setIm(m, n, vectorJ[i])
        return _field

    def indexMap(self, index: int) -> Tuple:
        assert 1 <= index <= 2*self.ntor+2*self.mpol*(2*self.ntor+1)
        if index <= self.ntor+self.mpol*(2*self.ntor+1):
            if 1 <= index <= self.ntor:
                m = 0
                n = index
            elif self.ntor+1 <= index <= self.ntor+self.mpol*(2*self.ntor+1):
                m = (index-self.ntor-1) // (2*self.ntor+1) + 1
                n = (index-self.ntor-1) % (2*self.ntor+1) - self.ntor 
            label = "re"
        elif self.ntor+self.mpol*(2*self.ntor+1)+1 <= index <= 2*self.ntor+2*self.mpol*(2*self.ntor+1):
            if self.ntor+self.mpol*(2*self.ntor+1)+1 <= index <= 2*self.ntor+self.mpol*(2*self.ntor+1):
                m = 0
                n = index - self.ntor - self.mpol*(2*self.ntor+1)
            elif 2*self.ntor+self.mpol*(2*self.ntor+1)+1 <= index <= 2*self.ntor+2*self.mpol*(2*self.ntor+1):
                m = (index-2*self.ntor-self.mpol*(2*self.ntor+1)-1) // (2*self.ntor+1) + 1
                n = (index-2*self.ntor-self.mpol*(2*self.ntor+1)-1) % (2*self.ntor+1) - self.ntor 
            label = "im"
        return m, n, label

    def getVectorB(self) -> np.ndarray:
        vectorB = np.zeros(2*self.ntor+2*self.mpol*(2*self.ntor+1))
        for i in range(2*self.ntor+2*self.mpol*(2*self.ntor+1)): 
            m, n, label = self.indexMap(i+1)
            if label == "re":
                vectorB[i] = self.getRe_CoefMN(m,n,0,0)
            elif label == "im":
                vectorB[i] = self.getIm_CoefMN(m,n,0,0)
        vectorB *= (-self.aveJacobian)
        return vectorB

    def getMatrixCoef(self) -> np.ndarray:
        matrixCoef = np.zeros([2*self.ntor+2*self.mpol*(2*self.ntor+1), 2*self.ntor+2*self.mpol*(2*self.ntor+1)])
        for i in range(2*self.ntor+2*self.mpol*(2*self.ntor+1)):
            m, n, equationLabel = self.indexMap(i+1) 
            for j in range(2*self.ntor+2*self.mpol*(2*self.ntor+1)):
                _m, _n, variableLabel = self.indexMap(j+1) 
                if equationLabel == "re":
                    if variableLabel == "re":
                        matrixCoef[i,j] = self.getRe_CoefMN(m,n,_m,_n) + self.getRe_CoefMN(m,n,-_m,-_n)
                    elif variableLabel == "im":
                        matrixCoef[i,j] = - self.getIm_CoefMN(m,n,_m,_n) + self.getIm_CoefMN(m,n,-_m,-_n)
                elif equationLabel == "im":
                    if variableLabel == "re":
                        matrixCoef[i,j] = self.getIm_CoefMN(m,n,_m,_n) + self.getIm_CoefMN(m,n,-_m,-_n)
                    elif variableLabel == "im":
                        matrixCoef[i,j] = self.getRe_CoefMN(m,n,_m,_n) - self.getRe_CoefMN(m,n,-_m,-_n)
        return matrixCoef

    def getRe_CoefMN(self, m: int, n : int, _m : int, _n: int) -> float:
        return (
            self.D.getRe(m-_m,n-_n) 
            - _m*self.P.getIm(m-_m,n-_n) 
            - _n*self.nfp*self.Q.getIm(m-_m,n-_n)
        )

    def getIm_CoefMN(self, m: int, n : int, _m : int, _n: int) -> float:
        return (
            self.D.getIm(m-_m,n-_n) 
            + _m*self.P.getRe(m-_m,n-_n)
            + _n*self.nfp*self.Q.getRe(m-_m,n-_n)
        )


if __name__ == "__main__":
    pass
