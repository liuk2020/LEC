#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# problem.py


import numpy as np
from scipy.linalg import solve
from ..geometry import Surface
from ..toroidalField import ToroidalField
from ..toroidalField import derivatePol, derivateTor
from typing import Tuple


class SurfaceEquilibrium:

    def __init__(self, surf: Surface, iota: float, aveJacobian: float=1.0) -> None:
        self.surf = surf
        self.iota = iota
        self.aveJacobian = aveJacobian
        self.nfp = self.surf.r.nfp
        self.mpol = self.surf.r.mpol
        self.ntor = self.surf.r.ntor
        self.P = surf.mertic[0][1]*iota + surf.mertic[1][1]
        self.Q = surf.mertic[0][0]*iota + surf.mertic[0][1]
        self.D = derivatePol(self.P) - derivateTor(self.Q)
    
    def run(self):
        self.Jacobian = self.getJacobian()

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
        JacobianGrid = self.Jacobian.getValue(thetaGrid, zetaGrid)
        metric = self.surf.mertic
        g_thetathetaGrid = metric[0][0].getValue(thetaGrid, zetaGrid)
        g_thetazetaGrid = metric[0][1].getValue(thetaGrid, zetaGrid)
        g_zetazetaGrid = metric[1][1].getValue(thetaGrid, zetaGrid)
        B2Grid = self.iota*self.iota*g_thetathetaGrid + 2*self.iota*g_thetazetaGrid + g_zetazetaGrid / np.power(JacobianGrid,2)
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
        for i in range(self.ntor+self.mpol*(2*self.ntor+1)):
            m, n, label = self.indexMap(i+1)
            if label == "re":
                reArr[self.ntor+(2*self.ntor+1)*(m-1)+(n+self.ntor+1)] = vectorJ[i]
            elif label == "im":
                imArr[self.ntor+(2*self.ntor+1)*(m-1)+(n+self.ntor+1)] = vectorJ[i]
        return ToroidalField(
            nfp = self.nfp, 
            mpol = self.mpol, 
            ntor = self.ntor, 
            reArr = reArr,
            imArr = imArr
        )

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
                        matrixCoef[i,j] = self.getIm_CoefMN(m,n,_m,_n) - self.getIm_CoefMN(m,n,-_m,-_n)
                    elif variableLabel == "im":
                        matrixCoef[i,j] = self.getRe_CoefMN(m,n,_m,_n) + self.getRe_CoefMN(m,n,-_m,-_n)
        return matrixCoef

    def getRe_CoefMN(self, m: int, n : int, _m : int, _n: int) -> float:
        return (
            self.D.getRe(m-_m,n-_n) 
            + _m*self.P.getIm(m-_m,n-_n) 
            + _n*self.nfp*self.Q.getIm(m-_m,n-_n)
        )

    def getIm_CoefMN(self, m: int, n : int, _m : int, _n: int) -> float:
        return (
            self.D.getIm(m-_m,n-_n) 
            - _m*self.P.getRe(m-_m,n-_n)
            - _n*self.nfp*self.Q.getRe(m-_m,n-_n)
        )


if __name__ == "__main__":
    pass
