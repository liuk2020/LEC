#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# problem.py


import numpy as np
from scipy.linalg import solve
from ..geometry import Surface
from ..toroidalField import ToroidalField
from ..toroidalField import derivatePol, derivateTor


class SurfaceEquilibrium:

    def __init__(self, surf: Surface, iota: float, aveJacobian: float) -> None:
        self.surf = surf
        self.iota = iota
        self.aveJacobian = aveJacobian
        self.nfp = self.surf.r.nfp
        self.mpol = self.surf.r.mpol
        self.ntor = self.surf.r.ntor
        self.P = iota*surf.mertic[0][1] + surf.mertic[1][1]
        self.Q = iota*surf.mertic[0][0] + surf.mertic[0][1]
        self.D = derivatePol(self.P) - derivateTor(self.Q)

    def getJacobian(self) -> ToroidalField:
        pass

    def getVectorB(self) -> np.ndarray:
        vectorB = np.zeros(2*self.ntor+2*self.mpol*(2*self.ntor+1))
        for i in range(2*self.ntor+2*self.mpol*(2*self.ntor+1)): 
            index = i + 1
            if index <= self.ntor+self.mpol*(2*self.ntor+1):
                if 1 <= index <= self.ntor:
                    m = 0
                    n = index
                elif self.ntor+1 <= index <= self.ntor+self.mpol*(2*self.ntor+1):
                    m = (index-self.ntor-1) // (2*self.ntor+1) + 1
                    n = (index-self.ntor-1) % (2*self.ntor+1) - self.ntor 
                vectorB[i] = self.getRe_CoefMN(m, n, 0, 0)
            elif self.ntor+self.mpol*(2*self.ntor+1)+1 <= index <= 2*self.ntor+2*self.mpol*(2*self.ntor+1):
                if self.ntor+self.mpol*(2*self.ntor+1)+1 <= index <= 2*self.ntor+self.mpol*(2*self.ntor+1):
                    m = 0
                    n = index - self.ntor - self.mpol*(2*self.ntor+1)
                elif 2*self.ntor+self.mpol*(2*self.ntor+1)+1 <= index <= 2*self.ntor+2*self.mpol*(2*self.ntor+1):
                    m = (index-2*self.ntor-self.mpol*(2*self.ntor+1)-1) // (2*self.ntor+1) + 1
                    n = (index-2*self.ntor-self.mpol*(2*self.ntor+1)-1) % (2*self.ntor+1) - self.ntor 
                vectorB[i] = self.getIm_CoefMN(m, n, 0, 0)
        vectorB *= self.aveJacobian
        return vectorB

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
