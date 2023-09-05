#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# problem.py


from scipy.linalg import solve
from ..geometry import Surface
from ..toroidalField import ToroidalField
from ..toroidalField import derivatePol, derivateTor


class SurfaceEquilibrium:

    def __init__(self, surf: Surface, iota: float) -> None:
        self.surf = surf
        self.iota = iota
        self.nfp = self.surf.r.nfp
        self.P = iota*surf.mertic[0][1] + surf.mertic[1][1]
        self.Q = iota*surf.mertic[0][0] + surf.mertic[0][1]
        self.D = derivatePol(self.P) - derivateTor(self.Q)

    def getJacobian(self) -> ToroidalField:
        pass

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
