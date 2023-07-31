#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# field.py


import numpy as np
from typing import Tuple


class ToroidalField:
    r"""
    ## The Fourier representation of the field f defined on the toroidal surface. 
    $$ f(\theta, \varphi) = \sum_{m,n} F_{m,n}\exp(i(m\theta-nN_{fp}\varphi)) $$
    """

    def __init__(self, nfp: int, mpol: int, ntor: int, reArr: np.ndarray, imArr: np.ndarray) -> None:
        """
        ### Initialization with Fourier harmonics. 
        Args:
            nfp: the number of field periods. 
            mpol, ntor: the resolution in the poloidal/toroidal direction. 
            reArr, imArr: the real/imaginary part of the Fourier coefficients. 
        """
        assert reArr.shape == imArr.shape
        assert (2*ntor+1)*mpol+ntor+1 == reArr.size
        self.nfp = nfp
        self.mpol = mpol
        self.ntor = ntor
        self.reArr = reArr
        self.imArr = imArr

    @property
    def xm(self) -> np.ndarray:
        return np.array([
            self.indexReverseMap(i)[0] for i in range(self.mpol*(2*self.ntor+1)+self.ntor+1)
        ])

    @property
    def xn(self) -> np.ndarray:
        return np.array([
            self.indexReverseMap(i)[1] for i in range(self.mpol*(2*self.ntor+1)+self.ntor+1)
        ])

    def indexMap(self, m: int, n: int) -> int:
        assert abs(m) <= self.mpol and abs(n) <= self.ntor
        return self.ntor + (2*self.ntor+1)*(m-1) + (n+self.ntor+1)

    def indexReverseMap(self, index: int) -> Tuple[int]: 
        assert index < (self.mpol*(2*self.ntor+1)+self.ntor+1)
        if index <= self.ntor:
            return 0, index
        else:
            return (index-self.ntor-1)//(2*self.ntor+1)+1, (index-self.ntor-1)%(2*self.ntor+1)-self.ntor

    def getValue(self, thetaArr: np.ndarray, zetaArr: np.ndarray) -> np.ndarray:
        assert type(thetaArr) == type(zetaArr)
        if not isinstance(thetaArr, np.ndarray):
            try:
                thetaArr, zetaArr = np.array(thetaArr), np.array(zetaArr)
            except:
                thetaArr, zetaArr = np.array([thetaArr]), np.array([zetaArr])
        angleMat = (
            np.dot(self.xm.reshape(-1,1), thetaArr.reshape(1,-1)) - 
            self.nfp * np.dot(self.xn.reshape(-1,1), zetaArr.reshape(1,-1))
        )
        valueArr = 2 * (
            np.dot(self.reArr.reshape(1,-1), np.cos(angleMat)) - 
            np.dot(self.imArr.reshape(1,-1), np.sin(angleMat))
        )
        valueArr -= self.reArr[0]
        try:
            m, n = thetaArr.shape
            return valueArr.reshape(m, n)
        except:
            return valueArr

    def getRe(self, m: int=0, n: int=0) -> float: 
        if abs(m) > self.mpol or abs(n) > self.ntor:
            return 0
        else:
            return self.reArr[self.indexMap(m, n)] 

    def getIm(self, m: int, n: int) -> float:
        if abs(m) > self.mpol or abs(n) > self.ntor:
            return 0
        else:
            return self.imArr[self.indexMap(m, n)]


if __name__ == "__main__":
    pass
    