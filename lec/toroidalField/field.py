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

    def setRe(self, m: int=0, n: int=0, value: float=0): 
        assert 0 <= m <= self.mpol and -self.ntor <= n <= self.ntor
        self.reArr[self.indexMap(m, n)] = value 

    def setIm(self, m: int=0, n: int=0, value: float=0): 
        assert 0 <= m <= self.mpol and -self.ntor <= n <= self.ntor
        self.imArr[self.indexMap(m, n)] = value

    # operator overloading ####################################################
    def __add__(self, other):
        assert self.nfp == other.nfp
        assert self.mpol == other.mpol
        assert self.ntor == other.ntor
        return ToroidalField(
            nfp = self.nfp, 
            mpol = self.mpol, 
            ntor = self.ntor,
            reArr = self.reArr + other.reArr, 
            imArr = self.imArr + other.imArr
        )

    def __sub__(self, other):
        assert self.nfp == other.nfp
        assert self.mpol == other.mpol
        assert self.ntor == other.ntor
        return ToroidalField(
            nfp = self.nfp, 
            mpol = self.mpol, 
            ntor = self.ntor,
            reArr = self.reArr - other.reArr, 
            imArr = self.imArr - other.imArr
        )

    def __mul__(self, other):
        if isinstance(other, ToroidalField):
            assert self.nfp == other.nfp
            mpol, ntor = self.mpol, self.ntor
            nums = (2*ntor+1)*mpol+ntor+1
            reArr, imArr = np.zeros(nums), np.zeros(nums)
            for i in range(nums):
                m, n = self.indexReverseMap(i)
                for _m in range(-mpol, mpol+1):
                    for _n in range(-ntor, ntor+1):
                        reArr[i] += (
                            self.getRe(_m,_n)*other.getRe(m-_m,n-_n) - 
                            self.getIm(_m,_n)*other.getIm(m-_m,n-_n)
                        )
                        imArr[i] += (
                            self.getRe(_m,_n)*other.getIm(m-_m,n-_n) + 
                            self.getIm(_m,_n)*other.getRe(m-_m,n-_n)
                        )
            return ToroidalField(
                nfp = self.nfp, 
                mpol = mpol, 
                ntor = ntor,
                reArr = reArr,
                imArr = imArr
            )
        else:
            return ToroidalField(
                nfp = self.nfp, 
                mpol = self.mpol, 
                ntor = self.ntor, 
                reArr = other * self.reArr,
                imArr = other * self.imArr
            )

    def __eq__(self, other) -> bool:
        try:
            assert self.nfp == other.nfp
            assert self.mpol == other.mpol
            assert self.ntor == other.ntor
            assert self.reArr.any() == other.reArr.any()
            assert self.imArr.any() == other.imArr.any()
            return True
        except:
            return False


if __name__ == "__main__":
    pass
    