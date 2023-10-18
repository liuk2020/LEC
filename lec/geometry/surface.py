#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# surface.py


import numpy as np
from typing import List
from ..toroidalField import ToroidalField
from ..toroidalField import derivatePol, derivateTor 


class Surface:

    def __init__(self, rField: ToroidalField, zField: ToroidalField) -> None:
        self.r = rField
        self.z = zField

    @property
    def dRdTheta(self) -> ToroidalField:
        return derivatePol(self.r)

    @property
    def dRdPhi(self) -> ToroidalField:
        return derivateTor(self.r)

    @property
    def dZdTheta(self) -> ToroidalField:
        return derivatePol(self.z)

    @property
    def dZdPhi(self) -> ToroidalField:
        return derivateTor(self.z)

    @property
    def mertic(self) -> List[List[ToroidalField]]:
        g_thetatheta = self.dRdTheta*self.dRdTheta + self.dZdTheta*self.dZdTheta
        g_thetaphi = self.dRdTheta*self.dRdPhi + self.dZdTheta*self.dZdPhi
        g_phiphi = self.dRdPhi*self.dRdPhi + self.r*self.r + self.dZdPhi*self.dZdPhi
        return [[g_thetatheta, g_thetaphi], [g_thetaphi, g_phiphi]]

    # fileio ##################################################################
    # TODO: read surface form booz_xform
    @classmethod
    def readBoozXform(cls, boozFile: str, surfaceIndex: int=-1):
        """
        This function cannot be trusted... 
        """
        import booz_xform as bx
        b = bx.Booz_xform()
        try:
            b.read_boozmn(boozFile)
        except:
            print("Cannot open " + boozFile + "...")
        nfp = b.nfp
        mpol = b.mpol-1
        ntor = b.ntor
        rbc = b.rmnc_b[surfaceIndex, :]
        zbs = b.zmns_b[surfaceIndex, :]
        print(mpol)
        print(ntor)
        print((2*ntor+1)*mpol+ntor+1)
        print(rbc.size)
        try:
            rbs = b.rmns_b[surfaceIndex, :]
            zbc = b.zmnc_b[surfaceIndex, :]
        except:
            rbs = np.zeros_like(rbc)
            zbc = np.zeros_like(rbc)
        rField = ToroidalField(
            nfp = nfp, mpol = mpol, ntor = ntor,
            reArr = rbc, imArr = -rbs
        )
        zField = ToroidalField(
            nfp = nfp, mpol = mpol, ntor = ntor,
            reArr = zbc, imArr = -zbs
        )
        return cls(
            rField, zField
        )

    @classmethod
    def readVMECInput(cls, vmecFile: str):
        import f90nml
        try:
            inData = f90nml.read(vmecFile)["indata"]
        except:
            raise FileNotFoundError(
                "Cannot open " + vmecFile + "..."
            )
        nfp = inData["nfp"]
        mpol = inData["mpol"]-1
        ntor = inData["ntor"]
        _rField = ToroidalField(
            nfp=nfp, mpol=mpol, ntor=ntor,
            reArr=np.zeros((2*ntor+1)*mpol+ntor+1), imArr=np.zeros((2*ntor+1)*mpol+ntor+1)
        )
        _zField = ToroidalField(
            nfp=nfp, mpol=mpol, ntor=ntor,
            reArr=np.zeros((2*ntor+1)*mpol+ntor+1), imArr=np.zeros((2*ntor+1)*mpol+ntor+1)
        )
        arrRbc = np.array(inData["rbc"]) 
        arrZbs = np.array(inData["zbs"]) 
        arrRbc[arrRbc==None] = 0
        arrZbs[arrZbs==None] = 0
        try:
            arrRbs = np.array(inData["rbs"]) 
            arrZbc = np.array(inData["zbc"]) 
            arrRbs[arrRbs==None] = 0
            arrZbc[arrZbc==None] = 0
        except KeyError:
            arrRbs = np.zeros_like(arrRbc)
            arrZbc = np.zeros_like(arrRbc)
        nmin, mmin = inData.start_index["rbc"]
        mlen, nlen = np.shape(inData["rbc"])
        for i in range(mlen):
            for j in range(nlen):
                m = i + mmin
                n = j + nmin
                if m==0 and n==0:
                    _rField.setRe(m, n, value=arrRbc[i,j])
                    _rField.setIm(m, n, value=-arrRbs[i,j])
                    _zField.setRe(m, n, value=arrZbc[i,j])
                    _zField.setIm(m, n, value=-arrZbs[i,j])
                else:
                    _rField.setRe(m, n, value=arrRbc[i,j]/2)
                    _rField.setIm(m, n, value=-arrRbs[i,j]/2)
                    _zField.setRe(m, n, value=arrZbc[i,j]/2)
                    _zField.setIm(m, n, value=-arrZbs[i,j]/2)
        return cls(
            _rField, _zField
        )

    @classmethod
    def readVMECOutput(cls, vmecFile: str, surfaceIndex: int=-1):
        import xarray
        try:
            inData = xarray.open_dataset(vmecFile)
        except:
            raise FileNotFoundError(
                "Cannot open " + vmecFile + "..."
            )
        nfp = int(inData["nfp"].values)
        mpol = int(inData["mpol"].values)-1
        ntor = int(inData["ntor"].values)
        rbc = inData["rmnc"].values[surfaceIndex,:]
        zbs = inData["zmns"].values[surfaceIndex,:]
        try:
            rbs = inData["rmns"].values[surfaceIndex,:]
            zbc = inData["zmnc"].values[surfaceIndex,:]
        except:
            rbs = np.zeros_like(rbc)
            zbc = np.zeros_like(rbc)
        rbc[1:-1] = rbc[1:-1] / 2
        zbs[1:-1] = zbs[1:-1] / 2
        rbs[1:-1] = rbs[1:-1] / 2
        zbs[1:-1] = zbs[1:-1] / 2
        rField = ToroidalField(
            nfp = nfp, mpol = mpol, ntor = ntor,
            reArr = rbc, imArr = -rbs
        )
        zField = ToroidalField(
            nfp = nfp, mpol = mpol, ntor = ntor,
            reArr = zbc, imArr = -zbs
        )
        return cls(
            rField, zField
        )
        

if __name__ == "__main__":
    pass
