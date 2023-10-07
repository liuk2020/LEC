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
    @classmethod
    def readBoozXform(cls, boozFile: str, surfaceIndex: int=-1):
        import booz_xform as bx
        b = bx.Booz_xform()
        try:
            b.read_boozmn(boozFile)
        except:
            print("Cannot open " + boozFile + "...")
        nfp = b.nfp
        mpol = b.mpol
        ntor = b.ntor
        rbc = b.rmnc_b[surfaceIndex, :]
        zbs = b.zmns_b[surfaceIndex, :]
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
            print("Cannot open " + vmecFile + "...")
        nfp = inData["nfp"]
        mpol = inData["mpol"]
        ntor = inData["ntor"]
        rbc = np.array(inData["rbc"])
        zbs = np.array(inData["zbs"])
        try:
            rbs = np.array(inData["rbs"])
            zbc = np.array(inData["zbc"])
        except KeyError:
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
    def readVMECOutput(cls, vmecFile: str, surfaceIndex: int=-1):
        import xarray
        try:
            inData = xarray.open_dataset(vmecFile)
        except:
            print("Cannot open " + vmecFile + "...")
        nfp = int(inData["nfp"].values)
        mpol = int(inData["mpol"].values)
        ntor = int(inData["ntor"].values)
        rbc = inData["rmnc"].values[surfaceIndex,:]
        zbs = inData["zmns"].values[surfaceIndex,:]
        try:
            rbs = inData["rmns"].values[surfaceIndex,:]
            zbc = inData["zmnc"].values[surfaceIndex,:]
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
        

if __name__ == "__main__":
    pass
