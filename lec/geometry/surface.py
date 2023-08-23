#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# surface.py


from typing import List
from ..toroidalField import ToroidalField
from ..toroidalField import derivatePol, derivateTor 
from ..toroidalField import multiply


class Surface:

    def __init__(self, rField: ToroidalField, zField: ToroidalField) -> None:
        self.r = rField
        self.z = zField
        self.dRdTheta = derivatePol(rField)
        self.dRdPhi = derivateTor(rField)
        self.dZdTheta = derivatePol(zField)
        self.dZdPhi = derivateTor(zField)

    @property
    def mertic(self) -> List[List[ToroidalField]]:
        g_thetatheta = multiply(self.dRdTheta, self.dRdTheta) + multiply(self.dZdTheta, self.dZdTheta)
        g_thetaphi = multiply(self.dRdTheta, self.dRdPhi) + multiply(self.dZdTheta, self.dZdPhi)
        g_phiphi = multiply(self.dRdPhi, self.dRdPhi) + multiply(self.r, self.r) + multiply(self.dZdPhi, self.dZdPhi)
        return [[g_thetatheta, g_thetaphi], [g_thetaphi, g_phiphi]]



if __name__ == "__main__":
    pass
