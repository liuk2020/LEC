#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# problem.py


from scipy.linalg import solve
from ..geometry import Surface
from ..toroidalField import ToroidalField


class SurfaceEquilibrium:

    def __init__(self, surf: Surface, iota: float) -> None:
        self.surf = surf
        self.iota = iota

    def getJacobian(self) -> ToroidalField:
        pass


if __name__ == "__main__":
    pass
