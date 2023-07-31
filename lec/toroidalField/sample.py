#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# sample.py


import numpy as np
from scipy import fft
from .field import ToroidalField


def fftToroidalField(sampleValue: np.ndarray, nfp: int=1) -> ToroidalField:
    """
    ### Get a toroidal field by fft. 
        `sampleValue.shape = (numsTheta, numsVarphi)` 
    Args:
        sampleValue: the samples. 
        nfp: the number of field periods. 
    Returns:
        (class)ToroidalField
    """
    mlen, nlen = sampleValue.shape
    mpol, ntor = (mlen-1)//2, (nlen-1)//2
    f = fft.fftshift(fft.fft2(sampleValue)/sampleValue.size)
    freal, fimag = f[:,:].real.flatten(), f[:,:].imag.flatten()
    return ToroidalField(
        nfp=nfp, mpol=mpol, ntor=ntor,
        reArr=freal[nlen*mpol+ntor :], imArr=fimag[nlen*mpol+ntor :]
    )


def fftToroidalField_toroidalReversed(sampleValue: np.ndarray, nfp: int=1) -> ToroidalField:
    return fftToroidalField(np.flip(sampleValue, 1), nfp)
    


if __name__ == "__main__":
    pass
