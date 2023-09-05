#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# operate.py


import numpy as np
from .field import ToroidalField


def conjugate(field: ToroidalField) -> ToroidalField: 
    return ToroidalField(
        nfp = field.nfp, 
        mpol = field.mpol, 
        ntor = field.ntor, 
        reArr = field.reArr,
        imArr = - field.imArr
    )


def add(fieldA: ToroidalField, fieldB: ToroidalField) -> ToroidalField: 
    assert fieldA.nfp == fieldB.nfp
    assert fieldA.mpol == fieldB.mpol
    assert fieldA.ntor == fieldB.ntor
    return ToroidalField(
        nfp = fieldA.nfp, 
        mpol = fieldA.mpol, 
        ntor = fieldA.ntor,
        reArr = fieldA.reArr + fieldB.reArr, 
        imArr = fieldA.imArr + fieldB.imArr
    )


def multiply(c: float, field: ToroidalField) -> ToroidalField:
    return ToroidalField(
        nfp = field.nfp, 
        mpol = field.mpol, 
        ntor = field.ntor, 
        reArr = c * field.reArr,
        imArr = c * field.imArr
    )


def multiply_im(field: ToroidalField) -> ToroidalField: 
    return ToroidalField(
        nfp = field.nfp, 
        mpol = field.mpol, 
        ntor = field.ntor, 
        reArr = - field.imArr,
        imArr = field.reArr
    )    


def conv(fieldA: ToroidalField, fieldB: ToroidalField) -> ToroidalField: 
    assert fieldA.nfp == fieldB.nfp
    mpol, ntor = fieldA.mpol, fieldA.ntor
    nums = (2*ntor+1)*mpol+ntor+1
    reArr, imArr = np.zeros(nums), np.zeros(nums)
    for i in range(nums):
        m, n = fieldA.indexReverseMap(i)
        for _m in range(-mpol, mpol+1):
            for _n in range(-ntor, ntor+1):
                reArr[i] += (
                    fieldA.getRe(_m,_n)*fieldB.getRe(m-_m,n-_n) - fieldA.getIm(_m,_n)*fieldB.getIm(m-_m,n-_n)
                )
                imArr[i] += (
                    fieldA.getRe(_m,_n)*fieldB.getIm(m-_m,n-_n) + fieldA.getIm(_m,_n)*fieldB.getRe(m-_m,n-_n)
                )
    return ToroidalField(
        nfp = fieldA.nfp, 
        mpol = mpol, 
        ntor = ntor,
        reArr = reArr,
        imArr = imArr
    )


if __name__ == "__main__": 
    pass
