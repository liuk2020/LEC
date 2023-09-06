#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# testIndexMap.py


def testX(mpol: int=3, ntor: int=4) -> None:
    
    def Re(m, n):
        return "Re("+str(m)+","+str(n)+")"
    def Im(m, n):
        return "Im("+str(m)+","+str(n)+")"
    
    for i in range(2*ntor+2*mpol*(2*ntor+1)):
        index = i + 1
        if 1 <= index <= ntor:
            print(Re(0, index))
        elif ntor+1 <= index <= ntor+mpol*(2*ntor+1):
            print(Re((index-ntor-1)//(2*ntor+1)+1, (index-ntor-1)%(2*ntor+1)-ntor))
        elif ntor+mpol*(2*ntor+1)+1 <= index <= 2*ntor+mpol*(2*ntor+1):
            print(Im(0, index-ntor-mpol*(2*ntor+1)))
        elif 2*ntor+mpol*(2*ntor+1)+1 <= index <= 2*ntor+2*mpol*(2*ntor+1):
            print(Im((index-2*ntor-mpol*(2*ntor+1)-1)//(2*ntor+1)+1, (index-2*ntor-mpol*(2*ntor+1)-1)%(2*ntor+1)-ntor))


if __name__ == "__main__":
    testX()
