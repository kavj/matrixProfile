import numba as nb
from numba.decorators import jit
import scipy
import numpy

@jit(nopython=True,nogil=True,cache=True)
def updateMP(ts,mp,mpI,sLen):
    return


@jit(nopython=True,nogil=True,cache=True)
def bivshift(C,dX,dG,dY,dF):
    return


@jit(nopython=True,nogil=True,cache=True)
def uvshift(S,dX,dF):
    return


@jit(nopython=True,nogil=True,cache=True)
def movmean(T,m):
    return


@jit(nopython=True,nogil=True,cache=True)
def compStats(T,m):
    return

@jit(nopython=True,nogil=True,cache=True)
def initBlock():
    return

@jit(nopython=True,nogil=True,cache=True)
def solveBlock():
    return

def probInst():
    n = x.size-slen+1
    for i in range(n):
    return


