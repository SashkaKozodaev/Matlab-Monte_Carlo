import numpy as np
import scipy
from math import pi, exp, sin, cos, sqrt, tan
import matplotlib.pyplot as plt
from scipy import integrate
import time
import multiprocessing
from progress.bar import IncrementalBar

def convert_to(number, nmax, m):
    base = nmax+1
    psi = '0'*m
    digits = '0123456789'
    result = ''
    while number > 0:
        result = digits[number % base] + result
        number //= base
    psi = psi+result
    psi = psi[-m:]
    psi = list(map(int, psi))
    return psi


def a_plus(psi_f, j, nmax):
    psi_f[j] = psi_f[j] + 1
    if psi_f[j]>nmax : psi_f[j] = 0
    koef = sqrt(psi_f[j])

    return koef, psi_f

def a_minus(psi_f, j):

    koef = sqrt(psi_f[j])
    psi_f[j] = psi_f[j]-1
    if psi_f[j]<0:
        koef = 0
        psi_f[j]=0


    return koef, psi_f



