import numpy as np
import scipy
from math import pi, exp, sin, cos, sqrt, tan
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.linalg import expm
import time
import random
import multiprocessing
from Functions import convert_to, a_plus, a_minus


def He(nmax, m, t, U, V, mu, dt):


    R = (nmax+1)**2

    Hamiltonian = np.zeros((R, R))
    psi = np.zeros((R, 2))

    for i in range(R):
        psi[i][:] = convert_to(i, nmax, 2)


    for i in range(0, R):
        psi0 = psi[i][:].copy()

        for j in range(0, 1):
            jj = j+1
            #if j == m-1: jj=0 #отключить в случае m=2

            koef1, psi1 = a_plus(psi0.copy(), j, nmax)
            koef2, psi2 = a_minus(psi1, jj)
            matr_el = koef1 * koef2
            index = np.where((psi==psi2).all(axis=1))[0][0]
            Hamiltonian[i][index] -= t * matr_el

            koef1, psi1 = a_plus(psi0.copy(), jj, nmax)
            koef2, psi2 = a_minus(psi1, j)
            matr_el = koef1 * koef2
            index = np.where((psi==psi2).all(axis=1))[0][0]
            Hamiltonian[i][index] -= t * matr_el

            matr_el = psi0[j]*psi0[j] + psi0[jj]*psi0[jj]
            #matr_el = psi0[j]*psi0[j]
            Hamiltonian[i][i] += U/2 * matr_el

            matr_el = psi0[j]*psi0[jj]
            Hamiltonian[i][i] += V * matr_el

    return expm(-dt*Hamiltonian), np.dot(Hamiltonian, expm(-dt*Hamiltonian))

