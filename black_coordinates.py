import numpy as np
import scipy
from math import pi, exp, sin, cos, sqrt, tan
import matplotlib.pyplot as plt
from scipy import integrate
import time
import random
import multiprocessing
from progress.bar import IncrementalBar
from bound_cond import bound_cond


def black_coordinates(m, M, index1, index2):


    w = np.array([index1, index2])  #нижний левый угол белого квадрата
    d_l = np.zeros((4, 2))   #матрица с нижними левыми углами чёрных квадратов


    d_l[0][:] = np.array([index1, index2+1])
    d_l[1][:] = np.array([index1-1, index2])
    d_l[2][:] = np.array([index1, index2-1])
    d_l[3][:] = np.array([index1+1, index2])

    for i in range(0, 4):
        d_l[i][:] =  bound_cond(d_l[i][:].copy(), m, M) #проверка граничных условий

    d_r = d_l.copy()    #матрица с нижними правыми углами чёрных квадратов
    d_r[:, 0] += 1
    for i in range(0, 4): d_r[i][:] =  bound_cond(d_r[i][:].copy(), m, M) #проверка граничных условий

    u_l= d_l.copy()    #матрица с верхними левыми углами чёрных квадратов
    u_l[:, 1] += 1
    for i in range(0, 4): u_l[i][:] =  bound_cond(u_l[i][:].copy(), m, M) #проверка граничных условий

    u_r = d_l.copy()    #матрица с верхними правыми углами чёрных квадратов
    u_r[:, :] += 1
    for i in range(0, 4): u_r[i][:] =  bound_cond(u_r[i][:].copy(), m, M) #проверка граничных условий

    sq_u = [d_l[0], d_r[0], u_l[0], u_r[0]]
    sq_l = [d_l[1], d_r[1], u_l[1], u_r[1]]
    sq_d = [d_l[2], d_r[2], u_l[2], u_r[2]]
    sq_r = [d_l[3], d_r[3], u_l[3], u_r[3]]
    #return d_l, d_r, u_l, u_r
    return sq_u, sq_l, sq_d, sq_r


