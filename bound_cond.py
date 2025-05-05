import numpy as np
import scipy
from math import pi, exp, sin, cos, sqrt, tan
import matplotlib.pyplot as plt
from scipy import integrate
import time
import random
import multiprocessing
from progress.bar import IncrementalBar

def bound_cond(b, m, M):

    if b[0] == -1:   b[0] = m-1
    if b[0] == m:    b[0] = 0

    if b[1] == -1:   b[1] = 2*M-1    #!!!!!!!
    if b[1] == 2*M:    b[1] = 0  #!!!!!!!

    return b
