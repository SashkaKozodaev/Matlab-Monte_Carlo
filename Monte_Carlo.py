import numpy as np
import scipy
from math import pi, exp, sin, cos, sqrt, tan
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.linalg import expm
import time
import random
import multiprocessing
from progress.bar import IncrementalBar
from bound_cond import bound_cond
from black_coordinates import black_coordinates
from weight import He

#задаём параметры системы
nmax = 5
m = 6
R = (nmax+1)**2
t = 1
U = 4
V = 3
M = 10
M2 = 2*M
mu = 0
beta = 5
dt = beta/(2*M)
ex, H_e = He(nmax, m, t, U, V, mu, dt)

#Задание начальной населённости
n = np.zeros((2*M, m))


n[:, 1] = 1
n[:, 2] = 2
n[:, 3] = 2
n[:, 4] = 1
#Рассчёт начальной энергии

E_ini = 0
for i in range(M2):
    for j in range(m):
        sq_u, sq_l, sq_d, sq_r = black_coordinates(m, M, i, j)
        c1 = sq_u.copy()[2]
        c2 = sq_u.copy()[3]
        c3 = sq_u.copy()[0]
        c4 = sq_u.copy()[1]
        c7 = sq_l.copy()[1]
        c8 = sq_d.copy()[3]
        c9 = sq_d.copy()[0]
        c10 = sq_d.copy()[1]

        n1 = int(n.copy()[(2*M-int(c1[1]))%M2, int(c1[0])%m])
        n2 = int(n.copy()[(2*M-int(c2[1]))%M2, int(c2[0])%m])
        n3 = int(n.copy()[(2*M-int(c3[1]))%M2, int(c3[0])%m])
        n4 = int(n.copy()[(2*M-int(c4[1]))%M2, int(c4[0])%m])
        n7 = int(n.copy()[(2*M-int(c7[1]))%M2, int(c7[0])%m])
        n8 = int(n.copy()[(2*M-int(c8[1]))%M2, int(c8[0])%m])
        n9 = int(n.copy()[(2*M-int(c9[1]))%M2, int(c9[0])%m])
        n10 = int(n.copy()[(2*M-int(c10[1]))%M2, int(c10[0])%m])


        if i%2 == 0:
            if j%2 == 0:
                E_ini += H_e.copy()[n2 + n1*(nmax+1)][n4 + n3*(nmax+1)]/ex.copy()[n2 + n1*(nmax+1)][n4 + n3*(nmax+1)]
        else:
            if j%2 == 1:
                E_ini += H_e.copy()[n8 + n7*(nmax+1)][n10 + n9*(nmax+1)]/ex.copy()[n8 + n7*(nmax+1)][n10 + n9*(nmax+1)]

E_ini = E_ini/M


x_ax_pogr_mas = []
y_ax_pogr_mas = []
delta_mas = []
dE_mas = []
delta34_mas = []
N_steps = 10**4

N_term = N_steps//20
Error = 0.01

E_mas = np.zeros(N_steps)
E_av_mas = np.zeros(N_steps)
E_mas[0] = E_ini
count = 1
finish_flag = 0
finish_step = 0

for i in range(1, N_steps, 1):

    #1 шаг, случайный выбор одного из белых квадратов
    def step1(m, M):
        white_sq_1 = random.randint(0, m/2-1) * 2
        white_sq_2 = random.randint(0, 2*M/2-1) * 2 #!!!!!!!!

        return [white_sq_1, white_sq_2]

    white_sq_1, white_sq_2 = step1(m, M)

    #2 шаг, Определение координат всех его соседей черных квадратов с учетом периодических граничных условий
    #функция выдаёт 4 вектора для квадратов u, l, d, r
    #координаты квадратов записаны в порядке d_l, d_r, u_l, u_r
    sq_u, sq_l, sq_d, sq_r = black_coordinates(m, M, white_sq_1, white_sq_2)

    c1 = sq_u.copy()[2]    #    1___2
    c2 = sq_u.copy()[3]    #    | ч |
    c3 = sq_u.copy()[0]    #5___3___4___11
    c4 = sq_u.copy()[1]    #| ч | б | ч |
    c5 = sq_l.copy()[2]    #6___7___8___12
    c6 = sq_l.copy()[0]    #    | ч |
    c7 = sq_l.copy()[1]    #    9___10
    c8 = sq_d.copy()[3]
    c9 = sq_d.copy()[0]
    c10 = sq_d.copy()[1]
    c11 = sq_r.copy()[3]
    c12 = sq_r.copy()[1]

    #3 шаг, выбор процедуры
    ksi = random.random()
    direct_proc = False
    if ksi < 0.5: direct_proc = True


    #4 шаг, проверка возможности процедуры
    if direct_proc:
        opportunity_proc = False
        if (n[int((2*M-sq_u.copy()[0][1]))%M2, int(sq_u[0][0])]>0 and n[(2*M-int(sq_d.copy()[2][1]))%M2, int(sq_d[2][0])]>0 and \
            n[int(2*M-sq_r.copy()[0][1])%M2, int(sq_r[0][0])]<nmax and n[(2*M-int(sq_r.copy()[2][1]))%M2, int(sq_r[2][0])]<nmax): opportunity_proc = True

    else:
        opportunity_proc = False
        if (n[(2*M-int(sq_u.copy()[0][1]))%M2, int(sq_u[0][0])]<nmax and n[(2*M-int(sq_d.copy()[2][1]))%M2, int(sq_d[2][0])]<nmax and \
            n[(2*M-int(sq_r.copy()[0][1]))%M2, int(sq_r[0][0])]>0 and n[(2*M-int(sq_r.copy()[2][1]))%M2, int(sq_r[2][0])]>0): opportunity_proc = True


    #5 шаг, Если процедура возможна, изменение траектории происходит в соответствии с алгоритмом Метрополиса

    dE = 0

    if opportunity_proc:

        n1 = int(n.copy()[(2*M-int(c1[1]))%M2, int(c1[0])])    #    1___2
        n2 = int(n.copy()[(2*M-int(c2[1]))%M2, int(c2[0])])    #    | ч |
        n3 = int(n.copy()[(2*M-int(c3[1]))%M2, int(c3[0])])    #5___3___4___11
        n4 = int(n.copy()[(2*M-int(c4[1]))%M2, int(c4[0])])    #| ч | б | ч |
        n5 = int(n.copy()[(2*M-int(c5[1]))%M2, int(c5[0])])    #6___7___8___12
        n6 = int(n.copy()[(2*M-int(c6[1]))%M2, int(c6[0])])    #    | ч |
        n7 = int(n.copy()[(2*M-int(c7[1]))%M2, int(c7[0])])    #    9___10
        n8 = int(n.copy()[(2*M-int(c8[1]))%M2, int(c8[0])])
        n9 = int(n.copy()[(2*M-int(c9[1]))%M2, int(c9[0])])
        n10 = int(n.copy()[(2*M-int(c10[1]))%M2, int(c10[0])])
        n11 = int(n.copy()[(2*M-int(c11[1]))%M2, int(c11[0])])
        n12 = int(n.copy()[(2*M-int(c12[1]))%M2, int(c12[0])])

        if direct_proc:
            w_old_1 = ex.copy()[n2 + n1*(nmax+1)][n4 + n3*(nmax+1)]
            w_old_2 = ex.copy()[n3 + n5*(nmax+1)][n7 + n6*(nmax+1)]
            w_old_3 = ex.copy()[n8 + n7*(nmax+1)][n10 + n9*(nmax+1)]
            w_old_4 = ex.copy()[n11 + n4*(nmax+1)][n12 + n8*(nmax+1)]
            w_old = w_old_1*w_old_2*w_old_3*w_old_4

            E_old_1 = H_e.copy()[n2 + n1*(nmax+1)][n4 + n3*(nmax+1)]
            E_old_2 = H_e.copy()[n3 + n5*(nmax+1)][n7 + n6*(nmax+1)]
            E_old_3 = H_e.copy()[n8 + n7*(nmax+1)][n10 + n9*(nmax+1)]
            E_old_4 = H_e.copy()[n11 + n4*(nmax+1)][n12 + n8*(nmax+1)]

            n3 = n3 - 1
            n4 = n4 + 1
            n7 = n7 - 1
            n8 = n8 + 1


            w_new_1 = ex.copy()[n2 + n1*(nmax+1)][n4 + n3*(nmax+1)]
            w_new_2 = ex.copy()[n3 + n5*(nmax+1)][n7 + n6*(nmax+1)]
            w_new_3 = ex.copy()[n8 + n7*(nmax+1)][n10 + n9*(nmax+1)]
            w_new_4 = ex.copy()[n11 + n4*(nmax+1)][n12 + n8*(nmax+1)]
            w_new = w_new_1*w_new_2*w_new_3*w_new_4

            E_new_1 = H_e.copy()[n2 + n1*(nmax+1)][n4 + n3*(nmax+1)]
            E_new_2 = H_e.copy()[n3 + n5*(nmax+1)][n7 + n6*(nmax+1)]
            E_new_3 = H_e.copy()[n8 + n7*(nmax+1)][n10 + n9*(nmax+1)]
            E_new_4 = H_e.copy()[n11 + n4*(nmax+1)][n12 + n8*(nmax+1)]

            dE = (E_new_1/w_new_1 + E_new_2/w_new_2 + E_new_3/w_new_3 + E_new_4/w_new_4 - E_old_1/w_old_1 - E_old_2/w_old_2 - E_old_3/w_old_3 - E_old_4/w_old_4)/M
        else:

            w_old_1 = ex.copy()[n2 + n1*(nmax+1)][n4 + n3*(nmax+1)]
            w_old_2 = ex.copy()[n3 + n5*(nmax+1)][n7 + n6*(nmax+1)]
            w_old_3 = ex.copy()[n8 + n7*(nmax+1)][n10 + n9*(nmax+1)]
            w_old_4 = ex.copy()[n11 + n4*(nmax+1)][n12 + n8*(nmax+1)]
            w_old = w_old_1*w_old_2*w_old_3*w_old_4

            E_old_1 = H_e.copy()[n2 + n1*(nmax+1)][n4 + n3*(nmax+1)]
            E_old_2 = H_e.copy()[n3 + n5*(nmax+1)][n7 + n6*(nmax+1)]
            E_old_3 = H_e.copy()[n8 + n7*(nmax+1)][n10 + n9*(nmax+1)]
            E_old_4 = H_e.copy()[n11 + n4*(nmax+1)][n12 + n8*(nmax+1)]


            n3 = n3 + 1
            n4 = n4 - 1
            n7 = n7 + 1
            n8 = n8 - 1


            w_new_1 = ex.copy()[n2 + n1*(nmax+1)][n4 + n3*(nmax+1)]
            w_new_2 = ex.copy()[n3 + n5*(nmax+1)][n7 + n6*(nmax+1)]
            w_new_3 = ex.copy()[n8 + n7*(nmax+1)][n10 + n9*(nmax+1)]
            w_new_4 = ex.copy()[n11 + n4*(nmax+1)][n12 + n8*(nmax+1)]
            w_new = w_new_1*w_new_2*w_new_3*w_new_4

            E_new_1 = H_e.copy()[n2 + n1*(nmax+1)][n4 + n3*(nmax+1)]
            E_new_2 = H_e.copy()[n3 + n5*(nmax+1)][n7 + n6*(nmax+1)]
            E_new_3 = H_e.copy()[n8 + n7*(nmax+1)][n10 + n9*(nmax+1)]
            E_new_4 = H_e.copy()[n11 + n4*(nmax+1)][n12 + n8*(nmax+1)]


            dE = (E_new_1/w_new_1 + E_new_2/w_new_2 + E_new_3/w_new_3 + E_new_4/w_new_4 - E_old_1/w_old_1 - E_old_2/w_old_2 - E_old_3/w_old_3 - E_old_4/w_old_4)/M

        R = w_new/w_old

        accept = False
        if R>=1:
            accept = True
        else:
            csi = random.random()
            if csi < R: accept = True
        if accept==False:
            dE = 0
        else:
            if direct_proc:
                n[(2*M-int(c3[1]))%M2, int(c3[0])] -= 1
                n[(2*M-int(c4[1]))%M2, int(c4[0])] += 1
                n[(2*M-int(c7[1]))%M2, int(c7[0])] -= 1
                n[(2*M-int(c8[1]))%M2, int(c8[0])] += 1
            else:
                n[(2*M-int(c3[1]))%M2, int(c3[0])] += 1
                n[(2*M-int(c4[1]))%M2, int(c4[0])] -= 1
                n[(2*M-int(c7[1]))%M2, int(c7[0])] += 1
                n[(2*M-int(c8[1]))%M2, int(c8[0])] -= 1
    dE_mas.append(dE)
    E_mas[i] = E_mas[i-1]+dE

    N_term = 0 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if i>N_term:
        E_av_mas[i] = 1/(i-N_term)*((i-N_term-1)*E_av_mas[i-1] + E_mas[i])


    count+=1
    if count % 10**5 ==0: print(count)

    #6 шаг, оценка погрешности
    #if (i>N_term+500 and i%10**4==0):   #@@@1
    if i%10**3==0:                       #@@@1
        C = 1
        k  = 1
        while C>0.1:
            #X = 1/(i-N_term)*np.sum(E_mas[N_term:]) #@@@2
            X = 1/i*np.sum(E_mas[:i])                #@@@2
            #Y = 1/(i-N_term)*np.sum(E_mas[N_term:]**2)  #@@@3
            Y = 1/i*np.sum(E_mas[:i]**2)                 #@@@3
            Z = 0
            #for j in range(i - N_term - k):  #@@@4
            for j in range(i - k):   #@@@4
                #Z += E_mas[N_term+j]*E_mas[N_term+j+k]  #@@@4
                Z += E_mas[j]*E_mas[j+k]                 #@@@4

            #Z = Z/(i-N_term -k) #@@@5
            Z = Z/(i - k)         #@@@5
            C = (Z - X**2)/(Y - X**2)

            k +=1

        #E_autocorr = E_mas.copy()[N_term : i:k] #@@@5
        E_autocorr = E_mas.copy()[0 : i : k]     #@@@5

        delta = np.std(E_autocorr)
        #mean = np.mean(E_mas.copy()[N_term:i]) #@@@6
        #mean = np.mean(E_mas.copy()[:i])  #@@@6
        #print("delta = ", delta)
        if delta/E_av_mas[i] < Error and finish_flag == 0 :
            #print("!!!", " Step = ",i ," delta = ", delta, " E_av = ", E_av_mas[i], "delta/E_av = ", delta/E_av_mas[i])
            #print("!!!", " Step = ", i)
            finish_flag = 1
            finish_step = i

        delta_mas.append(delta)

        delta34 = (max(E_mas.copy()[i*3//4:i]) - min(E_mas.copy()[i*3//4:i]))/2
        delta34_mas.append(delta34)

        #print("delta = ", delta, " delta/E = ",delta/E_av_mas[i] )
        x_ax_pogr_mas.append(i)
        y_ax_pogr_mas.append(delta/E_av_mas[i])


print("U = ", U, " ,Finish step = ", finish_step)
plt.plot(np.linspace(0, N_steps, N_steps), E_mas)


plt.show()

plt.plot(x_ax_pogr_mas, y_ax_pogr_mas)
#plt.plot(x_ax_pogr_mas, np.ones(len(x_ax_pogr_mas))*0.05)
plt.plot(x_ax_pogr_mas, np.ones(len(x_ax_pogr_mas))*Error)
plt.xlabel('Monte-Carlo steps')
plt.ylabel('Relative error')

plt.show()

