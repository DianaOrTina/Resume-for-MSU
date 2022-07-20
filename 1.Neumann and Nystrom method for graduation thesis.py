import scipy as sp
import scipy.integrate as spi
import numpy as np
import numpy.linalg as LA
import math
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from scipy.interpolate import make_interp_spline
from sympy import *

size = 500
R1 = -30
R2 = 30
ite = 100



def Fx(x):
    return 1 + (24 / 71) / (x ** 2 + 1) + (40 / 71) / (x ** 2 + 4)


b = 1
s = 1  # d'
p_mw = 1
A_w = 1
n_w = 2
N = A_w * math.pi * (A_w + 5 * p_mw ** 2) * (A_w + 8 * p_mw ** 2) / (
        p_mw * (A_w ** 2 + 21 * p_mw ** 2 * A_w + 120 * p_mw ** 4))
# N=1
wj = (R2 - R1) / size


def m(x):
    return (p_mw / math.pi) / (x ** 2 + p_mw ** 2)


def w(x):
    return A_w / (x ** 2 + ((n_w + 1) ** 2) * (p_mw ** 2))


'''



def m(x):
    return math.exp(-2*math.fabs(x))

A=0
B=1
b=1
s=1
N=(2/3)*B+(52/27)*A
wj = (R2 - R1) / size

def Q(x):
    return (1/3)*A*x**2-(16/9)*A*math.fabs(x)+(56/27)*A+B/3

def R(x):
    return A*x**2 + B

def w(x):
    return math.exp(-math.fabs(x))*Q(x)/(1+math.exp(-math.fabs(x))*R(x))

def Fx(x):
    return 1+math.exp(-math.fabs(x))*R(x)


'''


def K(x, y):
    return b * m(x - y) / (1 + (s / b) * w(x))


def f(x):
    return (N * m(x) - s * w(x)) / (b + s * w(x))


def Multiply(K, ff):  # получить резурьтат (матрица K * вектор f = вектор ff)
    c = [0 for i in range(size)]
    for i in range(size):
        for j in range(size):
            c[i] += K[i][j] * ff[j]
    return c


def MatrixMultiply(K1, K2):
    K3 = [[0] * size for i in range(size)]
    for i in range(size):
        for j in range(size):
            for k in range(size):
                K3[i][j] += K1[i][k] * K2[k][j];
    return K3


def InverseMatrix(K1):
    k1 = np.array(K1)
    k2 = np.linalg.inv(k1)
    K2 = np.asarray(k2).reshape(size, size)
    return K2


def solve1(x, ff):  # 直接暴算的 Neumann 法...
    c = [0 for i in range(size)]  # result
    tj = [0 for i in range(size)]  # 网格点
    Kij = [[0] * size for i in range(size)]

    c0 = [0 for i in range(size)]
    c1 = [0 for i in range(size)]
    c2 = [0 for i in range(size)]
    c3 = [0 for i in range(size)]
    c4 = [0 for i in range(size)]
    c5 = [0 for i in range(size)]
    c6 = [0 for i in range(size)]
    c7 = [0 for i in range(size)]
    c8 = [0 for i in range(size)]
    c9 = [0 for i in range(size)]


    for i in range(0, size):
        tj[i] = R1 + i * (R2 - R1) / size
    for i in range(0, size):
        for j in range(0, size):
            Kij[i][j] = wj * K(tj[i], tj[j])
    fg = ff
    for j in range(ite):  # iteration :ite = 100
        for i in range(size):
            c[i] += fg[i]
        if j == 0:
            for i in range(size):
                c0[i] += c[i]
        elif j == 2:
            for i in range(size):
                c1[i] += c[i]
        elif j == 5:
            for i in range(size):
                c2[i] += c[i]
        elif j == 10:
            for i in range(size):
                c3[i] += c[i]
        elif j == 20:
            for i in range(size):
                c4[i] += c[i]
        elif j == 50:
            for i in range(size):
                c5[i] += c[i]
        elif j == 90:
            for i in range(size):
                c6[i] += c[i]
        elif j == 50:
            for i in range(size):
                c7[i] += c[i]
        elif j == 70:
            for i in range(size):
                c8[i] += c[i]
        elif j == 90:
            for i in range(size):
                c9[i] += c[i]

        fg = Multiply(Kij, fg)
    return c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c


def solve2(x, ff):  # 用到了 (E-K) 矩阵, 不用那么麻烦的...
    cc = [0 for i in range(size)]  # result
    tj = [0 for i in range(size)]  # 网格点
    Kij = [[0] * size for i in range(size)]
    KK = [[0] * size for i in range(size)]
    for i in range(0, size):
        tj[i] = R1 + i * (R2 - R1) / size
    for i in range(0, size):
        for j in range(0, size):
            Kij[i][j] = wj * K(tj[i], tj[j])
    for i in range(size):
        for j in range(size):
            if (i == j):
                KK[i][j] = 1.0 - Kij[i][j]
            else:
                KK[i][j] = -Kij[i][j]
    KK = InverseMatrix(KK)
    cc = Multiply(KK, ff)
    return cc


if __name__ == '__main__':
    goal = [0 for i in range(size)]
    c = [0 for i in range(size)]
    x = [0 for i in range(size)]
    ff = [0 for i in range(size)]
    cc1 = [0 for i in range(size)]
    cc = [0 for i in range(size)]
    cc2 = [0 for i in range(size)]

    for i in range(0, size):
        x[i] = R1 + (R2 - R1) * i / size
        goal[i] = Fx(x[i])
        ff[i] = f(x[i])

    c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c= solve1(x, ff)
    #cc = solve2(x, ff)

    for i in range(0, size):
        c[i] += 1
        cc[i] += 1
        c0[i] += 1
        c1[i] += 1
        c2[i] += 1
        c3[i] += 1
        c4[i] += 1
        c5[i] += 1
        c6[i] += 1
        c7[i] += 1
        c8[i] += 1
        c9[i] += 1

    for i in range(0, size):
        cc1[i] = math.fabs(goal[i] - c[i]) / goal[i]

    for i in range(0, size):
        cc2[i] = math.fabs(goal[i] - cc[i]) / goal[i]

    # plt.title("")
    # plt.xlabel('x')
    # plt.ylabel('f(x)')


    def to_percent(temp, position):
        return '%1.2f' % (100 * temp) + '%'


    l=plt.plot(x[(int)(size / 2):size - 1], goal[(int)(size / 2):size - 1],label="exact solution")
    ll=plt.plot(x[(int)(size / 2):size - 1], c[(int)(size / 2):size - 1],label="Neumann series solution")

    l0=plt.plot(x[(int)(size / 2):size - 1], c0[(int)(size / 2):size - 1],label="0 iteration")
    l1 =plt.plot(x[(int)(size / 2):size - 1], c1[(int)(size / 2):size - 1],label="2 iteration")
    l2 =plt.plot(x[(int)(size / 2):size - 1], c2[(int)(size / 2):size - 1],label="5 iteration")
    l3= plt.plot(x[(int)(size / 2):size - 1], c3[(int)(size / 2):size - 1],label="10 iteration")
    l4 =plt.plot(x[(int)(size / 2):size - 1], c4[(int)(size / 2):size - 1],label="20 iteration")
    l5 =plt.plot(x[(int)(size / 2):size - 1], c5[(int)(size / 2):size - 1],label="50 iteration")
    l6 =plt.plot(x[(int)(size / 2):size - 1], c6[(int)(size / 2):size - 1],label="90 iteration")
    #l7= plt.plot(x[(int)(size / 2):size - 1], c7[(int)(size / 2):size - 1],label="50 iteration")
    #l8 = plt.plot(x[(int)(size / 2):size - 1], c8[(int)(size / 2):size - 1],label="70 iteration")
    #l9 = plt.plot(x[(int)(size / 2):size - 1], c9[(int)(size / 2):size - 1],label="90 iteration")
    plt.legend(loc='upper right')
    plt.show()


    plt.title("relative error in %")
    plt.axis([x[(int)(size/2)],x[size-1],0,0.01])
    #plt.plot(x[(int)(size / 2):size - 1], goal[(int)(size / 2):size - 1])
    plt.plot(x[(int)(size / 2):size - 1], cc1[(int)(size / 2):size - 1])

    plt.gca().yaxis.set_major_formatter(FuncFormatter(to_percent))
    plt.show()



    l1 = plt.plot(x[(int)(size / 2):size - 1], goal[(int)(size / 2):size - 1], label="exact solution")
    l3 = plt.plot(x[(int)(size / 2):size - 1], cc[(int)(size / 2):size - 1], label="Nystrom method solution")

    plt.legend(loc='upper right')
    plt.show()


    plt.title("relative error in %")
    plt.axis([x[(int)(size / 2)], x[size - 1], 0, 0.005])
    # plt.plot(x[(int)(size / 2):size - 1], goal[(int)(size / 2):size - 1])
    plt.plot(x[(int)(size / 2):size - 1], c2[(int)(size / 2):size - 1])

    def to_percent(temp, position):
        return '%1.2f' % (100 * temp) + '%'
    plt.gca().yaxis.set_major_formatter(FuncFormatter(to_percent))
    plt.show()

'''
    ax1=plt.subplot(121)
    ax2=plt.subplot(122)

    def to_percent(temp, position):
        return '%1.2f' % (100 * temp) + '%'

    ax1.set_title("numerical solution")
    l1=ax1.plot(x[(int)(size / 2):size - 1], goal[(int)(size / 2):size - 1],label="exact solution")
    l2=ax1.plot(x[(int)(size / 2):size - 1], c[(int)(size / 2):size - 1],label="Neumann series solution")





    ax1.legend(loc='upper right')


    ax2.set_title("relative error in %")
    ax2.axis([x[(int)(size/2)],x[size-1],0,0.05])
    ax2.plot(x[(int)(size / 2):size - 1], cc1[(int)(size / 2):size - 1])

    plt.gca().yaxis.set_major_formatter(FuncFormatter(to_percent))
    plt.show()




    ax1=plt.subplot(121)
    ax2=plt.subplot(122)

    def to_percent(temp, position):
        return '%1.2f' % (100 * temp) + '%'

    ax1.set_title("numerical solution")
    l1=ax1.plot(x[(int)(size / 2):size - 1], goal[(int)(size / 2):size - 1],label="exact solution")
    l2=ax1.plot(x[(int)(size / 2):size - 1], cc[(int)(size / 2):size - 1],label="Nystrom method solution")
    ax1.legend(loc='upper right')


    ax2.set_title("relative error in %")
    ax2.axis([x[(int)(size/2)],x[size-1],0,0.005])
    ax2.plot(x[(int)(size / 2):size - 1], cc2[(int)(size / 2):size - 1])

    plt.gca().yaxis.set_major_formatter(FuncFormatter(to_percent))
    plt.show()

'''