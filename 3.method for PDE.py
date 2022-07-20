import numpy as np
import matplotlib.pyplot as plt
from math import pi, sin, cos
from mpl_toolkits.mplot3d import Axes3D
from numpy.linalg import norm

sigma = 0

L = 1
T = 1
maxk = 50
dt = T/maxk
n = 50
dx = L/n
cond = 0.025


x = np.zeros(n + 1)
u = np.zeros((n + 1, maxk + 1))
u1 = np.zeros((n + 1, maxk + 1))
time = np.zeros(maxk + 1)
for i in range(0,n+1):
    x[i] = i * dx

for k in range(0,maxk+1):
    time[k] = k * dt

for i in range(0,n+1):
    for j in range(0,maxk+1):
        u1[i,j]=sin(2*pi*x[i])*cos(2*pi*time[j])


for i in range(0, n + 1):
    u[i, 0] = sin(2 * pi * x[i])
for k in range(0, maxk + 1):
    u[0, k] = 0
    u[n, k] = 0


for i in range(0,n+1):
    for j in range(0,maxk+1):
        u1[i,j]=sin(2*pi*x[i])*cos(2*pi*time[j])


b = 2 * 0.25 * dt/(dx ** 2)

for k in range(maxk):
    for i in range(1, n):
        u[i, k + 1] = u[i, k] + 0.5 * b * (u[i - 1, k] + u[i + 1, k] - 2 * u[i, k])




figure = plt.figure()
ax1 = Axes3D(figure)
x1 = np.linspace(0, 1, n+1)
x2 = np.linspace(0, 1, maxk+1)
X, Y = np.meshgrid(x1, x2)

Z = np.array(u1)

plt.xlabel('x')
plt.ylabel('t')

ax1.plot_surface(Y, X, Z, rstride=1, cstride=1, cmap='viridis')
plt.show()
