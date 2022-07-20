
import math
import matplotlib.pyplot as plt



def f(t,x,y):
    a = y
    return a


def g(t,x,y):

    a = -math.sin(x)
    return a


def RK4(t, x, y, h):
    """
    :param t: t 的初始值
    :param x: x 的初始值
    :param y: y 的初始值
    :param h: 时间步长
    :return: 迭代新解
    """
    tarray, xarray, yarray = [], [], []
    while t <= 4*math.pi:
        tarray.append(t)
        xarray.append(x)
        yarray.append(y)
        t += h

        K_1 = f(t, x, y)
        L_1 = g(t, x, y)
        K_2 = f(t + h / 2, x + h / 2 * K_1, y + h / 2 * L_1)
        L_2 = g(t + h / 2, x + h / 2 * K_1, y + h / 2 * L_1)
        K_3 = f(t + h / 2, x + h / 2 * K_2, y + h / 2 * L_2)
        L_3 = g(t + h / 2, x + h / 2 * K_2, y + h / 2 * L_2)
        K_4 = f(t + h, x + h * K_3, y + h * L_3)
        L_4 = g(t + h, x + h * K_3, y + h * L_3)

        x = x + (K_1 + 2 * K_2 + 2 * K_3 + K_4) * h / 6
        y = y + (L_1 + 2 * L_2 + 2 * L_3 + L_4) * h / 6
    return tarray, xarray, yarray


def main():
    tarray, xarray, yarray = RK4(0, 1, 0, 0.01)
    print("Ronge-kutta ==> numerical result".center(168))
    #print('-' * 178)
    #print("对象\\时刻\t", "t=0\t\t", "  t=0.5\t\t", "  t=1\t\t\t", "  t=1.5\t\t", "  t=2\t\t\t", " t=2.5\t\t\t",
    #      "  t=3\t\t\t", "  t=3.5\t\t", "  t=4\t\t\t", " t=4.5\t\t\t", "  t=5")
    #print('-' * 178)
    #print("x:", end='')
    #for i in range(len(xarray)):
    #    if i % 40 == 0:
    #        print("\t\t", "%8.4f" % xarray[i], end='')
    print('\n', '-' * 177)
    print("y:", end='')
    #for i in range(len(yarray)):
    #    if i % 40 == 0:
    #        print("\t\t", "%8.4f" % yarray[i], end='')
    #print('\n', '-' * 177)
    #plt.figure('Ronge-kutta ==> numerical result')
    #plt.subplot(221)

    '''
    plt.scatter(tarray, xarray, label='x_scatter', s=1, c='#DC143C', alpha=0.6)

    plt.legend()
    plt.subplot(222)

    plt.scatter(tarray, yarray, label='y_scatter', s=1, c='#DC143C', alpha=0.6)

    plt.legend()
    plt.subplot(212)
    '''

    plt.scatter(tarray, xarray,  s=1, c='#DC143C', alpha=0.6)
    plt.xlabel('t')
    plt.legend()
    plt.show()


if __name__ == "__main__":
    main()

