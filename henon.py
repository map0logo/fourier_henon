from __future__ import division
from numba import autojit
import matplotlib.pyplot as plt
import pylab


@autojit
def HenonMap(x, y, a, b, alpha, gamma):
    """
    Henon map iteration
    :param x: x_n coordinate
    :param y: y_n coordinate
    :param a: a parameter
    :param b: b parameter
    :param alpha: alpha parameter
    :param gamma: gamma parameter
    :return: X-n+1 and y_n+1 coordinates
    """
    return a - alpha * x**2 + b * y, gamma * x


def HenonIterations(n):
    """
    List of the first n iterations of Henon map given x_0=0 y_0=0
    :param n: number of iterations
    :return: list of iterations
    """
    iterations = []
    a = 1.4
    b = 0.3
    alpha = 1.
    gamma = 1.
    x = 0.
    y = 0.
    iterations.append((x, y))
    for i in xrange(n):
        x, y = HenonMap(x, y, a, b, alpha, gamma)
        iterations.append((x, y))
    return iterations


def HenonPlot():
    """
    Plot of Henon Attractor
    :return:
    """
    iterations = HenonIterations(1000)
    x = [i for i, j in iterations]
    y = [j for i, j in iterations]
    T = range(len(iterations))
    pylab.figure()
    plt.scatter(x, y, c=T, marker='+', linewidth='1')


@autojit
def H(x, y, gamma):
    """

    :param x: x_n coordinate
    :param y: y_n coordinate
    :param gamma: gamma parameter
    :return: Henon map for a = 1, b = 1, and alpha = 0.2
    """
    return HenonMap(x, y, 1., 1., 0.2, gamma)


@autojit
def HenonIterate(x0, y0, n, gamma):
    """
    Iterates on previous H map starting on x_0, y_0
    :param x0: initial x_ coordinate
    :param y0: initial y_coordinate
    :param n: number of iterations
    :param gamma: gamma parameter
    :return: number of iterations until x_i^2 + y_i^2 > 100 or n if never
             having x_i^2 + y_i^2 > 100
    """
    i = 0
    x = x0
    y = y0
    while (x**2 + y**2 <= 100) and (i < n):
        x, y = H(x, y, gamma)
        i += 1
    return i


def HenonPlot2(n, gamma):
    """
    Plot of HenonIterate inside [-5, 5] x [5, 5] square
    :param n: number of iterations
    :param: gamma: gamma parameter
    :return:
    """
    lower = -5.
    upper = 5.
    step = (float(upper) - float(lower)) / float(n)
    val = [lower + i * step for i in range(n + 1)]
    mat = []
    for i in xrange(len(val)):
        mat.append([])
        for j in xrange(len(val)):
            mat[i].append(HenonIterate(val[i], val[j], n, gamma))
    # x = len(val) * val
    # y = [j for j in val for i in xrange(len(val))]
    # T = [mat[i][j] for i in xrange(len(val)) for j in xrange(len(val))]
    pylab.figure(2)
    pylab.contourf(val, val, mat, 10)
    # pylab.scatter(x, y, c=T, alpha=0.10)
    pylab.show()


if __name__ == '__main__':
    HenonPlot()
    HenonPlot2(200, 1.03)