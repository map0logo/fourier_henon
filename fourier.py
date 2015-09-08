from __future__ import division
from numba import autojit
from math import sin, cos, pi
import matplotlib.pyplot as plt


def trapez_integration(f, xmin, xmax, n=1000):
    """
    Numerical integration using trapezoidal rule with an uniform grid
    :param f: function of x
    :param xmin: interval lower bound
    :param xmax: interval upper bound
    :param n: number of grid points
    :return: Estimated value of the definite integral of f in the interval [a, b]
    """

    h = (float(xmax) - float(xmin)) / float(n)
    xvals = [xmin + i * h for i in range(n + 1)]
    fvals = [f(x) for x in xvals]
    integral = 0.5 * h * sum([p + q for p, q in zip(fvals[1:], fvals[:-1])])
    return integral


def fourier_coefficients(f, P, n):
    """
    Calculates Fourier coefficients of f on the interval [0, P]
    :param f: function of x
    :param P: interval upper bound and interval length
    :param n: Number of Fourier coefficients to calculate
    :return: fourier coeffcientes lists a = [a_0, a_1, ..., a_n], and
             b = [0, b_1, ..., b_N]
    """

    def afun(x):
        """
        Function inside integral for a's coefficients
        :return: f x cos(2 pi k x / P)
        """

        arg = (2. * pi * k * x) / P
        return f(x) * cos(arg)

    def bfun(x):
        """
        Function inside integral for b's coefficients
        :return: f x sin(2 pi k x / P)
        """

        arg = (2. * pi * k * x) / P
        return f(x) * sin(arg)

    a = []
    b = []
    a.append((2. / P) * trapez_integration(f, 0, P))
    b.append(0.)
    for k in xrange(1, n + 1):
        a.append((2. / P) * trapez_integration(afun, 0, P))
        b.append((2. / P) * trapez_integration(bfun, 0, P))
    return a, b


@autojit
def fourier_series(x, a, b, P, n):
    """
    Fourier series for F(x) for [0, P] with n terms, using coefficients a, b
    :param x: input variable
    :param a: a's coefficients
    :param b: b's coefficients
    :param P: interval upper bound and interval length
    :param n: Number of Fourier coefficients
    :return: Fourier series evaluated on x
    """
    F = a[0] / 2.
    for k in xrange(1, n + 1):
        arg = (2. * pi * k * x) / P
        F += a[n] * cos(arg) + b[n] * sin(arg)
    return F


def fourier_plot():
    """
    Plot of Fourier transform for 1 - x
    :return:
    """
    def f(x):
        """

        :param x:
        :return:
        """
        return 1. - x

    P = 1
    xmin = -1
    xmax = 2
    n_steps = 5000
    step = (float(xmax) - float(xmin)) / float(n_steps)
    xvals = [xmin + i * step for i in range(n_steps + 1)]
    for n in [5, 10, 25, 100]:
        a, b = fourier_coefficients(f, P, n)
        Fvals = [fourier_series(x, a, b, P, n) for x in xvals]
        plt.plot(xvals, Fvals, label=str(n))
    plt.legend()
    plt.show()

if __name__ == "__main__":
    fourier_plot()